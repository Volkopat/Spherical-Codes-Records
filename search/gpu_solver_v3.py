"""Wave 3 GPU Engine (Math KB arXiv:2607.17933 CARDAL + OpenMathOS Constructor/Referee Synthesis):
  1. Prong T  : Dimensional Tunneling (Rank-Adaptive d -> d+r -> d PHR Augmented Lagrangian)
  2. Prong DR : Batched GPU Relaxed Douglas-Rachford Elliptope Reflection (C_tau <-> K_d)
  3. Prong AL : Pure-GPU Rockafellar PHR Primal-Dual Augmented Lagrangian Polisher
"""

from __future__ import annotations

import math
import time
from typing import Tuple

import numpy as np
import torch


def _project_tangent(X: torch.Tensor, G: torch.Tensor) -> torch.Tensor:
    return G - (X * G).sum(dim=-1, keepdim=True) * X


def _retract(X: torch.Tensor) -> torch.Tensor:
    return X / torch.linalg.norm(X, dim=-1, keepdim=True).clamp_min(1e-15)


@torch.no_grad()
def batch_max_cos(X: torch.Tensor) -> torch.Tensor:
    n = X.shape[1]
    G = torch.bmm(X, X.transpose(1, 2))
    eye = torch.eye(n, device=X.device, dtype=torch.bool).unsqueeze(0)
    G.masked_fill_(eye, -2.0)
    return G.amax(dim=(1, 2))


def _phr_grad(
    X: torch.Tensor,
    lam: torch.Tensor,
    tau: torch.Tensor,
    rho: float,
    eye: torch.Tensor,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Exact gradient of Rockafellar PHR Augmented Lagrangian for inequality constraints G_ij - tau <= 0:
      Phi_rho(X, lam) = (1/(2*rho)) sum_{i<j} (max(0, lam_ij + rho*(G_ij - tau))^2 - lam_ij^2)
    Returns (grad_X, cur_max_cos, shifted_multiplier W = max(0, lam + rho*(G - tau))).
    """
    G = torch.bmm(X, X.transpose(1, 2))
    G.masked_fill_(eye, -2.0)
    cur_c = G.amax(dim=(1, 2))
    W = torch.clamp(lam + rho * (G - tau.view(-1, 1, 1)), min=0.0)
    W.masked_fill_(eye, 0.0)
    # Sum over i < j: d/dx_i [ (1/(2*rho)) sum_{a<b} W_{ab}^2 ] = sum_{j != i} W_{ij} x_j
    grad = torch.bmm(W, X)
    return grad, cur_c, W


@torch.no_grad()
def gpu_dimensional_tunneling(
    X_rec: np.ndarray,
    *,
    batch_size: int = 96,
    extra_r: int = 3,
    steps: int = 2200,
    seed: int = 0,
    top_k: int = 4,
    device: str = "cuda",
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Prong T (from Math KB arXiv:2607.17933 CARDAL + OpenMathOS):
    Temporarily lift X in OB(n, d) -> [U, V] in OB(n, d+r) (joint row-normalization),
    escape jammed d-dimensional contact graphs via low-rho PHR Augmented Lagrangian in R^{d+r},
    then smoothly ramp quadratic squeeze (alpha_k/2)*||V||_F^2 to return to S^{d-1}.
    """
    t0 = time.perf_counter()
    X0 = np.asarray(X_rec, dtype=np.float64)
    X0 = X0 / np.linalg.norm(X0, axis=1, keepdims=True)
    n, d = X0.shape

    g = torch.Generator(device=device)
    g.manual_seed(seed)

    U0 = torch.as_tensor(X0, dtype=torch.float64, device=device).unsqueeze(0).expand(batch_size, n, d).clone()
    # Multi-scale lift perturbations in the extra r dimensions + mild tangent kick in U
    scales = torch.linspace(0.08, 0.45, batch_size, device=device, dtype=torch.float64).view(batch_size, 1, 1)
    V0 = scales * torch.randn((batch_size, n, extra_r), generator=g, device=device, dtype=torch.float64)
    U0 = U0 + (0.15 * scales) * torch.randn((batch_size, n, d), generator=g, device=device, dtype=torch.float64)
    X = _retract(torch.cat([U0, V0], dim=-1))  # (B, n, d + extra_r)

    eye = torch.eye(n, device=device, dtype=torch.bool).unsqueeze(0)
    lam = torch.zeros((batch_size, n, n), device=device, dtype=torch.float64)
    m = torch.zeros_like(X)
    v = torch.zeros_like(X)
    b1, b2 = 0.85, 0.99

    c_init = float(batch_max_cos(torch.as_tensor(X0, dtype=torch.float64, device=device).unsqueeze(0))[0].item())
    # Ladder of target caps slightly below live record
    deltas = torch.linspace(0.0005, 0.012, batch_size, device=device, dtype=torch.float64)
    tau = torch.full((batch_size,), c_init, device=device, dtype=torch.float64) - deltas

    best_U = U0.clone()
    best_c_d = torch.full((batch_size,), 2.0, device=device, dtype=torch.float64)

    for step in range(1, steps + 1):
        frac = (step - 1) / max(1, steps - 1)
        # Keep rho modest during lift phase (0..0.45) so negative curvature in V is not masked (CARDAL pitfall)
        rho = 15.0 * ((250.0 / 15.0) ** frac)
        # Smooth squeeze of extra dimensions V -> 0 after frac > 0.30
        if frac <= 0.30:
            alpha = 0.0
        else:
            s_frac = (frac - 0.30) / 0.70
            alpha = 4.0 * (s_frac ** 2.2)

        grad, _, W = _phr_grad(X, lam, tau, rho, eye)
        if alpha > 0.0:
            grad[:, :, d:] += alpha * X[:, :, d:]

        rgrad = _project_tangent(X, grad)
        gnorm = torch.linalg.norm(rgrad.reshape(batch_size, -1), dim=-1).view(batch_size, 1, 1).clamp_min(1e-15)
        rgrad = rgrad / gnorm

        lr = 0.020 * ((3e-4 / 0.020) ** frac)
        m = b1 * m + (1.0 - b1) * rgrad
        v = b2 * v + (1.0 - b2) * (rgrad * rgrad)
        m_hat = m / (1.0 - b1**step)
        v_hat = v / (1.0 - b2**step)
        X = _retract(X - lr * m_hat / (torch.sqrt(v_hat) + 1e-12))

        if step % 20 == 0:
            lam = torch.clamp(W, max=50.0)

        # Periodically project U back to S^{d-1} and track best d-dimensional configuration
        if step >= int(0.65 * steps) and step % 25 == 0:
            U_proj = _retract(X[:, :, :d])
            c_d = batch_max_cos(U_proj)
            imp = c_d < best_c_d
            if imp.any():
                best_c_d = torch.where(imp, c_d, best_c_d)
                best_U[imp] = U_proj[imp]

    # Final projection check
    U_proj = _retract(X[:, :, :d])
    c_d = batch_max_cos(U_proj)
    imp = c_d < best_c_d
    if imp.any():
        best_c_d = torch.where(imp, c_d, best_c_d)
        best_U[imp] = U_proj[imp]

    order = torch.argsort(best_c_d)[: min(top_k, batch_size)]
    return best_U[order].cpu().numpy(), best_c_d[order].cpu().numpy(), time.perf_counter() - t0


@torch.no_grad()
def gpu_dr_elliptope_jump(
    X_rec: np.ndarray,
    *,
    batch_size: int = 64,
    steps: int = 240,
    beta_raar: float = 0.85,
    seed: int = 0,
    top_k: int = 4,
    device: str = "cuda",
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Prong DR (Batched GPU Relaxed Douglas-Rachford / RAAR in Gram Space):
    Alternates/reflects between off-diagonal cap polytope C_tau = {G_ii=1, G_ij <= tau}
    and rank-d PSD cone K_d = {G >= 0, rank(G) <= d} via batched torch.linalg.eigh.
    """
    t0 = time.perf_counter()
    X0 = np.asarray(X_rec, dtype=np.float64)
    X0 = X0 / np.linalg.norm(X0, axis=1, keepdims=True)
    n, d = X0.shape

    g = torch.Generator(device=device)
    g.manual_seed(seed)

    X_t = torch.as_tensor(X0, dtype=torch.float64, device=device)
    G0 = X_t @ X_t.T
    eye = torch.eye(n, device=device, dtype=torch.bool)
    G0_off = G0.masked_fill(eye, -2.0)
    c_ref = float(G0_off.max().item())

    # Diverse Gram initializations + ladder of target caps tau < c_ref
    noise = torch.randn((batch_size, n, n), generator=g, device=device, dtype=torch.float64)
    noise = 0.5 * (noise + noise.transpose(1, 2))
    sigmas = torch.linspace(0.02, 0.25, batch_size, device=device, dtype=torch.float64).view(batch_size, 1, 1)
    Z = G0.unsqueeze(0) + sigmas * noise

    # Target caps ladder: bracket tight and aggressive improvements
    deltas = torch.logspace(-4.0, -1.6, batch_size, device=device, dtype=torch.float64)
    tau = (c_ref - deltas).view(batch_size, 1, 1)
    eye_b = eye.unsqueeze(0).expand(batch_size, n, n)

    best_X = X_t.unsqueeze(0).expand(batch_size, n, d).clone()
    best_c = torch.full((batch_size,), 2.0, device=device, dtype=torch.float64)

    for step in range(1, steps + 1):
        # P_{C_tau}(Z): clamp off-diagonals <= tau, set diagonal = 1
        Y = torch.minimum(Z, tau)
        Y = Y.masked_fill(eye_b, 1.0)
        Y = 0.5 * (Y + Y.transpose(1, 2))

        # Reflect: R_C = 2Y - Z
        R = 2.0 * Y - Z
        R = 0.5 * (R + R.transpose(1, 2))

        # P_{K_d}(R): top-d positive eigenpairs
        evals, evecs = torch.linalg.eigh(R)
        top_vals = torch.clamp(evals[:, -d:], min=0.0)  # (B, d)
        top_vecs = evecs[:, :, -d:]  # (B, n, d)
        F = top_vecs * torch.sqrt(top_vals).unsqueeze(1)  # (B, n, d)
        W = torch.bmm(F, F.transpose(1, 2))

        # RAAR relaxed Douglas-Rachford update (prevents nonconvex limit cycles)
        Z = beta_raar * (Z + W - Y) + (1.0 - beta_raar) * Y

        if step % 10 == 0 or step == steps:
            X_cand = _retract(F)
            c_cand = batch_max_cos(X_cand)
            imp = c_cand < best_c
            if imp.any():
                best_c = torch.where(imp, c_cand, best_c)
                best_X[imp] = X_cand[imp]

    order = torch.argsort(best_c)[: min(top_k, batch_size)]
    return best_X[order].cpu().numpy(), best_c[order].cpu().numpy(), time.perf_counter() - t0


@torch.no_grad()
def gpu_phr_al_polish(
    X_init: np.ndarray,
    *,
    steps: int = 1800,
    top_k: int = 4,
    device: str = "cuda",
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Prong AL (Pure-GPU Rockafellar PHR Primal-Dual Augmented Lagrangian Polisher):
    Maintains explicit dual multipliers lam_ij >= 0 on active constraints so GPU converges
    much closer to the exact KKT minimax vertex before final CPU HiGHS SLP certification.
    """
    t0 = time.perf_counter()
    X = _retract(torch.as_tensor(X_init, dtype=torch.float64, device=device))
    B, n, _ = X.shape
    eye = torch.eye(n, device=device, dtype=torch.bool).unsqueeze(0)

    lam = torch.zeros((B, n, n), device=device, dtype=torch.float64)
    m = torch.zeros_like(X)
    v = torch.zeros_like(X)
    b1, b2 = 0.85, 0.99

    best_X = X.clone()
    best_c = batch_max_cos(X)

    for step in range(1, steps + 1):
        frac = (step - 1) / max(1, steps - 1)
        margin = 0.012 * ((5e-5 / 0.012) ** frac)
        tau = best_c - margin
        rho = 50.0 * ((2000.0 / 50.0) ** frac)
        lr = 0.010 * ((5e-5 / 0.010) ** frac)

        grad, cur_c, W = _phr_grad(X, lam, tau, rho, eye)
        imp = cur_c < best_c
        if imp.any():
            best_c = torch.where(imp, cur_c, best_c)
            best_X[imp] = X[imp]

        rgrad = _project_tangent(X, grad)
        gnorm = torch.linalg.norm(rgrad.reshape(B, -1), dim=-1).view(B, 1, 1).clamp_min(1e-15)
        rgrad = rgrad / gnorm

        m = b1 * m + (1.0 - b1) * rgrad
        v = b2 * v + (1.0 - b2) * (rgrad * rgrad)
        m_hat = m / (1.0 - b1**step)
        v_hat = v / (1.0 - b2**step)
        X = _retract(X - lr * m_hat / (torch.sqrt(v_hat) + 1e-12))

        if step % 15 == 0:
            # Normalize dual multiplier scale so it stays well-conditioned as tau tightens
            lam = torch.clamp(W, max=100.0)

    final_c = batch_max_cos(best_X)
    order = torch.argsort(final_c)[: min(top_k, B)]
    return best_X[order].cpu().numpy(), final_c[order].cpu().numpy(), time.perf_counter() - t0
