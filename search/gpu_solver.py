"""GPU-batched Riemannian Homotopy (LogSumExp + Dynamic Hinge-Square) for Spherical Codes."""

from __future__ import annotations

import math
import time
from typing import Tuple

import numpy as np
import torch


def _project_tangent(X: torch.Tensor, G: torch.Tensor) -> torch.Tensor:
    """Project Euclidean gradient G onto tangent space of unit spheres: G - <X, G> X."""
    return G - (X * G).sum(dim=-1, keepdim=True) * X


def _retract(X: torch.Tensor) -> torch.Tensor:
    return X / torch.linalg.norm(X, dim=-1, keepdim=True).clamp_min(1e-15)


def make_initial_batch(
    n: int, d: int, batch_size: int, seed: int, device: str, dtype: torch.dtype
) -> torch.Tensor:
    """Diverse ensemble: 50% pure i.i.d. sphere, 35% antipodal-paired, 15% noisy-antipodal."""
    g = torch.Generator(device=device)
    g.manual_seed(seed)
    X = torch.randn((batch_size, n, d), generator=g, device=device, dtype=dtype)

    n_anti = int(0.35 * batch_size)
    n_noisy = int(0.15 * batch_size)
    half = n // 2
    if n_anti > 0 and half > 0:
        U = torch.randn((n_anti, half, d), generator=g, device=device, dtype=dtype)
        X[:n_anti, :half, :] = U
        X[:n_anti, half : 2 * half, :] = -U
    if n_noisy > 0 and half > 0:
        s0, s1 = n_anti, n_anti + n_noisy
        U = torch.randn((n_noisy, half, d), generator=g, device=device, dtype=dtype)
        noise = 0.15 * torch.randn((n_noisy, half, d), generator=g, device=device, dtype=dtype)
        X[s0:s1, :half, :] = U
        X[s0:s1, half : 2 * half, :] = -U + noise

    return _retract(X)


@torch.no_grad()
def batch_max_cos(X: torch.Tensor) -> torch.Tensor:
    """Compute max_{i != j} <x_i, x_j> for each configuration in batch (B, n, d)."""
    B, n, _ = X.shape
    G = torch.bmm(X, X.transpose(1, 2))
    eye = torch.eye(n, device=X.device, dtype=torch.bool).unsqueeze(0)
    G.masked_fill_(eye, -2.0)
    return G.amax(dim=(1, 2))


def _lse_grad(X: torch.Tensor, beta: float, eye: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """Exact gradient of (1/beta) log sum_{i!=j} exp(beta <x_i, x_j>) w.r.t. X."""
    G = torch.bmm(X, X.transpose(1, 2))
    G.masked_fill_(eye, -1e9)
    max_g = G.amax(dim=(1, 2), keepdim=True)
    W = torch.exp((G - max_g) * beta)
    W.masked_fill_(eye, 0.0)
    Z = W.sum(dim=(1, 2), keepdim=True).clamp_min(1e-30)
    P = W / Z  # symmetric normalized weights
    grad = 2.0 * torch.bmm(P, X)
    return grad, max_g.squeeze(-1).squeeze(-1)


def _hinge_grad(X: torch.Tensor, tau: torch.Tensor, eye: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """Exact gradient of 0.5 * sum_{i!=j} max(0, <x_i, x_j> - tau)^2 w.r.t. X."""
    G = torch.bmm(X, X.transpose(1, 2))
    G.masked_fill_(eye, -2.0)
    max_g = G.amax(dim=(1, 2))
    E = torch.clamp(G - tau.view(-1, 1, 1), min=0.0)
    grad = 2.0 * torch.bmm(E, X)
    return grad, max_g


@torch.no_grad()
def gpu_homotopy_search(
    d: int,
    n: int,
    *,
    batch_size: int = 256,
    seed: int = 0,
    lse_steps: int = 2400,
    hinge_steps: int = 2500,
    top_k: int = 16,
    device: str = "cuda",
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Two-phase GPU Riemannian Homotopy:
      Phase 1 (float32): LogSumExp Riemannian Adam annealing beta = 6 -> 80 -> 300.
      Phase 2 (float64): Dynamic Hinge-Square Riemannian Adam on top_k*2 survivors,
                         progressively tightening tau = max_cos - delta.
    Returns (top_k_X_f64, top_k_cos_f64, elapsed_s).
    """
    t0 = time.perf_counter()
    X = make_initial_batch(n, d, batch_size, seed, device, torch.float32)
    eye = torch.eye(n, device=device, dtype=torch.bool).unsqueeze(0)

    # Phase 1: Riemannian Adam on LogSumExp in float32
    m = torch.zeros_like(X)
    v = torch.zeros_like(X)
    b1, b2, eps = 0.85, 0.99, 1e-8
    beta_start, beta_end = 6.0, 220.0
    lr_start, lr_end = 0.08, 0.004

    best_X = X.clone()
    best_c = torch.full((batch_size,), 2.0, device=device, dtype=torch.float32)

    for step in range(1, lse_steps + 1):
        frac = (step - 1) / max(1, lse_steps - 1)
        beta = beta_start * ((beta_end / beta_start) ** frac)
        lr = lr_start * (0.5 * (1.0 + math.cos(math.pi * frac))) + lr_end * frac

        grad, cur_c = _lse_grad(X, beta, eye)
        improved = cur_c < best_c
        if improved.any():
            best_c = torch.where(improved, cur_c, best_c)
            best_X[improved] = X[improved]

        rgrad = _project_tangent(X, grad)
        m = b1 * m + (1.0 - b1) * rgrad
        v = b2 * v + (1.0 - b2) * (rgrad * rgrad)
        m_hat = m / (1.0 - b1**step)
        v_hat = v / (1.0 - b2**step)
        step_dir = m_hat / (torch.sqrt(v_hat) + eps)
        # Clip per-point tangent step for stability
        snorm = torch.linalg.norm(step_dir, dim=-1, keepdim=True).clamp_min(1e-12)
        step_dir = step_dir * torch.clamp(0.25 / snorm, max=1.0)
        X = _retract(X - lr * step_dir)

    # Select top survivors and upgrade to float64
    n_surv = min(batch_size, max(top_k * 2, 32))
    idx = torch.argsort(best_c)[:n_surv]
    X64 = _retract(best_X[idx].to(torch.float64))

    # Phase 2: Dynamic Hinge-Square Riemannian Adam in float64
    m64 = torch.zeros_like(X64)
    v64 = torch.zeros_like(X64)
    best_X64 = X64.clone()
    best_c64 = batch_max_cos(X64)

    for step in range(1, hinge_steps + 1):
        frac = (step - 1) / max(1, hinge_steps - 1)
        # Tightening margin below current best cosine
        margin = 0.015 * ((2e-4 / 0.015) ** frac)
        tau = best_c64 - margin
        lr = 0.012 * ((1.5e-4 / 0.012) ** frac)

        grad, cur_c = _hinge_grad(X64, tau, eye)
        improved = cur_c < best_c64
        if improved.any():
            best_c64 = torch.where(improved, cur_c, best_c64)
            best_X64[improved] = X64[improved]

        rgrad = _project_tangent(X64, grad)
        # Normalize per-batch gradient scale so Adam sees consistent magnitudes
        gnorm = torch.linalg.norm(rgrad.reshape(n_surv, -1), dim=-1).view(n_surv, 1, 1).clamp_min(1e-15)
        rgrad = rgrad / gnorm

        m64 = b1 * m64 + (1.0 - b1) * rgrad
        v64 = b2 * v64 + (1.0 - b2) * (rgrad * rgrad)
        m_hat = m64 / (1.0 - b1**step)
        v_hat = v64 / (1.0 - b2**step)
        X64 = _retract(X64 - lr * m_hat / (torch.sqrt(v_hat) + 1e-12))

    # Final ordering
    final_c = batch_max_cos(best_X64)
    order = torch.argsort(final_c)[:top_k]
    elapsed = time.perf_counter() - t0
    return (
        best_X64[order].cpu().numpy(),
        final_c[order].cpu().numpy(),
        elapsed,
    )


@torch.no_grad()
def gpu_hinge_refine(
    X_init: np.ndarray,
    *,
    hinge_steps: int = 2500,
    init_margin: float = 0.02,
    final_margin: float = 1e-4,
    lr_start: float = 0.015,
    lr_end: float = 1e-4,
    top_k: int = 8,
    device: str = "cuda",
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Refine a batch of starting configurations (B, n, d) in float64 on GPU via Dynamic Hinge-Square."""
    t0 = time.perf_counter()
    X64 = _retract(torch.as_tensor(X_init, dtype=torch.float64, device=device))
    B, n, _ = X64.shape
    eye = torch.eye(n, device=device, dtype=torch.bool).unsqueeze(0)
    b1, b2 = 0.85, 0.99
    m64 = torch.zeros_like(X64)
    v64 = torch.zeros_like(X64)
    best_X64 = X64.clone()
    best_c64 = batch_max_cos(X64)

    for step in range(1, hinge_steps + 1):
        frac = (step - 1) / max(1, hinge_steps - 1)
        margin = init_margin * ((final_margin / init_margin) ** frac)
        tau = best_c64 - margin
        lr = lr_start * ((lr_end / lr_start) ** frac)

        grad, cur_c = _hinge_grad(X64, tau, eye)
        improved = cur_c < best_c64
        if improved.any():
            best_c64 = torch.where(improved, cur_c, best_c64)
            best_X64[improved] = X64[improved]

        rgrad = _project_tangent(X64, grad)
        gnorm = torch.linalg.norm(rgrad.reshape(B, -1), dim=-1).view(B, 1, 1).clamp_min(1e-15)
        rgrad = rgrad / gnorm

        m64 = b1 * m64 + (1.0 - b1) * rgrad
        v64 = b2 * v64 + (1.0 - b2) * (rgrad * rgrad)
        m_hat = m64 / (1.0 - b1**step)
        v_hat = v64 / (1.0 - b2**step)
        X64 = _retract(X64 - lr * m_hat / (torch.sqrt(v_hat) + 1e-12))

    final_c = batch_max_cos(best_X64)
    k = min(top_k, B)
    order = torch.argsort(final_c)[:k]
    return best_X64[order].cpu().numpy(), final_c[order].cpu().numpy(), time.perf_counter() - t0


def make_basin_hop_batch(X_rec: np.ndarray, batch_size: int = 128, seed: int = 0) -> np.ndarray:
    """Create diverse basin-hopping perturbations + Voronoi-hole teleportations from a live record."""
    rng = np.random.default_rng(seed)
    X = np.asarray(X_rec, dtype=np.float64)
    X = X / np.linalg.norm(X, axis=1, keepdims=True)
    n, d = X.shape
    G = X @ X.T
    np.fill_diagonal(G, -2.0)
    c_max = float(G.max())
    # Contact degree of each point
    deg = (G >= c_max - 0.01).sum(axis=1)
    crowded = np.argsort(-deg)

    # Pre-sample Voronoi holes
    probes = rng.standard_normal((4096, d))
    probes /= np.linalg.norm(probes, axis=1, keepdims=True)
    hole_scores = (probes @ X.T).max(axis=1)  # lower = deeper hole
    best_holes = probes[np.argsort(hole_scores)[:64]]

    batch = np.repeat(X[None, :, :], batch_size, axis=0)
    sigmas = [0.02, 0.05, 0.10, 0.18, 0.30, 0.45]
    for b in range(batch_size):
        if b < batch_size // 2:
            # Multi-scale Gaussian tangent kicks
            sig = sigmas[b % len(sigmas)]
            batch[b] += sig * rng.standard_normal((n, d))
        else:
            # Voronoi-hole teleportation of 1..3 crowded points + mild kick
            n_tel = 1 + (b % 3)
            for t in range(n_tel):
                idx = int(crowded[rng.integers(0, min(12, n))])
                hole = best_holes[rng.integers(0, len(best_holes))]
                batch[b, idx] = hole + 0.05 * rng.standard_normal(d)
            batch[b] += 0.03 * rng.standard_normal((n, d))
        batch[b] /= np.linalg.norm(batch[b], axis=1, keepdims=True)
    return batch

