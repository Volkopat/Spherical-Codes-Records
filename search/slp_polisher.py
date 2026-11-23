"""Exact Active-Set Sequential Linear Programming (Trust-Region SLP) polisher via HiGHS."""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy import sparse
from scipy.optimize import linprog


def slp_polish(
    X_init: np.ndarray,
    *,
    max_iters: int = 45,
    init_radius: float = 0.02,
    min_radius: float = 1e-13,
    active_margin: float = 0.008,
    max_pairs: int = 900,
    lp_time_limit: float = 3.0,
    tol: float = 1e-13,
) -> Tuple[np.ndarray, float, int]:
    """
    Polish a spherical code X (n x d) to exact KKT minimax precision using Trust-Region
    Sequential Linear Programming (SLP) solved with SciPy HiGHS.

    Linearization on tangent space <x_i, dx_i> = 0:
      <x_i', x_j'> = <x_i, x_j> + <x_j, dx_i> + <x_i, dx_j> + O(||dx||^2) <= t
    """
    X = np.asarray(X_init, dtype=np.float64).copy()
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    n, d = X.shape
    nd = n * d

    def _eval(U: np.ndarray) -> Tuple[float, np.ndarray]:
        G = U @ U.T
        np.fill_diagonal(G, -2.0)
        return float(G.max()), G

    cur_c, G = _eval(X)
    radius = init_radius
    margin = active_margin
    n_accepted = 0
    n_timeouts = 0
    base_d = np.arange(d, dtype=np.int32)
    rows_eq = np.repeat(np.arange(n, dtype=np.int32), d)
    cols_eq = np.arange(nd, dtype=np.int32)
    b_eq = np.zeros(n, dtype=np.float64)
    c_obj = np.zeros(nd + 1, dtype=np.float64)
    c_obj[nd] = 1.0

    for it in range(max_iters):
        if radius < min_radius or n_timeouts >= 2:
            break

        # Gather active + near-active upper-triangular pairs
        ii, jj = np.where(np.triu(G, k=1) >= cur_c - margin)
        if len(ii) > max_pairs:
            vals = G[ii, jj]
            keep = np.argpartition(vals, -max_pairs)[-max_pairs:]
            ii, jj = ii[keep], jj[keep]
        m = len(ii)
        if m == 0:
            break

        # Vectorized sparse construction for inequality constraints:
        # <x_j, dx_i> + <x_i, dx_j> - t <= -<x_i, x_j>
        rows_ub = np.repeat(np.arange(m, dtype=np.int32), 2 * d + 1)
        cols_i = ii[:, None] * d + base_d[None, :]
        cols_j = jj[:, None] * d + base_d[None, :]
        cols_t = np.full((m, 1), nd, dtype=np.int32)
        cols_ub = np.concatenate([cols_i, cols_j, cols_t], axis=1).ravel()

        data_ub = np.concatenate([X[jj], X[ii], np.full((m, 1), -1.0, dtype=np.float64)], axis=1).ravel()
        A_ub = sparse.csr_matrix((data_ub, (rows_ub, cols_ub)), shape=(m, nd + 1))
        b_ub = -G[ii, jj]

        A_eq = sparse.csr_matrix((X.ravel(), (rows_eq, cols_eq)), shape=(n, nd + 1))
        bounds = [(-radius, radius)] * nd + [(None, None)]
        res = linprog(
            c_obj,
            A_ub=A_ub,
            b_ub=b_ub,
            A_eq=A_eq,
            b_eq=b_eq,
            bounds=bounds,
            method="highs",
            options={"time_limit": lp_time_limit},
        )
        if not res.success or res.x is None:
            if res.status == 1:  # time/iteration limit
                n_timeouts += 1
            radius *= 0.5
            continue
        n_timeouts = 0

        pred_t = float(res.x[nd])
        pred_improv = cur_c - pred_t
        if pred_improv <= tol:
            # Check if we're already at tight active set
            if margin > 1e-5:
                margin = max(margin * 0.25, 1e-5)
                continue
            break

        dx = res.x[:nd].reshape(n, d)
        X_new = X + dx
        X_new /= np.linalg.norm(X_new, axis=1, keepdims=True)
        new_c, G_new = _eval(X_new)

        if new_c < cur_c - 1e-15:
            actual_improv = cur_c - new_c
            X, cur_c, G = X_new, new_c, G_new
            n_accepted += 1
            ratio = actual_improv / max(pred_improv, 1e-16)
            if ratio > 0.75:
                radius = min(radius * 1.5, 0.05)
            elif ratio < 0.25:
                radius *= 0.5
            # Shrink active margin as we converge
            margin = max(min(margin, max(20.0 * radius, 5e-4)), 5e-5)
        else:
            # Overshot or an inactive pair became active -> shrink trust radius
            radius *= 0.4
            margin = min(margin * 1.5, 0.05)

    return X, cur_c, n_accepted
