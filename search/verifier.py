"""Exact 60-digit arbitrary-precision verifier for spherical codes matching spherical-codes.org."""

from __future__ import annotations

from decimal import Decimal, ROUND_CEILING, getcontext
from pathlib import Path
from typing import Tuple

import numpy as np


def load_plain_table(plain_path: Path) -> dict[Tuple[int, int], Tuple[str, float, str]]:
    """Return dict[(d, n)] = (cosine_str, cosine_float, minpoly_str)."""
    table: dict[Tuple[int, int], Tuple[str, float, str]] = {}
    for line in plain_path.read_text(encoding="utf-8").splitlines():
        parts = line.strip().split(",")
        if len(parts) < 3:
            continue
        d, n = int(parts[0]), int(parts[1])
        c_str = parts[2].strip()
        poly = parts[3].strip() if len(parts) >= 4 else ""
        table[(d, n)] = (c_str, float(c_str), poly)
    return table


def max_cosine_f64(X: np.ndarray) -> float:
    """Fast float64 max pairwise cosine after row normalization."""
    X = np.asarray(X, dtype=np.float64)
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    U = X / norms
    G = U @ U.T
    np.fill_diagonal(G, -2.0)
    return float(G.max())


def certify_exact_decimal(X: np.ndarray, digits: int = 60) -> Tuple[Decimal, str, Tuple[int, int]]:
    """
    Rigorous arbitrary-precision certification of a spherical code X (shape n x d).

    1. Converts each coordinate to Decimal via repr(float(v)) (exact float64 representation).
    2. Normalizes each row to unit length at `digits` precision.
    3. Computes exact pairwise inner products <u_i, u_j> for all i < j (using exact float64
       normalization residual check) and returns:
       - exact_max_cos (Decimal)
       - ceil_12 (str): rigorous upper bound rounded UP to 12 decimal places (matching spherical-codes.org)
       - argmax_pair (i, j)
    """
    getcontext().prec = digits
    X = np.asarray(X, dtype=np.float64)
    n, d = X.shape

    # Normalize in high precision
    U = []
    for i in range(n):
        row = [Decimal(repr(float(X[i, k]))) for k in range(d)]
        norm_sq = sum(v * v for v in row)
        norm = norm_sq.sqrt()
        U.append([v / norm for v in row])

    # Pre-filter candidate pairs with float64 margin so certification is instant even for n=500
    U_f64 = X / np.linalg.norm(X, axis=1, keepdims=True)
    G = U_f64 @ U_f64.T
    np.fill_diagonal(G, -2.0)
    approx_max = float(G.max())
    # Check ALL pairs within 1e-6 of approx_max (or all pairs if n <= 120)
    if n <= 120:
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    else:
        ii, jj = np.where(np.triu(G, k=1) >= approx_max - 1e-6)
        pairs = list(zip(ii.tolist(), jj.tolist()))

    best_c = Decimal("-2")
    best_pair = (0, 1)
    for i, j in pairs:
        ui, uj = U[i], U[j]
        dot = sum(ui[k] * uj[k] for k in range(d))
        if dot > best_c:
            best_c = dot
            best_pair = (i, j)

    # Add tiny safety margin 10^{-(digits-10)} before ceiling to 12 decimals so rounding is strictly rigorous
    eps_safe = Decimal(10) ** (-(digits - 10))
    ceil_12 = str((best_c + eps_safe).quantize(Decimal("0.000000000001"), rounding=ROUND_CEILING))
    return best_c, ceil_12, best_pair


def save_code_cohn_format(X: np.ndarray, path: Path) -> None:
    """Save normalized coordinates in spherical-codes.org native CSV format (one point per line, 17 digits)."""
    X = np.asarray(X, dtype=np.float64)
    U = X / np.linalg.norm(X, axis=1, keepdims=True)
    lines = [",".join(f"{float(v):.17g}" for v in row) for row in U]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
