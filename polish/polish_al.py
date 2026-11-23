"""Polish a spherical code to a local optimum of the minimal-angle problem.

One code:   python -I polish_al.py <in_file> <out_file> [--rho0 R] [--rho-max R] [--outer N] [--budget S]
A folder:   python -I polish_al.py --batch <in_dir> <out_dir> [--workers W] [--rho0s 1e2,1e3,1e4,1e5,1e6]

Problem: minimise t subject to <x_i, x_j> <= t for all pairs i < j, with |x_i| = 1.
Method: the classical method of multipliers (Powell-Hestenes-Rockafellar augmented Lagrangian) over all pairs,

    L(X, t; lam) = t + (1 / (2 rho)) * sum_{i<j} ( max(0, lam_ij + rho (<x_i, x_j> - t))^2 - lam_ij^2 ),

minimised over X (rows normalised inside the objective) and t with L-BFGS, followed by the multiplier update
lam <- max(0, lam + rho (G - t)). rho is raised when the infeasibility stops shrinking.
Different starting values of rho relax the code into different nearby local optima; batch mode tries several
and keeps the best. Everything is float64 on the CPU; results are to be certified separately in exact arithmetic.
"""
from __future__ import annotations

import argparse
import os
import re
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


def load(path):
    X = np.array([[float(t) for t in ln.split(",")] for ln in Path(path).read_text(encoding="utf-8").splitlines() if ln.strip()])
    return X / np.linalg.norm(X, axis=1, keepdims=True)


def save(path, X):
    tmp = str(path) + ".tmp"
    Path(tmp).write_text("\n".join(",".join(repr(float(v)) for v in row) for row in X) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def tight(g):
    m = g.max()
    return [int((g > m - t).sum()) for t in (1e-12, 1e-9, 1e-6)]


def equalise(X, ii, jj, iters=8):
    """Minimum-norm corrections that bring <x_i, x_j> to one common value on the given pairs, staying on the spheres.
    The common value is itself an unknown, so the system stays consistent when the pairs over-determine the code."""
    from scipy.sparse import csr_matrix, hstack, vstack
    from scipy.sparse.linalg import lsmr

    n, d = X.shape
    m = len(ii)
    base = np.arange(d)
    rows = np.repeat(np.arange(m), 2 * d)
    cols = np.concatenate([ii[:, None] * d + base[None, :], jj[:, None] * d + base[None, :]], axis=1).ravel()
    erow, ecol = np.repeat(np.arange(n), d), np.arange(n * d)
    minus1, zero = csr_matrix(-np.ones((m, 1))), csr_matrix((n, 1))
    best_X, best_s = X, np.inf
    for _ in range(iters):
        g = np.einsum("ij,ij->i", X[ii], X[jj])
        s = float(g.max() - g.min())
        if s < best_s:
            best_X, best_s = X, s
        if best_s < 5e-16:
            break
        J = csr_matrix((np.concatenate([X[jj], X[ii]], axis=1).ravel(), (rows, cols)), shape=(m, n * d))
        E = csr_matrix((X.ravel(), (erow, ecol)), shape=(n, n * d))
        A = vstack([hstack([J, minus1]), hstack([E, zero])], format="csr")
        sol = lsmr(A, np.concatenate([g.mean() - g, np.zeros(n)]), atol=1e-16, btol=1e-16, conlim=1e14, maxiter=50000)[0]
        X = X + sol[:-1].reshape(n, d)
        X = X / np.linalg.norm(X, axis=1, keepdims=True)
    return best_X, best_s


def polish(X, rho0=1e3, rho_max=1e9, outer_max=300, inner=3000, budget=600.0, log=None):
    """Returns (best_X, best_value, info). best_value is the float64 maximal cosine of best_X."""
    X = X / np.linalg.norm(X, axis=1, keepdims=True)
    n, d = X.shape
    iu = np.triu_indices(n, 1)
    start = float((X @ X.T)[iu].max())
    best, best_X = start, X.copy()
    lam = np.zeros(len(iu[0]))
    rho, t = rho0, start
    W = np.zeros((n, n))

    def fg(z, lam, rho):
        U = z[:-1].reshape(n, d)
        nrm = np.linalg.norm(U, axis=1, keepdims=True)
        Xn = U / nrm
        w = np.maximum(0.0, lam + rho * ((Xn @ Xn.T)[iu] - z[-1]))
        f = z[-1] + (w @ w - lam @ lam) / (2.0 * rho)
        W[:] = 0.0
        W[iu] = w
        gX = (W + W.T) @ Xn
        gU = (gX - (gX * Xn).sum(axis=1, keepdims=True) * Xn) / nrm
        return f, np.concatenate([gU.ravel(), [1.0 - w.sum()]])

    t0 = time.perf_counter()
    prev_feas, calm, reason, resid = np.inf, 0, "outer limit", np.inf
    for outer in range(1, outer_max + 1):
        if time.perf_counter() - t0 > budget:
            reason = "budget"
            break
        res = minimize(fg, np.concatenate([X.ravel(), [t]]), args=(lam, rho), jac=True, method="L-BFGS-B",
                       options={"maxiter": inner, "maxfun": 4 * inner, "maxcor": 40, "ftol": 0.0, "gtol": 1e-14})
        U = res.x[:-1].reshape(n, d)
        X = U / np.linalg.norm(U, axis=1, keepdims=True)
        t = float(res.x[-1])
        g = (X @ X.T)[iu]
        cur = float(g.max())
        feas = cur - t                                   # how far the largest cosine sits above t
        lam = np.maximum(0.0, lam + rho * (g - t))
        act = lam > 0
        spread = float(np.abs(g[act] - t).max()) if act.any() else 0.0   # among pairs that carry a multiplier
        resid = max(feas, spread)
        # stationarity: tangential gradient of sum_ij lam_ij <x_i, x_j> (the multipliers sum to one)
        W[:] = 0.0
        W[iu] = lam
        gX = (W + W.T) @ X
        stat = float(np.abs(gX - (gX * X).sum(axis=1, keepdims=True) * X).max()) / max(lam.sum(), 1e-300)
        if cur < best:
            best, best_X = cur, X.copy()
        if log and (outer <= 3 or outer % 20 == 0):
            log(f"  outer {outer:3d} rho {rho:.1e} inner {res.nit:5d} t {t:.15f} max {cur:.15f} resid {resid:.1e} stat {stat:.1e} active {int(act.sum())} tight {tight(g)} {time.perf_counter() - t0:.0f}s")
        calm = calm + 1 if resid < 3e-14 else 0
        if calm >= 3:
            reason = "converged"
            break
        if feas > 0.25 * prev_feas and rho < rho_max:
            rho = min(rho * 4.0, rho_max)
        prev_feas = feas
    # Final step: move the pairs that carry a multiplier onto one common value (minimum-norm Gauss-Newton),
    # which removes the last ~1e-13 of scatter left by the quasi-Newton inner solves.
    act = lam > 0
    eq = ""
    if act.sum() > 0:
        Xe, se = equalise(X, iu[0][act], iu[1][act])
        ce = float((Xe @ Xe.T)[iu].max())
        eq = f" equalise: spread {se:.0e} max {ce:.15f}"
        if ce < best:
            best, best_X = ce, Xe
            reason += "+equalised"
    if log:
        log(f"  final stat {stat:.1e} resid {resid:.1e}{eq}")
    gb = (best_X @ best_X.T)[iu]
    return best_X, best, {"start": start, "outer": outer, "rho": rho, "resid": resid, "stat": stat, "reason": reason, "tight": tight(gb),
                          "time": time.perf_counter() - t0, "iso": n * (d - 1) - d * (d - 1) // 2 + 1}


def _task(args):
    path, rho0, rho_max, budget = args
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    X = load(path)
    bX, b, info = polish(X, rho0=rho0, rho_max=rho_max, budget=budget)
    return path.name, rho0, b, info, bX


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--batch", action="store_true")
    ap.add_argument("--budget", type=float, default=600.0)
    ap.add_argument("--rho0", type=float, default=1e3)
    ap.add_argument("--rho0s", default="1e2,1e3,1e4,1e5,1e6")
    ap.add_argument("--rho-max", type=float, default=1e9)
    ap.add_argument("--outer", type=int, default=300)
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()

    if not a.batch:
        X = load(a.inp)
        n, d = X.shape
        bX, b, info = polish(X, rho0=a.rho0, rho_max=a.rho_max, outer_max=a.outer, budget=a.budget, log=lambda s: print(s, flush=True))
        save(a.out, bX)
        print(f"({d}, {n}) start {info['start']:.15f} end {b:.15f} drop {info['start'] - b:.3e} tight {info['tight']} (isostatic {info['iso']}) "
              f"outer {info['outer']} rho {info['rho']:.0e} resid {info['resid']:.1e} time {info['time']:.0f}s stop: {info['reason']}", flush=True)
        return

    from multiprocessing import Pool

    inp, out = Path(a.inp), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    files = sorted(inp.glob("d*_n*.txt"), key=lambda p: -int(re.match(r"d(\d+)_n(\d+)", p.name)[1]) * int(re.match(r"d(\d+)_n(\d+)", p.name)[2]))
    rho0s = [float(x) for x in a.rho0s.split(",")]
    tasks = [(p, r, a.rho_max, a.budget) for p in files for r in rho0s]
    best = {}
    t0 = time.perf_counter()
    with Pool(a.workers) as pool:
        for k, (name, rho0, b, info, bX) in enumerate(pool.imap_unordered(_task, tasks), 1):
            print(f"[{k}/{len(tasks)} {time.perf_counter() - t0:.0f}s] {name} rho0 {rho0:.0e}: {info['start']:.13f} -> {b:.13f} drop {info['start'] - b:.2e} "
                  f"tight {info['tight']} resid {info['resid']:.0e} {info['reason']} {info['time']:.0f}s", flush=True)
            if name not in best or b < best[name][0]:
                best[name] = (b, rho0, info)
                save(out / name, bX)
    print("\nbest per code")
    for name in sorted(best):
        b, rho0, info = best[name]
        print(f"{name:28s} {info['start']:.13f} -> {b:.13f} drop {info['start'] - b:.2e} rho0 {rho0:.0e} tight {info['tight']} iso {info['iso']} resid {info['resid']:.0e} {info['reason']}")


if __name__ == "__main__":
    main()
