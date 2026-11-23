"""Exact check of spherical codes against a copy of the Table of Spherical Codes (rational arithmetic).

usage: python verify_exact.py <codes_dir> <live_dir>

<codes_dir> holds files named d{d}_n{n}.txt (one point per line, comma-separated coordinates).
<live_dir>/plain/plain.txt is a copy of https://spherical-codes.org/plain; optionally
<live_dir>/coords/{d}_{n}.txt holds the listed coordinates for a side-by-side comparison.

For each (d, n):
  * parse the saved file and the live reference file as exact rationals (fractions.Fraction of the decimal strings);
  * c12 = the max cosine rounded up to 12 places, proved exactly: every pair has cos <= c12 and some pair has
    cos > c12 - 1e-12 (cos = <x,y> / (|x||y|), tested without square roots);
  * compare with the live table value and with values derivable from other live entries;
  * contact structure against the listed code.
"""
from __future__ import annotations

import hashlib
import math
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 80
REC, LIVE = Path(sys.argv[1]), Path(sys.argv[2])
import re  # noqa: E402

FILES = {tuple(int(x) for x in m.groups()): p for p in REC.glob("d*_n*.txt") if (m := re.match(r"d(\d+)_n(\d+)(?:_NEW_RECORD)?\.txt$", p.name))}
TARGETS = sorted(FILES)
E12 = Fraction(1, 10**12)
SUMMARY = []


def load(path: Path, d: int, n: int):
    rows = []
    for ln in path.read_text(encoding="utf-8").splitlines():
        ln = ln.strip().strip("[](){}").rstrip(",")
        if not ln:
            continue
        toks = [t.strip() for t in ln.split(",")]
        try:
            rows.append([Fraction(t) for t in toks])
        except ValueError:
            continue
    assert len(rows) == n, (path.name, len(rows), n)
    assert all(len(r) == d for r in rows), path.name
    return rows


def cos_leq(dot, nx, ny, B):  # cos <= B, exact (B > 0)
    return dot <= 0 or dot * dot <= B * B * nx * ny


def cos_lt(dot, nx, ny, B):  # cos < B, exact (B > 0)
    return dot <= 0 or dot * dot < B * B * nx * ny


def analyse(rows):
    n = len(rows)
    nsq = [sum(v * v for v in r) for r in rows]
    X = np.array([[float(v) for v in r] for r in rows])
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    G = X @ X.T
    iu = np.triu_indices(n, 1)
    g = G[iu]
    # exact dot products only where needed: all pairs (n <= 149, cheap enough)
    dots = {}
    best = None
    for i in range(n):
        ri = rows[i]
        for j in range(i + 1, n):
            dt = sum(a * b for a, b in zip(ri, rows[j]))
            dots[(i, j)] = dt
    # 80-digit value of the max cosine
    cand = np.argsort(g)[-50:]
    mx, arg = mp.mpf(-2), None
    for k in cand:
        i, j = int(iu[0][k]), int(iu[1][k])
        dt = dots[(i, j)]
        c = mp.mpf(dt.numerator) / mp.mpf(dt.denominator) / mp.sqrt(mp.mpf(nsq[i].numerator) / mp.mpf(nsq[i].denominator) * mp.mpf(nsq[j].numerator) / mp.mpf(nsq[j].denominator))
        if c > mx:
            mx, arg = c, (i, j)
    c12 = Fraction(int(mp.ceil(mx * 10**12)), 10**12)
    all_leq = all(cos_leq(dt, nsq[i], nsq[j], c12) for (i, j), dt in dots.items())
    some_gt = any(not cos_leq(dt, nsq[i], nsq[j], c12 - E12) for (i, j), dt in dots.items())
    norm_err = max(abs(math.sqrt(float(v)) - 1.0) for v in nsq)
    return {"n": n, "nsq": nsq, "dots": dots, "max": mx, "arg": arg, "c12": c12, "c12_proved": all_leq and some_gt, "norm_err": norm_err, "G": G, "g": np.sort(g), "X": X}


def contacts(a, tol):
    G = a["G"].copy()
    np.fill_diagonal(G, -2)
    m = float(a["max"])
    C = G > m - tol
    return int(C.sum() // 2), np.sort(C.sum(1))


def fmt(fr):
    return f"{float(fr):.12f}"


# ---- live table
table = {}
for ln in (LIVE / "plain" / "plain.txt").read_text(encoding="utf-8").splitlines():
    p = ln.strip().split(",")
    if len(p) >= 3:
        try:
            table[(int(p[0]), int(p[1]))] = Fraction(p[2].strip())
        except ValueError:
            pass
dims = sorted({d for d, _ in table})
print(f"live table: {len(table)} entries, dimensions {dims[0]}..{dims[-1]}")


def derived(d, n):
    """Upper bounds for (d, n) that follow from OTHER live entries by standard constructions."""
    out = []
    for (dd, nn), c in table.items():
        if (dd, nn) == (d, n):
            continue
        if dd <= d and nn >= n:  # embed in more dimensions, delete points
            out.append((c, f"({dd},{nn}) by deletion/embedding"))
    # orthogonal union: codes in R^a and R^b, a + b = d, n1 + n2 >= n, cosine max(c1, c2, 0)
    for a in range(1, d):
        b = d - a
        for n1 in range(1, n):
            n2 = n - n1
            c1 = best_known(a, n1)
            c2 = best_known(b, n2)
            if c1 is None or c2 is None:
                continue
            out.append((max(c1, c2, Fraction(0)), f"({a},{n1}) + ({b},{n2}) orthogonal union"))
    return min(out, key=lambda t: t[0]) if out else None


def best_known(d, n):
    """Smallest live value for at least n points in at most d dimensions; trivial cases filled in."""
    if n <= 1:
        return Fraction(-1)
    if d == 1:
        return Fraction(-1) if n == 2 else None
    if n <= d + 1:
        return Fraction(-1, n - 1)  # simplex
    if n <= 2 * d:
        return Fraction(0)  # cross polytope subset
    if d == 2:
        return None  # cos(2 pi / n): irrational, and never competitive here; skipped
    key = (d, n)
    if key not in _BK:
        vals = []
        for dd in range(3, d + 1):
            ns, sm = _SUF.get(dd, ([], []))
            k = bisect.bisect_left(ns, n)
            if k < len(ns):
                vals.append(sm[k])
        _BK[key] = min(vals) if vals else None
    return _BK[key]


import bisect  # noqa: E402

_BK: dict = {}
_SUF: dict = {}
for _dd in dims:
    _ns = sorted(nn for (x, nn) in table if x == _dd)
    _sm, _run = [], None
    for _nn in reversed(_ns):
        _run = table[(_dd, _nn)] if _run is None else min(_run, table[(_dd, _nn)])
        _sm.append(_run)
    _SUF[_dd] = (_ns, _sm[::-1])


for d, n in TARGETS:
    print(f"\n================ ({d}, {n}) ================")
    new_path = FILES[(d, n)]
    raw = new_path.read_bytes()
    print(f"file {new_path.name}: {len(raw)} bytes, sha256 {hashlib.sha256(raw).hexdigest().upper()[:16]}…")
    A = analyse(load(new_path, d, n))
    T = table.get((d, n))
    print(f"new : max cos {mp.nstr(A['max'], 22)} at pair {A['arg']}; rounded up {fmt(A['c12'])} (exact proof: {A['c12_proved']}); norm error {A['norm_err']:.2e}")
    dv0 = derived(d, n)
    listed = fmt(T) if T is not None else "not listed"
    bound = min(x for x in (T, dv0[0] if dv0 else None) if x is not None)
    SUMMARY.append((d, n, listed, fmt(A["c12"]), A["c12_proved"], A["c12"] < bound, float(bound - A["c12"]), (dv0[1] if dv0 and (T is None or dv0[0] < T) else "own entry")))
    ref_path = LIVE / "coords" / f"{d}_{n}.txt"
    if T is None or not ref_path.exists():
        print(f"table value {listed}; best derivable {fmt(dv0[0]) if dv0 else 'n/a'} via {dv0[1] if dv0 else ''}; beats both: {A['c12'] < bound}")
        continue
    R = analyse(load(ref_path, d, n))
    print(f"live: max cos {mp.nstr(R['max'], 22)}; rounded up {fmt(R['c12'])} (exact proof: {R['c12_proved']}); table says {fmt(T)}; match: {R['c12'] == T}")
    strictly = all(cos_lt(dt, A["nsq"][i], A["nsq"][j], T) for (i, j), dt in A["dots"].items())
    print(f"every pair of the new code strictly below the table value (exact): {strictly};  new c12 < table: {A['c12'] < T};  drop {float(T - A['c12']):.12f}")
    print(f"angle: live {math.degrees(math.acos(float(T))):.6f} deg -> new {math.degrees(math.acos(float(A['c12']))):.6f} deg")
    dv = derived(d, n)
    if dv:
        print(f"best value derivable from other live entries: {fmt(dv[0])} via {dv[1]};  new beats it: {A['c12'] < dv[0]}")
    nb = {k: table.get(k) for k in ((d, n - 1), (d, n + 1), (d - 1, n), (d + 1, n))}
    print("neighbours: " + ", ".join(f"{k}={fmt(v) if v is not None else 'n/a'}" for k, v in nb.items()))
    # same configuration as the listed code?
    for tol in (1e-9, 1e-6):
        ca, da = contacts(A, tol)
        cr, dr = contacts(R, tol)
        print(f"contact pairs within {tol:g} of the max: new {ca} (degrees {da.min()}..{da.max()}), listed {cr} (degrees {dr.min()}..{dr.max()})")
    diff = np.abs(A["g"] - R["g"])
    print(f"sorted pairwise cosines, new against listed: largest gap {diff.max():.3e}, mean gap {diff.mean():.3e}")
    # distinct points and smallest non-contact gap
    print(f"second-level structure: {int((A['g'] > float(A['max']) - 1e-3).sum())} pairs within 1e-3 of the max (listed: {int((R['g'] > float(R['max']) - 1e-3).sum())})")

print("\n================ summary against the live table ================")
print(f"{'(d, n)':>10s} {'live value':>15s} {'saved file':>15s} {'proved':>6s} {'beats':>6s} {'margin':>10s}  tightest live bound")
for d, n, listed, c12, ok, beats, margin, via in SUMMARY:
    print(f"{f'({d}, {n})':>10s} {listed:>15s} {c12:>15s} {str(ok):>6s} {str(beats):>6s} {margin:10.3e}  {via}")
