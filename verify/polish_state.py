"""python polish_state.py <codes_dir>: how tightly each code sits at a local optimum (float64 is enough here).

Counts the pairs whose cosine is within 1e-12, 1e-9 and 1e-6 of the largest one. At a converged local optimum
the pairs that determine the code all sit at the largest cosine, so the first two counts are large and close."""
import re
import sys
from pathlib import Path

import numpy as np

rows = []
for p in sorted(Path(sys.argv[1]).glob("d*_n*.txt")):
    d, n = (int(x) for x in re.match(r"d(\d+)_n(\d+)", p.name).groups())
    X = np.array([[float(t) for t in ln.split(",")] for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()])
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    g = (X @ X.T)[np.triu_indices(n, 1)]
    m = g.max()
    rows.append((d, n, int((g > m - 1e-12).sum()), int((g > m - 1e-9).sum()), int((g > m - 1e-6).sum())))
print(f"{'(d, n)':>10s} {'within 1e-12':>13s} {'within 1e-9':>12s} {'within 1e-6':>12s}")
for d, n, a, b, c in sorted(rows):
    print(f"{f'({d}, {n})':>10s} {a:13d} {b:12d} {c:12d}")
loose = sorted((d, n) for d, n, a, b, c in rows if b < 0.5 * c)
print(f"codes where fewer than half of the near-contacts are tight to 1e-9: {len(loose)} of {len(rows)}")
print(loose)
