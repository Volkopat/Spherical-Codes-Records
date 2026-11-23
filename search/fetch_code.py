"""Fetch live spherical-codes.org coordinate files and verify exact 60-digit match."""

from __future__ import annotations

import re
import urllib.request
from pathlib import Path
from typing import Optional

import numpy as np

from verifier import certify_exact_decimal, load_plain_table

CACHE_DIR = Path(__file__).resolve().parent / "data" / "codes"
CACHE_DIR.mkdir(parents=True, exist_ok=True)


def fetch_cohn_code(d: int, n: int) -> Optional[np.ndarray]:
    """Download and cache the live spherical-codes.org coordinates for (d, n)."""
    cache_path = CACHE_DIR / f"{d}_{n}.txt"
    if cache_path.exists():
        raw = cache_path.read_text(encoding="utf-8")
        if raw.lstrip().startswith("<"):
            cache_path.unlink()
            raw = ""
    else:
        raw = ""

    if not raw:
        url = f"https://spherical-codes.org/coordinates/{d}/{n}"
        try:
            raw = urllib.request.urlopen(url, timeout=15).read().decode("utf-8")
            if raw.lstrip().startswith("<"):
                return None
            cache_path.write_text(raw, encoding="utf-8")
        except Exception:
            return None

    try:
        vals = [float(x) for x in raw.replace(",", " ").split()]
    except ValueError:
        return None
    if len(vals) != n * d:
        return None
    return np.array(vals, dtype=np.float64).reshape(n, d)


if __name__ == "__main__":
    html = urllib.request.urlopen("https://spherical-codes.org/", timeout=15).read().decode("utf-8", "ignore")
    links = re.findall(r'href="([^"]*)"', html)
    code_links = [u for u in links if "12" in u and "3" in u][:8]
    print("Sample links on homepage:", code_links)

    tbl = load_plain_table(Path(__file__).resolve().parent / "data" / "plain.txt")
    # Try first matching link format
    for d, n in [(3, 12), (6, 29)]:
        X = fetch_cohn_code(d, n)
        if X is not None:
            _, c12, _ = certify_exact_decimal(X)
            print(f"Verified ({d},{n}): table={tbl[(d,n)][0]} vs our_60digit_ceil12={c12}")
        else:
            print(f"Could not fetch ({d},{n}) with default URL")
