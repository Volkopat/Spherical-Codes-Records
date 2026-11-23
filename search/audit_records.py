"""Completely independent, zero-shortcut, brute-force ALL-PAIRS 80-digit audit of new spherical code records."""

from __future__ import annotations

import hashlib
import json
import math
import urllib.request
from decimal import Decimal, ROUND_CEILING, getcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REC_DIR = ROOT / "records"


def brute_force_certify(csv_path: Path, d: int, n: int, prec: int = 80) -> dict:
    """Read raw saved CSV file and brute-force check ALL n*(n-1)/2 pairs at `prec`-digit Decimal precision."""
    getcontext().prec = prec
    raw = csv_path.read_text(encoding="utf-8")
    sha256 = hashlib.sha256(raw.encode("utf-8")).hexdigest().upper()

    lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]
    assert len(lines) == n, f"Expected {n} points, got {len(lines)}"

    U = []
    max_norm_err = Decimal("0")
    for idx, ln in enumerate(lines):
        toks = [t.strip() for t in ln.split(",") if t.strip()]
        assert len(toks) == d, f"Point {idx}: expected {d} coords, got {len(toks)}"
        row = [Decimal(t) for t in toks]
        norm_sq = sum(v * v for v in row)
        norm = norm_sq.sqrt()
        err = abs(norm - Decimal("1"))
        if err > max_norm_err:
            max_norm_err = err
        U.append([v / norm for v in row])

    best_c = Decimal("-2")
    best_pair = (-1, -1)
    n_pairs = 0
    for i in range(n):
        ui = U[i]
        for j in range(i + 1, n):
            uj = U[j]
            dot = sum(ui[k] * uj[k] for k in range(d))
            if dot > best_c:
                best_c = dot
                best_pair = (i, j)
            n_pairs += 1

    assert n_pairs == n * (n - 1) // 2
    eps_safe = Decimal(10) ** (-(prec - 15))
    ceil_12 = str((best_c + eps_safe).quantize(Decimal("0.000000000001"), rounding=ROUND_CEILING))
    angle_deg = math.degrees(math.acos(float(best_c)))
    return {
        "d": d,
        "n": n,
        "file": str(csv_path),
        "sha256": sha256,
        "n_pairs_checked": n_pairs,
        "max_saved_norm_err": f"{max_norm_err:.3e}",
        "exact_max_cos_30d": f"{best_c:.30f}",
        "rigorous_ceil_12": ceil_12,
        "min_angle_deg": round(angle_deg, 8),
        "argmax_pair": list(best_pair),
    }


def main():
    # Fetch fresh live table directly from spherical-codes.org
    live_txt = urllib.request.urlopen("https://spherical-codes.org/plain", timeout=20).read().decode("utf-8")
    live_tbl = {}
    for ln in live_txt.splitlines():
        parts = ln.strip().split(",")
        if len(parts) >= 3:
            live_tbl[(int(parts[0]), int(parts[1]))] = parts[2].strip()

    targets = [(11, 47), (9, 97), (10, 133), (7, 149)]
    report = []
    print("================ INDEPENDENT BRUTE-FORCE ALL-PAIRS 80-DIGIT AUDIT ================", flush=True)
    for d, n in targets:
        csv_path = REC_DIR / f"d{d}_n{n}_NEW_RECORD.txt"
        cert = brute_force_certify(csv_path, d, n, prec=80)
        live_str = live_tbl[(d, n)]
        live_np1 = live_tbl.get((d, n + 1), "N/A")
        live_dec = Decimal(live_str)
        our_dec = Decimal(cert["rigorous_ceil_12"])
        improv = live_dec - our_dec
        old_angle = math.degrees(math.acos(float(live_dec)))
        angle_gain = cert["min_angle_deg"] - old_angle

        assert our_dec < live_dec, f"AUDIT FAILED for ({d},{n}): {our_dec} >= {live_dec}"
        if live_np1 != "N/A":
            assert our_dec < Decimal(live_np1), f"Not better than n+1 for ({d},{n})"

        cert["live_table_ceil_12"] = live_str
        cert["live_table_np1_ceil_12"] = live_np1
        cert["cosine_improvement"] = f"{improv:.12f}"
        cert["old_angle_deg"] = round(old_angle, 8)
        cert["angle_gain_deg"] = round(angle_gain, 8)
        cert["status"] = "CERTIFIED_WORLD_RECORD"
        report.append(cert)

        print(
            f"  PASS (d={d:2d}, n={n:3d}): ALL {cert['n_pairs_checked']:5d} pairs checked at 80 digits | "
            f"live={live_str} -> OUR={cert['rigorous_ceil_12']} "
            f"(delta_cos=-{improv:.12f}, delta_angle=+{angle_gain:.6f} deg, norm_err={cert['max_saved_norm_err']})",
            flush=True,
        )

    out = REC_DIR / "audit_report_80digit.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nALL {len(report)}/4 WORLD RECORDS INDEPENDENTLY CERTIFIED AT 80-DIGIT PRECISION -> {out}", flush=True)


if __name__ == "__main__":
    main()
