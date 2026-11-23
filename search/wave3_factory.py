"""Wave 3 Compounded World-Record Factory (Math KB CARDAL arXiv:2607.17933 + OpenMathOS + Wave 2 Cascade):
  - Prong B+AL : GPU Voronoi-Hole Basin-Hopping + PHR Primal-Dual Augmented Lagrangian
  - Prong T+AL : Dimensional Tunneling (d -> d+3 -> d, CARDAL low-rho + smooth squeeze) + PHR AL
  - Prong C+AL : (d, n+1) -> (d, n) Deletion Transfer (preferring our own Wave 2 records as parents!)
  - Prong E+AL : (d, n-1) -> (d, n) 64-Hole Voronoi Insertion Transfer (preferring our own Wave 2 records!)
  - Prong F+AL : Deep Second-Round Tunneling + Basin-Hop on every newly broken record
  - Immediate 80-digit brute-force all-pairs certification against live spherical-codes.org/plain
"""

from __future__ import annotations

import json
import math
import time
import urllib.request
from decimal import Decimal
from pathlib import Path
from typing import Optional

import numpy as np
import torch

from audit_records import brute_force_certify
from fetch_code import fetch_cohn_code
from gpu_solver import make_basin_hop_batch
from gpu_solver_v3 import gpu_dimensional_tunneling, gpu_phr_al_polish
from slp_polisher import slp_polish
from turbo_factory import make_insertion_batch
from verifier import certify_exact_decimal, load_plain_table, save_code_cohn_format

ROOT = Path(__file__).resolve().parent
REC_DIR = ROOT / "records"
REC_DIR.mkdir(parents=True, exist_ok=True)


def fetch_best_available_code(d: int, n: int) -> Optional[np.ndarray]:
    """Use our own previously broken NEW_RECORD CSV if available; otherwise fetch live spherical-codes.org."""
    local_rec = REC_DIR / f"d{d}_n{n}_NEW_RECORD.txt"
    if local_rec.exists():
        vals = [float(x) for x in local_rec.read_text(encoding="utf-8").replace(",", " ").split()]
        if len(vals) == n * d:
            return np.array(vals, dtype=np.float64).reshape(n, d)
    return fetch_cohn_code(d, n)


def refresh_live_plain() -> dict:
    txt = urllib.request.urlopen("https://spherical-codes.org/plain", timeout=20).read().decode("utf-8")
    p = ROOT / "data" / "plain_wave3.txt"
    p.write_text(txt, encoding="utf-8")
    return load_plain_table(p)


def select_wave3_targets(tbl: dict, max_targets: int = 80) -> list[tuple[int, int]]:
    """
    Tier 1: Neighbors (d, n-1) and (d, n+1) of our top Wave 2 breakthroughs (cascade!).
    Tier 2: High-bump unexplored open targets across d=7..32, n=45..230.
    """
    already_broken = set()
    for f in REC_DIR.glob("d*_n*_NEW_RECORD.txt"):
        parts = f.stem.split("_")
        if len(parts) >= 2 and parts[0].startswith("d") and parts[1].startswith("n"):
            already_broken.add((int(parts[0][1:]), int(parts[1][1:])))

    done_wave12 = set(already_broken)
    for fname in ("world_records.json", "factory_records.json"):
        p = REC_DIR / fname
        if p.exists():
            for r in json.loads(p.read_text(encoding="utf-8")):
                done_wave12.add((r["d"], r["n"]))

    # Tier 1: Cascade around our largest Wave 2 breakthroughs
    picked: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    cert_path = REC_DIR / "ALL_CERTIFIED_WORLD_RECORDS_80D.json"
    if cert_path.exists():
        data = json.loads(cert_path.read_text(encoding="utf-8"))
        recs = data.get("live_world_records", data) if isinstance(data, dict) else data
        for x in recs[:18]:
            d, n = x["d"], x["n"]
            for nn in (n - 1, n + 1):
                if (d, nn) in tbl and not tbl[(d, nn)][2] and (d, nn) not in already_broken and (d, nn) not in seen:
                    if 40 <= nn <= 240:
                        picked.append((d, nn))
                        seen.add((d, nn))

    # Tier 2: Unexplored high-bump targets in d=7..32, n=45..230 (especially d=25..32 and n=181..230!)
    by_d: dict[int, dict[int, float]] = {}
    for (d, n), (_, cf, _) in tbl.items():
        by_d.setdefault(d, {})[n] = cf

    cands = []
    for d in range(7, 33):
        if d not in by_d:
            continue
        for n, c in by_d[d].items():
            if not (45 <= n <= 230):
                continue
            if tbl[(d, n)][2]:  # skip algebraic minimal polynomials
                continue
            if (d, n) in done_wave12 or (d, n) in seen:
                continue
            if (n - 1) in by_d[d] and (n + 1) in by_d[d]:
                bump = c - 0.5 * (by_d[d][n - 1] + by_d[d][n + 1])
                # Slight priority boost for completely unexplored dimensions d=25..32 and n>=181
                bonus = 0.002 if (d >= 25 or n >= 181) else 0.0
                cands.append((bump + bonus, d, n))

    cands.sort(reverse=True)
    per_d: dict[int, int] = {}
    for _, d, n in cands:
        if per_d.get(d, 0) >= 6:
            continue
        picked.append((d, n))
        seen.add((d, n))
        per_d[d] = per_d.get(d, 0) + 1
        if len(picked) >= max_targets:
            break
    return picked


def attack_wave3(d: int, n: int, tbl: dict) -> dict:
    rec_str, rec_f, _ = tbl[(d, n)]
    rec_dec = Decimal(rec_str)
    print(f"\n=== WAVE3 TARGET ({d:2d}, {n:3d}) | LIVE={rec_str} ===", flush=True)

    best_X = None
    best_dec = rec_dec
    best_c12 = rec_str
    best_prong = "live_record"
    t0 = time.perf_counter()

    def _try_candidate(X_cand: np.ndarray, prong: str, force: bool = False) -> bool:
        nonlocal best_X, best_dec, best_c12, best_prong
        U = X_cand / np.linalg.norm(X_cand, axis=1, keepdims=True)
        G = U @ U.T
        np.fill_diagonal(G, -2.0)
        if not force and float(G.max()) > float(best_dec) + 0.004:
            return False
        Xp, cp, nacc = slp_polish(X_cand, max_iters=25, max_pairs=700, lp_time_limit=1.5)
        if cp < float(best_dec) + 1e-8:
            _, c12, _ = certify_exact_decimal(Xp, digits=60)
            c12_dec = Decimal(c12)
            if c12_dec < best_dec:
                delta = float(rec_dec - c12_dec)
                best_X = Xp.copy()
                best_dec = c12_dec
                best_c12 = c12
                best_prong = prong
                save_code_cohn_format(best_X, REC_DIR / f"d{d}_n{n}_NEW_RECORD.txt")
                print(
                    f"  *** NEW WORLD RECORD ({prong})! {rec_str} -> {c12} "
                    f"(delta=-{delta:.12f}, slp_acc={nacc}) ***",
                    flush=True,
                )
                return True
        return False

    X_rec = fetch_cohn_code(d, n)
    if X_rec is not None:
        # Prong A: Direct SLP check
        _try_candidate(X_rec, "A_slp_polish", force=True)

        # Prong T+AL (Math KB CARDAL + OpenMathOS Dimensional Tunneling d -> d+3 -> d)
        Xt, ct, dt_t = gpu_dimensional_tunneling(X_rec, batch_size=96, extra_r=3, steps=2000, seed=d * 777 + n, top_k=6)
        Xtal, ctal, dt_tal = gpu_phr_al_polish(Xt, steps=1600, top_k=3)
        print(f"  Prong T+AL (Tunnel {dt_t:.1f}s + AL {dt_tal:.1f}s): gpu={ctal[0]:.10f}", flush=True)
        for k in range(len(Xtal)):
            _try_candidate(Xtal[k], f"T_tunnel_al_{k}")

        # Prong B+AL (Voronoi Basin-Hop + PHR AL)
        bh_batch = make_basin_hop_batch(X_rec, batch_size=160, seed=d * 1000 + n)
        Xs_b, cs_b, dt_b = gpu_phr_al_polish(bh_batch, steps=1800, top_k=3)
        print(f"  Prong B+AL (BasinHop+AL {dt_b:.1f}s): gpu={cs_b[0]:.10f}", flush=True)
        for k in range(len(Xs_b)):
            _try_candidate(Xs_b[k], f"B_basinhop_al_{k}")

    # Prong C+AL: (d, n+1) -> (d, n) Deletion Transfer (using OUR Wave 2/3 record if available!)
    if (d, n + 1) in tbl:
        X_np1 = fetch_best_available_code(d, n + 1)
        if X_np1 is not None:
            # Cap deletion batch to 128 most-crowded deletions if n > 128 for speed
            if n + 1 <= 128:
                del_idx = list(range(n + 1))
            else:
                G1 = X_np1 @ X_np1.T
                np.fill_diagonal(G1, -2.0)
                deg1 = (G1 >= G1.max() - 0.015).sum(axis=1)
                del_idx = np.argsort(-deg1)[:128].tolist()
            del_batch = np.stack([np.delete(X_np1, i, axis=0) for i in del_idx], axis=0)
            Xs_c, cs_c, dt_c = gpu_phr_al_polish(del_batch, steps=1800, top_k=3)
            print(f"  Prong C+AL (Del n+1 {dt_c:.1f}s): gpu={cs_c[0]:.10f}", flush=True)
            for k in range(len(Xs_c)):
                _try_candidate(Xs_c[k], f"C_del_al_{k}")

    # Prong E+AL: (d, n-1) -> (d, n) 64-Hole Voronoi Insertion Transfer (using OUR Wave 2/3 record if available!)
    if (d, n - 1) in tbl:
        X_nm1 = fetch_best_available_code(d, n - 1)
        if X_nm1 is not None:
            ins_batch = make_insertion_batch(X_nm1, batch_size=64, seed=d * 2000 + n)
            Xs_e, cs_e, dt_e = gpu_phr_al_polish(ins_batch, steps=1800, top_k=3)
            print(f"  Prong E+AL (Ins n-1 {dt_e:.1f}s): gpu={cs_e[0]:.10f}", flush=True)
            for k in range(len(Xs_e)):
                _try_candidate(Xs_e[k], f"E_ins_al_{k}")

    # Prong F+AL: Deep Second-Round Tunneling + Basin-Hop on OUR new record!
    if best_X is not None and float(rec_dec - best_dec) >= 1e-6:
        Xt2, _, _ = gpu_dimensional_tunneling(best_X, batch_size=64, extra_r=2, steps=1600, seed=d * 4000 + n, top_k=8)
        bh2 = make_basin_hop_batch(best_X, batch_size=96, seed=d * 3000 + n)
        comb = np.concatenate([Xt2, bh2], axis=0)
        Xs_f, cs_f, dt_f = gpu_phr_al_polish(comb, steps=1800, top_k=3)
        print(f"  Prong F+AL (DeepTunnel+Hop {dt_f:.1f}s): gpu={cs_f[0]:.10f}", flush=True)
        for k in range(len(Xs_f)):
            _try_candidate(Xs_f[k], f"F_deeptunnel_al_{k}")

    torch.cuda.empty_cache()
    elapsed = time.perf_counter() - t0
    broken = best_dec < rec_dec

    cert = None
    if broken:
        cert = brute_force_certify(REC_DIR / f"d{d}_n{n}_NEW_RECORD.txt", d, n, prec=80)
        assert Decimal(cert["rigorous_ceil_12"]) < rec_dec, "80-digit audit mismatch!"
        old_ang = math.degrees(math.acos(float(rec_dec)))
        cert["live_table_ceil_12"] = rec_str
        cert["cosine_improvement"] = f"{rec_dec - Decimal(cert['rigorous_ceil_12']):.12f}"
        cert["angle_gain_deg"] = round(cert["min_angle_deg"] - old_ang, 8)
        print(
            f"  [80-DIGIT AUDIT PASS] ({d},{n}): {rec_str} -> {cert['rigorous_ceil_12']} "
            f"(delta_cos=-{cert['cosine_improvement']}, +{cert['angle_gain_deg']:.6f} deg, {elapsed:.1f}s)",
            flush=True,
        )
    else:
        print(f"  No break ({d},{n}) in {elapsed:.1f}s", flush=True)

    return {
        "d": d,
        "n": n,
        "old_record_ceil12": rec_str,
        "new_record_ceil12": best_c12,
        "broken": broken,
        "improvement": float(rec_dec - best_dec),
        "winning_prong": best_prong,
        "elapsed_s": round(elapsed, 2),
        "audit_80d": cert,
    }


def main():
    tbl = refresh_live_plain()
    if torch.cuda.is_available():
        _ = torch.randn((512, 512), device="cuda") @ torch.randn((512, 512), device="cuda")
        torch.cuda.synchronize()

    out_path = REC_DIR / "wave3_records.json"
    all_res = json.loads(out_path.read_text(encoding="utf-8")) if out_path.exists() else []
    done = {(r["d"], r["n"]) for r in all_res}

    targets = [t for t in select_wave3_targets(tbl, max_targets=80) if t not in done]
    print(f"WAVE 3 (Math KB CARDAL + OpenMathOS + Cascade): {len(targets)} fresh targets queued!", flush=True)
    print("First 16 targets (Tier 1 Cascades):", targets[:16], flush=True)

    n_broken = sum(1 for r in all_res if r["broken"])
    for idx, (d, n) in enumerate(targets, 1):
        r = attack_wave3(d, n, tbl)
        all_res.append(r)
        if r["broken"]:
            n_broken += 1
        out_path.write_text(json.dumps(all_res, indent=2), encoding="utf-8")
        print(f"  [Wave 3 Progress {idx}/{len(targets)} | Wave 3 broken: {n_broken} (+44 Wave 1-2 live = {n_broken + 44} total)]", flush=True)

    print(f"\n================ WAVE 3 COMPLETE: {n_broken}/{len(targets)} NEW WORLD RECORDS BROKEN ================", flush=True)


if __name__ == "__main__":
    main()
