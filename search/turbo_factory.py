"""Continuous Turbo World-Record Factory for Henry Cohn's spherical-codes.org (d=6..24)."""

from __future__ import annotations

import json
import math
import time
from decimal import Decimal
from pathlib import Path

import numpy as np
import torch

from audit_records import brute_force_certify
from fetch_code import fetch_cohn_code
from gpu_solver import gpu_hinge_refine, make_basin_hop_batch
from slp_polisher import slp_polish
from verifier import certify_exact_decimal, load_plain_table, save_code_cohn_format

ROOT = Path(__file__).resolve().parent
REC_DIR = ROOT / "records"
REC_DIR.mkdir(parents=True, exist_ok=True)


def make_insertion_batch(X_nm1: np.ndarray, batch_size: int = 64, seed: int = 0) -> np.ndarray:
    """Insert 1 point into the deepest Voronoi holes of (d, n-1) to form (batch_size, n, d)."""
    rng = np.random.default_rng(seed)
    X = np.asarray(X_nm1, dtype=np.float64)
    X = X / np.linalg.norm(X, axis=1, keepdims=True)
    nm1, d = X.shape
    probes = rng.standard_normal((8192, d))
    probes /= np.linalg.norm(probes, axis=1, keepdims=True)
    scores = (probes @ X.T).max(axis=1)
    holes = probes[np.argsort(scores)[:batch_size]]
    batch = np.empty((batch_size, nm1 + 1, d), dtype=np.float64)
    for b in range(batch_size):
        batch[b, :nm1] = X + (0.015 if b >= 16 else 0.0) * rng.standard_normal((nm1, d))
        batch[b, nm1] = holes[b]
        batch[b] /= np.linalg.norm(batch[b], axis=1, keepdims=True)
    return batch


def select_factory_targets(tbl: dict, max_targets: int = 160) -> list[tuple[int, int]]:
    """Select top 160 highest-bump non-algebraic open targets across d=6..24, n=30..180."""
    by_d: dict[int, dict[int, float]] = {}
    for (d, n), (_, cf, _) in tbl.items():
        by_d.setdefault(d, {})[n] = cf

    cands = []
    for d in range(6, 25):
        if d not in by_d:
            continue
        for n, c in by_d[d].items():
            if not (30 <= n <= 180):
                continue
            if tbl[(d, n)][2]:  # skip algebraic minimal polynomial entries
                continue
            if (n - 1) in by_d[d] and (n + 1) in by_d[d]:
                bump = c - 0.5 * (by_d[d][n - 1] + by_d[d][n + 1])
                cands.append((bump, d, n))

    cands.sort(reverse=True)
    picked = []
    per_d: dict[int, int] = {}
    for _, d, n in cands:
        if per_d.get(d, 0) >= 12:
            continue
        picked.append((d, n))
        per_d[d] = per_d.get(d, 0) + 1
        if len(picked) >= max_targets:
            break
    return picked


def attack_target_turbo(d: int, n: int, tbl: dict) -> dict:
    rec_str, rec_f, _ = tbl[(d, n)]
    rec_dec = Decimal(rec_str)
    print(f"\n=== TARGET ({d:2d}, {n:3d}) | LIVE={rec_str} ===", flush=True)

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
        if not force and float(G.max()) > float(best_dec) + 0.006:
            return False
        Xp, cp, nacc = slp_polish(X_cand, max_iters=30, max_pairs=700, lp_time_limit=1.5)
        if cp < float(best_dec) + 1e-8:
            _, c12, _ = certify_exact_decimal(Xp, digits=60)
            c12_dec = Decimal(c12)
            if c12_dec < best_dec:
                delta = float(rec_dec - c12_dec)
                best_X = Xp.copy()
                best_dec = c12_dec
                best_c12 = c12
                best_prong = prong
                csv_path = REC_DIR / f"d{d}_n{n}_NEW_RECORD.txt"
                save_code_cohn_format(best_X, csv_path)
                print(
                    f"  *** NEW WORLD RECORD ({prong})! {rec_str} -> {c12} "
                    f"(delta=-{delta:.12f}, slp_acc={nacc}) ***",
                    flush=True,
                )
                return True
        return False

    # Prong A + B: Live record polish + GPU Voronoi Basin-Hopping
    X_rec = fetch_cohn_code(d, n)
    if X_rec is not None:
        _try_candidate(X_rec, "A_slp_polish", force=True)

        bh_batch = make_basin_hop_batch(X_rec, batch_size=192, seed=d * 1000 + n)
        Xs_b, cs_b, dt_b = gpu_hinge_refine(bh_batch, hinge_steps=2000, top_k=4)
        print(f"  Prong B (BasinHop {dt_b:.1f}s): gpu={cs_b[0]:.10f}", flush=True)
        for k in range(len(Xs_b)):
            _try_candidate(Xs_b[k], f"B_basinhop_{k}")

    # Prong C: (d, n+1) -> (d, n) Single-Point Deletion Transfer
    if (d, n + 1) in tbl:
        X_np1 = fetch_cohn_code(d, n + 1)
        if X_np1 is not None:
            del_batch = np.stack([np.delete(X_np1, i, axis=0) for i in range(n + 1)], axis=0)
            Xs_c, cs_c, dt_c = gpu_hinge_refine(del_batch, hinge_steps=2000, top_k=4)
            print(f"  Prong C (Del n+1  {dt_c:.1f}s): gpu={cs_c[0]:.10f}", flush=True)
            for k in range(len(Xs_c)):
                _try_candidate(Xs_c[k], f"C_del_np1_{k}")

    # Prong E: (d, n-1) -> (d, n) Voronoi-Hole Insertion Transfer
    if (d, n - 1) in tbl:
        X_nm1 = fetch_cohn_code(d, n - 1)
        if X_nm1 is not None:
            ins_batch = make_insertion_batch(X_nm1, batch_size=64, seed=d * 2000 + n)
            Xs_e, cs_e, dt_e = gpu_hinge_refine(ins_batch, hinge_steps=2000, top_k=3)
            print(f"  Prong E (Ins n-1  {dt_e:.1f}s): gpu={cs_e[0]:.10f}", flush=True)
            for k in range(len(Xs_e)):
                _try_candidate(Xs_e[k], f"E_ins_nm1_{k}")

    # Prong F: Second-round deep basin-hop on OUR new record if we broke it!
    if best_X is not None and float(rec_dec - best_dec) >= 1e-6:
        bh2 = make_basin_hop_batch(best_X, batch_size=128, seed=d * 3000 + n)
        Xs_f, cs_f, dt_f = gpu_hinge_refine(bh2, hinge_steps=2200, top_k=3)
        print(f"  Prong F (DeepHop  {dt_f:.1f}s): gpu={cs_f[0]:.10f}", flush=True)
        for k in range(len(Xs_f)):
            _try_candidate(Xs_f[k], f"F_deephop_{k}")

    torch.cuda.empty_cache()
    elapsed = time.perf_counter() - t0
    broken = best_dec < rec_dec

    cert = None
    if broken:
        # Immediate brute-force all-pairs 80-digit audit on saved CSV!
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
    tbl = load_plain_table(ROOT / "data" / "plain.txt")
    if torch.cuda.is_available():
        _ = torch.randn((512, 512), device="cuda") @ torch.randn((512, 512), device="cuda")
        torch.cuda.synchronize()

    # Load previous results so we don't repeat already-run targets
    prev_path = REC_DIR / "world_records.json"
    fac_path = REC_DIR / "factory_records.json"
    done_pairs = set()
    if prev_path.exists():
        for r in json.loads(prev_path.read_text(encoding="utf-8")):
            done_pairs.add((r["d"], r["n"]))
    all_res = json.loads(fac_path.read_text(encoding="utf-8")) if fac_path.exists() else []
    for r in all_res:
        done_pairs.add((r["d"], r["n"]))

    targets = [t for t in select_factory_targets(tbl, max_targets=160) if t not in done_pairs]
    print(f"TURBO WORLD-RECORD FACTORY: {len(targets)} fresh open targets queued (d=6..24, n=30..180)", flush=True)

    n_broken_session = sum(1 for r in all_res if r["broken"])
    for idx, (d, n) in enumerate(targets, 1):
        r = attack_target_turbo(d, n, tbl)
        all_res.append(r)
        if r["broken"]:
            n_broken_session += 1
        fac_path.write_text(json.dumps(all_res, indent=2), encoding="utf-8")
        print(f"  [Progress {idx}/{len(targets)} | Factory records broken so far: {n_broken_session} (+4 from Wave 1 = {n_broken_session + 4} total)]", flush=True)

    print(f"\n================ FACTORY WAVE COMPLETE: {n_broken_session} NEW + 4 WAVE-1 = {n_broken_session + 4} TOTAL WORLD RECORDS ================", flush=True)


if __name__ == "__main__":
    main()
