"""4-Pronged RTX 5090 GPU + Exact HiGHS SLP World-Record Hunter for spherical-codes.org."""

from __future__ import annotations

import json
import time
from decimal import Decimal
from pathlib import Path

import numpy as np
import torch

from fetch_code import fetch_cohn_code
from gpu_solver import gpu_hinge_refine, gpu_homotopy_search, make_basin_hop_batch
from slp_polisher import slp_polish
from verifier import certify_exact_decimal, load_plain_table, save_code_cohn_format

ROOT = Path(__file__).resolve().parent
REC_DIR = ROOT / "records"
REC_DIR.mkdir(parents=True, exist_ok=True)


def select_targets(tbl: dict) -> list[tuple[int, int]]:
    """Select high-yield non-algebraic open targets across d=6..16, n=40..160."""
    by_d: dict[int, dict[int, float]] = {}
    for (d, n), (_, cf, _) in tbl.items():
        by_d.setdefault(d, {})[n] = cf

    cands = []
    for d in range(6, 17):
        if d not in by_d:
            continue
        for n, c in by_d[d].items():
            if not (40 <= n <= 160):
                continue
            if tbl[(d, n)][2]:  # skip algebraic minimal polynomial entries
                continue
            if (n - 1) in by_d[d] and (n + 1) in by_d[d]:
                bump = c - 0.5 * (by_d[d][n - 1] + by_d[d][n + 1])
                cands.append((bump, d, n))

    cands.sort(reverse=True)
    # Pick top 18 diverse targets (at most 3 per dimension so we explore broadly!)
    picked = []
    per_d: dict[int, int] = {}
    for _, d, n in cands:
        if per_d.get(d, 0) >= 3:
            continue
        picked.append((d, n))
        per_d[d] = per_d.get(d, 0) + 1
        if len(picked) >= 18:
            break
    return picked


def attack_target(d: int, n: int, tbl: dict) -> dict:
    rec_str, rec_f, _ = tbl[(d, n)]
    rec_dec = Decimal(rec_str)
    print(f"\n==================== TARGET d={d}, n={n} | LIVE RECORD = {rec_str} ====================", flush=True)

    best_X = None
    best_dec = rec_dec
    best_c12 = rec_str
    best_prong = "live_record"
    t0 = time.perf_counter()

    def _try_candidate(X_cand: np.ndarray, prong: str, force: bool = False):
        nonlocal best_X, best_dec, best_c12, best_prong
        # Quick float64 pre-screen so we only SLP-polish promising candidates
        U = X_cand / np.linalg.norm(X_cand, axis=1, keepdims=True)
        G = U @ U.T
        np.fill_diagonal(G, -2.0)
        if not force and float(G.max()) > float(best_dec) + 0.005:
            return False
        Xp, cp, nacc = slp_polish(X_cand, max_iters=40, max_pairs=900, lp_time_limit=3.0)
        if cp < float(best_dec) + 1e-8:
            c_ex, c12, pair = certify_exact_decimal(Xp, digits=60)
            c12_dec = Decimal(c12)
            if c12_dec < best_dec:
                delta = float(rec_dec - c12_dec)
                best_X = Xp.copy()
                best_dec = c12_dec
                best_c12 = c12
                best_prong = prong
                print(
                    f"  *** NEW WORLD RECORD ({prong})! ceil12={c12} < {rec_str} "
                    f"(IMPROVEMENT = -{delta:.12f}, slp_acc={nacc}) ***",
                    flush=True,
                )
                save_code_cohn_format(best_X, REC_DIR / f"d{d}_n{n}_NEW_RECORD.txt")
                return True
        return False

    # If we already have a saved new record on disk for (d, n), load & try to polish it even further!
    prev_rec_path = REC_DIR / f"d{d}_n{n}_NEW_RECORD.txt"
    if prev_rec_path.exists():
        vals = [float(x) for x in prev_rec_path.read_text(encoding="utf-8").replace(",", " ").split()]
        if len(vals) == n * d:
            X_prev = np.array(vals, dtype=np.float64).reshape(n, d)
            _, p12, _ = certify_exact_decimal(X_prev, digits=60)
            if Decimal(p12) < best_dec:
                best_X = X_prev.copy()
                best_dec = Decimal(p12)
                best_c12 = p12
                best_prong = "saved_prev_record"
                print(f"  Loaded existing NEW RECORD on disk: ceil12={p12}", flush=True)
            _try_candidate(X_prev, "A_repolish_saved", force=True)

    # Download live record
    X_rec = fetch_cohn_code(d, n)
    if X_rec is not None:
        _, ref12, _ = certify_exact_decimal(X_rec, digits=60)
        print(f"  Live record downloaded & verified: ceil12={ref12} (table={rec_str})", flush=True)

        # Prong A: Direct SLP polish of live record
        _try_candidate(X_rec, "A_slp_polish")

        # Prong B: GPU Basin-Hopping + Voronoi-Hole Teleportation
        bh_batch = make_basin_hop_batch(X_rec, batch_size=128, seed=d * 1000 + n)
        Xs_b, cs_b, dt_b = gpu_hinge_refine(bh_batch, hinge_steps=2000, top_k=4)
        print(f"  Prong B (Basin-Hop GPU {dt_b:.1f}s): best_gpu={cs_b[0]:.10f} (rec={rec_f:.10f})", flush=True)
        for k in range(len(Xs_b)):
            _try_candidate(Xs_b[k], f"B_basinhop_{k}")

    # Prong C: Point-Deletion Transfer from (d, n+1)
    if (d, n + 1) in tbl:
        X_np1 = fetch_cohn_code(d, n + 1)
        if X_np1 is not None:
            del_batch = np.stack([np.delete(X_np1, i, axis=0) for i in range(n + 1)], axis=0)
            Xs_c, cs_c, dt_c = gpu_hinge_refine(del_batch, hinge_steps=2000, top_k=4)
            print(f"  Prong C (Del from n+1 GPU {dt_c:.1f}s): best_gpu={cs_c[0]:.10f} (rec={rec_f:.10f})", flush=True)
            for k in range(len(Xs_c)):
                _try_candidate(Xs_c[k], f"C_del_np1_{k}")

    # Prong D: Pure From-Scratch Global GPU Homotopy
    Xs_d, cs_d, dt_d = gpu_homotopy_search(d, n, batch_size=256, seed=d * 100 + n, lse_steps=1800, hinge_steps=1800, top_k=4)
    print(f"  Prong D (From-Scratch GPU {dt_d:.1f}s): best_gpu={cs_d[0]:.10f} (rec={rec_f:.10f})", flush=True)
    for k in range(len(Xs_d)):
        _try_candidate(Xs_d[k], f"D_scratch_{k}")

    torch.cuda.empty_cache()
    elapsed = time.perf_counter() - t0
    broken = best_dec < rec_dec
    res = {
        "d": d,
        "n": n,
        "old_record_ceil12": rec_str,
        "new_record_ceil12": best_c12,
        "broken": broken,
        "improvement": float(rec_dec - best_dec),
        "winning_prong": best_prong,
        "elapsed_s": round(elapsed, 2),
    }
    print(
        f"  Summary ({d},{n}) in {elapsed:.1f}s: broken={broken}, old={rec_str}, best={best_c12}, "
        f"improv={res['improvement']:.12f} ({best_prong})",
        flush=True,
    )
    return res


def main():
    tbl = load_plain_table(ROOT / "data" / "plain.txt")
    # Warm up CUDA
    if torch.cuda.is_available():
        _ = torch.randn((512, 512), device="cuda") @ torch.randn((512, 512), device="cuda")
        torch.cuda.synchronize()

    targets = select_targets(tbl)
    print(f"Selected {len(targets)} high-yield open targets on spherical-codes.org:", targets, flush=True)

    out_json = REC_DIR / "world_records.json"
    all_res = json.loads(out_json.read_text(encoding="utf-8")) if out_json.exists() else []
    done_pairs = {(r["d"], r["n"]) for r in all_res}
    for d, n in targets:
        if (d, n) in done_pairs:
            print(f"Skipping already completed ({d},{n})", flush=True)
            continue
        r = attack_target(d, n, tbl)
        all_res.append(r)
        out_json.write_text(json.dumps(all_res, indent=2), encoding="utf-8")

    broken = [r for r in all_res if r["broken"]]
    print(f"\n================ HUNT COMPLETE: {len(broken)}/{len(all_res)} WORLD RECORDS BROKEN ================", flush=True)
    for r in broken:
        print(
            f"  d={r['d']:2d}, n={r['n']:3d}: {r['old_record_ceil12']} -> {r['new_record_ceil12']} "
            f"(improv={r['improvement']:.12f}, prong={r['winning_prong']})",
            flush=True,
        )


if __name__ == "__main__":
    main()
