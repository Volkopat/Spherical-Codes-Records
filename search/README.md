# Search code

These are the scripts as they were run, kept as a record of what produced the first upload. They are not a
polished package.

| File | What it does |
|---|---|
| `gpu_solver.py` | Batched GPU optimisation: soft-maximum annealing, squared-hinge refinement, and the perturbation batches |
| `gpu_solver_v3.py` | Augmented-Lagrangian refinement, lifting to extra dimensions, and the Douglas-Rachford variant that was dropped |
| `slp_polisher.py` | Trust-region sequential linear programming polish (HiGHS) |
| `hunt_records.py` | First run: 18 targets |
| `turbo_factory.py` | Second run: 142 targets |
| `wave3_factory.py` | Third run: 80 targets planned, 50 run |
| `fetch_code.py` | Downloads and caches the listed code for a table entry |
| `verifier.py`, `audit_records.py` | Decimal-arithmetic checks used inside the runs |

Things to know before reusing them:

- The run scripts download the table and the listed codes from spherical-codes.org into `data/` and write their
  results into `records/`, both next to the scripts.
- They compare against the copy of the table saved when the run starts. The table moves during a run, so a
  code they report as an improvement may already have been overtaken; check against a fresh copy with
  `verify/verify_exact.py`.
- `fetch_code.py` treats any download error as "no code available", so a target can be skipped silently.
- The polish in `slp_polisher.py` is called with tight caps (25 to 30 iterations, 700 pairs, 1.5 seconds per
  linear programme) and stops short of a local optimum on most codes. Use `polish/polish_al.py` afterwards.
- Messages printed by these scripts call an improvement a "world record". That means only that the value was
  below the cached table entry when the script ran.
