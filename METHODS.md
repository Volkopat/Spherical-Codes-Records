# Methods

This document describes how the codes in this repository were found, polished and checked, and what is
and is not claimed about them.

## Problem

A spherical code is a set of `n` unit vectors in `R^d`. The task is to make the smallest angle between two
distinct points as large as possible, that is, to minimise the largest cosine

    c(X) = max over i < j of <x_i, x_j>.

The [Table of Spherical Codes](https://spherical-codes.org/) lists, for each dimension and number of points, the
best code known, with the cosine given as a rigorous upper bound rounded up to twelve decimal places.

## 1. Search

The search code is in `search/`. It runs batches of candidate codes on one GPU (PyTorch) and is followed by a
polishing step on the CPU.

**Targets.** Entries of the table whose value sits well above the average of its two neighbours in the number of
points, `c(d, n) - (c(d, n-1) + c(d, n+1)) / 2`, were tried first, since such a bump suggests a weak entry.
The first two runs covered dimensions 6 to 24 and 30 to 180 points; the third added the neighbours `(d, n-1)` and
`(d, n+1)` of earlier improvements and went up to dimension 32 and 230 points. Entries listed with a minimal
polynomial were skipped.

**Starting points.** Almost every improvement starts from a code that is already in the table:

- the listed code for `(d, n)`, perturbed;
- the listed code for `(d, n+1)` with one point deleted;
- the listed code for `(d, n-1)` with one point inserted into a deep hole;
- in the third run, the same three constructions applied to our own earlier improvements;
- in the first run only, random starts.

The codes here are therefore local improvements of the work of the table's earlier contributors, not
independent constructions.

**Perturbation.** Each start is copied into a batch of 64 to 192 candidates. Half receive Gaussian kicks at six
scales (0.02 to 0.45). In the other half, one to three of the most crowded points are moved to deep holes of the
code, found by random probing, followed by a small kick.

**Optimisation on the GPU.**

- *Soft maximum* (random starts only): Riemannian Adam on `(1/b) log sum exp(b <x_i, x_j>)`, with `b` raised from 6
  to 220, in single precision.
- *Squared hinge*: Riemannian Adam on `sum max(0, <x_i, x_j> - s)^2` in double precision, where the level `s` is
  kept slightly below the best largest cosine seen so far and the gap is tightened from about 1e-2 to 1e-4.
- *Augmented Lagrangian* (third run): the same loop with a multiplier for every pair
  (the Powell-Hestenes-Rockafellar form), which approaches the optimum more closely than the plain penalty.
- *Lifting* (third run): the code is embedded in `R^(d+2)` or `R^(d+3)`, optimised there, and the extra coordinates
  are then squeezed back to zero by a growing quadratic penalty. Raising the rank to escape a local optimum is a
  standard device in low-rank semidefinite programming; the immediate prompt for trying it here was a
  rank-adaptive augmented-Lagrangian solver for semidefinite programs (arXiv:2607.17933).
- A Douglas-Rachford iteration on the Gram matrix was also tried and dropped: in a test on three targets it was
  about fifty times slower than the augmented-Lagrangian loop and gave worse codes.

**First polish.** The best few candidates of each batch were refined by sequential linear programming (HiGHS) with
a trust region. This step was capped at 25 to 30 iterations, 700 pairs and 1.5 seconds per linear programme, which
turned out to be too little: most of the codes it produced were still short of a local optimum (see below).

## 2. Final polish

`polish/polish_al.py` takes a code to a local optimum of the problem

    minimise t   subject to   <x_i, x_j> <= t for all i < j,   |x_i| = 1.

It is the classical method of multipliers over all pairs, in double precision on the CPU:

    L(X, t; lam) = t + (1 / (2 rho)) * sum over i < j of ( max(0, lam_ij + rho (<x_i, x_j> - t))^2 - lam_ij^2 )

is minimised over `X` and `t` by L-BFGS, the multipliers are updated by `lam <- max(0, lam + rho (G - t))`, and
`rho` is raised when the infeasibility stops shrinking. A last Gauss-Newton step moves the pairs that carry a
multiplier onto one common value. Codes of this size take seconds to about two minutes each.

Different starting values of `rho` relax a code into different nearby local optima, whose values differ in the
sixth to eighth decimal place, so each code is polished from several starting values and the best is kept.
The result agrees with the local optimum to about 1e-13; the last unit of the twelfth decimal place is not
guaranteed.

The first upload was made before this step existed. Polishing lowered 65 of those 71 codes at twelve places; the
comparison is in the README.

## 3. Verification

`verify/verify_exact.py` checks each code file in exact rational arithmetic:

1. the decimal strings in the file are read as exact fractions;
2. `c12`, the largest cosine rounded up to twelve places, is computed at 80 digits and then proved: every pair
   satisfies `cos <= c12` and some pair satisfies `cos > c12 - 1e-12`, tested without square roots as
   `<x, y>^2 <= c12^2 |x|^2 |y|^2`;
3. `c12` is compared with a fresh copy of the table, including the values that follow from other entries by
   deleting points, embedding in a higher dimension, or taking a union of codes in orthogonal subspaces.

A code is counted as an improvement only if `c12` is strictly below all of these. The table's own upload check
applies the same rule; every file in `codes/first_upload` was accepted by it on 23 November 2026.

## What is not claimed

- **Optimality.** These are the best codes known for their parameters at the time of upload, nothing more.
- **New techniques.** Soft-maximum annealing, penalty and augmented-Lagrangian methods, rank lifting,
  perturb-and-reoptimise search and sequential linear programming are all standard.
- **Which component matters.** The three runs used different targets, so the higher success rate of the third run
  cannot be attributed to its solver; the solver variants were never compared on the same targets.
- **Independence of the entries.** Many improvements are neighbours in the number of points within one
  dimension and were found from one another.

## Known gaps

- The table is changing quickly: other contributors lowered more than nine thousand of its entries during the
  two days of this work. Of the 86 codes the search saved as improvements, 14 had been overtaken or tied before
  upload and one was made redundant by another of our own codes; 71 were uploaded. Any entry here may since
  have been improved, by us or by others.
- The polish finds the local optimum nearest to a code, not the best one nearby. Within hours of upload, 16 of
  the 71 entries were improved by another contributor by between 2e-8 and 1.4e-4, so better optima exist close
  to those codes that this search did not reach.
- Thirty targets of the third run were skipped because their starting codes could not be downloaded at the time,
  and were never tried.
- The search scripts compare against a copy of the table saved when a run starts, which goes stale during the
  run. The exact check in section 3 against a fresh copy is what decides whether a code counts.
