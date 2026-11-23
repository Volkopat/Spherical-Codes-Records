# Spherical Codes Records

Improved spherical codes for the [Table of Spherical Codes](https://spherical-codes.org/) maintained by Henry Cohn,
together with the code that found them, an exact verification script and a description of the method.

A spherical code is a set of `n` unit vectors in `R^d`; the aim is to make the smallest angle between two points
as large as possible, that is, to make the largest cosine as small as possible.

## Status

- **23 November 2026, first upload.** 71 codes in dimensions 6 to 31, with 36 to 215 points, were accepted into
  the table. Each lowers the listed cosine of the minimal angle, 20 of them by more than 1e-3.
- **Polished versions.** The same 71 codes after a final polish are being prepared and will be added in `codes/polished`.

The table changes quickly, so any of these entries may since have been improved, by us or by others. The live
table is the authority; this repository records what was uploaded and how it was obtained.

## The Codes

Cosines are rigorous upper bounds rounded up to twelve decimal places, as in the table. "Table before" is the
listed value just before the first upload. "Total drop" and "Angle gain" compare the best of our two versions
with it.

| d | n | Table before | First upload | Polished | Total drop | Angle gain |
|---:|---:|---|---|---|---:|---:|
| 6 | 63 | 0.498122491823 | 0.498115244407 |  | 7.25e-06 | +0.0005° |
| 6 | 74 | 0.521539568537 | 0.521520006171 |  | 1.96e-05 | +0.0013° |
| 6 | 115 | 0.596004493388 | 0.595701844045 |  | 3.03e-04 | +0.0216° |
| 6 | 121 | 0.603684460657 | 0.603684460653 |  | 4.00e-12 | +0.0000° |
| 6 | 127 | 0.616047073470 | 0.615810912727 |  | 2.36e-04 | +0.0172° |
| 7 | 42 | 0.329083455200 | 0.329046843051 |  | 3.66e-05 | +0.0022° |
| 7 | 183 | 0.586168985388 | 0.585983225851 |  | 1.86e-04 | +0.0131° |
| 8 | 79 | 0.394977978307 | 0.394977634281 |  | 3.44e-07 | +0.0000° |
| 8 | 119 | 0.458816566687 | 0.458464628615 |  | 3.52e-04 | +0.0227° |
| 9 | 102 | 0.378387942959 | 0.377672507745 |  | 7.15e-04 | +0.0443° |
| 10 | 85 | 0.314387838769 | 0.314381011297 |  | 6.83e-06 | +0.0004° |
| 10 | 88 | 0.322419246773 | 0.322012551456 |  | 4.07e-04 | +0.0246° |
| 10 | 89 | 0.322730324266 | 0.322508831910 |  | 2.21e-04 | +0.0134° |
| 10 | 115 | 0.352502539478 | 0.352490310886 |  | 1.22e-05 | +0.0007° |
| 10 | 133 | 0.384724603106 | 0.384722545132 |  | 2.06e-06 | +0.0001° |
| 10 | 134 | 0.387500256997 | 0.387391983482 |  | 1.08e-04 | +0.0067° |
| 11 | 126 | 0.328447661914 | 0.328429429201 |  | 1.82e-05 | +0.0011° |
| 11 | 127 | 0.332098604310 | 0.331436410570 |  | 6.62e-04 | +0.0402° |
| 11 | 128 | 0.332166503541 | 0.331535252262 |  | 6.31e-04 | +0.0383° |
| 11 | 167 | 0.378095864980 | 0.377933007828 |  | 1.63e-04 | +0.0101° |
| 11 | 215 | 0.409713439170 | 0.409073301313 |  | 6.40e-04 | +0.0402° |
| 12 | 79 | 0.241071142239 | 0.238883076505 |  | 2.19e-03 | +0.1291° |
| 12 | 80 | 0.243279277482 | 0.241900795717 |  | 1.38e-03 | +0.0814° |
| 12 | 92 | 0.258007696391 | 0.257948221348 |  | 5.95e-05 | +0.0035° |
| 12 | 114 | 0.290350338240 | 0.290219644029 |  | 1.31e-04 | +0.0078° |
| 12 | 141 | 0.310962949342 | 0.309348528084 |  | 1.61e-03 | +0.0973° |
| 12 | 171 | 0.322345337630 | 0.321811773354 |  | 5.34e-04 | +0.0323° |
| 12 | 172 | 0.332989693827 | 0.325359224605 |  | 7.63e-03 | +0.4630° |
| 12 | 173 | 0.333022278252 | 0.332944812555 |  | 7.75e-05 | +0.0047° |
| 12 | 181 | 0.344130099215 | 0.343604852221 |  | 5.25e-04 | +0.0320° |
| 12 | 186 | 0.351891117772 | 0.351610649413 |  | 2.80e-04 | +0.0172° |
| 13 | 70 | 0.201811558513 | 0.201193114439 |  | 6.18e-04 | +0.0362° |
| 13 | 71 | 0.204333058956 | 0.203126915411 |  | 1.21e-03 | +0.0706° |
| 13 | 123 | 0.272084487026 | 0.272081111911 |  | 3.38e-06 | +0.0002° |
| 13 | 124 | 0.275373068352 | 0.274881882109 |  | 4.91e-04 | +0.0293° |
| 13 | 139 | 0.291298832564 | 0.291033633055 |  | 2.65e-04 | +0.0159° |
| 13 | 140 | 0.294790367807 | 0.292339708355 |  | 2.45e-03 | +0.1469° |
| 13 | 141 | 0.295463243126 | 0.295000930305 |  | 4.62e-04 | +0.0277° |
| 13 | 209 | 0.339851756282 | 0.338134930945 |  | 1.72e-03 | +0.1046° |
| 14 | 93 | 0.211614877846 | 0.211609073555 |  | 5.80e-06 | +0.0003° |
| 14 | 100 | 0.223656652929 | 0.223300703402 |  | 3.56e-04 | +0.0209° |
| 17 | 36 | 0.052701899358 | 0.052695903487 |  | 6.00e-06 | +0.0003° |
| 17 | 90 | 0.163935338020 | 0.163923217779 |  | 1.21e-05 | +0.0007° |
| 17 | 166 | 0.209186347459 | 0.209023681833 |  | 1.63e-04 | +0.0095° |
| 17 | 178 | 0.215890665682 | 0.215868215445 |  | 2.25e-05 | +0.0013° |
| 17 | 179 | 0.223180793284 | 0.221773458354 |  | 1.41e-03 | +0.0827° |
| 17 | 181 | 0.226940470105 | 0.226923780714 |  | 1.67e-05 | +0.0010° |
| 18 | 43 | 0.066413840190 | 0.065744588196 |  | 6.69e-04 | +0.0384° |
| 18 | 175 | 0.208896958003 | 0.208891532279 |  | 5.43e-06 | +0.0003° |
| 18 | 176 | 0.212563244199 | 0.209841970423 |  | 2.72e-03 | +0.1595° |
| 18 | 177 | 0.213961001219 | 0.212596689376 |  | 1.36e-03 | +0.0800° |
| 21 | 120 | 0.139344261657 | 0.133355223356 |  | 5.99e-03 | +0.3464° |
| 21 | 121 | 0.139344261810 | 0.133538570034 |  | 5.81e-03 | +0.3358° |
| 22 | 53 | 0.057403269990 | 0.057383937843 |  | 1.93e-05 | +0.0011° |
| 22 | 63 | 0.076899986863 | 0.075392736227 |  | 1.51e-03 | +0.0866° |
| 23 | 120 | 0.119777994765 | 0.119231406771 |  | 5.47e-04 | +0.0315° |
| 23 | 122 | 0.124520867516 | 0.123067950381 |  | 1.45e-03 | +0.0839° |
| 23 | 176 | 0.152489921745 | 0.148600843119 |  | 3.89e-03 | +0.2254° |
| 23 | 177 | 0.152871506347 | 0.150697001281 |  | 2.17e-03 | +0.1261° |
| 24 | 61 | 0.063646394908 | 0.063627102576 |  | 1.93e-05 | +0.0011° |
| 24 | 70 | 0.076518930179 | 0.074691184365 |  | 1.83e-03 | +0.1050° |
| 24 | 71 | 0.076910440718 | 0.076664178538 |  | 2.46e-04 | +0.0142° |
| 24 | 127 | 0.127010481768 | 0.126863836126 |  | 1.47e-04 | +0.0085° |
| 24 | 128 | (not listed) | 0.129361781681 |  | | |
| 24 | 129 | 0.131456819733 | 0.131008694242 |  | 4.48e-04 | +0.0259° |
| 25 | 65 | 0.064022037337 | 0.062376397304 |  | 1.65e-03 | +0.0945° |
| 26 | 67 | 0.062956992866 | 0.060546855401 |  | 2.41e-03 | +0.1384° |
| 28 | 184 | 0.141515334952 | 0.141342539545 |  | 1.73e-04 | +0.0100° |
| 29 | 64 | 0.043272075045 | 0.043152598075 |  | 1.19e-04 | +0.0069° |
| 30 | 63 | 0.037664488335 | 0.037653889484 |  | 1.06e-05 | +0.0006° |
| 31 | 65 | 0.037040514725 | 0.037026775144 |  | 1.37e-05 | +0.0008° |

The entry for 128 points in 24 dimensions was not listed separately before the upload; its value then followed
from the 129-point entry, 0.131456819733.

Files are named `d{d}_n{n}.txt`, with one point per line and comma-separated coordinates. The same numbers are in
`results/` as CSV, next to the table's response to the first upload.

## Checking A Code

```
pip install numpy scipy mpmath
mkdir -p live/plain
curl -o live/plain/plain.txt https://spherical-codes.org/plain
python verify/verify_exact.py codes/polished live
```

The script reads the coordinates as exact fractions, proves the twelve-place bound for every pair, and compares it
with the downloaded table, including values that follow from other entries. An entry that is still ours shows as
equal to the table, not below it.

## Method

Each code was found by local search from codes already in the table: the listed code is perturbed, or a point is
deleted from the listed code with one more point, or a point is inserted into the listed code with one fewer.
Batches of such starts are optimised on a GPU with penalty and augmented-Lagrangian methods, in some cases after
lifting the code into two or three extra dimensions and squeezing it back. A final polish by the method of
multipliers takes each code to a local optimum. Every value is then checked in exact rational arithmetic.

The full description is in [METHODS.md](METHODS.md).

## What Is Not Claimed

- No code here is claimed to be optimal.
- None of the techniques is new; they are standard tools applied at scale to one table.
- The entries are not independent: many are neighbours in the number of points within one dimension and were
  found from one another.
- The search was run in three stages with different targets, so nothing is claimed about which solver component
  is responsible for the improvements.

## How This Was Made

The search code was written and run by an AI coding agent under the author's direction. The audit against the
live table, the final polisher and this write-up were produced with a second AI agent. No value in this
repository rests on either agent's report: each was established by the exact check in `verify/` and by the
table's own acceptance test.

## Credit

The starting points are entries of the table contributed by others. The entries we looked up were credited to
Henry Cohn, Ilmari Vahteristo, Vadim Borisov, Antoine Jay and A. Canli; each entry's page on spherical-codes.org
gives its reference.

## Citation

```
Digvijay Patil, Spherical Codes Records: Improved Entries for the Table of Spherical Codes, 2026.
https://github.com/Volkopat/Spherical-Codes-Records
```

## Licence

MIT; see [LICENSE](LICENSE).
