# Spherical Codes Records

Improved spherical codes for the [Table of Spherical Codes](https://spherical-codes.org/) maintained by Henry Cohn,
together with the code that found them, an exact verification script and a description of the method.

A spherical code is a set of `n` unit vectors in `R^d`; the aim is to make the smallest angle between two points
as large as possible, that is, to make the largest cosine as small as possible.

## Status

- **23 November 2026, first upload.** 71 codes in dimensions 6 to 31, with 36 to 215 points, were accepted into
  the table. Each lowers the listed cosine of the minimal angle, 20 of them by more than 1e-3.
- **Polished versions.** A final polish lowered 65 of the 71 further, by between 1e-11 and 5e-5; the other six
  were already at a local optimum to twelve places. The polished codes are in `codes/polished` and were
  prepared for a second upload on the same day.

The table changes quickly, so any of these entries may since have been improved, by us or by others. The live
table is the authority; this repository records what was uploaded and how it was obtained.

## The Codes

Cosines are rigorous upper bounds rounded up to twelve decimal places, as in the table. "Table before" is the
listed value just before the first upload. "Total drop" and "Angle gain" compare the best of our two versions
with it.

| d | n | Table before | First upload | Polished | Total drop | Angle gain |
|---:|---:|---|---|---|---:|---:|
| 6 | 63 | 0.498122491823 | 0.498115244407 | (unchanged) | 7.25e-06 | +0.0005° |
| 6 | 74 | 0.521539568537 | 0.521520006171 | (unchanged) | 1.96e-05 | +0.0013° |
| 6 | 115 | 0.596004493388 | 0.595701844045 | 0.595701844033 | 3.03e-04 | +0.0216° |
| 6 | 121 | 0.603684460657 | 0.603684460653 | (unchanged) | 4.00e-12 | +0.0000° |
| 6 | 127 | 0.616047073470 | 0.615810912727 | 0.615810904536 | 2.36e-04 | +0.0172° |
| 7 | 42 | 0.329083455200 | 0.329046843051 | (unchanged) | 3.66e-05 | +0.0022° |
| 7 | 183 | 0.586168985388 | 0.585983225851 | 0.585975202564 | 1.94e-04 | +0.0137° |
| 8 | 79 | 0.394977978307 | 0.394977634281 | 0.394977619599 | 3.59e-07 | +0.0000° |
| 8 | 119 | 0.458816566687 | 0.458464628615 | 0.458463391892 | 3.53e-04 | +0.0228° |
| 9 | 102 | 0.378387942959 | 0.377672507745 | 0.377651721871 | 7.36e-04 | +0.0456° |
| 10 | 85 | 0.314387838769 | 0.314381011297 | 0.314380222714 | 7.62e-06 | +0.0005° |
| 10 | 88 | 0.322419246773 | 0.322012551456 | 0.322001446626 | 4.18e-04 | +0.0253° |
| 10 | 89 | 0.322730324266 | 0.322508831910 | 0.322507152264 | 2.23e-04 | +0.0135° |
| 10 | 115 | 0.352502539478 | 0.352490310886 | 0.352449817908 | 5.27e-05 | +0.0032° |
| 10 | 133 | 0.384724603106 | 0.384722545132 | 0.384716529631 | 8.07e-06 | +0.0005° |
| 10 | 134 | 0.387500256997 | 0.387391983482 | 0.387385140712 | 1.15e-04 | +0.0072° |
| 11 | 126 | 0.328447661914 | 0.328429429201 | 0.328424878277 | 2.28e-05 | +0.0014° |
| 11 | 127 | 0.332098604310 | 0.331436410570 | 0.331408056846 | 6.91e-04 | +0.0419° |
| 11 | 128 | 0.332166503541 | 0.331535252262 | 0.331529911757 | 6.37e-04 | +0.0387° |
| 11 | 167 | 0.378095864980 | 0.377933007828 | 0.377911406358 | 1.84e-04 | +0.0114° |
| 11 | 215 | 0.409713439170 | 0.409073301313 | 0.409062863475 | 6.51e-04 | +0.0409° |
| 12 | 79 | 0.241071142239 | 0.238883076505 | 0.238881236220 | 2.19e-03 | +0.1292° |
| 12 | 80 | 0.243279277482 | 0.241900795717 | 0.241898498338 | 1.38e-03 | +0.0815° |
| 12 | 92 | 0.258007696391 | 0.257948221348 | 0.257946208564 | 6.15e-05 | +0.0036° |
| 12 | 114 | 0.290350338240 | 0.290219644029 | 0.290213688468 | 1.37e-04 | +0.0082° |
| 12 | 141 | 0.310962949342 | 0.309348528084 | 0.309343419095 | 1.62e-03 | +0.0976° |
| 12 | 171 | 0.322345337630 | 0.321811773354 | 0.321807118337 | 5.38e-04 | +0.0326° |
| 12 | 172 | 0.332989693827 | 0.325359224605 | 0.325328389266 | 7.66e-03 | +0.4649° |
| 12 | 173 | 0.333022278252 | 0.332944812555 | 0.332939491114 | 8.28e-05 | +0.0050° |
| 12 | 181 | 0.344130099215 | 0.343604852221 | 0.343596114412 | 5.34e-04 | +0.0326° |
| 12 | 186 | 0.351891117772 | 0.351610649413 | 0.351605408297 | 2.86e-04 | +0.0175° |
| 13 | 70 | 0.201811558513 | 0.201193114439 | 0.201189656366 | 6.22e-04 | +0.0364° |
| 13 | 71 | 0.204333058956 | 0.203126915411 | 0.203125858215 | 1.21e-03 | +0.0706° |
| 13 | 123 | 0.272084487026 | 0.272081111911 | 0.272076663259 | 7.82e-06 | +0.0005° |
| 13 | 124 | 0.275373068352 | 0.274881882109 | 0.274874385982 | 4.99e-04 | +0.0297° |
| 13 | 139 | 0.291298832564 | 0.291033633055 | 0.291027386415 | 2.71e-04 | +0.0163° |
| 13 | 140 | 0.294790367807 | 0.292339708355 | 0.292332011713 | 2.46e-03 | +0.1473° |
| 13 | 141 | 0.295463243126 | 0.295000930305 | 0.294990407336 | 4.73e-04 | +0.0284° |
| 13 | 209 | 0.339851756282 | 0.338134930945 | 0.338128488770 | 1.72e-03 | +0.1050° |
| 14 | 93 | 0.211614877846 | 0.211609073555 | 0.211604236112 | 1.06e-05 | +0.0006° |
| 14 | 100 | 0.223656652929 | 0.223300703402 | 0.223293823763 | 3.63e-04 | +0.0213° |
| 17 | 36 | 0.052701899358 | 0.052695903487 | (unchanged) | 6.00e-06 | +0.0003° |
| 17 | 90 | 0.163935338020 | 0.163923217779 | 0.163906680653 | 2.87e-05 | +0.0017° |
| 17 | 166 | 0.209186347459 | 0.209023681833 | 0.209003551235 | 1.83e-04 | +0.0107° |
| 17 | 178 | 0.215890665682 | 0.215868215445 | 0.215859219271 | 3.14e-05 | +0.0018° |
| 17 | 179 | 0.223180793284 | 0.221773458354 | 0.221734671032 | 1.45e-03 | +0.0850° |
| 17 | 181 | 0.226940470105 | 0.226923780714 | 0.226919117570 | 2.14e-05 | +0.0013° |
| 18 | 43 | 0.066413840190 | 0.065744588196 | (unchanged) | 6.69e-04 | +0.0384° |
| 18 | 175 | 0.208896958003 | 0.208891532279 | 0.208883551458 | 1.34e-05 | +0.0008° |
| 18 | 176 | 0.212563244199 | 0.209841970423 | 0.209823018944 | 2.74e-03 | +0.1606° |
| 18 | 177 | 0.213961001219 | 0.212596689376 | 0.212586793964 | 1.37e-03 | +0.0806° |
| 21 | 120 | 0.139344261657 | 0.133355223356 | 0.133330364538 | 6.01e-03 | +0.3478° |
| 21 | 121 | 0.139344261810 | 0.133538570034 | 0.133533577074 | 5.81e-03 | +0.3361° |
| 22 | 53 | 0.057403269990 | 0.057383937843 | 0.057377446010 | 2.58e-05 | +0.0015° |
| 22 | 63 | 0.076899986863 | 0.075392736227 | 0.075389923530 | 1.51e-03 | +0.0868° |
| 23 | 120 | 0.119777994765 | 0.119231406771 | 0.119223271761 | 5.55e-04 | +0.0320° |
| 23 | 122 | 0.124520867516 | 0.123067950381 | 0.123058974673 | 1.46e-03 | +0.0844° |
| 23 | 176 | 0.152489921745 | 0.148600843119 | 0.148576442429 | 3.91e-03 | +0.2268° |
| 23 | 177 | 0.152871506347 | 0.150697001281 | 0.150689292236 | 2.18e-03 | +0.1265° |
| 24 | 61 | 0.063646394908 | 0.063627102576 | 0.063611940810 | 3.45e-05 | +0.0020° |
| 24 | 70 | 0.076518930179 | 0.074691184365 | 0.074681599642 | 1.84e-03 | +0.1056° |
| 24 | 71 | 0.076910440718 | 0.076664178538 | 0.076653551374 | 2.57e-04 | +0.0148° |
| 24 | 127 | 0.127010481768 | 0.126863836126 | 0.126849815947 | 1.61e-04 | +0.0093° |
| 24 | 128 | 0.131456819733* | 0.129361781681 | 0.129311200388 | 2.15e-03 | +0.1240° |
| 24 | 129 | 0.131456819733 | 0.131008694242 | 0.130994338315 | 4.62e-04 | +0.0267° |
| 25 | 65 | 0.064022037337 | 0.062376397304 | 0.062372842048 | 1.65e-03 | +0.0947° |
| 26 | 67 | 0.062956992866 | 0.060546855401 | 0.060539953455 | 2.42e-03 | +0.1388° |
| 28 | 184 | 0.141515334952 | 0.141342539545 | 0.141300683499 | 2.15e-04 | +0.0124° |
| 29 | 64 | 0.043272075045 | 0.043152598075 | 0.043139261587 | 1.33e-04 | +0.0076° |
| 30 | 63 | 0.037664488335 | 0.037653889484 | 0.037645388296 | 1.91e-05 | +0.0011° |
| 31 | 65 | 0.037040514725 | 0.037026775144 | 0.037015237882 | 2.53e-05 | +0.0014° |

\* The entry for 128 points in 24 dimensions was not listed separately before the upload; the value shown
followed from the 129-point entry.

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
