# DemoKin (fast)

A faster drop-in version of [DemoKin](https://github.com/IvanWilli/DemoKin) (version 1.0.3): the same
functions, arguments and outputs, 10-100x faster. The kinship projections run in a C++ kernel (Rcpp) that
multiplies only the non-zero entries of the projection matrices, accumulating in the same order as R's
`%*%`, so results are **bit-identical** to DemoKin 1.0.3.

## Install

```r
# install.packages("remotes")
remotes::install_github("alburezg/DemoKin-fast")
```

The package is compiled from source, so Windows needs [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
and macOS needs the Xcode command line tools. The package is still called `DemoKin`, so installing it
replaces DemoKin 1.0.3 from CRAN; `install.packages("DemoKin")` brings the original back.

## Use

Exactly as before:

```r
library(DemoKin)
swe <- kin(p = swe_px[, "2015"], f = swe_asfr[, "2015"])
head(swe$kin_summary)
```

## Checks

- `tests/testthat/test-equivalence_with_v1.R`: outputs compared with `identical()` to fixtures produced
  by DemoKin 1.0.3, plus the package's original Caswell (2020) tests.
- `dev/equivalence/`: harness that runs 44 calls covering every exported function and option with both
  versions and compares them (44/44 identical).
- `dev/compare_DemoKin_original_vs_fast.R`: self-contained version of that comparison; run it from a
  folder that contains both `DemoKin/` (original source) and `DemoKin_fast/`.

See `NEWS.md` for the list of changes.
