# DemoKinFast 2.2.0

The same treatment as 2.1.0, extended to every function. Outputs are still
bit-identical to DemoKin 1.0.3.

* `kin()` re-groups kin types (`"s"`, `"c"`, `"a"`, `"n"`) with the exact C++
  group sum of 2.1.0, keeping the order of `summarise(.by = ...)` (groups in
  order of first appearance) and its output class (a plain data frame unless
  the input is a tibble), instead of `dplyr::summarise()`.
* `MASS::ginv()` of a diagonal matrix is computed directly (`ginv_diag()`):
  the SVD of a diagonal matrix only rescales it (by its largest entry, for
  more than 25 entries, as LAPACK's `dbdsdc` does), so the result is
  reproduced bit for bit without the SVD. Verified on 267 inputs from 1 to
  1010 entries, with zeros, negative and tied values. This was half the run
  time of `kin_multi_stage_time_variant_2sex()` (`get_G()` every year, and
  the parents' mixing) and is used in `kin_multi_stage()` too.
* `kin_multi_stage_time_variant_2sex()`: the boundary year decomposed the same
  matrix twice with `eigen()` (population structure and parents'
  distribution); the second decomposition is reused. `keep_parity_block()`
  writes the kept entries directly instead of 101 S4 column assignments. The
  output tables build only the kin types asked for instead of building all
  and filtering.
* `kin_multi_stage()`: the parity restriction of mothers no longer builds a
  606 x 606 Kronecker product per age (the entries it keeps are copied, as
  the product returned them), and the projection matrices are bound as dense
  matrices instead of through Matrix's `rbind()`/`cbind()`.
* Time-invariant two-sex models and both cause-of-death models build their
  block matrices without `Matrix::bdiag()` (`bdiag_dense()`, the same
  entries), as the time-varying two-sex model did in 2.1.0.
* The time-varying two-sex model passes its yearly projection matrices to the
  kernel as triplet lists of their non-zero entries, so the 404 x 404 dense
  matrices are neither allocated nor scanned; the one-sex models pass Focal's
  ageing operator once, as triplets.

# DemoKinFast 2.1.0

* `kin2sex()` is 8-30x faster than 2.0.0 on large outputs and when kin types
  are grouped (`"s"`, `"c"`, `"a"`, `"n"`), still bit-identical to DemoKin
  1.0.3 (46 scenarios in `dev/equivalence/`, three of them new, with grouped
  kin). On a two-sex time-varying model with 76 years of rates (101 ages,
  eight grouped kin types, both focal sexes), a run with output for 46 years
  and the summary takes 8 s (2.0.0: 262 s; DemoKin 1.0.3: 469 s), and a run
  with three output years 0.8 s (8.4 s; 65 s).
* Grouped kin types are summed by an exact C++ group sum (long double, as
  R's `sum()`) instead of `dplyr::group_by() |> summarise_at(funs(sum))`,
  which evaluated `sum()` once per group and took most of the run time. This
  also removes the `funs()` deprecation warning.
* The group boundaries of the regrouping and of the summary table are found
  in C++, reading the keys through the sort order instead of copying them.
* Output tables index the id columns directly (indexing a data frame by
  repeated rows spent its time making unique row names), build only the
  requested output years, and skip the issue-51 subset when it drops nothing.
* The time-varying two-sex projection builds its yearly block matrices as
  dense matrices directly (no round trip through `Matrix::bdiag()`), passes
  Focal's ageing operator to the kernel once as triplets, and the kernel
  converts each matrix in one pass and visits only the non-empty rows of the
  fertility matrices (adding 0 to the other rows changed nothing).

# DemoKin 1.0.0

* Added a `NEWS.md` file to track changes to the package.
* Change stable/non-stable references to time varying/non-varying rates.
* Add multi-state process.

# DemoKin 1.0.1

* Submitted to CRAN
* Death counts are placed in the age where Focal experience the death.
* Aggregated kin types are allowed (`s` for older and younger sisters, for example).

# DemoKin 1.0.2

* Added the Two sex multi-state time-varying framework. 

# DemoKin 1.0.3

* Website
* Archive CRAN
* Typos
# DemoKinFast 2.0.0

* Package renamed `DemoKinFast` so it can be installed next to `DemoKin`.
* Same functions, arguments and outputs as 1.0.3, 10-100x faster. The kin
  projections (`U %*% x + F %*% z` recursions) now run in a C++ kernel that
  multiplies only the non-zero entries, accumulating in the same order as R's
  `%*%`, so results are bit-identical to 1.0.3 (verified on 44 scenarios
  covering every function and option; see `dev/equivalence/`).
* Output data frames are assembled directly instead of through
  `pivot_longer()`/`pivot_wider()`/`melt()`/`dcast()`/`merge()`; time-varying
  runs keep only the requested years in memory.
* `kin()`/`kin2sex()` summary indicators computed in C++ with R's long-double
  accumulation (identical values).
* Dependencies dropped: `matrixcalc`, `purrr`, `reshape2`. Added: `Rcpp`.
