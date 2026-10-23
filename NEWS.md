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
# DemoKin 2.0.0

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
