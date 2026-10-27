# Does ginv_diag(s) reproduce MASS::ginv(diag(s)) bit for bit?
#   Rscript dev/regroup/test_ginv_diag.R
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))   # run from the package root
set.seed(1)
cases <- list()
for (n in c(1:30, 40, 60, 101, 202, 303, 404, 606, 1010)) {
  cases[[length(cases) + 1]] <- runif(n)                              # (0, 1)
  cases[[length(cases) + 1]] <- runif(n) * 10^runif(n, -8, 8)         # wide range
  s <- runif(n); s[sample(n, max(1, n %/% 5))] <- 0                   # zeros
  cases[[length(cases) + 1]] <- s
  s <- runif(n); s[sample(n, max(1, n %/% 7))] <- 1e-12               # below the tolerance
  cases[[length(cases) + 1]] <- s
  cases[[length(cases) + 1]] <- rep(c(0.3, 0.7), length.out = n)      # ties
  s <- runif(n) - 0.3                                                 # negative entries
  cases[[length(cases) + 1]] <- s
  cases[[length(cases) + 1]] <- 1 - runif(n)^8                         # survival-like, near 1
}
cases[[length(cases) + 1]] <- numeric(5)                               # all zero
ok <- vapply(cases, function(s) {
  ref <- MASS::ginv(Matrix::diag(s, length(s)))
  identical(ginv_diag(s), ref)
}, logical(1))
cat(sum(ok), "of", length(ok), "identical\n")
if (!all(ok)) {
  bad <- which(!ok)[1:min(5, sum(!ok))]
  for (b in bad) {
    s <- cases[[b]]; r <- MASS::ginv(Matrix::diag(s, length(s))); g <- ginv_diag(s)
    cat("n =", length(s), " max |diff| =", max(abs(r - g)), "\n")
  }
}
