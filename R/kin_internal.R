# Internal helpers shared by the kin functions: thin wrappers around the C++
# projection kernel and the builders that turn kin matrices into the long
# data frames the package returns.
#
# The data frame builders reproduce exactly (values, column order, row order,
# column types and class) what the original reshape (pivot / melt / dcast /
# merge) pipeline produced, without materialising the intermediate long tables.

kin_names_14 <- c("d", "gd", "ggd", "m", "gm", "ggm", "os", "ys",
                  "nos", "nys", "oa", "ya", "coa", "cya")

# Focal's ageing operator in age-only models: a shift within the living rows,
# so that phi[, i] is the unit vector of age i-1 and F %*% phi[, i] = F[, i].
focal_shift <- function(n, nlive) {
  G <- matrix(0, n, n)
  if (nlive > 1) G[cbind(2:nlive, 1:(nlive - 1))] <- 1
  G
}

# One-sex projection: Focal's whole trip (prev = NULL) or one time step.
project_kin_1sex <- function(Ut, ft, pit, ages, prev = NULL) {
  n <- 2 * ages
  pit <- pit[1:ages]
  k <- kin_project_cpp(Ut, ft, ft, focal_shift(n, ages),
                       phi0 = c(1, numeric(n - 1)), m0 = pit,
                       pi_all = pit, pi_os = pit,
                       ncol = ages, nlive = ages, prev = prev)
  k[kin_names_14]
}

# Layouts: which matrix rows are living kin (in the row order of the output)
# and which rows hold their deaths, plus the id columns that label them.
layout_1sex <- function(age) {
  ages <- length(age)
  list(ids = data.frame(age_kin = age),
       live = seq_len(ages),
       dead = list(dead = ages + seq_len(ages)))
}

# Sort kin names like merge()/data.table did: plain ASCII, locale-independent.
sort_kin <- function(x) sort(x, method = "radix")

# Values of one kin matrix block in output row order: id rows (lay$live or a
# death row set) by Focal age, Focal age fastest.
block_values <- function(x, rows, cols) as.vector(t(x[rows, cols, drop = FALSE]))

# Time-invariant long data frame: one block per kin (sorted), rows ordered by
# the id columns then Focal's age; deaths relocated one Focal age earlier
# (NA at the last age), as merge(all.x = TRUE) did.
build_kin_df_ti <- function(kin_list, age, lay) {
  ages <- length(age)
  nid <- nrow(lay$ids)
  kin_sorted <- sort_kin(names(kin_list))
  nk <- length(kin_sorted)
  nblock <- nid * ages
  out <- data.frame(kin = rep(kin_sorted, each = nblock), stringsAsFactors = FALSE)
  ids_rep <- lay$ids[rep(rep(seq_len(nid), each = ages), nk), , drop = FALSE]
  rownames(ids_rep) <- NULL
  out <- cbind(out, ids_rep)
  out$age_focal <- rep(as.integer(age), times = nid * nk)
  out$living <- unlist(lapply(kin_sorted, function(k) block_values(kin_list[[k]], lay$live, seq_len(ages))), use.names = FALSE)
  for (dn in names(lay$dead)) {
    out[[dn]] <- unlist(lapply(kin_sorted, function(k) {
      x <- kin_list[[k]]
      as.vector(t(cbind(x[lay$dead[[dn]], -1, drop = FALSE], NA_real_)))
    }), use.names = FALSE)
  }
  out
}

# Time-varying long data table: `store` is a list by year (character names,
# ascending) of list(ages = selected Focal ages, kin = list of kin matrices
# restricted to those ages). Deaths at (year + step, age + step) are relocated
# to (year, age) when that cell was selected, NA otherwise.
build_kin_df_tv <- function(store, year_step, lay) {
  years <- as.integer(names(store))
  nid <- nrow(lay$ids)
  ny <- length(store)
  year_col <- kin_col <- af_col <- liv_col <- id_idx <- vector("list", ny)
  dead_cols <- lapply(lay$dead, function(d) vector("list", ny))
  for (yi in seq_len(ny)) {
    Y <- years[yi]
    A <- store[[yi]]$ages
    nA <- length(A)
    nxt <- match(as.character(Y + year_step), names(store))
    # column (in next year's stored matrices) holding deaths for each selected age
    dcols <- if (is.na(nxt)) rep(NA_integer_, nA) else match(A + year_step, store[[nxt]]$ages)
    ok <- !is.na(dcols)
    kin_sorted <- sort_kin(names(store[[yi]]$kin))
    nk <- length(kin_sorted)
    nblock <- nid * nA
    year_col[[yi]] <- rep(Y, nblock * nk)
    kin_col[[yi]] <- rep(kin_sorted, each = nblock)
    id_idx[[yi]] <- rep(rep(seq_len(nid), each = nA), nk)
    af_col[[yi]] <- rep(A, times = nid * nk)
    liv_col[[yi]] <- unlist(lapply(kin_sorted, function(k) block_values(store[[yi]]$kin[[k]], lay$live, seq_len(nA))), use.names = FALSE)
    for (dn in names(lay$dead)) {
      dead_cols[[dn]][[yi]] <- unlist(lapply(kin_sorted, function(k) {
        dcol <- matrix(NA_real_, length(lay$dead[[dn]]), nA)
        if (any(ok)) dcol[, ok] <- store[[nxt]]$kin[[k]][lay$dead[[dn]], dcols[ok], drop = FALSE]
        as.vector(t(dcol))
      }), use.names = FALSE)
    }
  }
  out <- data.frame(year = unlist(year_col, use.names = FALSE), kin = unlist(kin_col, use.names = FALSE),
                    stringsAsFactors = FALSE)
  ids_rep <- lay$ids[unlist(id_idx, use.names = FALSE), , drop = FALSE]
  rownames(ids_rep) <- NULL
  out <- cbind(out, ids_rep)
  out$age_focal <- unlist(af_col, use.names = FALSE)
  out$cohort <- out$year - out$age_focal
  out$living <- unlist(liv_col, use.names = FALSE)
  for (dn in names(lay$dead)) out[[dn]] <- unlist(dead_cols[[dn]], use.names = FALSE)
  data.table::setDT(out)
  out
}

# Selected Focal ages (ascending) per year, from output_period_cohort_combination()
selected_ages_by_year <- function(out_selected) {
  sp <- split(out_selected$age, out_selected$year)
  lapply(sp, function(a) sort(unique(a)))
}
