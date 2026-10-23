# Two-sex helpers (with or without causes of death): projection wrapper and
# row layouts. The matrix block structure follows Caswell (2022) and
# Caswell, Margolis & Verdery (2023): living rows are [females, males], death
# rows follow, by sex, and (for causes of death) by age with cause fastest.

project_kin_2sex <- function(Ut, Ft, Ft_star, pit, sex_focal, ages, prev = NULL) {
  agess <- 2 * ages
  n <- nrow(Ut)
  pif <- pit[1:ages]
  pim <- pit[(ages + 1):agess]
  # G matrix moves focal by age
  G <- matrix(0, nrow = ages, ncol = ages)
  G[row(G) - 1 == col(G)] <- 1
  Gt <- matrix(0, n, n)
  Gt[1:agess, 1:agess] <- as.matrix(Matrix::bdiag(G, G))
  # locate focal at age 0 depending sex
  phi0 <- numeric(n)
  phi0[ifelse(sex_focal == "f", 1, ages + 1)] <- 1
  # initial conditions on os and nos depend on focal sex
  pios <- if (sex_focal == "f") pif else pim
  k <- kin_project_cpp(Ut, Ft, Ft_star, Gt, phi0, m0 = pit,
                       pi_all = pif + pim, pi_os = pios,
                       ncol = ages, nlive = agess, prev = prev)
  k[kin_names_14]
}

# time invariant: rows sorted by (age_kin, sex_kin)
layout_2sex_ti <- function(age, causes = 0) {
  ages <- length(age)
  agess <- 2 * ages
  live <- as.vector(rbind(seq_len(ages), ages + seq_len(ages)))
  ids <- data.frame(age_kin = rep(age, each = 2), sex_kin = rep(c("f", "m"), ages), stringsAsFactors = FALSE)
  if (causes == 0) {
    dead <- list(dead = agess + live)
  } else {
    sex0 <- rep(0:1, ages)      # 0 female, 1 male, in `live` order
    a0 <- rep(age, each = 2)
    dead <- lapply(seq_len(causes), function(cc) agess + sex0 * ages * causes + a0 * causes + cc)
    names(dead) <- paste0("deadcause", seq_len(causes))
  }
  list(ids = ids, live = live, dead = dead)
}

# time varying: rows sorted by (sex_kin, age_kin); death columns in C-locale order (dcast)
layout_2sex_tv <- function(age, causes = 0) {
  ages <- length(age)
  agess <- 2 * ages
  live <- seq_len(agess)
  ids <- data.frame(sex_kin = rep(c("f", "m"), each = ages), age_kin = rep(age, 2), stringsAsFactors = FALSE)
  if (causes == 0) {
    dead <- list(dead = agess + live)
  } else {
    sex0 <- rep(0:1, each = ages)
    a0 <- rep(age, 2)
    dead <- lapply(seq_len(causes), function(cc) agess + sex0 * ages * causes + a0 * causes + cc)
    names(dead) <- paste0("deadcause", seq_len(causes))
    dead <- dead[sort(names(dead), method = "radix")]
  }
  list(ids = ids, live = live, dead = dead)
}
