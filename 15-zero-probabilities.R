## Section 3.2, Remark 2: distributions with zero probabilities. Each product of Proposition 1
## is read from left to right, and the first factor that is not positive decides the entry:
## if it is zero the entry is zero, and if it is undefined the entry is free. Then P is in the
## model if and only if it satisfies the conditional independencies and p = B r holds on the
## entries that are not free, for some probability vector r.
source("common.R")
set.seed(20261102)
zero <- 1e-13

setup <- function(g) {
  cn   <- meraconstraints(g)
  obs  <- subset(g$nodes_data, latent == FALSE)$name
  ord  <- intersect(topological_sort(g$caugi_graph), obs)
  grid <- do.call(expand.grid, lapply(setNames(obs, obs), \(v) seq_len(subset(g$nodes_data, name == v)$nvals) - 1))
  list(g = g, cn = cn, obs = obs, ord = ord, grid = grid)
}

## The entries of p for district j under the rule of Remark 2, with NA for a free entry
pstar <- function(s, prob, j) {
  lev <- entry_levels(s$cn, j)
  pl  <- s$cn$functional_relations$district_parms$plist[[j]]
  mem <- s$ord[s$ord %in% names(pl)]
  vapply(seq_len(nrow(lev)), \(r) {
    val <- 1
    for (m in mem) {
      mask <- rep(TRUE, nrow(s$grid))
      for (w in pl[[m]]) mask <- mask & s$grid[[w]] == lev[r, w]
      den <- sum(prob[mask])
      if (den < zero) return(NA_real_)
      f <- sum(prob[mask & s$grid[[m]] == lev[r, m]]) / den
      if (f < zero) return(0)
      val <- val * f
    }
    val
  }, numeric(1))
}

## Smallest t such that |B r - p| <= t on the entries that are not free, by linear programming
gap <- function(s, p, j) {
  B  <- t(rcdd::q2d(s$cn$polytope$Vrep[[j]])[, -(1:2), drop = FALSE])
  ok <- !is.na(p)
  if (!any(ok)) return(0)
  Bs <- B[ok, , drop = FALSE]; ps <- p[ok]; nc <- ncol(B)
  H  <- rbind(cbind(0, ps, -Bs, 1), cbind(0, -ps, Bs, 1),   # b - A x >= 0, x = (r, t)
              cbind(0, 0, diag(nc + 1)), c(1, 1, rep(-1, nc), 0))
  rcdd::q2d(rcdd::lpcdd(rcdd::d2q(unname(H)), rcdd::d2q(c(rep(0, nc), 1)), minimize = TRUE)$optimal.value)
}
in_model <- function(s, prob, tol = 1e-8)
  all(vapply(seq_along(s$cn$polytope$Hrep), \(j) gap(s, pstar(s, prob, j), j), numeric(1)) <= tol)

sparse_r <- function(n) {
  k <- sample(1:min(4, n), 1); w <- rgamma(k, 1); r <- numeric(n)
  r[sample(n, k)] <- w / sum(w); r
}

graphs <- list(
  "IV"          = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y)),
  "Front door"  = dhvm(caugi(A %-->% M %-->% Y, U1 %-->% A + Y)),
  "Verma graph" = dhvm(caugi(A %-->% B %-->% C %-->% D, U %-->% B + D)),
  "Fig 4 right" = dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4)),
  "Bell"        = dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B)),
  "Sequential"  = dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4 %-->% V5, U2 %-->% V2 + V3, U3 %-->% V4 + V5))
)
sets <- lapply(graphs, setup)

## ---- 1. Distributions in the model, from sparse response function distributions ------
## The entries that are not free must equal the true district factor, and the linear
## program must accept.
section("Distributions in the model, with many zero probabilities")
necessity <- function(s, nsim = 50) {
  mism <- 0; rej <- 0; nfree <- 0; nzero <- 0
  poterms <- s$cn$functional_relations$district_parms$poterms
  roots <- setdiff(s$obs, unlist(lapply(poterms, `[[`, "out")))
  for (it in seq_len(nsim)) {
    prob <- rep(1, nrow(s$grid)); truth <- list()
    for (j in seq_along(s$cn$polytope$Hrep)) {
      B <- t(rcdd::q2d(s$cn$polytope$Vrep[[j]])[, -(1:2), drop = FALSE])
      truth[[j]] <- as.vector(B %*% sparse_r(ncol(B)))
      lev  <- entry_levels(s$cn, j)
      key  <- do.call(paste, lev[names(lev)])
      gkey <- do.call(paste, s$grid[names(lev)])
      prob <- prob * truth[[j]][match(gkey, key)]
    }
    for (v in roots) {
      nv <- subset(s$g$nodes_data, name == v)$nvals
      pv <- rgamma(nv, 1); if (runif(1) < .3) pv[sample(nv, 1)] <- 0
      if (sum(pv) == 0) pv[1] <- 1
      prob <- prob * (pv / sum(pv))[s$grid[[v]] + 1]
    }
    stopifnot(abs(sum(prob) - 1) < 1e-10)
    stopifnot(independences_hold(s$cn, cbind(s$grid, joint = prob)))
    for (j in seq_along(s$cn$polytope$Hrep)) {
      p <- pstar(s, prob, j); ok <- !is.na(p)
      nfree <- nfree + sum(!ok); nzero <- nzero + sum(p[ok] == 0)
      if (any(abs(p[ok] - truth[[j]][ok]) > 1e-8)) mism <- mism + 1
      if (gap(s, p, j) > 1e-8) rej <- rej + 1
    }
  }
  c(wrong_entries = mism, rejected = rej, free_entries = nfree, zero_entries = nzero)
}
print(t(sapply(sets, necessity)))

## ---- 2. Verma graph with C = B ----------------------------------------------------------
## Then D = g(B, U), so (A, B, D) is an instrumental variable model with instrument A, and P
## is in the model if and only if P(b, d | a) satisfies the instrumental inequalities.
section("Verma graph with C = B, compared with the instrumental inequalities")
sv <- sets[["Verma graph"]]; si <- sets[["IV"]]
Hiv <- si$cn$polytope$Hrep[[1]]; levi <- entry_levels(si$cn, 1)
iv_excess <- function(pabd) {
  p <- mapply(\(x, y, z) pabd[z + 1, x + 1, y + 1], levi$X, levi$Y, levi$Z)
  max(-(Hiv[Hiv[, 1] == 0, 2] + Hiv[Hiv[, 1] == 0, -(1:2), drop = FALSE] %*% p))
}
tab <- matrix(0, 2, 2, dimnames = list(remark = c("rejects", "accepts"), iv = c("violated", "satisfied")))
pa <- c(.4, .6)
for (it in 1:500) {
  pabd <- array(0, c(2, 2, 2))
  for (a in 1:2) {
    x <- rgamma(4, .4); x[runif(4) < .35] <- 0
    if (sum(x) == 0) x[sample(4, 1)] <- 1
    pabd[a, , ] <- matrix(x / sum(x), 2, 2)
  }
  ex <- iv_excess(pabd)
  if (abs(ex) < 1e-6 && ex != 0) next
  prob <- with(sv$grid, ifelse(C == B, pa[A + 1] * mapply(\(a, b, d) pabd[a + 1, b + 1, d + 1], A, B, D), 0))
  i <- if (in_model(sv, prob)) "accepts" else "rejects"
  j <- if (ex > 1e-9) "violated" else "satisfied"
  tab[i, j] <- tab[i, j] + 1
}
print(tab)

## ---- 3. Settings that are never observed ----------------------------------------------
section("Settings that are never observed")
acc <- 0
for (it in 1:100) {
  x <- rgamma(4, .5); x <- x / sum(x)
  acc <- acc + in_model(si, with(si$grid, ifelse(Z == 0, x[X + 2 * Y + 1], 0)))
}
cat("IV with P(Z = 1) = 0, any P(x, y | Z = 0): accepted", acc, "of 100\n")

sf <- sets[["Front door"]]; acc <- 0
for (it in 1:100) {
  x <- rgamma(4, .5); x <- x / sum(x)
  acc <- acc + in_model(sf, with(sf$grid, ifelse(M == A, x[A + 2 * Y + 1], 0)))
}
cat("Front door with M = A, any P(a, y): accepted", acc, "of 100\n")

sb <- sets[["Bell"]]
res <- matrix(0, 2, 2, dimnames = list(remark = c("rejects", "accepts"), truth = c("A, Y dependent", "A, Y independent")))
for (it in 1:200) {
  if (it %% 2) {
    pA <- runif(1); pB <- matrix(runif(4), 2, 2)
    f <- function(a, b, y) (if (a) pA else 1 - pA) * (if (b) pB[a + 1, y + 1] else 1 - pB[a + 1, y + 1])
    tr <- "A, Y independent"
  } else {
    arr <- array(rgamma(8, 1), c(2, 2, 2)); arr <- arr / rep(apply(arr, 3, sum), each = 4)
    f <- function(a, b, y) arr[a + 1, b + 1, y + 1]
    tr <- "A, Y dependent"
  }
  prob <- with(sb$grid, ifelse(X == 0, .5 * mapply(f, A, B, Y), 0))
  i <- if (in_model(sb, prob)) "accepts" else "rejects"
  res[i, tr] <- res[i, tr] + 1
}
cat("Bell graph with P(X = 1) = 0:\n")
print(res)
