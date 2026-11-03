## Shared setup for the reproducibility scripts.
## Requires the R packages meraconstraints (https://github.com/sachsmc/meraconstraints),
## caugi, rcdd, data.table and MASS.
suppressPackageStartupMessages({
  library(meraconstraints)
  library(data.table)
})
if (!"nontrivial_screen" %in% getNamespaceExports("meraconstraints"))
  stop("These scripts need a version of meraconstraints that has nontrivial_screen(). ",
       "Install the current master branch of https://github.com/sachsmc/meraconstraints.")

## Text of the constraints of district j, optionally only the nontrivial ones
## (Proposition 4), in terms of observed (star = FALSE) or interventional (star = TRUE)
## probabilities.
constraint_text <- function(cn, j = 1, nontrivial_only = TRUE, star = FALSE) {
  H    <- cn$polytope$Hrep[[j]]
  prob <- if (star) cn$functional_relations$probstarform[[j]] else cn$functional_relations$probform[[j]]
  txt  <- hrep_to_text(H, prob)
  if (nontrivial_only) txt[cn$nontrivial[[j]]] else txt
}

## Values of the variables of district j, one row per entry of p, in the order of the
## columns of the H-representation.
entry_levels <- function(cn, j = 1) {
  rel <- cn$functional_relations
  lev <- as.data.frame(rel$problevels[[j]])
  lev[match(names(unlist(rel$probform[[j]])), do.call(paste, c(as.list(lev), sep = "."))), , drop = FALSE]
}

## Column `value` of the table tab, which has one row per entry of p, in the order of the
## columns of the H-representation.
as_entries <- function(cn, j, tab, value) {
  lev  <- entry_levels(cn, j)
  tab  <- as.data.frame(tab)
  key  <- function(d) do.call(paste, c(as.list(d[, names(lev), drop = FALSE]), sep = "."))
  x    <- tab[[value]][match(key(lev), key(tab))]
  stopifnot(!anyNA(x))
  x
}

## For each inequality row of district j, the entry p_k that is zero at every vertex on its
## facet, or NA if there is none. By Proposition 4(ii), the row is then the nonnegativity
## constraint p_k >= 0 in a different form, and NA means that the row is nontrivial.
nonnegativity_entries <- function(cn, j = 1) {
  H  <- cn$polytope$Hrep[[j]]
  H  <- H[H[, 1] == 0, , drop = FALSE]
  V  <- rcdd::q2d(cn$polytope$Vrep[[j]])[, -(1:2), drop = FALSE]
  on <- abs(H[, 2] + H[, -(1:2), drop = FALSE] %*% t(V)) < 1e-9
  z  <- (on %*% V) == 0
  apply(z, 1, \(r) if (any(r)) which(r)[1] else NA_integer_)
}

## Do the conditional independencies of cn hold for the joint distribution in tab, a table
## of the observed variables with the probabilities in column `joint`?
independences_hold <- function(cn, tab, tol = 1e-10) {
  tab <- as.data.frame(tab)
  sum_over <- function(cols) ave(tab$joint, tab[cols], FUN = sum)
  all(vapply(cn$independences, \(ii) {
    a <- ii$indep[1]; b <- ii$indep[2]; z <- setdiff(ii$given, "")
    pz <- if (length(z)) sum_over(z) else rep(sum(tab$joint), nrow(tab))
    all(abs(sum_over(c(a, b, z)) * pz - sum_over(c(a, z)) * sum_over(c(b, z))) < tol)
  }, logical(1)))
}

section <- function(title) cat("\n==== ", title, " ====\n", sep = "")
