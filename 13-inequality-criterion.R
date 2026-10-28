## Discussion, Proposition 8 and the conjecture after it: exhaustive check of the graphical
## criterion for when a c-degree 1 district has inequality constraints. Districts: members M1..Mk sharing one latent U, with exogenous
## binary external parents E1..Em. For each graph, compare
##   (a) whether the district is "saturated" after adding HLP edges, with
##   (b) whether the district polytope has inequality constraints beyond the probability
##       constraints and its own equality constraints (exact polytope comparison).
source("common.R")

## ---- graphical criterion ----------------------------------------------------------
## pa: named list of parent sets (members and externals; every member also has "U").
is_ancestor <- function(pa, a, b) {          # is a an ancestor of b?
  front <- pa[[b]]; seen <- character()
  while (length(front)) {
    x <- front[1]; front <- front[-1]
    if (x == a) return(TRUE)
    if (!(x %in% seen) && !is.null(pa[[x]])) { seen <- c(seen, x); front <- c(front, pa[[x]]) }
  }
  FALSE
}
saturated_after_hlp <- function(pa, members) {
  ## depth-first search over HLP edge additions; TRUE if some sequence reaches a saturated district
  key <- function(pa) paste(sapply(members, \(m) paste(sort(pa[[m]]), collapse = "")), collapse = "|")
  visited <- new.env()
  is_saturated <- function(pa) {
    ## all pairs adjacent and each member's parents contain all earlier members and their parents
    ord <- members[order(sapply(members, \(m) sum(sapply(members, \(o) o != m && is_ancestor(pa, o, m)))))]
    for (j in seq_along(ord)) {
      earlier <- ord[seq_len(j - 1)]
      need <- unique(c(earlier, unlist(pa[earlier])))
      if (!all(need %in% pa[[ord[j]]])) return(FALSE)
    }
    TRUE
  }
  dfs <- function(pa) {
    k <- key(pa); if (!is.null(visited[[k]])) return(FALSE); visited[[k]] <- TRUE
    if (is_saturated(pa)) return(TRUE)
    for (a in members) for (b in members) {
      if (a == b || a %in% pa[[b]] || b %in% pa[[a]]) next
      if (all(pa[[a]] %in% pa[[b]]) && !is_ancestor(pa, b, a)) {
        pa2 <- pa; pa2[[b]] <- c(pa2[[b]], a)
        if (dfs(pa2)) return(TRUE)
      }
    }
    FALSE
  }
  dfs(pa)
}

## ---- polytope check -------------------------------------------------------------
has_inequalities <- function(cn, j) {
  H <- cn$polytope$Hrep[[j]]
  lev <- cn$functional_relations$problevels[[j]]
  dov <- cn$functional_relations$district_parms$poterms[[j]]$do
  n <- ncol(H) - 2
  blocks <- if (length(dov)) interaction(lev[, dov, drop = FALSE], drop = TRUE) else factor(rep(1, n))
  simplex <- t(sapply(levels(blocks), \(b) as.numeric(blocks == b)))
  if (is.null(dim(simplex))) simplex <- matrix(simplex, nrow = 1)
  Heq <- H[H[, 1] == 1, , drop = FALSE]
  Hrel <- rbind(Heq, cbind(0, 0, diag(n)), cbind(1, -1, simplex))      # equalities + simplices
  Vrel <- rcdd::scdd(rcdd::d2q(Hrel), representation = "H")$output
  Vrel <- rcdd::q2d(Vrel[Vrel[, 1] == "0", -(1:2), drop = FALSE])
  ## the relaxed polytope equals ours iff all its vertices satisfy our inequalities
  !all(apply(Vrel, 1, \(p) hrep_test(H, p)))
}

## ---- enumeration ----------------------------------------------------------------
results <- list()
run_case <- function(k, ext_targets, internal) {
  members <- paste0("M", 1:k); exts <- if (length(ext_targets)) paste0("E", seq_along(ext_targets)) else character()
  pa <- setNames(lapply(members, \(m) "U"), members)
  for (e in seq_along(ext_targets)) for (t in ext_targets[[e]]) pa[[members[t]]] <- c(pa[[members[t]]], exts[e])
  for (ed in internal) pa[[members[ed[2]]]] <- c(pa[[members[ed[2]]]], members[ed[1]])
  edges <- unlist(lapply(members, \(m) {
    p <- setdiff(pa[[m]], "U")
    if (length(p)) paste0(p, " %-->% ", m) else character()
  }))
  edges <- c(edges, paste0("U %-->% ", paste(members, collapse = " + ")))
  g <- dhvm(eval(parse(text = paste0("caugi(", paste(edges, collapse = ", "), ")"))))
  cn <- meraconstraints(g)
  j <- which(sapply(cn$functional_relations$district_parms$poterms, \(x) setequal(x$out, members)))
  if (!length(j)) return(NULL)
  data.frame(graph = paste(edges, collapse = "; "), saturated = saturated_after_hlp(pa, members),
             inequalities = has_inequalities(cn, j))
}

subsets <- function(k) unlist(lapply(1:k, \(s) combn(k, s, simplify = FALSE)), recursive = FALSE)
cases <- list()
for (k in 2:3) {
  pairs <- combn(k, 2, simplify = FALSE)
  internals <- unlist(lapply(0:length(pairs), \(s) combn(length(pairs), s, \(i) pairs[i], simplify = FALSE)), recursive = FALSE)
  for (m in 0:(if (k == 2) 2 else 1)) {
    ext_choices <- if (m == 0) list(list()) else apply(expand.grid(rep(list(seq_along(subsets(k))), m)), 1, \(r) subsets(k)[r], simplify = FALSE)
    for (et in ext_choices) for (int in internals) cases[[length(cases) + 1]] <- list(k = k, et = et, int = int)
  }
}
cat("Number of districts:", length(cases), "\n")
for (i in seq_along(cases)) {
  cs <- cases[[i]]
  r <- tryCatch(run_case(cs$k, cs$et, cs$int), error = function(e) data.frame(graph = conditionMessage(e), saturated = NA, inequalities = NA))
  if (!is.null(r)) { results[[length(results) + 1]] <- r; cat(sprintf("%3d %-5s %-5s %s\n", i, r$saturated, r$inequalities, r$graph)); flush.console() }
}
res <- do.call(rbind, results)
cat("\nSaturated after adding HLP edges vs. has inequality constraints:\n")
print(table(saturated = res$saturated, inequalities = res$inequalities, useNA = "ifany"))
write.csv(res, "inequality-criterion-results.csv", row.names = FALSE)
