## Section 5.1: Propositions 5 (HLP) and 6 (Strong Face Splitting). The Figure 6 example, the
## Figure 7 district that neither rule can reduce, and the half Bell graph (not in the paper).
source("common.R")

latents <- function(g) subset(g$nodes_data, latent == TRUE)$name
observed <- function(g) subset(g$nodes_data, latent == FALSE)$name
pa_all <- function(g, v) parents(g$caugi_graph, nodes = v)
pa_obs <- function(g, v) intersect(unlist(lapply(v, \(x) pa_all(g, x))), observed(g))
ch_obs <- function(g, u) intersect(children(g$caugi_graph, nodes = u), observed(g))

## Proposition 5: may the edge a -> b be added (absent) or removed (present)?
hlp_applies <- function(g, a, b) {
  all(pa_all(g, a) %in% setdiff(pa_all(g, b), a)) &&             # Pa(a) is a subset of Pa(b)
    length(intersect(pa_all(g, a), latents(g))) > 0 &&           # so they share a latent parent
    !(b %in% ancestors(g$caugi_graph, nodes = a))                # and the edge makes no cycle
}

## Proposition 6: the two conditions for splitting the unobserved vertices Us, with D a
## nonempty subset of their common children
split_conditions <- function(g, Us, D) {
  C <- lapply(Us, \(u) setdiff(ch_obs(g, u), D))
  cond1 <- all(vapply(C, \(Ci) all(vapply(D, \(d) all(c(pa_obs(g, Ci), Ci) %in% pa_obs(g, d)), logical(1))), logical(1)))
  allC <- unique(unlist(C))
  cond2 <- all(vapply(setdiff(latents(g), Us), \(u) {
    ch <- ch_obs(g, u)
    !any(ch %in% allC) || any(vapply(C, \(Ci) all(ch %in% Ci), logical(1))) || all(D %in% ch)
  }, logical(1)))
  c(cond1 = cond1, cond2 = cond2)
}

section("Figure 6")
left   <- dhvm(caugi(X %-->% Y, U1 %-->% X + Z, U2 %-->% X + Y))
center <- dhvm(caugi(Z %-->% X, X %-->% Y, U1 %-->% X + Z, U2 %-->% X + Y))
right  <- dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y))
cat("c-degree of the left, center and right DAGs:", max(dhvm_c_degree(left)), max(dhvm_c_degree(center)), max(dhvm_c_degree(right)), "\n")
cat("Proposition 5 adds Z -> X to the left DAG:", hlp_applies(left, "Z", "X"), "\n")
cat("Proposition 6 applied directly to the left DAG, U1, D = {X}:", split_conditions(left, "U1", "X"), "\n")
cat("Proposition 6 applied to the center DAG, U1, D = {X}:", split_conditions(center, "U1", "X"), "\n")
cat("After splitting, both new unobserved vertices have a single child, which leaves the right DAG.\n")
cat("Conditional independencies of the three DAGs:",
    length(dhvm_list_independences(left)), length(dhvm_list_independences(center)), length(dhvm_list_independences(right)), "\n")

section("Figure 7: the district {V2, V4, V5}")
g <- dhvm(caugi(V1 %-->% V2 + V6, V2 %-->% V3, V3 %-->% V4, V4 %-->% V5, V6 %-->% V3,
                U1 %-->% V2 + V4, U2 %-->% V1 + V3 + V6, U3 %-->% V4 + V5))
mem <- c("V2", "V4", "V5")
pairs <- expand.grid(a = mem, b = mem, stringsAsFactors = FALSE)
pairs <- pairs[pairs$a != pairs$b, ]
pairs$hlp <- mapply(\(a, b) hlp_applies(g, a, b), pairs$a, pairs$b)
print(pairs, row.names = FALSE)
cat("\nEvery way of splitting U1, U3 or both, with C_i not all empty:\n")
for (Us in list("U1", "U3", c("U1", "U3"))) {
  common <- Reduce(intersect, lapply(Us, \(u) ch_obs(g, u)))
  subs <- unlist(lapply(seq_along(common), \(k) combn(common, k, simplify = FALSE)), recursive = FALSE)
  for (D in subs) {
    if (all(vapply(Us, \(u) all(ch_obs(g, u) %in% D), logical(1)))) next
    cat(sprintf("  %-8s D = {%s}: ", paste(Us, collapse = ", "), paste(D, collapse = ", ")))
    print(split_conditions(g, Us, D))
  }
}

section("Half of a Bell scenario (not in the paper), W -> A <- U -> B")
hb <- dhvm(caugi(W %-->% A, U %-->% A + B))
cat("Proposition 5 adds B -> A:", hlp_applies(hb, "B", "A"), "\n")
cn <- meraconstraints(hb)
print(cn)
