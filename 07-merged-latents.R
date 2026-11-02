## Section 5: the graph of Figure 8, with c-degree 2, after merging U1 and U3.
## First district {V1, V3, V6}: no nontrivial constraints.
## Second district {V2, V4, V5}: 1096 nontrivial inequality constraints and 6 linearly
## independent nontrivial equality constraints, and a strictly positive distribution that
## satisfies the conditional independencies and all equality constraints but violates 13
## inequality constraints.
source("common.R")

g  <- dhvm(caugi(V1 %-->% V2 + V6, V2 %-->% V3, V3 %-->% V4, V4 %-->% V5, V6 %-->% V3,
                 U1 %-->% V2 + V4, U2 %-->% V1 + V3 + V6, U3 %-->% V4 + V5))
print(dhvm_c_degree(g))
gm <- dhvm_merge_latents(g)
cn <- meraconstraints(gm)
print(cn)

counts <- sapply(names(cn$polytope$Hrep), \(j) { H <- cn$polytope$Hrep[[j]]; f <- cn$nontrivial[[j]]
  c(ineq = sum(H[, 1] == 0), nontrivial_ineq = sum(f & H[, 1] == 0),
    eq = sum(H[, 1] == 1), nontrivial_eq = sum(f & H[, 1] == 1)) })
section("Constraint counts (Proposition 6)")
print(counts)

d1 <- "V1_V3_V6"
H1 <- cn$polytope$Hrep[[d1]]
section("First district: equality constraints, in terms of interventional probabilities")
cat(constraint_text(cn, d1, nontrivial_only = FALSE, star = TRUE)[H1[, 1] == 1], sep = "\n")
section("First district: inequality constraints that are nonnegativity constraints")
k <- nonnegativity_entries(cn, d1)
cat(sum(!is.na(k)), "of", length(k), "inequality constraints\n")

section("First district: the inequality displayed in Section 5")
lev <- entry_levels(cn, d1)
a <- with(lev, as.numeric(V1 == 0 & V3 == 0 & V6 == 0 & V2 == 1) - as.numeric(V1 == 0 & V6 == 0 & V2 == 0))
V <- rcdd::q2d(cn$polytope$Vrep[[d1]])[, -(1:2)]
E <- H1[H1[, 1] == 1, -1, drop = FALSE]
b <- with(lev, as.numeric(V1 == 0 & V3 == 1 & V6 == 0 & V2 == 1))
cat("Largest value over the vertices of the left side minus the right side:", max(V %*% a), "\n")
cat("Adding equality constraints turns it into P*(V1 = 0, V3 = 1, V6 = 0 | V2 = 1) >= 0:",
    qr(rbind(E, c(0, a + b)))$rank == qr(E)$rank, "\n")

## Example distribution: every variable except V4 takes the value 1 with probability 1/2
## given its parents, and P(V4 = v4 | V1, V2, V3) = 1/8 + (3/4) a(V1, V2, V3, v4), where
## a = (1 - V1*V3) * 1{V2 == v4} + V1*V3 * 1{V2 != v4}.
pt <- data.table(expand.grid(V1 = 0:1, V2 = 0:1, V3 = 0:1, V4 = 0:1, V5 = 0:1, V6 = 0:1))
pt[, a := (1 - V1 * V3) * (V2 == V4) + V1 * V3 * (V2 != V4)]
pt[, joint := (1/2)^5 * (1/8 + (3/4) * a)]
section("Example distribution")
cat("Sum of probabilities:", sum(pt$joint), "  smallest probability:", as.character(MASS::fractions(min(pt$joint))), "\n")
cat("Conditional independencies hold:", independences_hold(cn, pt[, .(V1, V2, V3, V4, V5, V6, joint)]), "\n")

## District factor of {V2, V4, V5}: P(v2 | v1) P(v4 | v1, v2, v3) P(v5 | v1, v2, v3, v4)
pt[, p1 := sum(joint), by = V1]
pt[, p12 := sum(joint), by = .(V1, V2)]
pt[, p123 := sum(joint), by = .(V1, V2, V3)]
pt[, p1234 := sum(joint), by = .(V1, V2, V3, V4)]
pt[, p12345 := sum(joint), by = .(V1, V2, V3, V4, V5)]
pt[, pstar := (p12 / p1) * (p1234 / p123) * (p12345 / p1234)]
pu <- unique(pt[, .(V1, V2, V3, V4, V5, a, pstar)])
cat("Inner sum for each (v1, v3):\n")
print(pu[, .(inner = sum(a * pstar)), by = .(V1, V3)])
cat("Left side of the inequality of Clauser et al.:", pu[, sum(a * pstar)], "(bound 3)\n")
cat("P*(V4 = 1 | v1, v3) for each (v1, v3):", pu[V4 == 1, .(p = sum(pstar)), by = .(V1, V3)]$p, "\n")

## District factor of {V1, V3, V6}: P(v1) P(v6 | v1) P(v3 | v1, v6, v2)
pt[, p6 := sum(joint), by = .(V1, V6)]
pt[, p126 := sum(joint), by = .(V1, V2, V6)]
pt[, p1236 := sum(joint), by = .(V1, V2, V3, V6)]
pt[, pstar1 := p6 * p1236 / p126]
pv1 <- as_entries(cn, d1, unique(pt[, .(V1, V2, V3, V6, pstar1)]), "pstar1")
ok1 <- hrep_test(H1, pv1, which = TRUE)
cat("First district: all constraints hold:", all(ok1), "\n")

d2  <- "V2_V4_V5"
H   <- cn$polytope$Hrep[[d2]]
pv  <- as_entries(cn, d2, pu, "pstar")
ok  <- hrep_test(H, pv, which = TRUE)
cat("Second district: equality constraints satisfied:", all(ok[H[, 1] == 1]), "\n")
cat("Second district: inequality constraints violated:", sum(!ok), "\n")
cat("All violated constraints are nontrivial:", all(cn$nontrivial[[d2]][!ok]), "\n")
txt <- hrep_to_text(H, cn$functional_relations$probstarform[[d2]])[!ok]
section("Shortest violated constraint (as printed in Section 5)")
cat(txt[which.min(nchar(txt))], "\n")

## The condition displayed in Section 5: the margin of (V4, V5) under do(V1, V3) does not
## depend on V1. Together with (N) and (O) it spans all equality constraints of the district.
section("Second district: the Verma-type constraint")
rel <- cn$functional_relations
lv2 <- meraconstraints:::district_levels(rel, 2)
mem <- rel$district_parms$poterms[[2]]$out
ext <- rel$district_parms$poterms[[2]]$do
ext_pa <- lapply(setNames(mem, mem), \(m) intersect(parents(gm$caugi_graph, nodes = m), ext))
obs <- subset(gm$nodes_data, latent == FALSE)$name
anc <- lapply(setNames(obs, obs), \(v) intersect(ancestors(gm$caugi_graph, nodes = v), obs))
N <- meraconstraints:::indicator_rows(lv2, ext)
O <- do.call(rbind, lapply(meraconstraints:::closed_member_sets(mem, anc), meraconstraints:::margin_rows,
                           lev = lv2, ext = ext, ext_pa = ext_pa))
M <- meraconstraints:::margin_rows(lv2, c("V4", "V5"), ext, ext_pa)
Ee <- -H[H[, 1] == 1, -(1:2), drop = FALSE]
cat("rank(N, O):", qr(rbind(N, O))$rank, "  rank(N, O, margin):", qr(rbind(N, O, M))$rank,
    "  rank(N, O, H_e):", qr(rbind(N, O, Ee))$rank, "\n")
