## Section 3 and Appendix C: the instrumental variable model (Figure 3).
source("common.R")

iv <- dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y))
cn <- meraconstraints(iv)
print(cn)

section("Matrix B (rows: P*(x, y | z); columns: response function variable values)")
rel <- dhvm_linear_relations(iv)
print(t(rcdd::q2d(dhvm_relations_to_polytope(rel)$Vrep[[1]])[, -(1:2)]))

section("H-representation (rcdd format: columns are [equality flag, b, -A])")
print(cn$polytope$Hrep[[1]])

section("All constraints, in terms of observed probabilities")
cat(constraint_text(cn, 1, flagged_only = FALSE), sep = "\n")

section("Rows flagged as potentially nontrivial (Proposition 4)")
print(which(cn$trivial_if_false[[1]]))
cat(constraint_text(cn, 1), sep = "\n")
