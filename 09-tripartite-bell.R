## Section 4.3: the tripartite Bell graph (Figure 6). This takes up to about half an hour.
source("common.R")

bell <- dhvm(caugi(X %-->% V1, Y %-->% V2, Z %-->% V3, U %-->% V1 + V2 + V3))
t <- system.time(cn <- meraconstraints(bell))
print(t)
H <- cn$polytope$Hrep[[1]]
cat("Total constraints:", nrow(H), "  inequalities:", sum(H[, 1] == 0),
    "  equalities:", sum(H[, 1] == 1), "\n")
fl <- cn$trivial_if_false[[1]]
cat("Flagged as potentially nontrivial:", sum(fl), "  of which equalities:", sum(fl & H[, 1] == 1), "\n")
section("Flagged equality constraints")
cat(constraint_text(cn, 1, star = TRUE)[H[fl, 1] == 1], sep = "\n")
