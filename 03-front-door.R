## Section 3.3: the front door model. Proposition 6 finds no nontrivial
## constraints: the equality constraints are normalization and (O), and each inequality
## constraint is a nonnegativity constraint in a different form.
source("common.R")

fd <- dhvm(caugi(A %-->% M %-->% Y, U1 %-->% A + Y))
cn <- meraconstraints(fd)
print(cn)
j  <- "A_Y"
H  <- cn$polytope$Hrep[[j]]
cat("\nDistrict", j, ":", sum(H[, 1] == 0), "inequality and", sum(H[, 1] == 1), "equality constraints;",
    sum(cn$nontrivial[[j]]), "nontrivial\n")

section("Equality constraints, in terms of interventional probabilities")
eq <- H[, 1] == 1
cat(constraint_text(cn, j, nontrivial_only = FALSE, star = TRUE)[eq], sep = "\n")
section("... and in terms of observed probabilities")
cat(constraint_text(cn, j, nontrivial_only = FALSE)[eq], sep = "\n")

section("The entry p_k that makes each inequality a nonnegativity constraint")
lev <- entry_levels(cn, j)
k   <- nonnegativity_entries(cn, j)
print(data.frame(constraint = constraint_text(cn, j, nontrivial_only = FALSE, star = TRUE)[!eq],
                 lev[k, c("A", "Y", "M")], row.names = NULL), right = FALSE)
cat("Distinct entries:", length(unique(k)), "of", length(k), "\n")

section("The screen that ignores the equality constraints wrongly flags these")
old <- prob_simplex_screen(cn$functional_relations, cn$polytope$Hrep)
cat(constraint_text(cn, j, nontrivial_only = FALSE, star = TRUE)[old[[j]]], sep = "\n")
