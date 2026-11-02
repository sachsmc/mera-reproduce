## Section 4.1: sequential instrumental variables (Figure 3).
source("common.R")

sq <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4 %-->% V5, U2 %-->% V2 + V3, U3 %-->% V4 + V5))
cn <- meraconstraints(sq)
print(cn)
section("Conditional independencies")
print(dhvm_list_independences(sq))
for (j in seq_along(cn$polytope$Hrep)) {
  section(paste("Nontrivial constraints, district", names(cn$polytope$Hrep)[j]))
  cat(constraint_text(cn, j), sep = "\n")
}
