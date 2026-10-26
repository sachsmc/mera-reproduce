## Section 4.1: sequential instrumental variables (Figure 2).
source("common.R")

sq <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4 %-->% V5, U1 %-->% V2 + V3, U2 %-->% V4 + V5))
cn <- meraconstraints(sq)
print(cn)
section("Conditional independencies")
print(dhvm_list_independences(sq))
for (j in seq_along(cn$polytope$Hrep)) {
  section(paste("Flagged constraints, district", names(cn$polytope$Hrep)[j]))
  cat(constraint_text(cn, j), sep = "\n")
}
