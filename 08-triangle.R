## Section 5: the triangle graph (Figure 8). After merging the three unobserved variables,
## the algorithm finds no nontrivial constraints.
source("common.R")

tri <- dhvm(caugi(U1 %-->% V1 + V3, U2 %-->% V1 + V2, U3 %-->% V2 + V3))
print(dhvm_c_degree(tri))
cn <- meraconstraints(dhvm_merge_latents(tri))
print(cn)
cat("\nNumber of constraints flagged as potentially nontrivial:",
    sum(unlist(cn$trivial_if_false)), "\n")
