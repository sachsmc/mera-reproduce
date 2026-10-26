## Section 3.3: the front door model (Figure 4). The constraints flagged as potentially
## nontrivial are trivial once written in terms of observed probabilities.
source("common.R")

fd <- dhvm(caugi(A %-->% M %-->% Y, U1 %-->% A + Y))
cn <- meraconstraints(fd)
print(cn)
j <- which(sapply(cn$trivial_if_false, any))
section("Flagged constraints, in terms of observed probabilities")
cat(constraint_text(cn, j), sep = "\n")
section("The same constraints, in terms of interventional probabilities")
cat(constraint_text(cn, j, star = TRUE), sep = "\n")
