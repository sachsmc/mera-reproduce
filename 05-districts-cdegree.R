## Section 2.2: districts and c-degree of the graph in Figure 1.
source("common.R")

g1 <- dhvm(caugi(A %-->% B + E, B %-->% C, C %-->% D, D %-->% F, E %-->% C, G %-->% F,
                 U1 %-->% B + D, U2 %-->% A + E + C, U3 %-->% F + D))
print(g1)
section("Districts and their c-degrees")
print(dhvm_c_degree(g1))
