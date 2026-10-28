# Reproducibility scripts

Scripts that reproduce the examples, numbers and tables in
"Deriving Complete Constraints in Hidden Variable Models".

## Requirements

- R (tested with R 4.4.2)
- R packages: `meraconstraints` (>= 0.3.0, https://github.com/sachsmc/meraconstraints),
  `caugi`, `rcdd`, `data.table`, `MASS`

## Usage

Run the scripts from this directory, either one at a time,

    Rscript 01-iv.R

or all together, saving each script's output to `output/`:

    Rscript run_all.R          # everything except the two slow scripts
    Rscript run_all.R slow     # everything

## Contents

| Script | Paper section | What it reproduces |
|---|---|---|
| `01-iv.R` | 3, Appendix C | Matrix **B**, H-representation and the flagged constraints for the IV model (Figure 3) |
| `02-referee-example.R` | 2.4 | The example after Proposition 1: dropping the external parent Z gives 13/40 and 7/20 instead of 1/3 |
| `03-front-door.R` | 3.3 | Flagged but trivial constraints of the front door model (Figure 4) |
| `04-sequential-iv.R` | 4.1 | Conditional independencies and the two sets of instrumental inequalities (Figure 2) |
| `05-districts-cdegree.R` | 2.2 | Districts and c-degrees of the graph in Figure 1 |
| `06-nested-equivalent.R` | 4.2 | The 8 inequality constraints that distinguish the two graphs of Figure 5 |
| `07-figure7.R` | 5 | Constraints after merging latents in Figure 7, and the example distribution that violates 13 inequality constraints |
| `08-triangle.R` | 5 | The merged triangle graph (Figure 8) has no nontrivial constraints |
| `09-tripartite-bell.R` | 4.3 | Constraint counts for the tripartite Bell graph (Figure 6); takes up to about 30 minutes |
| `11-kedagni-mourifie.R` | 8 | Exact comparison with the inequalities of Kedagni and Mourifie (2020) for binary treatment and outcome |
| `12-pearl-bonet-insufficient.R` | 8 | A distribution satisfying Pearl's inequalities that is incompatible with the IV model (ternary instrument) |
| `13-inequality-criterion.R` | 8 | Proposition 8 and the conjecture after it: for 90 small districts, inequality constraints exist exactly when the district is not saturated after adding HLP edges (a few minutes) |
| `10-benchmarks.R` | 4.4, Table 1 | Matrix sizes, constraint counts and times; with the redundancy check it takes several hours (set `check_redundancy <- FALSE` to skip it) |

`common.R` holds shared setup. The closed-form screening of Proposition 4 is
`prob_simplex_screen()` in `meraconstraints`. `table1-benchmarks.csv` contains the results
reported in Table 1 (all cases except the tripartite Bell graph).

Variables in the scripts have the same names as in the paper's figures (for
example V1, ..., V6 and U1, U2, U3 in Figure 7), so the printed constraints can be
compared directly with the paper. 
