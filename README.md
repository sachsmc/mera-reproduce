# Reproducibility scripts

Scripts that reproduce the examples, numbers and tables in
"Deriving Complete Constraints in Hidden Variable Models".

## Requirements

- R (tested with R 4.4.2 and 4.5.3)
- R packages: `caugi`, `rcdd`, `data.table`, `MASS`, and `meraconstraints`
  (https://github.com/sachsmc/meraconstraints), in a version that has `nontrivial_screen()`.
  `common.R` stops with a message if the installed version does not.
- `caugi` needs a Rust toolchain (`cargo`, `rustc` >= 1.80) when installed from source, and
  `rcdd` needs the GMP library.

```r
install.packages(c("caugi", "rcdd", "data.table", "MASS", "remotes"))
remotes::install_github("sachsmc/meraconstraints")
```

## Usage

Run the scripts from this directory, either one at a time,

    Rscript 01-iv.R

or all together, saving each script's output to `output/`:

    Rscript run_all.R          # everything except the two slow scripts
    Rscript run_all.R slow     # everything

`run_all.R` ends with an error if a script fails.

## Rerunning the benchmarks (Table 1)

    Rscript 10-benchmarks.R > output/10-benchmarks.txt 2>&1

This takes about half an hour, nearly all of it for the tripartite Bell graph. To time only
the other cases, run `Rscript 10-benchmarks.R "Tripartite Bell"`, which writes to files ending
in `-partial` and leaves the full results alone. The script writes

| File | Contents |
|---|---|
| `table1-benchmarks.csv` | Results for all districts, saved after every case |
| `table1-rows.tex` | The rows of Table 1, ready to paste into the paper |
| `table1-machine.txt` | The computer and package versions used |

The counts do not depend on the computer, and the script stops if they differ from Table 1.
Only the times change. After a rerun on a different computer, also update the description of
the computer in the caption of Table 1 and the times in Section 4.2 and in the response letter.

`Rscript 10-benchmarks.R --redundancy` also runs `rcdd::redundant`, as reported in the response
letter (Reviewer 2, comment 2). This takes several hours. It does not run `rcdd::redundant` for
the tripartite Bell graph, but it still does that graph's V-to-H conversion, about 20 minutes,
and it overwrites the full `table1-*` files. `redundancy-check.csv` holds the results of an earlier
run, in the column `n_redundant_removed`.

## Contents

| Script | Paper section | What it reproduces |
|---|---|---|
| `01-iv.R` | 3, Appendix B | Matrix **B**, H-representation and the nontrivial constraints of the IV model (Figure 4) |
| `02-referee-example.R` | 2.5 (supporting check) | The example in Reviewer 1's first comment, for Proposition 3: dropping the external parent Z gives 13/40 and 7/20 instead of 1/3 |
| `03-front-door.R` | 3.3 | The front door model has no nontrivial constraints |
| `04-sequential-iv.R` | 3.3 | Conditional independencies and the two sets of instrumental inequalities (Figure 3) |
| `05-districts-cdegree.R` | 2.2 | Districts and c-degrees of the graph in Figure 1 |
| `06-nested-equivalent.R` | 4.1 | The 8 inequality constraints that distinguish the two graphs of Figure 5, their form in the paper, and an example distribution (about a minute) |
| `07-merged-latents.R` | 5 | Constraints after merging latents in Figure 7, the Verma-type constraint, the example distribution that violates 13 inequality constraints, and the classes of the 1096 nontrivial inequalities |
| `08-triangle.R` | 5 | The merged triangle scenario has no nontrivial constraints |
| `09-tripartite-bell.R` | 4.2 | Constraint counts and the displayed inequality for the tripartite Bell graph (Figure 6); takes up to about 30 minutes |
| `10-benchmarks.R` | 6, Table 1 | Matrix sizes, constraint counts and times; see above |
| `11-kedagni-mourifie.R` | 7 (supporting check) | Exact comparison with the inequalities of Kedagni and Mourifie (2020) for binary treatment and outcome |
| `12-pearl-bonet-insufficient.R` | 1 (supporting check) | A distribution satisfying Pearl's inequalities that is incompatible with the IV model (ternary instrument) |
| `13-inequality-criterion.R` | 4.3 | Remark 3, and a supporting check of its converse, which the paper does not claim: for 90 small districts, inequality constraints exist exactly when the district contains one of the two structures. Also classifies the examples of Section 4.3 |
| `14-nontrivial-criterion.R` | 3.3 | Proposition 6: counts, the relation to the nested Markov model, the example row, invariance to how rows are written, and a linear programming check |
| `15-zero-probabilities.R` | 3.2, Remark 2 | The rule for zero probabilities, including the Verma graph with C = B |
| `16-transformations.R` | 2.3 | Propositions 1 and 2: the example of Figure 2 and the district of Figure 7 that neither can reduce |

`common.R` holds shared setup. The criterion of Proposition 6 is `nontrivial_screen()` in
`meraconstraints`, and `meraconstraints()` stores its result in the element `nontrivial`.
`inequality-criterion-results.csv` holds the results of `13-inequality-criterion.R`.

Variables in the scripts mostly have the same names as in the paper's figures (for
example V1, ..., V6 and U1, U2, U3 in Figure 7), so the printed constraints can be
compared directly with the paper. The exceptions are the triangle, which uses V1, V2, V3 where
the paper has A, B, C; the front door, which uses U1 where the paper has U; and the merged
unobserved variable of Figure 7, which prints as U1 where the paper has U13.
