## Shared setup for the reproducibility scripts.
## Requires the R packages meraconstraints (https://github.com/sachsmc/meraconstraints),
## caugi, rcdd and data.table.
suppressPackageStartupMessages({
  library(meraconstraints)
  library(data.table)
})

## Text of the constraints of district j, optionally only those flagged as potentially
## nontrivial, in terms of observed (star = FALSE) or interventional (star = TRUE)
## probabilities.
constraint_text <- function(cn, j = 1, flagged_only = TRUE, star = FALSE) {
  H    <- cn$polytope$Hrep[[j]]
  prob <- if (star) cn$functional_relations$probstarform[[j]] else cn$functional_relations$probform[[j]]
  txt  <- hrep_to_text(H, prob)
  if (flagged_only) txt[cn$trivial_if_false[[j]]] else txt
}

section <- function(title) cat("\n==== ", title, " ====\n", sep = "")
