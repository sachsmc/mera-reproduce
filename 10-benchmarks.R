## Section 6, Table 1: size, constraint counts and computation times.
## Usage: Rscript 10-benchmarks.R [--redundancy] [names of cases to skip]
## e.g.   Rscript 10-benchmarks.R "Tripartite Bell"
##
## Output, written after every case so that a long run can be interrupted:
##   table1-benchmarks.csv   the results for all districts of all cases
##   table1-rows.tex         the rows of Table 1, ready to paste into the paper
##   table1-machine.txt      the computer and software used
## If cases are skipped, the output goes to files ending in -partial instead, so that
## the full results are not overwritten.
##
## The counts do not depend on the computer, and the script stops if they differ from the
## counts in the paper. The times do.
##
## --redundancy also runs rcdd::redundant, as reported in the response letter. It is
## slow, about 4 hours for Figure 7 and 7.5 hours for the ternary IV model, and is not run
## for the tripartite Bell graph.
source("common.R")

args <- commandArgs(trailingOnly = TRUE)
check_redundancy <- "--redundancy" %in% args
skip <- setdiff(args, "--redundancy")
stem <- if (length(skip)) "-partial" else ""

graphs <- list(
  "IV, |Z|=|X|=|Y|=2"      = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y)),
  "IV, |Z|=3, |X|=|Y|=2"   = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = 3)),
  "IV, |Z|=4, |X|=|Y|=2"   = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = 4)),
  "IV, |Z|=|X|=|Y|=3"      = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = 3, X = 3, Y = 3)),
  "Fig 5 right"            = dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4)),
  "X -> A <- U -> B <- Y"  = dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B)),
  "Fig 7 merged"           = dhvm_merge_latents(dhvm(caugi(V1 %-->% V2 + V6, V2 %-->% V3, V3 %-->% V4, V4 %-->% V5, V6 %-->% V3,
                                                           U1 %-->% V2 + V4, U2 %-->% V1 + V3 + V6, U3 %-->% V4 + V5))),
  "Tripartite Bell"        = dhvm(caugi(X %-->% V1, Y %-->% V2, Z %-->% V3, U %-->% V1 + V2 + V3))
)
stopifnot(all(skip %in% names(graphs)))
graphs <- graphs[setdiff(names(graphs), skip)]

## Table 1: label in the paper, and the counts for the largest district of each case:
## rows of B, columns of B, distinct columns, inequalities, equalities, nontrivial
## inequalities, nontrivial equalities
paper <- list(
  "IV, |Z|=|X|=|Y|=2"     = list("IV, $|Z|=|X|=|Y|=2$", c(8, 16, 12, 12, 2, 4, 0)),
  "IV, |Z|=3, |X|=|Y|=2"  = list("IV, $|Z|=3$, $|X|=|Y|=2$", c(12, 32, 28, 48, 3, 36, 0)),
  "IV, |Z|=4, |X|=|Y|=2"  = list("IV, $|Z|=4$, $|X|=|Y|=2$", c(16, 64, 60, 160, 4, 144, 0)),
  "IV, |Z|=|X|=|Y|=3"     = list("IV, $|Z|=|X|=|Y|=3$", c(27, 729, 333, 1395, 3, 1368, 0)),
  "Fig 5 right"           = list("Fig.~\\ref{fig:nested-ineq} (right), $\\{V_2,V_4\\}$", c(16, 64, 40, 24, 6, 8, 0)),
  "X -> A <- U -> B <- Y" = list("$X \\rightarrow A \\leftarrow U \\rightarrow B \\leftarrow Y$", c(16, 16, 16, 24, 8, 8, 0)),
  "Fig 7 merged"          = list("Fig.~\\ref{fig:incomp}, merged $\\{V_2,V_4,V_5\\}$", c(32, 64, 48, 1128, 12, 1096, 6)),
  "Tripartite Bell"       = list("Tripartite Bell, Fig.~\\ref{fig:tripart}", c(64, 64, 64, 53856, 38, 53792, 0))
)

## A best effort description of the computer
cpu_ram <- function() {
  out <- function(cmd, args) tryCatch(system2(cmd, args, stdout = TRUE, stderr = FALSE), error = function(e) character(), warning = function(w) character())
  sys <- Sys.info()[["sysname"]]
  if (sys == "Linux") {
    c(grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1], grep("MemTotal", readLines("/proc/meminfo"), value = TRUE))
  } else if (sys == "Darwin") {
    c(out("sysctl", c("-n", "machdep.cpu.brand_string")), out("sysctl", c("-n", "hw.memsize")))
  } else {
    out("powershell", c("-NoProfile", "-Command",
      "(Get-CimInstance Win32_Processor).Name; (Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"))
  }
}
loadNamespace("rcdd")                        # so that sessionInfo() records its version
writeLines(c(format(Sys.time(), "%Y-%m-%d %H:%M %Z"), unname(Sys.info()[c("sysname", "release", "machine")]),
             cpu_ram(), paste("Cores:", parallel::detectCores()), "", capture.output(sessionInfo())),
           paste0("table1-machine", stem, ".txt"))

## So that the first case does not pay for compiling the code
invisible(local({
  g <- graphs[[1]]; r <- dhvm_linear_relations(g); nontrivial_screen(g, r, dhvm_relations_to_polytope(r))
}))

format_time <- function(s) {
  if (s < 0.01) "$<0.01$" else if (s < 1) sprintf("%.2f", s) else if (s < 10) sprintf("%.1f", s) else sprintf("%.0f", s)
}
big <- function(n) if (n >= 10000) format(n, big.mark = ",", trim = TRUE) else as.character(n)

res <- list(); tex <- character()
for (nm in names(graphs)) {
  cat("Running:", nm, "\n"); flush.console()
  g   <- graphs[[nm]]
  rel <- dhvm_linear_relations(g)
  t_vh <- unname(system.time(tope <- dhvm_relations_to_polytope(rel))["elapsed"])
  t_sc <- unname(system.time(flag <- nontrivial_screen(g, rel, tope))["elapsed"])
  run_red <- check_redundancy && nm != "Tripartite Bell"
  t_ir <- NA_real_
  if (run_red) t_ir <- unname(system.time(Hirr <- lapply(tope$HrepQ, \(h) rcdd::redundant(h, representation = "H")$output))["elapsed"])
  cat("  V-to-H:", t_vh, "s  nontrivial screen:", t_sc, "s\n"); flush.console()
  case <- list()
  for (j in seq_along(tope$Hrep)) {
    H <- tope$Hrep[[j]]; V <- tope$Vrep[[j]]
    case[[j]] <- data.frame(
      graph = nm, district = names(tope$Hrep)[j], members = length(rel$district_parms$poterms[[j]]$out),
      rows_B = ncol(V) - 2, cols_B_all = length(unique(unlist(rel$relation[[j]]))), cols_B_distinct = nrow(V),
      n_ineq = sum(H[, 1] == 0), n_eq = sum(H[, 1] == 1),
      n_nontrivial_ineq = sum(flag[[j]] & H[, 1] == 0), n_nontrivial_eq = sum(flag[[j]] & H[, 1] == 1),
      n_redundant_removed = if (run_red) nrow(H) - nrow(Hirr[[j]]) else NA,
      sec_VtoH = t_vh, sec_screen = t_sc, sec_irredundant = t_ir)
  }
  case <- do.call(rbind, case)
  res[[nm]] <- case
  r <- case[which.max(case$n_ineq), ]
  counts <- unlist(r[c("rows_B", "cols_B_all", "cols_B_distinct", "n_ineq", "n_eq", "n_nontrivial_ineq", "n_nontrivial_eq")])
  if (!all(counts == paper[[nm]][[2]]))
    stop("The counts for ", nm, " differ from Table 1: ", paste(counts, collapse = " "))
  tex <- c(tex, paste0(paste(c(paper[[nm]][[1]], r$rows_B, r$cols_B_all, r$cols_B_distinct, big(r$n_ineq), r$n_eq,
                               big(r$n_nontrivial_ineq), r$n_nontrivial_eq, format_time(t_vh), format_time(t_sc)), collapse = " & "), " \\\\"))
  all_res <- do.call(rbind, res)
  write.csv(all_res, paste0("table1-benchmarks", stem, ".csv"), row.names = FALSE)
  writeLines(tex, paste0("table1-rows", stem, ".tex"))
}
print(all_res, row.names = FALSE)
cat("\nThe counts agree with Table 1.\n\nRows of Table 1:\n")
cat(tex, sep = "\n")
