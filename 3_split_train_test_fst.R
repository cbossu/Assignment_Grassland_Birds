#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Stratified 90/10 train/test splits of BREEDING birds for nested FST-outlier
# selection (so FST outliers are never chosen using the birds they're tested on).
#
# Input: a folder of individual files, one per breeding population, one sample ID
#        per line (e.g. vcftools --weir-fst-pop lists: pops/AZ.txt, pops/CA.txt ...).
#        The population name is taken from the file name.
#
# Output (in out_dir):
#   all_breeding_assignments.csv   sample_id, pop, fold, set (train/test) per replicate
#   rep01/train_<pop>.txt          90% of each pop  -> --weir-fst-pop lists
#   rep01/test.txt                 the held-out 10%  -> assign these with FEEMS
#   rep01/run_fst.sh               PAIRWISE vcftools FST (every pop pair, training
#                                  birds only) -> rep01/pairwise/fst_<A>__<B>.weir.fst
#   Then: 00b_select_top_fst.R picks top variants per pair and unions them.
#
# scheme = "kfold": 10 stratified folds, every bird is in exactly one test set
#                   (lets you assess assignment accuracy for all breeding birds).
# scheme = "random": n_rep independent random 90% draws per population.
# -----------------------------------------------------------------------------

## ---- settings ---------------------------------------------------------------
pop_dir     <- "input/pops"            # folder of per-population ID files
pop_pattern <- "\\.ind.txt$"         # which files to read
out_dir     <- "output/fst_splits"
vcf         <- "input/pairwiseFST/losh.impute3.filtered.vcf.gz"  # full imputed VCF (NOT the outlier subset)
scheme      <- "kfold"           # "kfold" or "random"
k           <- 10                # folds (kfold) -> 90% train
n_rep       <- 2                 # replicates (random)
train_frac  <- 0.9               # (random)
min_train   <- 3                 # warn if a pop has fewer training birds than this
seed        <- 20261021
## ----------------------------------------------------------------------------

set.seed(seed)

## ---- read individual files ---------------------------------------------------
files <- list.files(pop_dir, pattern = pop_pattern, full.names = TRUE)
if (length(files) == 0) stop("no files matching ", pop_pattern, " in ", pop_dir)

ind <- do.call(rbind, lapply(files, function(f) {
  ids <- trimws(readLines(f, warn = FALSE))
  ids <- ids[ids != "" & !startsWith(ids, "#")]
  data.frame(sample_id = ids,
             pop = sub(pop_pattern, "", basename(f)),
             stringsAsFactors = FALSE)
}))

dups <- ind$sample_id[duplicated(ind$sample_id)]
if (length(dups)) stop("samples in more than one pop file: ", paste(head(dups), collapse = ", "))

cat("Read", nrow(ind), "birds from", length(files), "populations\n")
print(table(ind$pop))

## ---- assign folds within each population --------------------------------------
# kfold: shuffle within pop, then deal 1..k like cards, so each fold holds ~10%
# of every pop. Pops with < k birds simply have no test bird in some folds.
if (scheme == "kfold") {
  ind$fold <- NA_integer_
  for (p in unique(ind$pop)) {
    i <- which(ind$pop == p)
    ind$fold[i] <- sample(rep_len(sample(k), length(i)))
  }
  reps <- seq_len(k)
  is_test <- function(r) ind$fold == r
} else if (scheme == "random") {
  test_mat <- sapply(seq_len(n_rep), function(r) {
    out <- logical(nrow(ind))
    for (p in unique(ind$pop)) {
      i <- which(ind$pop == p)
      n_test <- max(1, round(length(i) * (1 - train_frac)))
      if (length(i) - n_test < 1) n_test <- 0      # pop of 1: keep it in training
      out[sample(i, n_test)] <- TRUE
    }
    out
  })
  reps <- seq_len(n_rep)
  is_test <- function(r) test_mat[, r]
} else stop("scheme must be 'kfold' or 'random'")

## ---- write lists ------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
long <- list()

for (r in reps) {
  rdir <- file.path(out_dir, sprintf("rep%02d", r))
  dir.create(rdir, showWarnings = FALSE)
  test <- is_test(r)
  train <- ind[!test, ]

  ntr <- table(factor(train$pop, levels = unique(ind$pop)))
  if (any(ntr < min_train))
    warning(sprintf("rep %d: pops with < %d training birds: %s", r, min_train,
                    paste(names(ntr)[ntr < min_train], collapse = ", ")))

  pop_files <- character()
  for (p in unique(train$pop)) {
    f <- file.path(rdir, paste0("train_", p, ".txt"))
    writeLines(train$sample_id[train$pop == p], f)
    pop_files <- c(pop_files, f)
  }
  writeLines(ind$sample_id[test], file.path(rdir, "test.txt"))
  writeLines(train$sample_id, file.path(rdir, "train_all.txt"))

  # vcftools: PAIRWISE Weir & Cockerham FST, training birds only.
  # One run per pop pair -> pairwise/fst_<A>__<B>.weir.fst
  pops_r <- sort(unique(train$pop))
  pairs  <- combn(pops_r, 2, simplify = FALSE)
  # absolute path, since run_fst.sh cd's into the replicate folder
  vcf_abs <- shQuote(if (startsWith(vcf, "/")) vcf else file.path(getwd(), vcf))
  fst_cmd <- vapply(pairs, function(pp) paste0(
    "vcftools --gzvcf ", vcf_abs,
    " --weir-fst-pop train_", pp[1], ".txt",
    " --weir-fst-pop train_", pp[2], ".txt",
    " --out pairwise/fst_", pp[1], "__", pp[2],
    " 2> pairwise/fst_", pp[1], "__", pp[2], ".stderr"), "")
  # pairs run in parallel with xargs; set NCPU (e.g. $SLURM_NTASKS_PER_NODE)
  writeLines(c("#!/usr/bin/env bash", "set -euo pipefail",
               paste("cd", shQuote(normalizePath(rdir))),
               "mkdir -p pairwise",
               "NCPU=${NCPU:-1}",
               "cat > pairwise/commands.txt <<'EOF'", fst_cmd, "EOF",
               "xargs -P \"$NCPU\" -I{} bash -c '{}' < pairwise/commands.txt",
               sprintf("echo \"done: %d pairwise FST runs\"", length(pairs))),
             file.path(rdir, "run_fst.sh"))

  long[[length(long) + 1]] <- data.frame(
    replicate = r, sample_id = ind$sample_id, pop = ind$pop,
    set = ifelse(test, "test", "train"))
}

long <- do.call(rbind, long)
write.csv(long, file.path(out_dir, "all_breeding_assignments.csv"), row.names = FALSE)

## ---- sanity check -----------------------------------------------------------
smry <- aggregate(set ~ replicate, data = long,
                  FUN = function(x) round(mean(x == "train"), 3))
names(smry)[2] <- "prop_train"
print(smry)
if (scheme == "kfold") {
  n_test <- table(long$sample_id[long$set == "test"])
  stopifnot(all(n_test == 1), length(n_test) == nrow(ind))
  cat("Every breeding bird is held out exactly once across", k, "folds.\n")
}
cat("Wrote splits to", normalizePath(out_dir), "\n")

