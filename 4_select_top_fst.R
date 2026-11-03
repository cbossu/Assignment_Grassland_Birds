#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Pick top FST variants from PAIRWISE vcftools output and union them across pairs.
#
# Usage:  Rscript 4_select_top_fst.R output/fst_splits/rep01
#         (run once per replicate, or once on the all-breeding-birds run)
#
# Reads   <rep_dir>/pairwise/fst_<A>__<B>.weir.fst
# Writes  <rep_dir>/top_fst_positions.txt   CHROM<TAB>POS  -> vcftools --positions /
#                                                               bcftools view -T
#         <rep_dir>/top_fst_by_pair.csv     which pair(s) selected each SNP
#         <rep_dir>/top_fst_summary.csv     per-pair threshold and counts
# -----------------------------------------------------------------------------

## ---- settings ---------------------------------------------------------------
select_by   <- "n"      # "n" = top_n SNPs per pair; "quantile" = above top_q per pair
top_n       <- 2000
top_q       <- 0.995    # 99.5th percentile
min_fst     <- 0        # also require FST >= this (0 = off)
vcf         <- "input/pairwiseFST/losh.impute3.filtered.vcf.gz"  # full imputed VCF (NOT the outlier subset)
## ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
rep_dir <- if (length(args)) args[1] else "output/fst_splits/rep01"
files <- list.files(file.path(rep_dir, "pairwise"), pattern = "\\.weir\\.fst$",
                    full.names = TRUE)
if (!length(files)) stop("no .weir.fst files in ", file.path(rep_dir, "pairwise"))

read_fst <- function(f) {
  x <- read.table(f, header = TRUE, sep = "\t", na.strings = c("-nan", "nan", "NA"),
                  colClasses = c("character", "integer", "numeric"))
  names(x) <- c("CHROM", "POS", "FST")
  x <- x[is.finite(x$FST), ]   # drop sites invariant within the pair
  x$pair <- sub("^fst_(.*)\\.weir\\.fst$", "\\1", basename(f))
  x
}

sel <- list(); smry <- list()
for (f in files) {
  x <- read_fst(f)
  x <- x[order(-x$FST), ]
  thr <- if (select_by == "n") {
    if (nrow(x) >= top_n) x$FST[top_n] else min(x$FST)
  } else quantile(x$FST, top_q, names = FALSE)
  thr <- max(thr, min_fst)
  # ties at the threshold are all kept
  top <- x[x$FST >= thr, ]
  sel[[f]] <- top
  smry[[f]] <- data.frame(pair = x$pair[1], n_sites = nrow(x),
                          fst_threshold = signif(thr, 4), n_selected = nrow(top),
                          mean_fst_all = signif(mean(x$FST), 4))
}

top_all <- do.call(rbind, sel)
smry    <- do.call(rbind, smry)

# union across pairs, recording which pairs flagged each SNP
by_snp <- aggregate(cbind(pairs = pair) ~ CHROM + POS, data = top_all,
                    FUN = function(p) paste(sort(unique(p)), collapse = ";"))
by_snp$n_pairs <- lengths(strsplit(by_snp$pairs, ";"))
max_fst <- aggregate(FST ~ CHROM + POS, data = top_all, FUN = max)
names(max_fst)[3] <- "max_pairwise_fst"
by_snp <- merge(by_snp, max_fst, by = c("CHROM", "POS"))
by_snp <- by_snp[order(by_snp$CHROM, by_snp$POS), ]

write.table(by_snp[, c("CHROM", "POS")], file.path(rep_dir, "top_fst_positions.train.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(by_snp, file.path(rep_dir, "top_fst_by_pair.csv"), row.names = FALSE)
write.csv(smry,   file.path(rep_dir, "top_fst_summary.csv"), row.names = FALSE)

print(smry, row.names = FALSE)
cat(sprintf("\n%d pairs; %d SNPs selected (union); %d flagged by >1 pair\n",
            nrow(smry), nrow(by_snp), sum(by_snp$n_pairs > 1)))
cat("Subset with:\n  bcftools view -T", file.path(rep_dir, "top_fst_positions.train.txt"),
    " $vcf -Oz -o", file.path(rep_dir, "LOSH_topfst.vcf.gz"), "\n")

