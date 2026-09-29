#!/usr/bin/env Rscript
#BiocManager::install("DESeq2")
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage:\n")
  cat("  Rscript deseqDEG.R <featurecounts.txt> <GenoA_vs_GenoB> <CondA_vs_CondB> <Out_name>\n\n")
  cat("Examples:\n")
  cat("  Rscript deseqDEG.R Gyula_fcounts.txt WTvsWT TuMVvsMock\n")
  cat("  Rscript deseqDEG.R Gyula_fcounts.txt WTvsWT TVCVvsMock\n")
  cat("  Rscript deseqDEG.R Liu_CMV_fcounts.txt 2bvsWT CMVvsMock\n")
  quit(status = 1)
}

 #counts_file <- "FeatureCounts/Clavel_TCV_TyMV_fcounts.txt"
 #geno_comp   <- "WTvsWT"
 #cond_comp   <- "TCVvsMock"
 #
counts_file <- args[1]
geno_comp   <- args[2]
cond_comp   <- args[3]

project <- gsub(".*/|_.*","",counts_file)
out_file <- paste("FeatureCounts/",project,"_",cond_comp, "_",geno_comp, sep = "") 
clean_colname <- function(x) {
  x <- basename(x)
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x
}

parse_vs <- function(x, label) {
  parts <- strsplit(x, "vs", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop(sprintf("Bad %s comparison '%s' (expected like A vs B: e.g. WTvsatg2).", label, x))
  parts
}

make_tpm <- function(counts_mat, gene_length_bp) {
  len_kb <- gene_length_bp / 1000
  len_kb[len_kb <= 0] <- NA_real_
  rpk <- sweep(counts_mat, 1, len_kb, "/")
  scale <- colSums(rpk, na.rm = TRUE) / 1e6
  tpm <- sweep(rpk, 2, scale, "/")
  tpm[is.na(tpm)] <- 0
  tpm
}

# Parse sample name:
#   Prefix(-GENO)?_COND(-GENO)?_RepN
# Genotype is only detected if it matches the user-requested genotype labels (non-WT).
parse_sample <- function(sample, ...) {
  sample <- sub("\\.bam$", "", sample, ignore.case = TRUE)
  
  # replicate: _Rep1 / -Rep1 / _R1 / -R1
  rep <- NA_integer_
  core <- sample
  mrep <- regexec("([_-])(Rep|R)([0-9]+)$", core, ignore.case = TRUE)
  mr <- regmatches(core, mrep)[[1]]
  if (length(mr) == 4) {
    rep <- as.integer(mr[4])
    core <- sub("([_-])(Rep|R)[0-9]+$", "", core, ignore.case = TRUE)
  }
  
  toks <- strsplit(core, "_", fixed = TRUE)[[1]]
  if (length(toks) < 2) return(list(ok=FALSE, sample=sample))
  
  # assume last token is condition, everything before is prefix
  cond <- toks[length(toks)]
  prefix <- paste(toks[-length(toks)], collapse = "_")
  
  genotype <- "WT"
  strain <- prefix
  
  # 1) genotype encoded in PREFIX like Clavel-atg2_...
  # strain = part before first '-', genotype = the rest
  if (grepl("-", prefix, fixed = TRUE)) {
    parts <- strsplit(prefix, "-", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      strain <- parts[1]
      genotype <- paste(parts[-1], collapse = "-")
    }
  }
  
  # 2) genotype encoded in CONDITION like CMV-2b_...
  # but DO NOT treat Mock-TV / Mock-Tu as genotype
  if (genotype == "WT" && grepl("-", cond, fixed = TRUE)) {
    parts <- strsplit(cond, "-", fixed = TRUE)[[1]]
    if (length(parts) >= 2 && parts[1] != "Mock") {
      cond <- parts[1]
      genotype <- paste(parts[-1], collapse = "-")
    }
  }
  
  list(ok=TRUE, sample=sample, prefix=strain, genotype=genotype, condition=cond, replicate=rep)
}

# Auto-resolve virus-specific mock labels:
# If user asks X vs Mock, but "Mock" doesn't exist and "Mock-*" exists,
# try Mock-<first2letters(X)> (TuMV -> Mock-Tu, TVCV -> Mock-TV).
resolve_mock <- function(condA, condB, observed_conditions) {
  if (condB == "Mock" && !("Mock" %in% observed_conditions)) {
    mock_variants <- observed_conditions[grepl("^Mock-", observed_conditions)]
    if (length(mock_variants) > 0) {
      tag <- substr(condA, 1, 2)
      candidate <- paste0("Mock-", tag)
      if (candidate %in% mock_variants) {
        message(sprintf("Auto-resolved baseline: '%s' -> '%s' (because '%s' not present).", condB, candidate, condB))
        condB <- candidate
      } else if (length(mock_variants) == 1) {
        message(sprintf("Auto-resolved baseline: '%s' -> '%s' (only mock variant found).", condB, mock_variants[1]))
        condB <- mock_variants[1]
      } else {
        warning("Multiple Mock-* variants exist and none matched by tag; using ALL Mock-* variants.")
        # return special marker handled downstream
        condB <- "Mock-*"
      }
    }
  }
  list(condA = condA, condB = condB)
}

# ---------- read featureCounts ----------
dt <- fread(counts_file, header = TRUE, sep = "\t", data.table = FALSE,  check.names = FALSE)
required <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
miss <- setdiff(required, colnames(dt))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

len_idx <- match("Length", colnames(dt))
count_cols <- colnames(dt)[(len_idx + 1):ncol(dt)]
if (length(count_cols) < 2) stop("Found <2 count columns; check featureCounts output.")

cleaned <- make.unique(vapply(count_cols, clean_colname, character(1)))
colnames(dt)[(len_idx + 1):ncol(dt)] <- cleaned

counts_df <- dt[, cleaned, drop = FALSE]
rownames(counts_df) <- dt$Geneid

# fractional counts handling
nonint <- any(abs(as.matrix(counts_df) - round(as.matrix(counts_df))) > 1e-6, na.rm = TRUE)
if (nonint) {
  warning("Non-integer counts detected (likely featureCounts --fraction). Rounding for DESeq2.\n",
          "Best practice: rerun featureCounts without --fraction for DESeq2.")
  counts_df <- round(counts_df)
}
counts_mat <- as.matrix(counts_df)
storage.mode(counts_mat) <- "integer"

gene_len <- as.numeric(dt$Length)
names(gene_len) <- dt$Geneid

# ---------- parse comparisons ----------
geno_parts <- parse_vs(geno_comp, "genotype")
cond_parts <- parse_vs(cond_comp, "condition")
genoA <- geno_parts[1]; genoB <- geno_parts[2]
condA <- cond_parts[1]; condB <- cond_parts[2]
geno_levels <- unique(c(genoA, genoB))

# ---------- parse samples ----------
sample_info <- lapply(colnames(counts_mat), parse_sample)
ok <- vapply(sample_info, function(x) isTRUE(x$ok), logical(1))
if (!all(ok)) {
  bad <- vapply(sample_info[!ok], function(x) x$sample, character(1))
  stop("Could not parse these sample names. Expected something like Prefix_COND_RepN:\n  ",
       paste(bad, collapse = "\n  "))
}

coldata <- data.frame(
  sample    = vapply(sample_info, `[[`, "", "sample"),
  prefix    = vapply(sample_info, `[[`, "", "prefix"),
  genotype  = vapply(sample_info, `[[`, "", "genotype"),
  condition = vapply(sample_info, `[[`, "", "condition"),
  replicate = vapply(sample_info, `[[`, NA_integer_, "replicate"),
  stringsAsFactors = FALSE
)
rownames(coldata) <- coldata$sample

# Auto-resolve Mock if needed
resolved <- resolve_mock(condA, condB, unique(coldata$condition))
condA <- resolved$condA
condB <- resolved$condB

# Determine condition set to keep
if (condB == "Mock-*") {
  cond_keep <- unique(c(condA, unique(coldata$condition[grepl("^Mock-", coldata$condition)])))
} else {
  cond_keep <- unique(c(condA, condB))
}

# subset to requested geno/cond
keep <- (coldata$genotype %in% geno_levels) & (coldata$condition %in% cond_keep)
coldata_sub <- coldata[keep, , drop = FALSE]
counts_sub  <- counts_mat[, rownames(coldata_sub), drop = FALSE]
if (ncol(counts_sub) < 2) stop("After subsetting, <2 samples remain. Check your 'vs' arguments vs column names.")

# set reference levels = right-hand side of 'vs'
# (Only keep the two focal condition labels for DE results; if condB==Mock-* we keep those samples but compare condA vs each mock-variant isn't what DESeq2 contrast means)
# In that case, we collapse all Mock-* into "Mock" for testing.
if (condB == "Mock-*") {
  coldata_sub$condition <- ifelse(grepl("^Mock-", coldata_sub$condition), "Mock", coldata_sub$condition)
  condB_test <- "Mock"
} else {
  condB_test <- condB
}

coldata_sub$genotype  <- factor(coldata_sub$genotype,  levels = unique(c(genoB, genoA)))
coldata_sub$condition <- factor(coldata_sub$condition, levels = unique(c(condB_test, condA)))

two_genotypes <- length(levels(droplevels(coldata_sub$genotype))) == 2
design_formula <- if (two_genotypes) ~ genotype + condition + genotype:condition else ~ condition

dds <- DESeqDataSetFromMatrix(countData = counts_sub[, rownames(coldata_sub), drop=FALSE],
                             colData   = coldata_sub,
                             design    = design_formula)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)

# outputs
out_prefix <- sub("\\.(txt|tsv|tab|csv)$", "", basename(counts_file))
tag <- paste0(geno_comp, "__", cond_comp)

# condition effect in reference genotype (genoB)
res_ref <- results(dds, contrast = c("condition", condA, condB_test))

# interaction (if two genotypes)
res_int <- NULL
res_genoA <- NULL
if (two_genotypes) {
  rn <- resultsNames(dds)
  cond_coef <- rn[grepl("^condition_", rn) & grepl(condA, rn)]
  if (length(cond_coef) == 0) cond_coef <- rn[grepl("^condition_", rn)][1]

  int_coef <- rn[grepl("genotype", rn) & grepl("condition", rn) & grepl(condA, rn)]
  if (length(int_coef) == 0) int_coef <- rn[grepl("genotype", rn) & grepl("condition", rn)][1]

  if (!is.na(int_coef) && length(int_coef) > 0) {
    res_int <- results(dds, name = int_coef[1])
    # condition effect in genotype A = main condition + interaction
    res_genoA <- results(dds, contrast = list(c(cond_coef[1], int_coef[1])))
  }
}

# TPM for reporting
tpm <- make_tpm(counts(dds), gene_len[rownames(dds)])
tpm_df <- data.frame(Geneid = rownames(tpm), tpm, check.names = FALSE)

# write
write.csv(coldata_sub, file = paste0(out_file, ".samples.csv"), row.names = FALSE)
write.table(data.frame(Geneid = rownames(counts(dds)), counts(dds), check.names = FALSE),
            file = paste0(out_file, ".counts_used.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_df,
            file = paste0(out_file, ".TPM.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

res_to_df <- function(res) {
  df <- as.data.frame(res)
  df$Geneid <- rownames(df)
  df[, c("Geneid", setdiff(colnames(df), "Geneid")), drop = FALSE]
}

write.table(res_to_df(res_ref),
            file = paste0(out_file, ".DESeq2_", condA, "_vs_", condB_test, "_in_", genoB, ".tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(res_genoA)) {
  write.table(res_to_df(res_genoA),
              file = paste0(out_file,  ".DESeq2_", condA, "_vs_", condB_test, "_in_", genoA, ".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
if (!is.null(res_int)) {
  write.table(res_to_df(res_int),
              file = paste0(oout_file, ".DESeq2_interaction.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("Done.\n")
cat("Design: ", deparse(design_formula), "\n")
cat("Samples used: ", nrow(coldata_sub), "\n")
cat("Genotypes: ", paste(levels(droplevels(coldata_sub$genotype)), collapse = ", "), "\n")
cat("Conditions: ", paste(levels(droplevels(coldata_sub$condition)), collapse = ", "), "\n")
