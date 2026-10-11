#!/usr/bin/env Rscript

# GSE183947 processed-expression analysis.
# The GEO supplementary matrix contains FPKM, not raw counts. Do not use DESeq2
# on these values. This script uses limma on log2(FPKM + 1).
suppressPackageStartupMessages({
  library(limma)
})

options(timeout = 300)
outdir <- "cancer-results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
datadir <- "data/GSE183947"
dir.create(datadir, showWarnings = FALSE, recursive = TRUE)

fpkm_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183947/suppl/GSE183947_fpkm.csv.gz"
matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183947/matrix/GSE183947_series_matrix.txt.gz"
fpkm_file <- file.path(datadir, "GSE183947_fpkm.csv.gz")
matrix_file <- file.path(datadir, "GSE183947_series_matrix.txt.gz")

download_if_missing <- function(url, path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    message("Downloading ", basename(path))
    download.file(url, path, mode = "wb", quiet = FALSE)
  }
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Download failed: ", url)
  }
}
download_if_missing(fpkm_url, fpkm_file)
download_if_missing(matrix_url, matrix_file)

# Read GEO sample metadata from the series matrix.
geo_lines <- readLines(gzfile(matrix_file), warn = FALSE)
extract_geo_row <- function(prefix) {
  line <- geo_lines[startsWith(geo_lines, prefix)]
  if (length(line) != 1) stop("Expected exactly one GEO metadata row: ", prefix)
  values <- strsplit(sub(paste0("^", prefix), "", line), "\t", fixed = TRUE)[[1]]
  values <- gsub('^"|"$', "", values)
  values
}
geo_ids <- extract_geo_row("!Sample_geo_accession")
geo_titles <- extract_geo_row("!Sample_title")
if (length(geo_ids) != 60 || length(geo_titles) != 60) {
  stop("Expected 60 GEO samples; found ", length(geo_ids), " IDs and ",
       length(geo_titles), " titles.")
}
metadata <- data.frame(
  geo_accession = geo_ids,
  title = geo_titles,
  stringsAsFactors = FALSE
)
metadata$condition <- ifelse(grepl("tumou?r", metadata$title, ignore.case = TRUE),
                             "tumor",
                             ifelse(grepl("normal", metadata$title, ignore.case = TRUE),
                                    "normal", NA_character_))
metadata$pair <- sub(".*rep\\s*([0-9]+).*", "\\1", metadata$title, ignore.case = TRUE)
if (anyNA(metadata$condition) || any(metadata$pair == metadata$title)) {
  stop("Could not reliably parse tumour/normal condition and pair number from GEO titles.")
}
metadata$pair <- paste0("pair", metadata$pair)
if (anyDuplicated(metadata$pair[metadata$condition == "tumor"]) ||
    anyDuplicated(metadata$pair[metadata$condition == "normal"])) {
  stop("Duplicate pair IDs found within a condition.")
}
if (!identical(sort(metadata$pair[metadata$condition == "tumor"]),
               sort(metadata$pair[metadata$condition == "normal"]))) {
  stop("Tumour and normal sample pair IDs do not match.")
}
if (any(table(metadata$condition) != 30)) {
  stop("Expected 30 tumour and 30 normal samples; observed: ",
       paste(names(table(metadata$condition)), table(metadata$condition), collapse = ", "))
}

# Matrix columns are in GEO sample order for this accession. Assert dimensions
# before attaching metadata; never silently recycle or truncate labels.
fpkm <- read.csv(gzfile(fpkm_file), check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(fpkm) != 61) {
  stop("Expected one gene-ID column plus 60 sample columns; observed ", ncol(fpkm), " columns.")
}
gene_ids <- as.character(fpkm[[1]])
expr <- as.matrix(fpkm[, -1, drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- make.unique(gene_ids)
if (ncol(expr) != nrow(metadata)) stop("Expression matrix and GEO metadata sample counts differ.")
if (any(!is.finite(expr)) || any(expr < 0)) stop("FPKM matrix contains invalid or negative values.")
if (anyDuplicated(rownames(expr))) stop("Gene identifiers are not unique after import.")

# Reorder metadata by the paired design labels; expression columns follow the
# GEO series order. The pairing labels are checked explicitly above.
metadata$condition <- factor(metadata$condition, levels = c("normal", "tumor"))
metadata$pair <- factor(metadata$pair, levels = paste0("pair", 1:30))
rownames(metadata) <- metadata$geo_accession
colnames(expr) <- metadata$geo_accession
write.csv(metadata, file.path(outdir, "sample_metadata.csv"), row.names = FALSE)

# This is processed FPKM data, not raw counts. Paired limma is used for an
# exploratory tumour-vs-normal comparison; results should not be described as
# count-based DESeq2 output.
log_expr <- log2(expr + 1)
design <- model.matrix(~ pair + condition, data = metadata)
fit <- lmFit(log_expr, design)
fit <- eBayes(fit, robust = TRUE)
coef_name <- "conditiontumor"
if (!(coef_name %in% colnames(design))) stop("Tumour-vs-normal coefficient not found.")
res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
res$gene_id <- rownames(res)
res <- res[, c("gene_id", setdiff(colnames(res), "gene_id"))]
write.csv(res, file.path(outdir, "limma_tumor_vs_normal_all_genes.csv"), row.names = FALSE)
write.csv(subset(res, adj.P.Val < 0.05),
          file.path(outdir, "limma_tumor_vs_normal_FDR05.csv"), row.names = FALSE)

# Compact QC plots, using base graphics to keep dependencies minimal.
pdf(file.path(outdir, "expression_pca.pdf"), width = 8, height = 6)
pca <- prcomp(t(log_expr), scale. = FALSE)
plot(pca$x[, 1], pca$x[, 2],
     col = ifelse(metadata$condition == "tumor", "firebrick", "steelblue"),
     pch = 19, xlab = "PC1", ylab = "PC2",
     main = "GSE183947: log2(FPKM + 1) PCA")
legend("topright", legend = c("Normal", "Tumor"),
       col = c("steelblue", "firebrick"), pch = 19, bty = "n")
dev.off()

pdf(file.path(outdir, "limma_volcano.pdf"), width = 8, height = 6)
with(res, plot(logFC, -log10(pmax(P.Value, .Machine$double.xmin)),
               pch = 16, cex = 0.45,
               xlab = "log2 fold-change (tumor vs normal)",
               ylab = "-log10 nominal P-value",
               main = "GSE183947 paired limma"))
abline(v = c(-1, 1), lty = 2, col = "grey50")
abline(h = -log10(0.05), lty = 2, col = "grey50")
dev.off()

summary <- data.frame(
  metric = c("dataset", "samples", "tumor_samples", "normal_samples",
             "matched_pairs", "genes_tested", "FDR_lt_0.05",
             "input_type", "method"),
  value = c("GSE183947", ncol(expr), sum(metadata$condition == "tumor"),
            sum(metadata$condition == "normal"), nlevels(metadata$pair),
            nrow(expr), sum(res$adj.P.Val < 0.05, na.rm = TRUE),
            "FPKM (processed; not raw counts)",
            "limma on log2(FPKM + 1), paired design ~ pair + condition")
)
write.csv(summary, file.path(outdir, "analysis_summary.csv"), row.names = FALSE)
writeLines(c(
  "Dataset: GSE183947 (breast cancer)",
  "Source: NCBI GEO supplementary FPKM matrix and GEO series matrix metadata",
  "Input values are FPKM, not raw gene counts. DESeq2 was not used.",
  "Analysis: limma on log2(FPKM + 1), paired design ~ pair + condition.",
  "This is a processed-expression reanalysis, not a raw-FASTQ pipeline run.",
  "Check original sample pairing and study methods before scientific interpretation."
), file.path(outdir, "README_results.txt"))
message("Analysis completed. Results written to ", outdir)
