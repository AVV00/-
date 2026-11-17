#!/usr/bin/env Rscript

# Figure 6C: average association of oxygen with malignant program scores.
# Run: Rscript Fig6C_GO_continuous_oxygen_mean_slope.R
# Dependencies: data.table, Matrix, ggplot2.

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(ggplot2)
})

root <- "/mnt/data1/linyanlin/scATAC"
gene_set_file <- file.path(root, "all_data/step5_pathway", "0407恶性表型基因集.txt")
cell_file <- file.path(root, "picture/check/261105", "All_malignantcell_mt_gene_expression.csv")
activity_file <- file.path(root, "all_data/step4_geneannotation",
                           "all_cells_gene_activity_gene_by_cell_sparse.rds")
pathway_file <- file.path(root, "all_data/step5_pathway", "pathway.csv")
output_stem <- file.path(root, "picture/pic6/261106", "Fig6C_GO_continuous_oxygen_mean_slope")

# Use the five Figure 6B programs and their deduplicated union.
programs <- c("Proliferation", "Stress or hypoxia", "Stemness", "Invasion", "Immune Escape")
gene_table <- fread(gene_set_file, sep = "\t", header = TRUE, fill = TRUE,
                    colClasses = "character", encoding = "UTF-8")
stopifnot(all(programs %in% names(gene_table)))
signatures <- lapply(as.list(gene_table[, ..programs]), function(symbols) {
    symbols <- toupper(trimws(symbols))
    symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
    symbols[symbols == "CD133"] <- "PROM1"
    symbols[symbols == "CD90"] <- "THY1"
    unique(symbols)
})
signatures[["Malignant program score"]] <- unique(unlist(signatures, use.names = FALSE))

# Decode the project sample and patient IDs from the canonical cell names.
cells <- fread(cell_file, select = c("cell", "candidate_mutation_status", "oxygen"))
cells[, cell := trimws(cell)]
cells[, candidate_mutation_status := trimws(candidate_mutation_status)]
cells[, oxygen := suppressWarnings(as.numeric(oxygen))]
cells <- cells[candidate_mutation_status %chin% c("G+", "WITHOUT") & is.finite(oxygen)]
stopifnot(nrow(cells) > 0L, !anyNA(cells$cell), !anyDuplicated(cells$cell))
cells[, G := as.integer(candidate_mutation_status == "G+")]
cells[, sample_id := sub("^.*_(HCC[0-9]+_.*)$", "\\1", cell)]
patient_map <- setNames(paste0("HCC", 1:10), paste0("HCC", c(3, 6:14)))
cells[, patient := unname(patient_map[sub("_.*$", "", sample_id)])]
stopifnot(!anyNA(cells$patient))

# Fit only within the shared observed oxygen range of G+ and WITHOUT.
min_cells <- 10L
coverage <- cells[, {
    positive <- oxygen[G == 1L]
    reference <- oxygen[G == 0L]
    paired <- length(positive) > 0L && length(reference) > 0L
    lo <- if (paired) max(min(positive), min(reference)) else NA_real_
    hi <- if (paired) min(max(positive), max(reference)) else NA_real_
    overlap <- paired && lo < hi
    positive_used <- if (overlap) positive[positive >= lo & positive <= hi] else numeric()
    reference_used <- if (overlap) reference[reference >= lo & reference <= hi] else numeric()
    .(
        Gplus_total = length(positive),
        WITHOUT_total = length(reference),
        lo = lo,
        hi = hi,
        Gplus_overlap = length(positive_used),
        WITHOUT_overlap = length(reference_used),
        eligible = min(length(positive_used), length(reference_used)) >= min_cells &&
            min(uniqueN(positive_used), uniqueN(reference_used)) >= 2L
    )
}, by = .(patient, sample_id)]
print(coverage, nrows = Inf)
keep <- coverage[eligible == TRUE, .(patient, sample_id, lo, hi)]
if (uniqueN(keep$patient) < 2L) {
    stop("Fewer than two patients with estimable within-sample slopes")
}
cells <- merge(cells, keep, by = c("patient", "sample_id"), sort = FALSE)
cells <- cells[oxygen >= lo & oxygen <= hi]
cat("Included:", nrow(cells), "cells,", nrow(keep), "samples,",
    uniqueN(keep$patient), "patients\n")
cells[, O_step := (oxygen - mean(oxygen)) / 0.01, by = sample_id]

# Match the gene activity matrix and exclude genes used to compute oxygen.
activity <- readRDS(activity_file)
gene_ids <- toupper(trimws(rownames(activity)))
cell_ids <- colnames(activity)
cell_cols <- match(cells$cell, cell_ids)
stopifnot(length(gene_ids) > 0L, !is.null(cell_ids), !anyDuplicated(cell_ids),
          !anyNA(cell_cols))

pathways <- fread(pathway_file)
ox_record <- pathways[gmt_geneset_name == "KEGG_OXIDATIVE_PHOSPHORYLATION", genes]
stopifnot(length(ox_record) == 1L, !is.na(ox_record))
ox_genes <- unique(toupper(trimws(strsplit(ox_record, ";", fixed = TRUE)[[1L]])))
ox_genes <- ox_genes[nzchar(ox_genes) & !grepl("^MT-", ox_genes)]
ox_aliases <- c(
    ATP5F1A = "ATP5A1", ATP5F1B = "ATP5B", ATP5F1C = "ATP5C1", ATP5F1D = "ATP5D",
    ATP5F1E = "ATP5E", ATP5MC1 = "ATP5G1", ATP5MC2 = "ATP5G2", ATP5MC3 = "ATP5G3",
    ATP5ME = "ATP5I", ATP5MF = "ATP5J2", ATP5MG = "ATP5L", ATP5PB = "ATP5F1",
    ATP5PD = "ATP5H", ATP5PF = "ATP5J", ATP5PO = "ATP5O"
)
fallback <- unname(ox_aliases[ox_genes])
use_alias <- !ox_genes %chin% gene_ids & fallback %chin% gene_ids
ox_genes[use_alias] <- fallback[use_alias]
ox_genes <- intersect(ox_genes, gene_ids)
stopifnot(length(ox_genes) == 116L)
scoring_genes <- lapply(signatures, function(symbols) {
    setdiff(intersect(symbols, gene_ids), ox_genes)
})
stopifnot(all(lengths(scoring_genes) > 0L),
          !any(unlist(scoring_genes, use.names = FALSE) %chin% gene_ids[duplicated(gene_ids)]))

# The score of a cell is the arithmetic mean of its matched gene activities.
chunk_size <- 20000L
for (program in names(scoring_genes)) {
    gene_rows <- match(scoring_genes[[program]], gene_ids)
    values <- numeric(length(cell_cols))
    for (first in seq.int(1L, length(cell_cols), by = chunk_size)) {
        index <- first:min(first + chunk_size - 1L, length(cell_cols))
        values[index] <- Matrix::colMeans(activity[gene_rows, cell_cols[index], drop = FALSE])
    }
    stopifnot(all(is.finite(values)))
    set(cells, j = program, value = values)
}

# Fit Y = intercept + G + O_step + G:O_step for all programs within each sample.
# WITHOUT slope = beta_O. G+ slope = beta_O + beta_G:O.
sample_effect <- cells[, {
    design <- cbind(1, G, O_step, G * O_step)
    fit <- lm.fit(design, as.matrix(.SD))
    stopifnot(fit$rank == 4L)
    coefficients <- fit$coefficients
    .(
        program = names(scoring_genes),
        slope_WITHOUT = as.numeric(coefficients[3L, ]),
        slope_Gplus = as.numeric(coefficients[3L, ] + coefficients[4L, ])
    )
}, by = .(patient, sample_id), .SDcols = names(scoring_genes)]

# Weight the two statuses equally within sample, then weight samples within
# patient and patients within the cohort equally. Use a patient-level t interval.
patient_effect <- sample_effect[, .(
    slope_WITHOUT = mean(slope_WITHOUT),
    slope_Gplus = mean(slope_Gplus),
    mean_slope = mean((slope_WITHOUT + slope_Gplus) / 2),
    n_samples = .N
), by = .(patient, program)]

result <- patient_effect[, {
    n <- .N
    estimate <- mean(mean_slope)
    se <- sd(mean_slope) / sqrt(n)
    half_width <- qt(0.975, df = n - 1L) * se
    .(
        n_patients = n,
        n_samples = sum(n_samples),
        slope_WITHOUT = mean(slope_WITHOUT),
        slope_Gplus = mean(slope_Gplus),
        mean_OXPHOS_slope = estimate,
        ci_low = estimate - half_width,
        ci_high = estimate + half_width,
        p = if (se > 0) 2 * pt(-abs(estimate / se), df = n - 1L) else NA_real_
    )
}, by = program]
result[, genes_used := lengths(scoring_genes)[program]]
result[, p_adj := p.adjust(p, method = "BH", n = length(scoring_genes))]
result <- result[match(names(scoring_genes), result$program)]

plot <- ggplot(result, aes(mean_OXPHOS_slope, factor(program, levels = rev(names(scoring_genes))))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#777777", linewidth = 0.45) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                  width = 0.16, linewidth = 0.55, colour = "#344452") +
    geom_point(size = 2.8, colour = "#27905F") +
    labs(x = "Mean score change per 0.01 increase in oxygen", y = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.text.y = element_text(colour = "#202020"), axis.ticks.y = element_blank())

dir.create(dirname(output_stem), recursive = TRUE, showWarnings = FALSE)
pdf_file <- paste0(output_stem, ".pdf")
csv_file <- paste0(output_stem, ".csv")
ggsave(pdf_file, plot, width = 6.5, height = 4.1, device = grDevices::pdf, useDingbats = FALSE)
fwrite(result, csv_file)
print(result)
cat(pdf_file, csv_file, sep = "\n")
