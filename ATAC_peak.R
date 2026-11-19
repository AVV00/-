#!/usr/bin/env Rscript

# PCBD2 accessibility across the full annotated gene body and flanking sequence.
# Four figures with two tracks each. Run with Rscript peak_full_locus.R.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(SeuratObject)
  library(Signac)
  library(GenomicRanges)
  library(Rsamtools)
  library(patchwork)
})

root <- "/mnt/data1/linyanlin/scATAC"
hdd <- "/mnt/hdd/linyanlin/sc-ATAC/singlecellfilter"
input_dir <- file.path(root, "picture/pic4/261117")
out_dir <- file.path(root, "picture/pic4/261119/peak")
check_dir <- file.path(root, "picture/check/261105")
stem <- "Fig4_PCBD2_ATAC_tracks_full_locus"
panels <- c("A. HCC1 | 12539", "B. HCC4 | 2690", "C. HCC7 | 16399", "D. Module C")
panel_files <- c("A_HCC1_12539", "B_HCC4_2690", "C_HCC7_16399", "D_ModuleC")
mutation_patients <- c("HCC1", "HCC4", "HCC7")
target_positions <- c(12539L, 2690L, 16399L)
# The mutation table and the ATAC samples use different patient IDs.
panel_atac_patients <- c("HCC3", "HCC8", "HCC11", NA_character_)
statuses <- c("Mutation-", "Mutation+")
tracks <- paste(rep(panels, each = 2L), rep(c("-", "+"), length(panels)))
bin_width <- 100L
flank_bp <- 15000L
n_tracks <- length(tracks)

object_name <- "Seuratobj_Patients_Harmony_gene_activity_manual_cellannotation"
objects <- c(
  old = file.path(hdd, "step4_picture/manual", paste0(object_name, "_updated.rds")),
  pv = file.path(hdd, "step4_HCC12", paste0(object_name, "_cluster_fixed.rds")),
  fix = file.path(hdd, "HCC3,8/step4", paste0(object_name, "_cluster_fixed.rds")),
  new = file.path(root, "new2sample_scATAC/step3_cellannotation/manual",
                  paste0(object_name, "_updated.rds"))
)

strip_prefix <- function(x) sub("^(old_|new_)", "", trimws(x))

short_cell_id <- function(x) {
  ok <- !is.na(x) & grepl("CELL[0-9]+_N[0-9]+$", x)
  ans <- rep(NA_character_, length(x))
  ans[ok] <- sub("^.*?(CELL[0-9]+_N[0-9]+)$", "\\1", x[ok], perl = TRUE)
  ans
}

sample_key <- function(x) {
  # The audit supports these HCC12 aliases; the other tissue names stay distinct.
  sub("^(HCC12_C[12])_CD(450|451)$", "\\1", strip_prefix(x))
}

load_groups <- function() {
  cohort <- fread(
    file.path(input_dir, "Fig4A_PCBD2_gene_activity_cells.csv"),
    select = c("panel", "cell", "sample", "tissue_cell", "has_linked_mutation_record")
  )[panel == "Module C"]
  stopifnot(nrow(cohort) == 168140L, !anyDuplicated(cohort$cell))

  bridge <- unique(fread(
    file.path(check_dir, "00_all_cells_cellname_match.csv"),
    select = c("cell", "sample", "mgatk_Cell")
  ))
  stopifnot(!anyDuplicated(bridge$cell))
  bridge_row <- match(cohort$cell, bridge$cell)
  stopifnot(!anyNA(bridge_row))
  rebuilt <- paste(trimws(bridge$sample[bridge_row]), trimws(bridge$mgatk_Cell[bridge_row]),
                   sep = "_")
  stopifnot(
    all(trimws(cohort$sample) == trimws(bridge$sample[bridge_row])),
    all(trimws(cohort$tissue_cell) == rebuilt)
  )
  cohort[, `:=`(sample = trimws(sample), tissue_cell = rebuilt)]
  bridge[, short_id := short_cell_id(trimws(mgatk_Cell))]
  bridge[, join_key := paste(sample_key(sample), short_id, sep = "|")]
  targets <- bridge[!is.na(short_id), .(n_targets = uniqueN(cell)), by = join_key]
  cohort[, join_key := bridge$join_key[bridge_row]]
  cohort[, n_targets := targets$n_targets[match(join_key, targets$join_key)]]
  cohort[, linked_before := as.logical(has_linked_mutation_record) %in% TRUE]

  mutations <- fread(
    file.path(check_dir, "绐佸彉姹囨€讳俊鎭痏鍚玂XPHOS鍒嗙粍.csv"),
    select = c("tissue_cell", "patient", "mutationtype", "Position"),
    colClasses = "character"
  )
  mutations[, (names(mutations)) := lapply(.SD, trimws)]
  mutations[, Position := suppressWarnings(as.numeric(gsub(",", "", Position, fixed = TRUE)))]
  mutations[, full_id := strip_prefix(tissue_cell)]
  mutations[grepl("_CELL[0-9]+_N[0-9]+$", full_id), join_key :=
    paste(sample_key(sub("_CELL[0-9]+_N[0-9]+$", "", full_id)), short_cell_id(full_id),
          sep = "|")]
  qualified <- unique(mutations[!is.na(join_key), .(join_key, tissue_cell, patient)])
  candidates <- cohort[!linked_before & join_key %chin% qualified$join_key]
  new_cells <- candidates[n_targets %in% 1L, .(join_key, mapped_key = tissue_cell)]
  stopifnot(!anyDuplicated(new_cells$join_key))
  new_links <- unique(
    merge(qualified, new_cells, by = "join_key")[, .(tissue_cell, patient, mapped_key)]
  )
  old_links <- fread(
    file.path(input_dir, "Fig4A_PCBD2_gene_activity_used_mutation_mapping.csv"),
    select = c("tissue_cell", "patient", "mapped_key")
  )
  old_links[, (names(old_links)) := lapply(.SD, trimws)]
  links <- unique(rbindlist(list(old_links, new_links), use.names = TRUE))
  links <- links[!is.na(mapped_key) & nzchar(mapped_key)]
  if (anyDuplicated(links[, .(tissue_cell, patient)])) {
    stop("A mutation row has conflicting links. Check the mutation annotation.")
  }
  linked <- merge(mutations, links, by = c("tissue_cell", "patient"))
  positives <- lapply(seq_along(target_positions), function(i) {
    unique(linked[patient == mutation_patients[i] & Position == target_positions[i], mapped_key])
  })
  positives[[4L]] <- unique(linked[
    toupper(gsub("[^[:alnum:]]", "", mutationtype)) == "MODULEC", mapped_key
  ])
  cohort[, atac_patient := sub("_.*$", "", strip_prefix(sample))]
  cohort[, (panels) := lapply(positives, function(ids) tissue_cell %chin% ids)]
  for (k in seq_along(target_positions)) {
    if (any(cohort[[panels[k]]] & cohort$atac_patient != panel_atac_patients[k])) {
      stop("Positive cells outside the expected ATAC patient for ", panels[k])
    }
  }
  cohort[, source := fcase(
    grepl("^HCC1[34]_", sample), "new",
    sample == "HCC12_PV_CD450", "pv",
    sample == "HCC8_C2" | grepl("^HCC3_", sample), "fix",
    grepl("^HCC(6|7|8|9|10|11|12)_", sample), "old",
    default = "unknown"
  )]
  if (any(cohort$source == "unknown")) {
    stop("Unassigned samples: ", paste(cohort[source == "unknown", unique(sample)],
                                       collapse = ", "))
  }
  cat("Input cells:", nrow(cohort), "\n")
  cat("New unique mutation links:", nrow(new_cells), "\n")
  cat("New ambiguous candidates:", nrow(candidates) - nrow(new_cells), "\n")
  input_counts <- rbindlist(lapply(seq_along(panels), function(k) {
    panel_name <- panels[k]
    patient_id <- panel_atac_patients[k]
    in_panel <- if (is.na(patient_id)) cohort else cohort[atac_patient == patient_id]
    positive <- sum(in_panel[[panel_name]])
    data.table(panel = panel_name, negative = nrow(in_panel) - positive, positive = positive)
  }))
  print(input_counts)
  if (any(input_counts$positive == 0L | input_counts$negative == 0L)) {
    stop("At least one panel has no positive or negative cells.")
  }
  list(cohort = cohort, new_links = new_links, input_counts = input_counts)
}

read_locus <- function(atac) {
  annotation <- Signac::Annotation(atac)
  names(annotation) <- NULL
  annotation <- as.data.table(as.data.frame(annotation))
  required <- c("gene_name", "seqnames", "start", "end", "strand")
  if (!all(required %in% names(annotation))) stop("ATAC gene annotation is incomplete.")
  if (!"tx_id" %in% names(annotation)) annotation[, tx_id := NA_character_]
  genes <- annotation[toupper(gene_name) == "PCBD2" & strand %in% c("+", "-")]
  tx <- genes[, .(lo = min(start), hi = max(end)), by = .(seqnames, strand, tx_id)]
  if (!nrow(tx)) stop("PCBD2 missing from the ATAC annotation.")
  tx[, span := hi - lo]
  setorder(tx, -span, seqnames, lo, hi, tx_id)
  chosen <- tx[1L]
  chr <- as.character(chosen$seqnames)
  strand <- as.character(chosen$strand)
  tss <- if (strand == "+") chosen$lo else chosen$hi

  # Use the full extent of the annotated PCBD2 locus on the chosen chromosome
  # and strand. Keep the longest transcript for its TSS and exon marks.
  locus_genes <- genes[seqnames == chosen$seqnames & strand == chosen$strand]
  gene_lo <- min(locus_genes$start)
  gene_hi <- max(locus_genes$end)
  lo <- max(1L, gene_lo - flank_bp)
  hi <- gene_hi + flank_bp
  stopifnot(is.finite(lo), is.finite(hi), lo <= gene_lo, gene_hi <= hi)

  exons <- data.table(start = integer(), end = integer())
  if ("type" %in% names(locus_genes) && !is.na(chosen$tx_id)) {
    exons <- unique(locus_genes[type == "exon" & tx_id == chosen$tx_id, .(start, end)])
  }
  list(
    chr = chr, strand = strand, tx_id = as.character(chosen$tx_id),
    tss = tss, lo = lo, hi = hi,
    gene_lo = gene_lo, gene_hi = gene_hi, exons = exons,
    promoter_lo = if (strand == "+") tss - 2000L else tss - 499L,
    promoter_hi = if (strand == "+") tss + 499L else tss + 2000L
  )
}

get_depth <- function(atac, metadata, sample_cells, ids, sample_name) {
  depth <- metadata[["nCount_ATAC"]][ids]
  if (length(depth) != length(sample_cells) || any(!is.finite(depth) | depth <= 0)) {
    counts <- tryCatch(
      SeuratObject::GetAssayData(atac, layer = "counts"),
      error = function(e) SeuratObject::GetAssayData(atac, slot = "counts")
    )
    cols <- match(sample_cells, colnames(counts))
    if (anyNA(cols)) stop("Missing cells in ATAC counts: ", sample_name)
    depth <- Matrix::colSums(counts[, cols, drop = FALSE])
  }
  if (any(!is.finite(depth) | depth <= 0)) stop("Invalid ATAC depth: ", sample_name)
  depth
}

read_fragment_ends <- function(path, locus, barcodes, all_barcodes) {
  if (!all(file.access(c(path, paste0(path, ".tbi")), 4L) == 0L)) {
    stop("Unreadable fragment file or index: ", path)
  }
  tabix <- TabixFile(path)
  indexed <- seqnamesTabix(tabix)
  chr <- locus$chr
  if (!chr %in% indexed) {
    chr <- if (startsWith(chr, "chr")) sub("^chr", "", chr) else paste0("chr", chr)
  }
  if (!chr %in% indexed) stop("PCBD2 chromosome missing in index: ", path)
  records <- scanTabix(tabix, param = GRanges(chr, IRanges::IRanges(locus$lo, locus$hi)))[[1L]]
  if (!length(records)) return(list(bin = integer(), cell = integer(), pos = integer()))

  fragments <- fread(text = paste0(paste(records, collapse = "\n"), "\n"),
                     sep = "\t", header = FALSE, select = 1:4, showProgress = FALSE)
  if (!any(fragments$V4 %chin% unname(all_barcodes))) {
    stop("Stored barcodes do not match fragments near PCBD2: ", path)
  }
  cell <- match(fragments$V4, barcodes)
  hit <- !is.na(cell)
  ends <- c(fragments$V2[hit] + 1L, fragments$V3[hit])
  keep <- ends >= locus$lo & ends <= locus$hi
  list(bin = 1L + (ends[keep] - locus$lo) %/% bin_width, cell = rep(cell[hit], 2L)[keep], pos = ends[keep])
}

count_tracks <- function(cohort, input_counts) {
  cells <- depth <- before_filter <- numeric(n_tracks)
  before_depth <- promoter_before <- promoter_after <- numeric(n_tracks)
  sample_qc <- list()
  locus <- NULL
  for (src in names(objects)[names(objects) %chin% cohort$source]) {
    path <- objects[[src]]
    cat("\nReading RDS:", path, "\n")
    flush.console()
    ob <- readRDS(path)
    if (!"ATAC" %in% Assays(ob)) stop("ATAC assay missing: ", path)
    atac <- ob[["ATAC"]]
    metadata <- ob[[]]
    fragments <- Signac::Fragments(atac)
    filenames <- vapply(fragments, function(f) {
      basename(Signac::GetFragmentData(f, slot = "path"))
    }, character(1L))
    current <- read_locus(atac)
    if (is.null(locus)) {
      locus <- current
      bin_starts <- seq.int(locus$lo, locus$hi, by = bin_width)
      bin_ends <- pmin(locus$hi, bin_starts + bin_width - 1L)
      n_bins <- length(bin_starts)
      raw <- matrix(0, n_bins, n_tracks)
      cat("PCBD2:", locus$chr, "gene", locus$gene_lo, "-", locus$gene_hi,
          "TSS", locus$tss, "strand", locus$strand, "transcript", locus$tx_id, "\n")
      cat("Counted window:", locus$lo, "-", locus$hi,
          "(", locus$hi - locus$lo + 1L, "bp;", n_bins, "bins )\n")
    } else if (sub("^chr", "", current$chr) != sub("^chr", "", locus$chr) ||
               current$strand != locus$strand || current$tss != locus$tss ||
               current$gene_lo != locus$gene_lo || current$gene_hi != locus$gene_hi) {
      stop("PCBD2 chromosome, TSS or gene bounds disagree across ATAC objects: ", path)
    }

    fragment_dir <- if (src == "new") {
      file.path(root, "new2sample_scATAC/step1_filter")
    } else {
      file.path(hdd, "step1_filter")
    }
    for (smp in sort(unique(cohort[source == src, sample]))) {
      cat("Counting:", smp, "\n")
      flush.console()
      u <- cohort[source == src & sample == smp]
      sample_patient <- unique(u$atac_patient)
      stopifnot(length(sample_patient) == 1L)
      sample_panels <- which(is.na(panel_atac_patients) | panel_atac_patients == sample_patient)
      sample_track_qc <- vector("list", length(panels))
      fn <- paste0(smp, "_fragments.tsv.gz")
      hit <- which(filenames == fn)
      if (length(hit) != 1L) stop("No unique fragment map for ", smp)
      map <- Signac::GetFragmentData(fragments[[hit]], slot = "cells")
      barcode <- unname(map[u$cell])
      ids <- match(u$cell, rownames(metadata))
      if (anyNA(ids) || length(barcode) != nrow(u) || anyNA(barcode) ||
          any(!nzchar(barcode)) || anyDuplicated(barcode)) {
        stop("Incomplete or duplicate ATAC cell/barcode matches in ", smp)
      }
      cell_depth <- get_depth(atac, metadata, u$cell, ids, smp)
      ends <- read_fragment_ends(file.path(fragment_dir, fn), locus, barcode, map)
      in_promoter <- ends$pos >= locus$promoter_lo & ends$pos <= locus$promoter_hi
      promoter_ends <- tabulate(ends$cell[in_promoter], nbins = nrow(u))
      promoter_score <- 1e6 * promoter_ends / cell_depth
      for (k in sample_panels) {
        grp <- 1L + as.integer(u[[panels[k]]])
        is_negative <- grp == 1L
        eligible <- which(is_negative & promoter_ends > 0L)
        n_drop <- min(length(eligible), floor(0.05 * sum(is_negative)))
        drop <- head(eligible[order(-promoter_score[eligible], u$cell[eligible], method = "radix")], n_drop)
        keep <- rep(TRUE, nrow(u))
        keep[drop] <- FALSE
        keep_ends <- keep[ends$cell]
        cols <- c(2L * k - 1L, 2L * k)
        input_cells <- tabulate(grp, nbins = 2L)
        kept_cells <- tabulate(grp[keep], nbins = 2L)
        before_filter[cols] <- before_filter[cols] + input_cells
        before_depth[cols] <- before_depth[cols] + vapply(1:2, function(j) sum(cell_depth[grp == j]), numeric(1L))
        promoter_before[cols] <- promoter_before[cols] + vapply(1:2, function(j) sum(promoter_ends[grp == j]), numeric(1L))
        promoter_after[cols] <- promoter_after[cols] + vapply(1:2, function(j) sum(promoter_ends[keep & grp == j]), numeric(1L))
        cells[cols] <- cells[cols] + kept_cells
        sample_track_qc[[k]] <- data.table(
          panel = panels[k], sample = smp, group = statuses, input_cells = input_cells,
          used_cells = kept_cells, excluded_high_PCBD2 = input_cells - kept_cells
        )
        depth[cols] <- depth[cols] + vapply(1:2, function(i) sum(cell_depth[keep & grp == i]),
                                          numeric(1L))
        raw[, cols] <- raw[, cols] + matrix(
          tabulate(ends$bin[keep_ends] + n_bins * (grp[ends$cell[keep_ends]] - 1L), nbins = 2L * n_bins),
          nrow = n_bins
        )
      }
      sample_qc[[smp]] <- rbindlist(sample_track_qc)
    }
    rm(ob, atac, metadata, fragments)
    invisible(gc())
  }
  stopifnot(!is.null(locus), all(cells > 0), all(depth > 0), sum(raw) > 0)
  expected_cells <- as.vector(rbind(input_counts$negative, input_counts$positive))
  stopifnot(all(before_filter == expected_cells))
  positive_cols <- seq(2L, n_tracks, by = 2L)
  stopifnot(all(cells[positive_cols] == before_filter[positive_cols]))
  centers <- (bin_starts + bin_ends) / 2
  coverage <- data.table(
    pos = rep(centers, n_tracks),
    bin_start = rep(bin_starts, n_tracks),
    bin_end = rep(bin_ends, n_tracks),
    panel = rep(panels, each = length(statuses) * n_bins),
    track = factor(rep(tracks, each = n_bins), levels = tracks),
    status = factor(rep(rep(statuses, length(panels)), each = n_bins), levels = statuses),
    fragment_ends = as.vector(raw),
    coverage = as.vector(sweep(raw, 2L, depth, "/") * 1e6)
  )
  qc <- data.table(panel = rep(panels, each = length(statuses)),
                   group = rep(statuses, length(panels)),
                   track = tracks, cells = cells, ATAC_peak_counts = depth,
                   fragment_ends = colSums(raw), input_cells = before_filter,
                   excluded_high_PCBD2 = before_filter - cells,
                   promoter_CPM_before = 1e6 * promoter_before / before_depth,
                   promoter_CPM_after = 1e6 * promoter_after / depth)
  list(coverage = coverage, qc = qc, sample_qc = rbindlist(sample_qc), locus = locus)
}

make_figure <- function(coverage, locus, panel_name) {
  palette <- c("Mutation-" = "#62A9CC", "Mutation+" = "#E98A50")
  window <- c(locus$lo, locus$hi)
  promoter <- function(top) annotate(
    "rect", xmin = locus$promoter_lo, xmax = locus$promoter_hi,
    ymin = 0, ymax = top, fill = "#F5DECF", alpha = 0.55
  )
  tss_line <- geom_vline(xintercept = locus$tss, linetype = "dashed",
                         colour = "#8B6049", linewidth = 0.35)

  coverage_plot <- ggplot(coverage, aes(pos, coverage)) +
    promoter(Inf) +
    geom_area(aes(fill = status), stat = "identity", position = "identity", alpha = 0.75) +
    geom_line(aes(colour = status), linewidth = 0.35) +
    tss_line +
    facet_grid(rows = vars(status), switch = "y") +
    scale_fill_manual(values = palette, guide = "none") +
    scale_colour_manual(values = palette, guide = "none") +
    coord_cartesian(xlim = window, expand = FALSE) +
    scale_y_sqrt() +
    labs(title = paste(panel_name, "- PCBD2 chromatin accessibility"), x = NULL,
         y = "Fragment ends per million ATAC peak counts\n(square-root scale)") +
    theme_minimal(base_size = 12) +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  gene_range <- c(locus$gene_lo, locus$gene_hi)
  arrow_range <- if (locus$strand == "+") gene_range else rev(gene_range)
  gene_plot <- ggplot() +
    promoter(1) +
    annotate("segment", x = arrow_range[1L], xend = arrow_range[2L], y = 0.45, yend = 0.45,
             linewidth = 0.45, arrow = grid::arrow(length = grid::unit(2, "mm"))) +
    geom_rect(data = locus$exons, aes(xmin = start, xmax = end), ymin = 0.32, ymax = 0.58,
              fill = "#313D4A") +
    tss_line +
    annotate("text", x = mean(gene_range), y = 0.8, label = "PCBD2",
             fontface = "bold", size = 3.5) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",",
                                                 scientific = FALSE, trim = TRUE)) +
    coord_cartesian(xlim = window, ylim = c(0, 1), expand = FALSE) +
    labs(x = paste(locus$chr, "position (bp)"), y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank())
  coverage_plot / gene_plot + plot_layout(heights = c(3, 0.8))
}

save_outputs <- function(results, new_links) {
  stage <- tempfile("peak-", tmpdir = out_dir)
  stopifnot(dir.create(stage))
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)

  for (k in seq_along(panels)) {
    panel_name <- panels[k]
    plot_stem <- file.path(stage, paste(stem, panel_files[k], sep = "_"))
    figure <- make_figure(results$coverage[panel == panel_name], results$locus, panel_name)
    ggsave(paste0(plot_stem, ".pdf"), figure, width = 12, height = 4.6,
           device = grDevices::pdf, useDingbats = FALSE, bg = "white")
    ggsave(paste0(plot_stem, ".png"), figure, width = 12, height = 4.6,
           dpi = 300, bg = "white")
  }

  table_stem <- file.path(stage, stem)
  fwrite(results$qc, paste0(table_stem, "_groups.csv"))
  fwrite(results$sample_qc, paste0(table_stem, "_sample_qc.csv"))
  fwrite(results$coverage, paste0(table_stem, "_coverage.csv"))
  fwrite(new_links, paste0(table_stem, "_new_mutation_links.csv"))

  staged <- list.files(stage, full.names = TRUE)
  if (length(staged) != 2L * length(panels) + 4L || any(file.size(staged) <= 0)) {
    stop("Failed to render the plots or result tables.")
  }
  outputs <- file.path(out_dir, basename(staged))
  if (!all(file.rename(staged, outputs))) stop("Could not save all outputs.")
  outputs
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create output directory: ", out_dir)
groups <- load_groups()
results <- count_tracks(groups$cohort, groups$input_counts)
print(results$qc[, .(panel, group, input_cells, cells, excluded_high_PCBD2,
                      promoter_CPM_before, promoter_CPM_after)])
outputs <- save_outputs(results, groups$new_links)
cat("\nSaved:\n", paste(outputs, collapse = "\n"), "\n", sep = "")


