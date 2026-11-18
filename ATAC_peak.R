suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(SeuratObject)
  library(Signac)
  library(GenomicRanges)
  library(Rsamtools)
})

gene = toupper(Sys.getenv("FIG6_GENE", "SLC25A21"))
stopifnot(grepl("^[A-Z0-9-]+$", gene))
r = "/mnt/data1/linyanlin/scATAC"
h = "/mnt/hdd/linyanlin/sc-ATAC/singlecellfilter"
o = file.path(r, "picture/pic6")
base = "Seuratobj_Patients_Harmony_gene_activity_manual_cellannotation"
files = c(
  old = file.path(h, "step4_picture/manual", paste0(base, "_updated.rds")),
  pv = file.path(h, "step4_HCC12", paste0(base, "_cluster_fixed.rds")),
  fix = file.path(h, "HCC3,8/step4", paste0(base, "_cluster_fixed.rds")),
  new = file.path(r, "new2sample_scATAC/step3_cellannotation/manual", paste0(base, "_updated.rds"))
)

a = fread(
  file.path(o, "Fig6_ModuleC_OXPHOS_cell_groups.csv"),
  select = c("cell", "ATAC_sample", "group"),
  colClasses = "character"
)
a[, (names(a)) := lapply(.SD, trimws)]
lev = c("Low-Low", "Low-High", "High-Low", "High-High")
stopifnot(nrow(a) > 0L, !anyNA(a), !anyDuplicated(a$cell), all(a$group %chin% lev))
a[, src := fcase(
  grepl("^HCC1[34]_", ATAC_sample), "new",
  ATAC_sample == "HCC12_PV_CD450", "pv",
  ATAC_sample == "HCC8_C2" | grepl("^HCC3_", ATAC_sample), "fix",
  grepl("^HCC(6|7|8|9|10|11|12)_", ATAC_sample), "old",
  default = "unknown"
)]
if (any(a$src == "unknown")) stop("Unassigned samples: ", paste(unique(a[src == "unknown", ATAC_sample]), collapse = ", "))
needed = unique(a$src)
if (any(!file.exists(files[needed]))) stop("Missing RDS: ", paste(files[needed][!file.exists(files[needed])], collapse = ", "))

bin = 100L
loc = NULL
nc = nd = numeric(4)
for (s0 in needed) {
  cat("Loading:", files[[s0]], "\n")
  ob = readRDS(files[[s0]])
  if (!"ATAC" %in% Assays(ob) || !inherits(ob[["ATAC"]], "ChromatinAssay")) stop("ATAC assay missing: ", files[[s0]])
  at = ob[["ATAC"]]
  md = ob@meta.data
  z = a[src == s0]
  ii = match(z$cell, rownames(md))
  if (anyNA(ii)) stop(sum(is.na(ii)), " grouped cells missing from ", s0)
  ct = intersect(c("manual_celltype", "celltype_final"), names(md))
  if (!length(ct)) stop("No manual cell type column in ", s0)
  typ = gsub("[^a-z]", "", tolower(trimws(as.character(md[[ct[1L]]][ii]))))
  ok = typ %chin% c("malignant", "malignantcell", "malignantcells")
  cat("Manual annotation:", s0, " | retained:", sum(ok), " | excluded:", sum(!ok), "\n")
  z = z[ok]
  if (!nrow(z)) {
    rm(ob, at, md)
    invisible(gc())
    next
  }

  if (is.null(loc)) {
    an0 = Signac::Annotation(at)
    names(an0) = NULL
    an = as.data.table(as.data.frame(an0))
    rm(an0)
    if (!all(c("gene_name", "seqnames", "start", "end", "strand") %in% names(an))) stop("Gene coordinates missing in ", s0)
    if (!"tx_id" %in% names(an)) an[, tx_id := NA_character_]
    an = an[toupper(gene_name) == gene & strand %in% c("+", "-"),
            .(lo = min(start), hi = max(end)),
            by = .(seqnames, strand, tx_id)]
    if (!nrow(an)) stop("Gene missing from annotation: ", gene)
    an[, span := hi - lo]
    setorder(an, -span, seqnames, lo)
    loc = an[1L]
    chr = as.character(loc$seqnames)
    tss = if (loc$strand == "+") loc$lo else loc$hi
    lo = max(1L, tss - 15000L)
    hi = tss + 15000L
    nb = ceiling((hi - lo + 1L) / bin)
    raw = matrix(0, nrow = nb, ncol = 4L)
  }

  dep = md$nCount_ATAC[match(z$cell, rownames(md))]
  if (length(dep) != nrow(z) || any(!is.finite(dep) | dep <= 0)) {
    pk = tryCatch(
      GetAssayData(at, layer = "counts"),
      error = function(e) GetAssayData(at, slot = "counts")
    )
    ic = match(z$cell, colnames(pk))
    if (anyNA(ic)) stop("Cells missing from ATAC counts: ", s0)
    dep = Matrix::colSums(pk[, ic, drop = FALSE])
    rm(pk)
  }
  stopifnot(all(is.finite(dep) & dep > 0))
  z[, depth := dep]
  fr = Signac::Fragments(at)
  fn = vapply(fr, function(f) basename(Signac::GetFragmentData(f, slot = "path")), character(1L))
  for (smp in unique(z$ATAC_sample)) {
    u = z[ATAC_sample == smp]
    j = match(u$group, lev)
    bn = paste0(smp, "_fragments.tsv.gz")
    jfr = which(fn == bn)
    if (length(jfr) != 1L) stop("No unique stored fragment map for ", smp)
    cm = Signac::GetFragmentData(fr[[jfr]], slot = "cells")
    bc = unname(cm[u$cell])
    if (length(bc) != nrow(u) || anyNA(bc) || any(!nzchar(bc)) || anyDuplicated(bc)) stop("Incomplete or duplicated fragment barcodes for ", smp)
    dir = if (grepl("^HCC1[34]_", smp)) file.path(r, "new2sample_scATAC/step1_filter") else file.path(h, "step1_filter")
    f = file.path(dir, bn)
    if (!all(file.access(c(f, paste0(f, ".tbi")), 4L) == 0L)) stop("Unreadable fragment or index: ", f)
    tf = TabixFile(f)
    chrom = chr
    if (!chrom %in% seqnamesTabix(tf)) stop("Chromosome not indexed in ", f, ": ", chr)
    v = scanTabix(tf, param = GRanges(chrom, IRanges::IRanges(lo, hi)))[[1L]]
    nc = nc + tabulate(j, 4L)
    nd = nd + vapply(1:4, function(k) sum(u$depth[j == k]), numeric(1L))
    if (length(v)) {
      d = fread(text = paste(v, collapse = "\n"), sep = "\t", header = FALSE, select = 1:4, showProgress = FALSE)
      m = match(d$V4, bc)
      if (!any(d$V4 %chin% unname(cm))) stop("Stored fragment barcodes do not match reads at this locus: ", smp)
      hit = !is.na(m)
      pos = c(d$V2[hit] + 1L, d$V3[hit])
      grp = rep(j[m[hit]], 2L)
      take = pos >= lo & pos <= hi
      bi = 1L + (pos[take] - lo) %/% bin
      raw = raw + matrix(tabulate(bi + (grp[take] - 1L) * nb, nbins = 4L * nb), nrow = nb, ncol = 4L)
    }
    cat("Counted:", smp, " | cells:", nrow(u), "\n")
  }
  rm(ob, at, md, fr)
  invisible(gc())
}

stopifnot(!is.null(loc), all(nc > 0), all(nd > 0), sum(raw) > 0)
p0 = max(lo, if (loc$strand == "+") tss - 2000L else tss - 499L)
p1 = min(hi, if (loc$strand == "+") tss + 499L else tss + 2000L)
z = data.table(
  pos = rep(seq(from = lo + bin / 2, by = bin, length.out = nb), 4L),
  group = factor(rep(lev, each = nb), levels = lev),
  coverage = as.vector(sweep(raw, 2, nd, "/") * 1e6)
)
pal = setNames(c("#244F73", "#77BBD6", "#E3A14D", "#D9685D"), lev)
lab = setNames(paste0(lev, "\nn = ", format(nc, big.mark = ",", trim = TRUE)), lev)
p = ggplot(z, aes(pos, coverage)) +
  annotate("rect", xmin = p0, xmax = p1, ymin = 0, ymax = Inf, fill = "#F5DECF", alpha = .55) +
  geom_area(aes(fill = group), alpha = .70) +
  geom_line(aes(colour = group), linewidth = .35) +
  geom_vline(xintercept = tss, linetype = 2, linewidth = .35, colour = "#8B6049") +
  facet_grid(rows = vars(group), switch = "y", labeller = labeller(group = as_labeller(lab))) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(labels = function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)) +
  labs(
    title = paste("Figure 6D. ATAC accessibility near", gene),
    subtitle = paste0(chr, ":", lo, "-", hi, " | 100-bp bins | Four cell states"),
    x = paste(chr, "position (bp)"),
    y = "Fragment ends per million ATAC peak counts",
    caption = "Dashed line = TSS. Shaded region = promoter. Tracks share one coverage scale."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )
out = file.path(o, paste0("Fig6D_ATAC_tracks_", gene, ".pdf"))
ggsave(out, p, width = 10, height = 7, device = grDevices::pdf, useDingbats = FALSE)
print(data.table(group = lev, cells = nc, peak_depth = nd))
cat("PDF:", out, "\n")
