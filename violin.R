
#!/usr/bin/env Rscript

# ============================================================
# Nuclear OXPHOS Complex I–V Violin Plot
# ============================================================
#
# Input:
#   OXPHOS_5complexes_61145cells_gene_activity.csv
#
# Analysis:
#   1. Rank cells by MSS
#   2. Select lowest 1000 cells (Low MSS)
#   3. Select highest 1000 cells (High MSS)
#   4. Calculate mean Gene Activity per cell per complex
#   5. Plot violin + boxplot + median point
#   6. Annotate median difference (High - Low)
#
# Colors:
#   High MSS: #D98C84
#   Low MSS:  #8FA5C9
#
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


# ------------------------------------------------------------
# 2. Set working directory and input
# ------------------------------------------------------------

setwd(
  "/mnt/data1/linyanlin/scATAC/picture/pic3/final/heatmap"
)

f <- "OXPHOS_5complexes_61145cells_gene_activity.csv"

d <- fread(f)

nm <- names(d)


# ------------------------------------------------------------
# 3. Identify relevant columns
# ------------------------------------------------------------

mc <- grep(
  "(^mss$|mss_score|mito.*stress|stress.*score)",
  nm,
  ignore.case = TRUE,
  value = TRUE
)[1]

cc <- grep(
  "(^cell$|cell_id|barcode)",
  nm,
  ignore.case = TRUE,
  value = TRUE
)[1]

xc <- grep(
  "^complex$|complex_id",
  nm,
  ignore.case = TRUE,
  value = TRUE
)[1]

sc <- grep(
  "^score$|^gene_activity$|^activity$|^mean_activity$|^avg_activity$|^value$",
  nm,
  ignore.case = TRUE,
  value = TRUE
)[1]

if (is.na(mc) || is.na(cc)) {
  stop("Cannot find MSS or Cell column")
}


# ------------------------------------------------------------
# 4. Extract cell ID and MSS
# ------------------------------------------------------------

d[, Cell := as.character(get(cc))]

d[, MSS_num := suppressWarnings(as.numeric(get(mc)))]

ord <- unique(
  d[is.finite(MSS_num), .(Cell, MSS_num)]
)

if (anyDuplicated(ord$Cell)) {
  stop("Duplicated cells with different MSS")
}

setorder(ord, MSS_num, Cell)

if (nrow(ord) < 2000) {
  stop("Fewer than 2000 MSS cells")
}


# ------------------------------------------------------------
# 5. Select extreme MSS cells
# ------------------------------------------------------------

# Lowest 1000 MSS cells
low_cells <- ord[
  1:1000,
  .(Cell, Group = "Low MSS")
]

# Highest 1000 MSS cells
high_cells <- ord[
  (.N - 999):.N,
  .(Cell, Group = "High MSS")
]

g <- rbind(
  low_cells,
  high_cells
)

# Keep only selected 2000 cells
d <- merge(
  d,
  g,
  by = "Cell"
)


# ------------------------------------------------------------
# 6. Calculate Gene Activity for each complex
# ------------------------------------------------------------

if (!is.na(xc) && !is.na(sc)) {

  # If complex and score columns already exist

  z <- d[
    ,
    .(
      Cell,
      Group,
      Complex = as.character(get(xc)),
      Score = as.numeric(get(sc))
    )
  ]

  z <- z[
    ,
    .(
      Score = mean(Score, na.rm = TRUE)
    ),
    by = .(Cell, Group, Complex)
  ]

} else {

  # Identify OXPHOS genes by gene name

  mp <- data.table(Gene = nm)

  mp[, Complex := NA_character_]


  # Complex I
  mp[
    grepl("^NDUF|^MT-ND", Gene, ignore.case = TRUE),
    Complex := "Complex I"
  ]


  # Complex II
  mp[
    grepl("^SDH", Gene, ignore.case = TRUE),
    Complex := "Complex II"
  ]


  # Complex III
  mp[
    grepl("^UQC|^UQCR|^CYC1|^MT-CYB",
          Gene, ignore.case = TRUE),
    Complex := "Complex III"
  ]


  # Complex IV
  mp[
    grepl("^COX|^MT-CO", Gene, ignore.case = TRUE),
    Complex := "Complex IV"
  ]


  # Complex V
  mp[
    grepl("^ATP5|^ATP6$|^ATP8$|^MT-ATP",
          Gene, ignore.case = TRUE),
    Complex := "Complex V"
  ]


  # Keep identified genes
  mp <- mp[!is.na(Complex)]


  # Print gene counts
  print(
    mp[
      ,
      .(Genes = .N),
      by = Complex
    ][order(Complex)]
  )


  if (uniqueN(mp$Complex) != 5) {
    stop("Unable to identify all five complexes")
  }


  # Save gene mapping
  fwrite(
    mp,
    "OXPHOS_nuclear_violin_gene_mapping.csv"
  )


  # Wide to long format
  v <- melt(
    d,
    id.vars = c("Cell", "Group"),
    measure.vars = mp$Gene,
    variable.name = "Gene",
    value.name = "Value"
  )


  # Add complex annotation
  v <- merge(
    v,
    mp,
    by = "Gene"
  )


  # Convert values to numeric
  v[, Value := as.numeric(Value)]


  # Mean Gene Activity per cell per complex
  z <- v[
    ,
    .(
      Score = mean(Value, na.rm = TRUE)
    ),
    by = .(Cell, Group, Complex)
  ]

}


# ------------------------------------------------------------
# 7. Prepare plotting data
# ------------------------------------------------------------

z <- z[is.finite(Score)]

z[
  ,
  Complex := factor(
    Complex,
    levels = c(
      "Complex I",
      "Complex II",
      "Complex III",
      "Complex IV",
      "Complex V"
    )
  )
]

z[
  ,
  Group := factor(
    Group,
    levels = c(
      "High MSS",
      "Low MSS"
    )
  )
]


# ------------------------------------------------------------
# 8. Calculate summary statistics
# ------------------------------------------------------------

sm <- z[
  ,
  .(
    N = .N,
    Mean = mean(Score),
    Median = median(Score)
  ),
  by = .(Complex, Group)
]


# ------------------------------------------------------------
# 9. Calculate median differences
# ------------------------------------------------------------

w <- dcast(
  sm,
  Complex ~ Group,
  value.var = "Median"
)

w[
  ,
  Delta_median :=
    get("High MSS") - get("Low MSS")
]


# ------------------------------------------------------------
# 10. Create panel labels
# ------------------------------------------------------------

w[
  ,
  Panel := paste0(
    as.character(Complex),
    "\nMedian difference ",
    sprintf("%+.3f", Delta_median)
  )
]

z <- merge(
  z,
  w[, .(Complex, Panel)],
  by = "Complex"
)

z[
  ,
  Panel := factor(
    Panel,
    levels = w$Panel
  )
]


# ------------------------------------------------------------
# 11. Save statistics
# ------------------------------------------------------------

fwrite(
  sm,
  "OXPHOS_nuclear_1000each_violin_enhanced_stats.csv"
)


# ------------------------------------------------------------
# 12. Plot violin
# ------------------------------------------------------------

p <- ggplot(
  z,
  aes(
    x = Group,
    y = Score,
    fill = Group
  )
) +

  # Violin
  geom_violin(
    width = 0.72,
    scale = "width",
    trim = TRUE,
    alpha = 0.88,
    color = "#424242",
    linewidth = 0.32
  ) +

  # Boxplot
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    color = "#202020",
    linewidth = 0.55
  ) +

  # Median point
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.5,
    fill = "#202020",
    color = "white",
    stroke = 0.45
  ) +

  # Five complexes
  facet_wrap(
    ~Panel,
    nrow = 1,
    scales = "free_y"
  ) +

  # MSS colors
  scale_fill_manual(
    values = c(
      "High MSS" = "#D98C84",
      "Low MSS" = "#8FA5C9"
    )
  ) +

  # Axis labels
  labs(
    x = NULL,
    y = "Nuclear OXPHOS gene activity"
  ) +

  # Main theme
  theme_classic(
    base_size = 13,
    base_family = "serif"
  ) +

  theme(

    # Facet titles
    strip.text = element_text(
      face = "bold",
      size = 11,
      lineheight = 1.1
    ),

    # Facet background
    strip.background = element_rect(
      fill = "#F7F7F7",
      color = "#404040",
      linewidth = 0.5
    ),

    # X-axis labels
    axis.text.x = element_text(
      face = "bold",
      size = 11
    ),

    # Y-axis title
    axis.title.y = element_text(
      size = 14
    ),

    # Panel borders
    panel.border = element_rect(
      color = "#404040",
      fill = NA,
      linewidth = 0.45
    ),

    # Panel spacing
    panel.spacing = grid::unit(
      0.9,
      "lines"
    ),

    # Hide legend
    legend.position = "none"
  )


# ------------------------------------------------------------
# 13. Save PDF
# ------------------------------------------------------------

ggsave(
  filename = "OXPHOS_nuclear_1000each_violin_enhanced.pdf",
  plot = p,
  width = 14,
  height = 4.8
)


# ------------------------------------------------------------
# 14. Save PNG
# ------------------------------------------------------------

if (requireNamespace("ragg", quietly = TRUE)) {

  ggsave(
    filename = "OXPHOS_nuclear_1000each_violin_enhanced.png",
    plot = p,
    width = 14,
    height = 4.8,
    dpi = 300,
    device = ragg::agg_png
  )

}


# ------------------------------------------------------------
# 15. Print results
# ------------------------------------------------------------

print(w)

cat("Finished nuclear enhanced violin\n")
