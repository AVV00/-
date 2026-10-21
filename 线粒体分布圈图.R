#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(circlize)
})

# ============================================================
# 1. Input / output
# ============================================================

workdir <- "/mnt/data1/linyanlin/scATAC/picture/pic1/circle"

infile <- file.path(
  workdir,
  "total-娣嬪反鐦ょ獊鍙樻眹鎬?xlsx"
)

outdir <- file.path(
  workdir,
  "output"
)

if (!dir.exists(outdir)) {
  dir.create(
    outdir,
    recursive = TRUE
  )
}

if (!file.exists(infile)) {
  stop(
    "Input file not found: ",
    infile
  )
}

cat("Input file:\n", infile, "\n\n")
cat("Output directory:\n", outdir, "\n\n")


# ============================================================
# 2. Six sheets
# ============================================================

sheet_names <- c(
  "T-soma-1%",
  "B-T",
  "T NK-T",
  "P-soma-1%",
  "B-P",
  "TNK-P"
)

real_sheets <- excel_sheets(infile)

missing_sheets <- setdiff(
  sheet_names,
  real_sheets
)

if (length(missing_sheets) > 0) {

  stop(
    "Missing sheets: ",
    paste(
      missing_sheets,
      collapse = ", "
    )
  )
}

cat("Sheets to plot:\n")
print(sheet_names)
cat("\n")


# ============================================================
# 3. Safe output name
# ============================================================

safe_name <- function(x) {

  y <- x

  y <- gsub(
    "%",
    "pct",
    y
  )

  y <- gsub(
    "[[:space:]]+",
    "_",
    y
  )

  y <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    y
  )

  y <- gsub(
    "_+",
    "_",
    y
  )

  y <- gsub(
    "^_|_$",
    "",
    y
  )

  return(y)
}


# ============================================================
# 4. Mutation classification
# ============================================================

class_mutation <- function(
  gene,
  region,
  variant_type,
  pos
) {

  gene2 <- toupper(
    trimws(
      as.character(gene)
    )
  )

  region2 <- tolower(
    trimws(
      as.character(region)
    )
  )

  type2 <- tolower(
    trimws(
      as.character(variant_type)
    )
  )


  gene2[
    gene2 %in% c(
      "",
      "NA",
      "NAN",
      "NULL"
    )
  ] <- NA_character_


  region2[
    region2 %in% c(
      "",
      "na",
      "nan",
      "null"
    )
  ] <- NA_character_


  type2[
    type2 %in% c(
      "",
      "na",
      "nan",
      "null"
    )
  ] <- NA_character_


  out <- rep(
    NA_character_,
    length(pos)
  )


  # ----------------------------------------------------------
  # Nonsynonymous
  # ----------------------------------------------------------

  idx <- (
    !is.na(type2) &
      grepl(
        "nonsynonymous|non_synonymous|non-synonymous|missense",
        type2
      )
  )

  out[idx] <- "Nonsynonymous"


  # ----------------------------------------------------------
  # Stopgain
  # ----------------------------------------------------------

  idx <- (
    is.na(out) &
      !is.na(type2) &
      grepl(
        "stopgain|stoploss|nonsense|stop_gained|stop_lost",
        type2
      )
  )

  out[idx] <- "Stopgain"


  # ----------------------------------------------------------
  # Synonymous
  # ----------------------------------------------------------

  idx <- (
    is.na(out) &
      !is.na(type2) &
      grepl(
        "synonymous",
        type2
      )
  )

  out[idx] <- "Synonymous"


  # ----------------------------------------------------------
  # rRNA
  # ----------------------------------------------------------

  idx <- (
    is.na(out) &
      (
        (
          !is.na(region2) &
            grepl(
              "rrna",
              region2
            )
        ) |
          (
            !is.na(gene2) &
              grepl(
                "^RNR1$|^RNR2$|RRNA",
                gene2
              )
          )
      )
  )

  out[idx] <- "rRNA"


  # ----------------------------------------------------------
  # tRNA
  # ----------------------------------------------------------

  idx <- (
    is.na(out) &
      (
        (
          !is.na(region2) &
            grepl(
              "trna",
              region2
            )
        ) |
          (
            !is.na(gene2) &
              grepl(
                "^TRN|^TRNA|TRNA|T-RNA",
                gene2
              )
          )
      )
  )

  out[idx] <- "tRNA"


  # ----------------------------------------------------------
  # D-loop by annotation
  # ----------------------------------------------------------

  idx <- (
    is.na(out) &
      (
        (
          !is.na(region2) &
            grepl(
              "d-loop|dloop|d_loop",
              region2
            )
        ) |
          (
            !is.na(gene2) &
              grepl(
                "D[_ -]?LOOP|CONTROL",
                gene2
              )
          )
      )
  )

  out[idx] <- "D-loop"


  # ----------------------------------------------------------
  # D-loop by mtDNA coordinate
  # ----------------------------------------------------------

  idx <- (
    is.na(out) &
      !is.na(pos) &
      (
        pos <= 576 |
          pos >= 16024
      )
  )

  out[idx] <- "D-loop"


  # ----------------------------------------------------------
  # Other
  # ----------------------------------------------------------

  out[
    is.na(out)
  ] <- "Other"

  return(out)
}


# ============================================================
# 5. Human mtDNA annotation
# ============================================================

mt_anno <- data.table(

  label = c(
    "D-loop",
    "tRNA-Phe",
    "RNR1",
    "tRNA-Val",
    "RNR2",
    "tRNA-Leu1",

    "ND1",
    "tRNA-Ile",
    "tRNA-Gln",
    "tRNA-Met",
    "ND2",

    "tRNA-Trp",
    "tRNA-Ala",
    "tRNA-Asn",
    "tRNA-Cys",
    "tRNA-Tyr",

    "COX1",
    "tRNA-Ser2",
    "tRNA-Asp",
    "COX2",
    "tRNA-Lys",

    "ATP8",
    "ATP6",
    "COX3",

    "tRNA-Gly",
    "ND3",
    "tRNA-Arg",

    "ND4L",
    "ND4",

    "tRNA-His",
    "tRNA-Ser1",
    "tRNA-Leu2",

    "ND5",
    "ND6",

    "tRNA-Glu",
    "CYTB",

    "tRNA-Thr",
    "tRNA-Pro",

    "D-loop"
  ),


  start = c(
    1,
    577,
    648,
    1602,
    1671,
    3230,

    3307,
    4263,
    4329,
    4402,
    4470,

    5512,
    5587,
    5657,
    5761,
    5826,

    5904,
    7446,
    7518,
    7586,
    8295,

    8366,
    8527,
    9207,

    9991,
    10059,
    10405,

    10470,
    10760,

    12138,
    12207,
    12266,

    12337,
    14149,

    14674,
    14747,

    15888,
    15956,

    16024
  ),


  end = c(
    576,
    647,
    1601,
    1670,
    3229,
    3304,

    4262,
    4331,
    4400,
    4469,
    5511,

    5579,
    5655,
    5729,
    5826,
    5891,

    7445,
    7514,
    7585,
    8269,
    8364,

    8572,
    9207,
    9990,

    10058,
    10404,
    10469,

    10766,
    12137,

    12206,
    12265,
    12336,

    14148,
    14673,

    14742,
    15887,

    15953,
    16023,

    16569
  )
)


# ============================================================
# 6. mtDNA annotation family
# ============================================================

anno_family <- function(x) {

  if (grepl("^ND", x)) {
    return("ND")
  }

  if (grepl("^COX", x)) {
    return("COX")
  }

  if (grepl("^ATP", x)) {
    return("ATP")
  }

  if (x == "CYTB") {
    return("CYTB")
  }

  if (grepl("^RNR", x)) {
    return("rRNA")
  }

  if (grepl("^tRNA", x)) {
    return("tRNA")
  }

  return("D-loop")
}


mt_anno[, family := vapply(
  label,
  anno_family,
  character(1)
)]


# ============================================================
# 7. Outer mtDNA annotation colors
# ============================================================

family_cols <- c(

  "ND" =
    "#F01C1C",

  "COX" =
    "#2738D6",

  "ATP" =
    "#E6D94C",

  "CYTB" =
    "#B931D6",

  "rRNA" =
    "#F0A52B",

  "tRNA" =
    "#43B5E8",

  "D-loop" =
    "#88D627"
)

mt_anno[, col := family_cols[family]]


# ============================================================
# 8. Mutation colors
#
# New requirement:
#
# Synonymous    = gray
# Nonsynonymous = pure red
# Stopgain      = dark red
# rRNA          = orange
# tRNA          = blue
# D-loop        = green
#
# ALL points:
# pch = 16
# cex = 0.60
# ============================================================

type_order <- c(
  "Synonymous",
  "Nonsynonymous",
  "Stopgain",
  "rRNA",
  "tRNA",
  "D-loop"
)


col_map <- c(

  "Synonymous" =
    "#9E9E9E",

  "Nonsynonymous" =
    "#FF0000",

  "Stopgain" =
    "#8B0000",

  "rRNA" =
    "#F0A52B",

  "tRNA" =
    "#43B5E8",

  "D-loop" =
    "#88D627"
)


# ============================================================
# All mutation point shapes are identical
# ============================================================

pch_map <- c(

  "Synonymous" =
    16,

  "Nonsynonymous" =
    16,

  "Stopgain" =
    16,

  "rRNA" =
    16,

  "tRNA" =
    16,

  "D-loop" =
    16
)


# ============================================================
# All mutation point sizes are identical: 0.60
# ============================================================

cex_map <- c(

  "Synonymous" =
    0.60,

  "Nonsynonymous" =
    0.60,

  "Stopgain" =
    0.60,

  "rRNA" =
    0.60,

  "tRNA" =
    0.60,

  "D-loop" =
    0.60
)


# ============================================================
# 9. Read one Excel sheet
# ============================================================

read_one_sheet <- function(sheet_name) {

  cat(
    "\n========================================\n"
  )

  cat(
    "Reading sheet: ",
    sheet_name,
    "\n",
    sep = ""
  )


  d <- read_excel(
    infile,
    sheet = sheet_name,
    .name_repair = "minimal"
  )

  d <- as.data.frame(d)


  # ----------------------------------------------------------
  # Only use first 20 columns.
  #
  # T-soma-1% and P-soma-1% contain extra summary tables
  # after column 20.
  # ----------------------------------------------------------

  if (ncol(d) < 20) {

    stop(
      "Sheet ",
      sheet_name,
      " has fewer than 20 columns."
    )
  }


  d <- d[
    ,
    1:20,
    drop = FALSE
  ]


  names(d) <- c(
    "sample",
    "patient",
    "disease_short",
    "disease_class",
    "position",
    "mutation_type",
    "count",
    "ref_base",
    "alt_base",
    "ref_context",
    "alt_context",
    "ref_L",
    "ref_H",
    "alt_L",
    "alt_H",
    "depth",
    "vaf",
    "gene",
    "region",
    "variant_type"
  )


  dt <- as.data.table(d)


  # ----------------------------------------------------------
  # Data type conversion
  # ----------------------------------------------------------

  dt[, sample :=
       as.character(sample)]

  dt[, patient :=
       as.character(patient)]

  dt[, position :=
       suppressWarnings(
         as.numeric(position)
       )]

  dt[, vaf :=
       suppressWarnings(
         as.numeric(vaf)
       )]

  dt[, gene :=
       as.character(gene)]

  dt[, region :=
       as.character(region)]

  dt[, variant_type :=
       as.character(variant_type)]


  # ----------------------------------------------------------
  # Valid mtDNA positions
  # ----------------------------------------------------------

  dt <- dt[
    !is.na(position) &
      position >= 1 &
      position <= 16569
  ]


  # ----------------------------------------------------------
  # VAF range
  # ----------------------------------------------------------

  dt[
    !is.na(vaf) &
      vaf < 0,
    vaf := 0
  ]

  dt[
    !is.na(vaf) &
      vaf > 1,
    vaf := 1
  ]


  # ----------------------------------------------------------
  # Missing VAF
  # ----------------------------------------------------------

  dt[
    is.na(vaf),
    vaf := 0.05
  ]


  # ----------------------------------------------------------
  # Mutation classification
  # ----------------------------------------------------------

  dt[, plot_type :=
       class_mutation(
         gene,
         region,
         variant_type,
         position
       )]


  # ----------------------------------------------------------
  # Radial position
  #
  # High VAF -> closer to outer ring
  # ----------------------------------------------------------

  set.seed(123)

  dt[, radius :=
       0.08 +
         0.84 * vaf +
         runif(
           .N,
           -0.015,
           0.015
         )]


  dt[
    radius < 0.03,
    radius := 0.03
  ]

  dt[
    radius > 0.97,
    radius := 0.97
  ]


  # ----------------------------------------------------------
  # Only six mutation classes are plotted
  # ----------------------------------------------------------

  dt_plot <- dt[
    plot_type %in%
      type_order
  ]


  cat(
    "Total valid mutation rows: ",
    nrow(dt),
    "\n",
    sep = ""
  )

  cat(
    "Plotted mutation rows: ",
    nrow(dt_plot),
    "\n",
    sep = ""
  )

  cat(
    "Other / gap rows not plotted: ",
    sum(
      dt$plot_type ==
        "Other"
    ),
    "\n\n",
    sep = ""
  )


  cat(
    "Mutation classes:\n"
  )


  print(

    dt_plot[
      ,
      .N,
      by =
        plot_type
    ][
      order(
        match(
          plot_type,
          type_order
        )
      )
    ]

  )


  return(dt_plot)
}


# ============================================================
# 10. Circle plotting function
# ============================================================

draw_circle <- function(sub) {

  circos.clear()


  circos.par(

    start.degree =
      90,

    gap.after =
      2,

    cell.padding =
      c(
        0,
        0,
        0,
        0
      ),

    track.margin =
      c(
        0.005,
        0.005
      ),

    points.overflow.warning =
      FALSE,

    canvas.xlim =
      c(
        -1.35,
        1.85
      ),

    canvas.ylim =
      c(
        -1.22,
        1.52
      )
  )


  circos.initialize(

    factors =
      "chrM",

    xlim =
      cbind(
        1,
        16569
      )
  )


  # ==========================================================
  # Outer mtDNA annotation ring
  # ==========================================================

  circos.trackPlotRegion(

    factors =
      "chrM",

    ylim =
      c(
        0,
        1
      ),

    track.height =
      0.13,

    bg.border =
      NA,

    panel.fun =
      function(x, y) {


        # ------------------------------------------------------
        # Colored mtDNA regions
        # ------------------------------------------------------

        for (
          i in seq_len(
            nrow(mt_anno)
          )
        ) {

          circos.rect(

            xleft =
              mt_anno$start[i],

            ybottom =
              0,

            xright =
              mt_anno$end[i],

            ytop =
              1,

            col =
              mt_anno$col[i],

            border =
              "white",

            lwd =
              0.7
          )
        }


        # ------------------------------------------------------
        # mtDNA gene labels
        # ------------------------------------------------------

        for (
          i in seq_len(
            nrow(mt_anno)
          )
        ) {

          mid <-
            (
              mt_anno$start[i] +
                mt_anno$end[i]
            ) / 2


          circos.text(

            x =
              mid,

            y =
              1.28,

            labels =
              mt_anno$label[i],

            facing =
              "clockwise",

            niceFacing =
              TRUE,

            adj =
              c(
                0,
                0.5
              ),

            cex =
              0.42
          )
        }
      }
  )


  # ==========================================================
  # Mutation ring
  # ==========================================================

  circos.trackPlotRegion(

    factors =
      "chrM",

    ylim =
      c(
        0,
        1
      ),

    track.height =
      0.52,

    bg.border =
      NA,

    panel.fun =
      function(x, y) {


        for (
          tp in type_order
        ) {


          dd <-
            sub[
              plot_type ==
                tp
            ]


          if (
            nrow(dd) ==
              0
          ) {

            next
          }


          circos.points(

            x =
              dd$position,

            y =
              dd$radius,

            pch =
              16,

            col =
              col_map[tp],

            cex =
              0.60
          )
        }
      }
  )


  # ==========================================================
  # Center text
  #
  # ONLY:
  #
  # xx mutations
  # ==========================================================

  text(

    x =
      0,

    y =
      0,

    labels =
      paste0(
        format(
          nrow(sub),
          big.mark = ","
        ),
        " mutations"
      ),

    cex =
      1.35,

    font =
      2
  )


  # ==========================================================
  # Legend
  # ==========================================================

  legend(

    x =
      0.82,

    y =
      1.40,

    bty =
      "n",

    cex =
      0.82,

    xpd =
      NA,

    y.intersp =
      1.10,

    x.intersp =
      0.60,

    pt.cex =
      rep(
        1.0,
        6
      ),

    legend =
      c(

        paste0(
          "Synonymous (n = ",
          format(
            sum(
              sub$plot_type ==
                "Synonymous"
            ),
            big.mark = ","
          ),
          ")"
        ),

        paste0(
          "Nonsynonymous (n = ",
          format(
            sum(
              sub$plot_type ==
                "Nonsynonymous"
            ),
            big.mark = ","
          ),
          ")"
        ),

        paste0(
          "Stopgain (n = ",
          format(
            sum(
              sub$plot_type ==
                "Stopgain"
            ),
            big.mark = ","
          ),
          ")"
        ),

        paste0(
          "rRNA (n = ",
          format(
            sum(
              sub$plot_type ==
                "rRNA"
            ),
            big.mark = ","
          ),
          ")"
        ),

        paste0(
          "tRNA (n = ",
          format(
            sum(
              sub$plot_type ==
                "tRNA"
            ),
            big.mark = ","
          ),
          ")"
        ),

        paste0(
          "D-loop (n = ",
          format(
            sum(
              sub$plot_type ==
                "D-loop"
            ),
            big.mark = ","
          ),
          ")"
        )
      ),

    pch =
      rep(
        16,
        6
      ),

    col =
      col_map[
        type_order
      ]
  )


  circos.clear()
}


# ============================================================
# 11. Generate six PDF files
# ============================================================

summary_list <- list()


for (
  sheet_name in sheet_names
) {


  sub <-
    read_one_sheet(
      sheet_name
    )


  if (
    nrow(sub) ==
      0
  ) {

    warning(
      "No mutations to plot in sheet: ",
      sheet_name
    )

    next
  }


  file_tag <-
    safe_name(
      sheet_name
    )


  out_pdf <-
    file.path(
      outdir,
      paste0(
        "mtDNA_circleplot_",
        file_tag,
        ".pdf"
      )
    )


  cat(
    "\nDrawing:\n",
    out_pdf,
    "\n"
  )


  pdf(

    out_pdf,

    width =
      9.5,

    height =
      7.4,

    useDingbats =
      FALSE
  )


  par(

    mar =
      c(
        0.5,
        0.5,
        1.0,
        0.5
      ),

    xpd =
      NA
  )


  draw_circle(
    sub
  )


  dev.off()


  summary_list[[sheet_name]] <-
    data.table(

      sheet =
        sheet_name,

      mutations_plotted =
        nrow(sub),

      synonymous =
        sum(
          sub$plot_type ==
            "Synonymous"
        ),

      nonsynonymous =
        sum(
          sub$plot_type ==
            "Nonsynonymous"
        ),

      stopgain =
        sum(
          sub$plot_type ==
            "Stopgain"
        ),

      rRNA =
        sum(
          sub$plot_type ==
            "rRNA"
        ),

      tRNA =
        sum(
          sub$plot_type ==
            "tRNA"
        ),

      D_loop =
        sum(
          sub$plot_type ==
            "D-loop"
        ),

      output_pdf =
        out_pdf
    )


  cat(
    "SUCCESS: ",
    out_pdf,
    "\n",
    sep = ""
  )
}


# ============================================================
# 12. Summary CSV
# ============================================================

summary_dt <-
  rbindlist(
    summary_list,
    fill = TRUE
  )


summary_file <-
  file.path(
    outdir,
    "mtDNA_circleplot_6sheets_summary.csv"
  )


fwrite(
  summary_dt,
  summary_file
)


cat(
  "\n\n========================================\n"
)

cat(
  "ALL DONE\n"
)

cat(
  "========================================\n\n"
)


print(
  summary_dt
)


cat(
  "\nOutput directory:\n",
  outdir,
  "\n"
)

cat(
  "\nSummary:\n",
  summary_file,
  "\n"
)
