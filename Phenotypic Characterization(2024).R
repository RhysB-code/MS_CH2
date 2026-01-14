# ================================================================
# GWAS phenotypes: rep-aware summaries, polished histograms, CSV (0–1,1–2,…)
# ================================================================

## 0) Setup ------------------------------------------------------
rm(list = ls())
setwd("C:/Users/RhysB/OneDrive/Desktop/Chill Spreadsheets")

packages <- c("readxl","dplyr","tidyr","stringr","ggplot2",
              "forcats","cowplot","tibble","purrr")
to_install <- packages[!packages %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(packages, library, character.only = TRUE))

## 1) Parameters -------------------------------------------------
xlsx_path         <- "GWAS_Data_Budbreak.xlsx"
date_levels       <- c("20-Dec","3-Jan","1-Jan","15-Jan")
expected_reps     <- 4
require_full_reps <- FALSE   # set TRUE to require all 4 reps for inclusion

## 2) Helpers ----------------------------------------------------
clean_names <- function(df){
  nm <- names(df)
  nm <- nm |>
    stringr::str_replace_all("#", "num ") |>
    stringr::str_replace_all("\\s+", " ") |>
    stringr::str_trim()
  names(df) <- nm
  df
}

extract_genotype <- function(x) sub("^(.*?)_.*$", "\\1", x)
extract_rep      <- function(x) ifelse(grepl("^.*?_(.*)$", x), sub("^.*?_(.*)$", "\\1", x), NA_character_)

find_measure_cols <- function(df, dates, term_regex){
  sapply(dates, function(d){
    idx <- which(stringr::str_detect(
      names(df),
      regex(paste0(term_regex, ".*\\b", d, "\\b"), ignore_case = TRUE)
    ))
    if (length(idx) == 0) NA_character_ else names(df)[idx[1]]
  }, USE.NAMES = TRUE)
}

prep_long <- function(df, dates_vec){
  stopifnot("Lug Order" %in% names(df))
  df_core <- df %>%
    dplyr::rename(LugOrder = `Lug Order`) %>%
    dplyr::mutate(
      Genotype  = extract_genotype(LugOrder),
      Replicate = extract_rep(LugOrder)
    )
  
  # Budbreak columns are exactly the date labels
  bb_cols <- intersect(dates_vec, names(df))
  if (length(bb_cols) == 0) stop("No budbreak columns for dates: ",
                                 paste(dates_vec, collapse=", "))
  
  floral_cols <- find_measure_cols(df, dates_vec, term_regex = "(^|\\s)floral")
  open_cols   <- find_measure_cols(df, dates_vec, term_regex = "open\\s*flowers?|\\bopen\\b")
  
  long <- df_core %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(bb_cols),
      names_to = "Date",
      values_to = "Budbreak",
      values_transform = list(Budbreak = as.numeric)
    )
  
  if (any(!is.na(floral_cols))) {
    keep <- c("LugOrder","Genotype","Replicate", na.omit(floral_cols))
    floral_long <- df_core[, keep, drop = FALSE] %>%
      tidyr::pivot_longer(
        cols = -c(LugOrder, Genotype, Replicate),
        names_to = "Date",
        values_to = "Floral",
        values_transform = list(Floral = as.numeric)
      ) %>%
      dplyr::mutate(Date = stringr::str_extract(Date, "(20-Dec|3-Jan|1-Jan|15-Jan)"))
    long <- dplyr::left_join(long, floral_long,
                             by = c("LugOrder","Genotype","Replicate","Date"))
  } else long$Floral <- NA_real_
  
  if (any(!is.na(open_cols))) {
    keep2 <- c("LugOrder","Genotype","Replicate", na.omit(open_cols))
    open_long <- df_core[, keep2, drop = FALSE] %>%
      tidyr::pivot_longer(
        cols = -c(LugOrder, Genotype, Replicate),
        names_to = "Date",
        values_to = "OpenFlowers",
        values_transform = list(OpenFlowers = as.numeric)
      ) %>%
      dplyr::mutate(Date = stringr::str_extract(Date, "(20-Dec|3-Jan|1-Jan|15-Jan)"))
    long <- dplyr::left_join(long, open_long,
                             by = c("LugOrder","Genotype","Replicate","Date"))
  } else long$OpenFlowers <- NA_real_
  
  long
}

## 3) Read sheets & build long table -----------------------------
s1 <- readxl::read_excel(xlsx_path, sheet = "GWAS1") |> clean_names()
s2 <- readxl::read_excel(xlsx_path, sheet = "GWAS2") |> clean_names()

dates_s1 <- c("20-Dec","3-Jan")
dates_s2 <- c("1-Jan","15-Jan")

long1 <- prep_long(s1, dates_s1)
long2 <- prep_long(s2, dates_s2)

dat_long <- dplyr::bind_rows(long1, long2) %>%
  dplyr::mutate(Date = factor(Date, levels = date_levels))

## 4) Aggregate to GENOTYPE × DATE (rep-aware) -------------------
geno_date <- dat_long %>%
  dplyr::group_by(Genotype, Date) %>%
  dplyr::summarise(
    n_reps            = sum(!is.na(Budbreak)),
    budbreak_mean     = mean(Budbreak, na.rm = TRUE),
    floral_mean       = mean(Floral, na.rm = TRUE),
    openflowers_mean  = mean(OpenFlowers, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    expected_reps = !!expected_reps,
    full_reps     = n_reps >= expected_reps
  )

if (require_full_reps) {
  message("Keeping only rows with full reps: n_reps >= ", expected_reps)
  geno_date <- dplyr::filter(geno_date, full_reps)
}

## 5) Polished histograms (0–1, 1–2, …) + counts on bars --------
# One panel: fixed integer bins; add labels; auto headroom for text
one_hist <- function(df, value_col, panel_title, xlab, xlim_max = NULL) {
  d <- df %>% dplyr::filter(!is.na(.data[[value_col]]))
  xmax <- if (is.null(xlim_max)) ceiling(max(d[[value_col]], na.rm = TRUE)) else xlim_max
  xmax <- max(1, xmax)                  # at least one bin
  breaks <- 0:xmax                      # [0,1), [1,2), ...
  
  # Bin once for headroom calculation
  bins <- cut(d[[value_col]], breaks = breaks, right = FALSE, include.lowest = TRUE)
  ymax <- if (length(bins)) max(table(bins)) else 0
  ypad <- if (ymax > 0) ymax * 1.18 else 1  # headroom for text labels
  
  ggplot(d, aes(x = .data[[value_col]])) +
    geom_histogram(breaks = breaks, color = "black", fill = "steelblue") +
    stat_bin(
      breaks = breaks,
      geom = "text",
      aes(label = after_stat(count)),
      vjust = -0.35, size = 5
    ) +
    labs(title = panel_title, x = xlab, y = "Count of genotypes") +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ypad)) +
    scale_x_continuous(breaks = 0:xmax) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text    = element_text(size = 11),
      plot.margin  = margin(t = 10, r = 12, b = 18, l = 12)
    )
}

# 2×2 grid with a separate title row (keeps the title far from panels)
make_hist_grid <- function(df, trait_col, grand_title, xlab) {
  xlim_max <- max(1, ceiling(max(df[[trait_col]], na.rm = TRUE)))
  
  d20 <- df %>% dplyr::filter(Date == "20-Dec")
  d3  <- df %>% dplyr::filter(Date == "3-Jan")
  d1  <- df %>% dplyr::filter(Date == "1-Jan")
  d15 <- df %>% dplyr::filter(Date == "15-Jan")
  
  p20 <- one_hist(d20, trait_col, "20-Dec", xlab, xlim_max)
  p3  <- one_hist(d3,  trait_col, "3-Jan",  xlab, xlim_max)
  p1  <- one_hist(d1,  trait_col, "1-Jan",  xlab, xlim_max)
  p15 <- one_hist(d15, trait_col, "15-Jan", xlab, xlim_max)
  
  # Chill-hour side labels
  row1_lab <- ggdraw() + draw_label("200 hrs", angle = 90, fontface = "bold", size = 13, vjust = 0.5)
  row2_lab <- ggdraw() + draw_label("263 hrs", angle = 90, fontface = "bold", size = 13, vjust = 0.5)
  
  row1 <- plot_grid(p20, p3, ncol = 2, rel_widths = c(1,1))
  row2 <- plot_grid(p1,  p15, ncol = 2, rel_widths = c(1,1))
  
  grid1 <- plot_grid(row1_lab, row1, ncol = 2, rel_widths = c(0.08, 1))
  grid2 <- plot_grid(row2_lab, row2, ncol = 2, rel_widths = c(0.08, 1))
  panels <- plot_grid(grid1, grid2, ncol = 1, rel_heights = c(1,1))
  
  # Title as its own row with generous padding
  title_grob <- ggdraw() +
    draw_label(grand_title, x = 0.5, y = 0.6, fontface = "bold", size = 18, hjust = 0.5, vjust = 0.5)
  
  plot_grid(title_grob, panels, ncol = 1, rel_heights = c(0.12, 1))  # bump first number to move title higher
}

## 6) Show plots -------------------------------------------------
p_bud <- make_hist_grid(
  geno_date, "budbreak_mean",
  "Budbreak per Genotype (mean across reps) — GWAS Intervals",
  "Mean buds broken per genotype"
)
print(p_bud)

# Optional extras:
# print(make_hist_grid(geno_date, "floral_mean", "Floral laterals per genotype", "Mean floral laterals per genotype"))
# print(make_hist_grid(geno_date, "openflowers_mean", "Open flowers per genotype", "Mean open flowers per genotype"))

## 7) CSV with ranges 0–1, 1–2, 2–3, … (zeros included in 0–1) ---


# ----- CSV with Excel-safe bin labels (0_to_1, 1_to_2, ...) -----
make_intbin_csv <- function(df, value_col, file_name) {
  v <- df[[value_col]]
  xmax <- max(1, suppressWarnings(ceiling(max(v, na.rm = TRUE))))  # at least 1
  
  breaks <- 0:xmax
  lefts  <- 0:(xmax - 1)
  rights <- 1:xmax
  
  # Choose ONE of these two label styles:
  labels <- sprintf("%d_to_%d", lefts, rights)        # safe: 0_to_1, 1_to_2, ...
  # labels <- sprintf("[%d,%d)", lefts, rights)       # also safe: [0,1), [1,2), ...
  
  out <- df %>%
    dplyr::filter(!is.na(.data[[value_col]])) %>%
    dplyr::mutate(
      Range = cut(
        .data[[value_col]],
        breaks = breaks,
        right  = FALSE,           # [a,b)
        include.lowest = TRUE,
        labels = labels
      )
    ) %>%
    dplyr::group_by(Date, Range) %>%
    dplyr::summarise(Genotypes = paste(sort(Genotype), collapse = ", "),
                     .groups = "drop") %>%
    # sort bins by left edge so they appear in numeric order
    dplyr::mutate(LeftEdge = as.numeric(gsub("^\\D*(\\d+).*", "\\1", Range))) %>%
    dplyr::arrange(factor(Date, levels = c("20-Dec","3-Jan","1-Jan","15-Jan")), LeftEdge) %>%
    dplyr::select(Date, Range, Genotypes)
  
  write.csv(out, file_name, row.names = FALSE, quote = TRUE)
  message("Wrote CSV: ", normalizePath(file_name))
}


# Write CSV(s)
make_intbin_csv(geno_date, "budbreak_mean",      "Budbreak_GenotypeBins_integer.csv")
# Optional:
# make_intbin_csv(geno_date, "floral_mean",        "Floral_GenotypeBins_integer.csv")
# make_intbin_csv(geno_date, "openflowers_mean",   "OpenFlowers_GenotypeBins_integer.csv")

