# =============================================================
# Blackberry budbreak GWAS (tetraploid) — robust end-to-end
#  - Traits: Binary (>0), Percent, BLUPs, BLUEs, PCA (PC1/PC2)
#  - Dates used: Jan 3, Jan 15 (ignores Dec 20 & Jan 1)
#  - Cohorts: ALL / PF-only (APF*) / FF-only (A-*)
#  - PF covariate: ONLY for ALL cohort (mixed PF/FF); skipped for homogeneous PF/FF cohorts
#  - Models: additive, 1-dom, general  (LOCO kinship)
#  - END: All Manhattan plots printed to Plots pane, labeled
# =============================================================

rm(list = ls()); options(stringsAsFactors = FALSE)

# ---- Packages ----
pkgs <- c("tidyverse","data.table","GWASpoly","stringr",
          "lme4","lmerTest","emmeans","parallel","cowplot")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Paths (your exact CSVs) ----
pheno_file <- "C:/Users/RhysB/OneDrive/Desktop/Chill Spreadsheets/GWAS_Data_Budbreak.csv"
geno_file  <- "C:/Users/RhysB/OneDrive/Desktop/Chill Spreadsheets/captureseq2025_GATK_postmean_70K_updog_filtered.csv"
out_dir    <- "gwaspoly_runsets"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("\n[PATH] PHENO:", pheno_file, "\n[PATH] GENO :", geno_file, "\n")
stopifnot(file.exists(pheno_file), file.exists(geno_file))

# ---- Plot registry (we'll print all Manhattans at the end) ----
RUN_REGISTRY <- list()

# =============================================================
# 0) Helpers
# =============================================================
.read_csv <- function(path, nrows = 0){
  data.table::fread(path, nrows = ifelse(nrows>0, nrows, -1),
                    data.table = TRUE, na.strings = c("", "NA"))
}
.normalize_names <- function(nm){ gsub("\\s+", " ", trimws(nm)) }
.strip_suffix <- function(nm) sub("\\.\\d+$", "", nm)
.find_by_base <- function(cands, names_vec){
  base <- .strip_suffix(names_vec)
  idx  <- match(TRUE, base %in% cands)
  if (is.na(idx)) NA_character_ else names_vec[idx]
}
.clamp_counts <- function(v, maxv){
  v <- suppressWarnings(as.numeric(v)); v <- pmax(0, v)
  if (!missing(maxv)) v <- pmin(v, maxv)
  v
}
get_name <- function(x) sub("(_.*)$","", x)   # used for *phenotype* IDs
get_rep  <- function(x) ifelse(grepl("_", x), sub("^[^_]+_","", x), NA_character_)

# =============================================================
# 1) Phenotypes: build Binary (>0), Percent
# =============================================================
ph_raw <- .read_csv(pheno_file)
names(ph_raw) <- make.unique(.normalize_names(names(ph_raw)), sep = ".")

name_col <- .find_by_base(c("Lug Order","LugOrder","Name","ID"), names(ph_raw))
if (is.na(name_col)) name_col <- names(ph_raw)[1]

base_names <- .strip_suffix(names(ph_raw))
tb_idx     <- base_names %in% c("Total Buds","TBuds","Buds","Total_Buds","Total Buds.1")
tbuds_cols <- names(ph_raw)[tb_idx]
if (!length(tbuds_cols)) stop("Couldn't find a 'Total Buds' family column in phenotype CSV.")

col_3jan  <- .find_by_base(c("3-Jan","3 Jan","Jan 3","3Jan"),     names(ph_raw))
col_15jan <- .find_by_base(c("15-Jan","15 Jan","Jan 15","15Jan"), names(ph_raw))
if (is.na(col_3jan) || is.na(col_15jan)) {
  stop("Missing date columns for 3-Jan and/or 15-Jan. Headers (first 30): ",
       paste(head(names(ph_raw), 30), collapse=" | "))
}

ph <- ph_raw %>%
  mutate(
    NAME  = get_name(.data[[name_col]]),
    REP   = get_rep(.data[[name_col]]),
    TBUDS = dplyr::coalesce(!!!rlang::syms(tbuds_cols)),
    VAL_3Jan  = .clamp_counts(.data[[col_3jan]],  TBUDS),
    VAL_15Jan = .clamp_counts(.data[[col_15jan]], TBUDS),
    BIN_3Jan  = as.integer(VAL_3Jan  > 0),    # responders = >0 buds broken
    BIN_15Jan = as.integer(VAL_15Jan > 0),
    PCT_3Jan  = ifelse(TBUDS > 0, 100*VAL_3Jan  / TBUDS, NA_real_),
    PCT_15Jan = ifelse(TBUDS > 0, 100*VAL_15Jan / TBUDS, NA_real_)
  )
ph$NAME <- as.factor(ph$NAME)

cat("\n[PHENO] Rows:", nrow(ph), "  Unique genotypes:", length(unique(ph$NAME)), "\n")
cat("[PHENO] Non-NA counts — 3-Jan:", sum(!is.na(ph$VAL_3Jan)),
    "  15-Jan:", sum(!is.na(ph$VAL_15Jan)), "\n")

# =============================================================
# 1A) PF/FF classification (with your overrides)
# =============================================================
is_pf_default <- function(nm) startsWith(nm, "APF")
is_ff_default <- function(nm) startsWith(nm, "A-")

pf_overrides_add <- c("PrimeArkFreedom","PrimeArkHorizon","APF-661TN")
ff_overrides_add <- c("Apache","Caddo","Immaculate","Lochness","Natchez",
                      "Navaho","Osage","Ouachita","Ponca","Stella",
                      "Superlicious","Tupy","Von")
pf_exclude <- c(); ff_exclude <- c()

all_names <- sort(unique(as.character(ph$NAME)))
PF_tbl <- tibble(NAME = all_names) %>%
  mutate(PF = dplyr::case_when(
    NAME %in% pf_exclude ~ NA_integer_,
    NAME %in% ff_exclude ~ NA_integer_,
    NAME %in% pf_overrides_add ~ 1L,
    NAME %in% ff_overrides_add ~ 0L,
    is_pf_default(NAME) ~ 1L,
    is_ff_default(NAME) ~ 0L,
    TRUE ~ NA_integer_
  ))

cat("\n[PF/FF] Assignment summary:\n")
PF_tbl %>%
  mutate(Group = dplyr::case_when(PF==1L ~ "PF", PF==0L ~ "FF", TRUE ~ "Unknown")) %>%
  count(Group) %>% print()

# =============================================================
# 2) BLUPs & BLUEs (robust)
# =============================================================
mk_blup_blue <- function(df, resp){
  d <- df %>% dplyr::filter(!is.na(.data[[resp]]), !is.na(NAME))
  d$NAME <- factor(d$NAME)
  
  nm_blup <- paste0("BLUP_", resp)
  nm_blue <- paste0("BLUE_", resp)
  
  if (dplyr::n_distinct(d$NAME) < 2 || nrow(d) < 5) {
    message("[BLUP/BLUE] Skipping ", resp, ": not enough data after NA filtering.")
    return(list(
      BLUP = tibble::tibble(NAME = character(), !!nm_blup := numeric()),
      BLUE = tibble::tibble(NAME = character(), !!nm_blue := numeric())
    ))
  }
  
  blups_tbl <- try({
    m <- lme4::lmer(stats::as.formula(paste0(resp, " ~ 1 + (1|NAME)")), data = d, REML = TRUE)
    re <- lme4::ranef(m, condVar = FALSE)$NAME
    tibble::tibble(NAME = rownames(re), !!nm_blup := as.numeric(re[["(Intercept)"]]))
  }, silent = TRUE)
  
  if (inherits(blups_tbl, "try-error")) {
    message("[BLUP] Falling back to centered genotype means for ", resp)
    mu <- mean(d[[resp]], na.rm = TRUE)
    blups_tbl <- d %>%
      dplyr::group_by(NAME) %>%
      dplyr::summarise(!!nm_blup := mean(.data[[resp]], na.rm = TRUE) - mu, .groups = "drop")
  }
  
  blues_tbl <- try({
    m <- stats::lm(stats::as.formula(paste0(resp, " ~ 0 + NAME")), data = d)
    est <- stats::coef(m)
    tibble::tibble(NAME = names(est), !!nm_blue := as.numeric(est))
  }, silent = TRUE)
  
  if (inherits(blues_tbl, "try-error")) {
    message("[BLUE] Falling back to simple genotype means for ", resp)
    blues_tbl <- d %>%
      dplyr::group_by(NAME) %>%
      dplyr::summarise(!!nm_blue := mean(.data[[resp]], na.rm = TRUE), .groups = "drop")
  }
  
  list(BLUP = blups_tbl, BLUE = blues_tbl)
}

bb3  <- mk_blup_blue(ph %>% dplyr::select(NAME, REP, PCT_3Jan)  %>% tidyr::drop_na(NAME), "PCT_3Jan")
bb15 <- mk_blup_blue(ph %>% dplyr::select(NAME, REP, PCT_15Jan) %>% tidyr::drop_na(NAME), "PCT_15Jan")

BLUP_tbl <- dplyr::full_join(bb3$BLUP,  bb15$BLUP,  by = "NAME")
BLUE_tbl <- dplyr::full_join(bb3$BLUE,  bb15$BLUE,  by = "NAME")

MEAN_tbl <- ph %>%
  group_by(NAME) %>%
  summarise(
    BIN_3Jan  = as.integer(any(BIN_3Jan  == 1, na.rm = TRUE)),
    BIN_15Jan = as.integer(any(BIN_15Jan == 1, na.rm = TRUE)),
    PCT_3Jan  = mean(PCT_3Jan,  na.rm = TRUE),
    PCT_15Jan = mean(PCT_15Jan, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================
# 3) PCA on BLUPs — safe
# =============================================================
pca_src <- BLUP_tbl %>% select(NAME, BLUP_PCT_3Jan, BLUP_PCT_15Jan)
PC_tbl  <- tibble(NAME = pca_src$NAME, PC1 = NA_real_, PC2 = NA_real_)
pc_ok   <- pca_src %>% filter(complete.cases(.))
if (nrow(pc_ok) >= 5) {
  X <- as.matrix(pc_ok[, -1, drop = FALSE])
  sds <- apply(X, 2, sd, na.rm = TRUE)
  keep_cols <- which(sds > 0 & is.finite(sds))
  if (length(keep_cols) >= 2) {
    X2 <- scale(X[, keep_cols, drop = FALSE])
    P  <- prcomp(X2, center = FALSE, scale. = FALSE)
    sc <- as.data.frame(P$x[, 1:2, drop = FALSE]); sc$NAME <- pc_ok$NAME
    PC_tbl <- PC_tbl %>% select(-PC1,-PC2) %>%
      left_join(sc %>% transmute(NAME, PC1, PC2 = ifelse(ncol(P$x) >= 2, PC2, NA_real_)), by = "NAME")
  } else if (length(keep_cols) == 1) {
    z <- as.numeric(scale(X[, keep_cols, drop = TRUE]))
    sc <- tibble(NAME = pc_ok$NAME, PC1 = z, PC2 = NA_real_)
    PC_tbl <- PC_tbl %>% select(-PC1,-PC2) %>% left_join(sc, by = "NAME")
    message("[PCA] Only one informative BLUP column; using its z-score as PC1 (PC2=NA).")
  } else {
    message("[PCA] No informative BLUP columns; skipping PC scores.")
  }
} else {
  message("[PCA] Not enough complete rows for PCA (need ≥5).")
}

# =============================================================
# 4) Master phenotype panel
# =============================================================
PH_MASTER <- MEAN_tbl %>%
  full_join(BLUP_tbl, by = "NAME") %>%
  full_join(BLUE_tbl, by = "NAME") %>%
  full_join(PC_tbl,  by = "NAME") %>%
  left_join(PF_tbl,  by = "NAME")

# =============================================================
# 5) Genotype matrix — harmonize sample names to phenotype NAMEs
#    (Main fix: replace "_" -> "-" in genotype headers)
# =============================================================
message("\nReading genotype matrix & harmonizing sample names...")
G <- .read_csv(geno_file)

# Identify structural columns (case-insensitive)
gn <- names(G)
idx_marker <- which(tolower(gn) %in% c("marker","markername","marker_id","markerid","id"))[1]
idx_chrom  <- which(tolower(gn) %in% c("chrom","chromosome","chr"))[1]
idx_pos    <- which(tolower(gn) %in% c("pos","position","bp","bp_position"))[1]
if (is.na(idx_marker) || is.na(idx_chrom) || is.na(idx_pos))
  stop("Genotype file must have marker/chrom/pos (any casing). Headers seen: ",
       paste(head(names(G), 30), collapse=" | "))

names(G)[idx_marker] <- "marker"
names(G)[idx_chrom]  <- "chrom"
names(G)[idx_pos]    <- "pos"

geno_samples <- setdiff(names(G), c("marker","chrom","pos"))
pheno_names  <- as.character(PH_MASTER$NAME)

# 1) Exact overlap
exact_overlap <- intersect(geno_samples, pheno_names)

# 2) Underscore→hyphen normalization for genotype headers
geno_norm <- tibble(
  orig = geno_samples,
  hyph = gsub("_","-", geno_samples, fixed = TRUE)
)

hyph_overlap <- intersect(geno_norm$hyph, pheno_names)

cat("[GENO] Exact overlap:", length(exact_overlap),
    " | Overlap after '_'→'-':", length(hyph_overlap), "\n")

# Build rename map for genotype columns where hyph matches a phenotype NAME
need_rename <- geno_norm %>%
  filter(hyph %in% pheno_names & orig != hyph)

# Avoid collisions: if a hyph name already exists as a column, drop that rename
need_rename <- need_rename %>%
  filter(!(hyph %in% names(G)))

# If multiple orig map to same hyph, keep the first
if (any(duplicated(need_rename$hyph))) {
  dups <- unique(need_rename$hyph[duplicated(need_rename$hyph)])
  message("[GENO] Multiple genotype columns map to the same NAME after '_'→'-'. Keeping first for: ",
          paste(dups, collapse=", "))
  need_rename <- need_rename %>% group_by(hyph) %>% slice(1) %>% ungroup()
}

# Apply renames
if (nrow(need_rename) > 0) {
  for (i in seq_len(nrow(need_rename))) {
    data.table::setnames(G, old = need_rename$orig[i], new = need_rename$hyph[i])
  }
}

# Final overlap
final_samples <- setdiff(names(G), c("marker","chrom","pos"))
samp_cols <- intersect(final_samples, pheno_names)
cat("[GENO] Final sample overlap with phenos:", length(samp_cols), "\n")

if (length(samp_cols) < 20) {
  miss_from_geno  <- setdiff(pheno_names, final_samples)
  miss_from_pheno <- setdiff(final_samples, pheno_names)
  message("[GENO] Low overlap; examples (first 20 each):\n  - In PHENO not in GENO: ",
          paste(head(miss_from_geno, 20), collapse=", "),
          "\n  - In GENO not in PHENO: ",
          paste(head(miss_from_pheno, 20), collapse=", "))
}

# Write subset for GWASpoly
keep_cols <- c("marker","chrom","pos", samp_cols)
G_sub <- as.data.table(G)[, ..keep_cols]
geno_out <- file.path(out_dir, "geno_numeric_tetra.csv")
data.table::fwrite(G_sub, geno_out)

# =============================================================
# 6) Run grid: cohorts × covariate × trait-groups
# =============================================================
cohorts <- list(
  ALL = PH_MASTER,                       # mixed PF/FF
  FF  = PH_MASTER %>% filter(PF == 0L)   # FF-only
)

cov_modes   <- c("noPFcov","withPFcov")
trait_lists <- list(
  BINARY = c("BIN_3Jan","BIN_15Jan"),
  PCT    = c("PCT_3Jan","PCT_15Jan"),
  BLUP   = c("BLUP_PCT_3Jan","BLUP_PCT_15Jan"),
  BLUE   = c("BLUE_PCT_3Jan","BLUE_PCT_15Jan"),
  PCA    = c("PC1","PC2")
)
models_use <- c("additive","1-dom","general")

# =============================================================
# =============================================================
# 7) Core runner
#    - clamps geno.freq into (0,1) to prevent errors when N is tiny
# =============================================================
run_one <- function(df, trait_vec, cov_mode, cohort_name, tag){
  df_use <- df %>% dplyr::select(NAME, dplyr::all_of(trait_vec), PF)
  
  keep <- trait_vec[sapply(trait_vec, function(t) {
    v <- df_use[[t]]; any(!is.na(v)) && (length(unique(na.omit(v))) > 1)
  })]
  if (!length(keep)) {
    message("[SKIP] No informative traits for ", tag, " | ", cohort_name, " | ", cov_mode)
    return(invisible(NULL))
  }
  
  pheno_out <- file.path(out_dir, paste0("pheno_", tag, "_", cohort_name, "_", cov_mode, ".csv"))
  if (cov_mode == "withPFcov") {
    if (!("PF" %in% names(df_use)) || length(unique(na.omit(df_use$PF))) < 2) {
      message("[INFO] PF covariate constant/absent in ", cohort_name, "; running WITHOUT PF covariate.")
      data.table::fwrite(df_use %>% dplyr::select(NAME, dplyr::all_of(keep)), pheno_out)
      fixed_formula <- NULL
    } else {
      data.table::fwrite(df_use %>% dplyr::select(NAME, dplyr::all_of(keep), PF), pheno_out)
      fixed_formula <- "PF"
    }
  } else {
    data.table::fwrite(df_use %>% dplyr::select(NAME, dplyr::all_of(keep)), pheno_out)
    fixed_formula <- NULL
  }
  
  data_gp <- read.GWASpoly(ploidy=4, pheno.file=pheno_out, geno.file=geno_out,
                           format="numeric", n.traits=length(keep), delim=",")
  data_gp <- set.K(data_gp, LOCO = TRUE, n.core = max(1, parallel::detectCores()-1))
  
  N <- nrow(data_gp@pheno)
  gf <- 1 - 5/N
  if (!is.finite(gf) || gf <= 0 || gf >= 1) gf <- 0.95  # safe fallback
  gf <- max(0.01, min(0.99, gf))
  
  params <- set.params(
    geno.freq  = gf,
    fixed      = if (is.null(fixed_formula)) NULL else fixed_formula,
    fixed.type = if (is.null(fixed_formula)) NULL else "factor"
  )
  
  fit <- GWASpoly(data=data_gp, models=models_use, traits=keep, params=params,
                  n.core=max(1, parallel::detectCores()-1), quiet=TRUE)
  fit_thr <- set.threshold(fit, method = "M.eff", level = 0.05)
  
  # Save tables (scores + additive QTL)
  for (tr in keep) {
    base <- file.path(out_dir, paste0("GWAS_", tag, "_", cohort_name, "_", cov_mode, "_", tr))
    write.GWASpoly(fit_thr, trait = tr, filename = paste0(base, "_scores.csv"), what = "scores")
    qtl <- get.QTL(fit_thr, traits = tr, models = "additive", bp.window = 5e6)
    data.table::fwrite(qtl, paste0(base, "_QTL.csv"))
  }
  
  RUN_REGISTRY[[length(RUN_REGISTRY) + 1]] <<- list(
    fit     = fit_thr,
    traits  = keep,
    tag     = tag,
    cohort  = cohort_name,
    covmode = cov_mode
  )
  invisible(TRUE)
}

# =============================================================
# 8) Execute full grid (PF covariate ONLY for ALL)
# =============================================================
for (cohort_name in names(cohorts)) {
  dfC <- cohorts[[cohort_name]]
  message("\n=== Cohort: ", cohort_name, " | N=", nrow(dfC), " ===")
  
  # Policy: ALL -> withPFcov only; FF -> noPFcov only
  cov_mode <- if (cohort_name == "ALL") "withPFcov" else "noPFcov"
  message("  Running covariate mode: ", cov_mode)
  
  for (tag in names(trait_lists)) {
    run_one(
      df        = dfC,
      trait_vec = trait_lists[[tag]],
      cov_mode  = cov_mode,
      cohort_name = cohort_name,
      tag = tag
    )
  }
}

# =============================================================
# 9) End-of-script: print ALL Manhattan plots to Plots pane
# =============================================================
plot_all_manhattans <- function(registry, arrange_grids = FALSE, ncol = 3){
  if (!is.list(registry) || length(registry) == 0) {
    message("No GWAS runs were registered; nothing to plot.")
    return(invisible(TRUE))
  }
  can_grid <- arrange_grids && requireNamespace("cowplot", quietly = TRUE)
  plot_list <- list()
  
  for (i in seq_along(registry)) {
    r <- registry[[i]]
    for (tr in r$traits) {
      lbl <- paste0(
        "Manhattan — ", r$tag,
        " | Cohort: ", r$cohort,
        " | ", ifelse(r$covmode=="withPFcov","with PF covariate","no PF covariate"),
        " | Trait: ", tr
      )
      p <- try(manhattan.plot(r$fit, traits = tr), silent = TRUE)
      if (inherits(p, "try-error") || !inherits(p, "ggplot")) {
        manhattan.plot(r$fit, traits = tr)
        try(title(lbl), silent = TRUE)
      } else {
        p <- p + ggplot2::ggtitle(lbl)
        if (can_grid) plot_list[[length(plot_list)+1]] <- p else print(p)
      }
    }
  }
  if (can_grid && length(plot_list)) {
    per_page <- ncol * 2
    chunks <- split(plot_list, ceiling(seq_along(plot_list)/per_page))
    for (pg in seq_along(chunks)) print(cowplot::plot_grid(plotlist = chunks[[pg]], ncol = ncol))
  }
  invisible(TRUE)
}

cat("\n--- Printing all Manhattan plots to the Plots pane ---\n")
plot_all_manhattans(RUN_REGISTRY, arrange_grids = FALSE, ncol = 3)
cat("\nDone. Tables saved in: ", normalizePath(out_dir), "\n")

# Optional: Summary of completed runs (add this if desired)
cat("Total GWAS runs completed:", length(RUN_REGISTRY), "\n")
if (length(RUN_REGISTRY) > 0) {
  summary_runs <- tibble(
    Cohort = sapply(RUN_REGISTRY, `[[`, "cohort"),
    Cov_Mode = sapply(RUN_REGISTRY, `[[`, "covmode"),
    Traits = sapply(RUN_REGISTRY, function(x) paste(x$traits, collapse = ", "))
  )
  print(summary_runs)
}
