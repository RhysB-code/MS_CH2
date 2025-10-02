# =============================================================
# GWASpoly pipeline for tetraploid blackberry budbreak (binary)
#  - Ignores Dec 20 & Jan 1
#  - Binary traits at Jan 3 and Jan 15
#  - LOCO kinship; additive + 1-dom models
#  - Optional PF covariate re-run
# =============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ---- Packages ----
pkgs <- c("tidyverse","data.table","GWASpoly","ggplot2")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- User paths (EDIT if needed) ----
pheno_file <- "GWAS_Data_Budbreak.csv"  # your uploaded phenos
geno_file  <- "captureseq2025_GATK_postmean_70K_updog_filtered.csv" # rec: 70K Updog-filtered
pf_file    <- "PF_covariate.csv"        # optional; NAME,PF (0/1). If missing, we skip PF run.

out_dir <- "gwaspoly_budbreak_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================
# 1) Build binary phenotypes (0/1) for Jan 3 and Jan 15
# =============================================================
ph <- fread(pheno_file)

# Try to derive an individual NAME that matches your genotype IDs.
# Many of your rows look like "APF-665T_3_19_B_1" in column "Lug Order".
# We'll use the token BEFORE the first underscore as the genotype name.
name_col <- if ("Lug Order" %in% names(ph)) "Lug Order" else names(ph)[1]
ph <- ph %>%
  mutate(NAME  = sub("(_.*)$","", .data[[name_col]]),
         TBUDS = coalesce(`Total Buds`, `Total Buds.1`))  # fall back if duplicated column

# Keep only evaluation columns we need
# NOTE: We deliberately ignore "20-Dec" and "1-Jan" as requested.
need_cols <- c("NAME","TBUDS","3-Jan","15-Jan")
missing_needed <- setdiff(need_cols, names(ph))
if (length(missing_needed)) {
  stop("Missing expected columns in phenotype file: ", paste(missing_needed, collapse=", "))
}

# Coerce bud counts to numeric and cap to TBUDS
ph <- ph %>%
  mutate(`3-Jan`  = suppressWarnings(as.numeric(`3-Jan`)),
         `15-Jan` = suppressWarnings(as.numeric(`15-Jan`))) %>%
  mutate(`3-Jan`  = pmax(0, pmin(`3-Jan`,  TBUDS)),
         `15-Jan` = pmax(0, pmin(`15-Jan`, TBUDS)))

# Make binary per cutting (1 if any bud broke)
ph <- ph %>%
  mutate(BB_3Jan_bin  = if_else(`3-Jan`  > 0, 1, 0),
         BB_15Jan_bin = if_else(`15-Jan` > 0, 1, 0))

# Aggregate to genotype level: if ANY cutting for a genotype broke bud → 1
ph_g <- ph %>%
  group_by(NAME) %>%
  summarize(
    BB_3Jan  = as.integer(any(BB_3Jan_bin  == 1, na.rm = TRUE)),
    BB_15Jan = as.integer(any(BB_15Jan_bin == 1, na.rm = TRUE)),
    n_obs    = n()
  ) %>% ungroup()

# Write GWASpoly-style phenotype file (NAME first, then traits; covariates can follow)
pheno_out_noPF <- file.path(out_dir, "phenotypes_binary_noPF.csv")
fwrite(ph_g[,c("NAME","BB_3Jan","BB_15Jan")], pheno_out_noPF)

# =============================================================
# 2) Prepare genotype (dosage) for GWASpoly
#    Expect columns: Marker, Chrom, Position, [REF, ALT optional], then one col per NAME
#    Format = "numeric" (dosage 0..4 for tetraploid).
# =============================================================
# Light reader that tolerates large files
ghead <- fread(geno_file, nrows = 5)
# Try to normalize column names
setnames(ghead, old = names(ghead),
         new = str_replace_all(tolower(names(ghead)), c("chromosome"="chrom","chr"="chrom","position"="pos")))

# Expect first three: marker, chrom, pos (any casing). If not exact, try to guess/rename.
required <- c("marker","chrom","pos")
cn <- tolower(colnames(ghead))
# Heuristics to rename
ren <- list()
if (!"marker" %in% cn) {
  # try SNP/ID column
  cand <- which(cn %in% c("id","snp","markername","marker_id","markerid"))[1]
  if (length(cand)==1 && !is.na(cand)) ren[[cand]] <- "marker"
}
if (!"chrom" %in% cn) {
  cand <- which(cn %in% c("chrom","chromosome","chr"))[1]
  if (length(cand)==1 && !is.na(cand)) ren[[cand]] <- "chrom"
}
if (!"pos" %in% cn) {
  cand <- which(cn %in% c("pos","position","bp","bp_position"))[1]
  if (length(cand)==1 && !is.na(cand)) ren[[cand]] <- "pos"
}

if (length(ren)) {
  g <- fread(geno_file)
  old <- names(g)
  for (i in seq_along(ren)) {
    names(g)[as.integer(names(ren)[i])] <- ren[[i]]
  }
} else {
  g <- fread(geno_file)
  names(g) <- tolower(names(g))
}
stopifnot(all(required %in% names(g)))

# Keep only samples present in phenotype
sample_cols <- intersect(names(g), ph_g$NAME)
if (length(sample_cols) == 0) {
  stop("No overlap between genotype sample names and phenotype NAME. ",
       "Double-check that NAME in phenotypes matches column headers in genotype file.")
}

g_sub <- g[, c(required, sample_cols), with = FALSE]

# Quick genotype sanity check: dosages should be numeric 0..4 for tetraploid
# (GWASpoly will also check.)
# fwrite for GWASpoly input
geno_out <- file.path(out_dir, "geno_numeric_tetra.csv")
fwrite(g_sub, geno_out)

# =============================================================
# 3) Run GWASpoly (LOCO), without PF covariate
# =============================================================
ploidy_level <- 4
data0 <- read.GWASpoly(ploidy = ploidy_level,
                       pheno.file = pheno_out_noPF,
                       geno.file  = geno_out,
                       format = "numeric", n.traits = 2, delim = ",")

# LOCO kinship to reduce proximal contamination
# (per GWASpoly vignette; see set.K with LOCO = TRUE)
data0 <- set.K(data0, LOCO = TRUE, n.core = max(1, parallel::detectCores()-1))

# Parameters: filter rare/low-information markers by max genotype frequency
N <- nrow(data0@pheno)
params0 <- set.params(geno.freq = 1 - 5/N)  # Jeff’s rule of thumb in vignette

# Test additive and 1-copy dominance in both directions
models <- c("additive","1-dom")

fit0 <- GWASpoly(data = data0,
                 models = models,
                 traits = c("BB_3Jan","BB_15Jan"),
                 params = params0,
                 n.core = max(1, parallel::detectCores()-1),
                 quiet = TRUE)

# Thresholds (M.eff is generally less over-conservative than Bonferroni)
fit0_thr <- set.threshold(fit0, method = "M.eff", level = 0.05)

# Save results (scores) & plots
write.GWASpoly(fit0_thr, trait = "BB_3Jan",
               filename = file.path(out_dir,"BB_3Jan_scores_noPF.csv"),
               what = "scores")
write.GWASpoly(fit0_thr, trait = "BB_15Jan",
               filename = file.path(out_dir,"BB_15Jan_scores_noPF.csv"),
               what = "scores")

png(file.path(out_dir,"QQ_noPF_BB_3Jan.png"), width=1200, height=900, res=150)
print(qq.plot(fit0_thr, trait="BB_3Jan")); dev.off()

png(file.path(out_dir,"QQ_noPF_BB_15Jan.png"), width=1200, height=900, res=150)
print(qq.plot(fit0_thr, trait="BB_15Jan")); dev.off()

png(file.path(out_dir,"Manhattan_noPF.png"), width=1600, height=900, res=150)
print(manhattan.plot(fit0_thr, traits=c("BB_3Jan","BB_15Jan"))); dev.off()

# Export top QTL (one lead per ~5 Mb window; tweak to your LD)
qtl_3 <- get.QTL(fit0_thr, traits="BB_3Jan",  models="additive", bp.window = 5e6)
qtl_15<- get.QTL(fit0_thr, traits="BB_15Jan", models="additive", bp.window = 5e6)
fwrite(qtl_3,  file.path(out_dir,"QTL_noPF_BB_3Jan.csv"))
fwrite(qtl_15, file.path(out_dir,"QTL_noPF_BB_15Jan.csv"))

# =============================================================
# 4) OPTIONAL: add PF covariate and re-run
#     PF_covariate.csv should be: NAME,PF  (PF = 0/1)
# =============================================================
if (file.exists(pf_file)) {
  pf <- fread(pf_file) %>% mutate(PF = as.factor(PF))
  ph_pf <- ph_g %>% left_join(pf, by = "NAME")
  if (any(is.na(ph_pf$PF))) {
    warning("PF missing for some genotypes; those entries will drop in covariate run.")
    ph_pf <- ph_pf %>% filter(!is.na(PF))
  }
  pheno_out_PF <- file.path(out_dir, "phenotypes_binary_withPF.csv")
  fwrite(ph_pf[,c("NAME","BB_3Jan","BB_15Jan","PF")], pheno_out_PF)
  
  data1 <- read.GWASpoly(ploidy = ploidy_level,
                         pheno.file = pheno_out_PF,
                         geno.file  = geno_out,
                         format = "numeric", n.traits = 2, delim = ",")
  # LOCO kinship
  data1 <- set.K(data1, LOCO = TRUE, n.core = max(1, parallel::detectCores()-1))
  
  # Include PF as fixed factor (Q+K style fixed covariate)
  params1 <- set.params(geno.freq = 1 - 5/nrow(data1@pheno),
                        fixed = "PF", fixed.type = "factor")
  
  fit1 <- GWASpoly(data = data1,
                   models = models,
                   traits = c("BB_3Jan","BB_15Jan"),
                   params = params1,
                   n.core = max(1, parallel::detectCores()-1),
                   quiet = TRUE)
  
  fit1_thr <- set.threshold(fit1, method = "M.eff", level = 0.05)
  
  write.GWASpoly(fit1_thr, trait = "BB_3Jan",
                 filename = file.path(out_dir,"BB_3Jan_scores_withPF.csv"),
                 what = "scores")
  write.GWASpoly(fit1_thr, trait = "BB_15Jan",
                 filename = file.path(out_dir,"BB_15Jan_scores_withPF.csv"),
                 what = "scores")
  
  png(file.path(out_dir,"Manhattan_withPF.png"), width=1600, height=900, res=150)
  print(manhattan.plot(fit1_thr, traits=c("BB_3Jan","BB_15Jan"))); dev.off()
  
  qtl_3_pf  <- get.QTL(fit1_thr, traits="BB_3Jan",  models="additive", bp.window = 5e6)
  qtl_15_pf <- get.QTL(fit1_thr, traits="BB_15Jan", models="additive", bp.window = 5e6)
  fwrite(qtl_3_pf,  file.path(out_dir,"QTL_withPF_BB_3Jan.csv"))
  fwrite(qtl_15_pf, file.path(out_dir,"QTL_withPF_BB_15Jan.csv"))
} else {
  message("PF_covariate.csv not found; skipped PF covariate run.")
}

# =============================================================
# Done. Outputs in: gwaspoly_budbreak_out/


# =============================================================
# Massive GWASpoly runset: Binary, Percent, BLUPs, BLUEs, PCA
# Cohorts: ALL / PF-only / FF-only; Re-run with PF covariate
# Models: additive, 1-dom, general; LOCO kinship
# =============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ---- Packages ----
pkgs <- c("tidyverse","data.table","GWASpoly","ggplot2","stringr",
          "lme4","lmerTest","emmeans","parallel")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- User paths ----
pheno_file <- "GWAS_Data_Budbreak.csv"  # your per-cutting phenos
geno_file  <- "captureseq2025_GATK_postmean_70K_updog_filtered.csv" # dosage 0..4
out_dir    <- "gwaspoly_runsets"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================
# 0) Helper: ID parsing, safe reading, and utilities
# =============================================================
.read_any <- function(path, nrows=0){
  # fread is fast; fall back to read.csv if needed
  if (!file.exists(path)) stop("File not found: ", path)
  tryCatch(
    fread(path, nrows = ifelse(nrows>0, nrows, -1)),
    error = function(e) read.csv(path, check.names = FALSE)
  )
}

get_name <- function(x){
  # Genotype name = token before first underscore
  sub("(_.*)$","", x)
}

get_rep  <- function(x){
  # Replicate token after first underscore, if present
  ifelse(grepl("_", x), sub("^[^_]+_","", x), NA_character_)
}

# ensure numeric and clamp to [0, TBUDS]
.clamp_counts <- function(v, maxv){
  v <- suppressWarnings(as.numeric(v))
  v <- pmax(0, v)
  v <- if (!missing(maxv)) pmin(v, maxv) else v
  v
}

# =============================================================
# 1) Ingest phenotypes; derive Binary, Percent; make tall + wide
#    (Ignores 20-Dec and 1-Jan, per your directive)
# =============================================================
ph_raw <- .read_any(pheno_file)

# Identify columns
name_col <- if ("Lug Order" %in% names(ph_raw)) "Lug Order" else names(ph_raw)[1]
buds_col <- c("Total Buds","Total Buds.1","TBuds","Buds","Total_Buds")
buds_col <- buds_col[buds_col %in% names(ph_raw)]
if (length(buds_col)==0) stop("Couldn't find a 'Total Buds' column.")

need_dates <- c("3-Jan","15-Jan")
miss <- setdiff(c(name_col, need_dates), names(ph_raw))
if (length(miss)) stop("Missing expected columns: ", paste(miss, collapse=", "))

ph <- ph_raw %>%
  mutate(NAME      = get_name(.data[[name_col]]),
         REP       = get_rep(.data[[name_col]]),
         TBUDS     = coalesce(!!!syms(buds_col))) %>%
  mutate(`3-Jan`   = .clamp_counts(`3-Jan`,  TBUDS),
         `15-Jan`  = .clamp_counts(`15-Jan`, TBUDS)) %>%
  mutate(BIN_3Jan  = as.integer(`3-Jan`  > 0),
         BIN_15Jan = as.integer(`15-Jan` > 0),
         PCT_3Jan  = ifelse(TBUDS>0, 100*`3-Jan`/TBUDS, NA_real_),
         PCT_15Jan = ifelse(TBUDS>0, 100*`15-Jan`/TBUDS, NA_real_))

# -------------------------------------------------------------
# 1A) PF/FF classification by name + override hooks
# -------------------------------------------------------------
# Default rule: APF* = PF; A-* = FF
is_pf_default <- function(nm) startsWith(nm, "APF")
is_ff_default <- function(nm) startsWith(nm, "A-")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TODO: Edit these override vectors to fix exceptions:
pf_overrides_add <- c(    # e.g., "A-2491T"
)
ff_overrides_add <- c(    # e.g., "APF-661TN"
)
# Optionally list specific IDs to remove from each class:
pf_exclude <- c()
ff_exclude <- c()
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

pf_df <- tibble(NAME = unique(ph$NAME)) %>%
  mutate(PF = case_when(
    NAME %in% pf_excludes ~ NA,       # not used
    NAME %in% pf_overrides_add ~ 1L,
    NAME %in% ff_overrides_add ~ 0L,
    is_pf_default(NAME) ~ 1L,
    is_ff_default(NAME) ~ 0L,
    TRUE ~ NA_integer_
  ))

# If you want unknowns to drop in subgroup runs but stay in ALL, keep NA here.

# =============================================================
# 2) BLUPs & BLUEs at Jan 3 and Jan 15 (rep-aware)
#     BLUPs: genotype random → conditional modes
#     BLUEs: genotype fixed → adjusted means via emmeans
# =============================================================
mk_blup_blue <- function(df, resp){
  # df: per-cutting rows with NAME, REP, plus covariates if desired
  out <- list()
  # ---- BLUPs (random genotype) ----
  # Minimal model (intercept only plus random G); add batch/block if you have them
  m_blup <- lmer(reformulate(termlabels = NULL, response = resp),
                 data = df %>% filter(!is.na(.data[[resp]])),
                 REML = TRUE)
  blups  <- ranef(m_blup, condVar = FALSE)$NAME
  blups  <- tibble(NAME = rownames(blups),
                   !!paste0("BLUP_", resp) := as.numeric(blups[["(Intercept)"]]))
  out$BLUP <- blups
  
  # ---- BLUEs (fixed genotype) ----
  m_blue <- lmer(as.formula(paste0(resp, " ~ 1 + NAME")),
                 data = df %>% filter(!is.na(.data[[resp]])),
                 REML = TRUE)
  emm <- emmeans(m_blue, ~ NAME)
  blues <- as_tibble(emm) %>%
    transmute(NAME, !!paste0("BLUE_", resp) := as.numeric(emmean))
  out$BLUE <- blues
  out
}

# Build BLUPs/BLUEs for Percent traits (more informative than counts)
bb_blup_blue_3  <- mk_blup_blue(ph %>% select(NAME, REP, PCT_3Jan)  %>% drop_na(NAME), "PCT_3Jan")
bb_blup_blue_15 <- mk_blup_blue(ph %>% select(NAME, REP, PCT_15Jan) %>% drop_na(NAME), "PCT_15Jan")

# Merge to one table per metric
BLUP_tbl <- full_join(bb_blup_blue_3$BLUP,  bb_blup_blue_15$BLUP,  by = "NAME")
BLUE_tbl <- full_join(bb_blup_blue_3$BLUE,  bb_blup_blue_15$BLUE,  by = "NAME")

# Also compute simple genotype-level means (quick continuous traits):
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
# 3) PCA from multi-date phenotypes (use BLUPs for robustness)
#     PC1 ~= overall responsiveness; PC2 ~= timing
# =============================================================
pca_src <- BLUP_tbl %>% select(NAME, BLUP_PCT_3Jan = BLUP_PCT_3Jan, BLUP_PCT_15Jan = BLUP_PCT_15Jan)
pc_complete <- pca_src %>% filter(complete.cases(.))
PC_tbl <- tibble(NAME = pca_src$NAME,
                 PC1 = NA_real_, PC2 = NA_real_)
if (nrow(pc_complete) >= 5) {
  X <- as.matrix(pc_complete[, -1])
  P <- prcomp(X, center = TRUE, scale. = TRUE)
  sc <- as.data.frame(P$x[, 1:2, drop = FALSE])
  sc$NAME <- pc_complete$NAME
  PC_tbl <- PC_tbl %>%
    select(-PC1, -PC2) %>%
    left_join(sc %>% transmute(NAME, PC1, PC2), by = "NAME")
}

# =============================================================
# 4) Assemble master phenotype panel with all requested traits
#     - Binary (>0): BIN_3Jan, BIN_15Jan
#     - Continuous (means): PCT_3Jan, PCT_15Jan
#     - BLUPs: BLUP_PCT_3Jan, BLUP_PCT_15Jan
#     - BLUEs: BLUE_PCT_3Jan, BLUE_PCT_15Jan
#     - PCA: PC1, PC2
# =============================================================
PH_MASTER <- MEAN_tbl %>%
  full_join(BLUP_tbl, by = "NAME") %>%
  full_join(BLUE_tbl, by = "NAME") %>%
  full_join(PC_tbl,  by = "NAME") %>%
  left_join(pf_df,   by = "NAME")    # adds PF as 0/1/NA

# =============================================================
# 5) Prepare genotype dosage for GWASpoly (0..4)
# =============================================================
message("Reading genotype matrix (this may take a bit)...")
G <- .read_any(geno_file)
names(G) <- tolower(names(G))
# Standardize first columns
ren_map <- c("markername"="marker","marker_id"="marker","markerid"="marker","id"="marker",
             "chromosome"="chrom","chr"="chrom","position"="pos","bp"="pos","bp_position"="pos")
for (k in names(ren_map)) {
  if (k %in% names(G)) names(G)[names(G)==k] <- ren_map[[k]]
}
req <- c("marker","chrom","pos")
stopifnot(all(req %in% names(G)))

# Sample intersection
samp_cols <- intersect(names(G), PH_MASTER$NAME)
if (!length(samp_cols)) stop("No overlap between genotype columns and phenotype NAMEs.")
G_sub <- G[, c(req, samp_cols), drop = FALSE]

# Save once (reused by all runs)
geno_out <- file.path(out_dir, "geno_numeric_tetra.csv")
fwrite(G_sub, geno_out)

# =============================================================
# 6) Define run grid: cohorts × covariate × traits
# =============================================================
# Cohorts
cohorts <- list(
  ALL = PH_MASTER,
  PF  = PH_MASTER %>% filter(PF == 1L),
  FF  = PH_MASTER %>% filter(PF == 0L)
)

# Covariate mode
cov_modes <- c("noPFcov","withPFcov")

# Trait sets to run (feel free to toggle)
trait_lists <- list(
  BINARY = c("BIN_3Jan","BIN_15Jan"),
  PCT    = c("PCT_3Jan","PCT_15Jan"),
  BLUP   = c("BLUP_PCT_3Jan","BLUP_PCT_15Jan"),
  BLUE   = c("BLUE_PCT_3Jan","BLUE_PCT_15Jan"),
  PCA    = c("PC1","PC2")         # PC2 optional; comment out if not needed
)

# GWASpoly models
models_use <- c("additive","1-dom","general")

# =============================================================
# 7) Core runner: builds a pheno CSV and runs GWASpoly
# =============================================================
run_one <- function(df, trait_vec, cov_mode, cohort_name, tag){
  df_use <- df %>% select(NAME, all_of(trait_vec), PF)
  # Drop traits with all NA (or constant)
  keep <- trait_vec[sapply(trait_vec, function(t) {
    v <- df_use[[t]]
    any(!is.na(v)) && (length(unique(na.omit(v))) > 1)
  })]
  if (length(keep) == 0) {
    message("[SKIP] No informative traits for ", tag, " ", cohort_name, " ", cov_mode)
    return(invisible(NULL))
  }
  
  # Build pheno table for GWASpoly
  pheno_out <- file.path(out_dir, paste0("pheno_", tag, "_", cohort_name, "_", cov_mode, ".csv"))
  if (cov_mode == "withPFcov") {
    # PF covariate requires variability; if all PF are same (e.g., PF-only cohort), drop covariate
    if (!("PF" %in% names(df_use)) || length(unique(na.omit(df_use$PF))) < 2) {
      message("[INFO] PF covariate constant/absent in ", cohort_name, "; running without PF covariate.")
      fwrite(df_use %>% select(NAME, all_of(keep)), pheno_out)
      fixed_formula <- NULL
    } else {
      fwrite(df_use %>% select(NAME, all_of(keep), PF), pheno_out)
      fixed_formula <- "PF"  # factor fixed in set.params
    }
  } else {
    fwrite(df_use %>% select(NAME, all_of(keep)), pheno_out)
    fixed_formula <- NULL
  }
  
  # ---- Read GWASpoly project ----
  data_gp <- read.GWASpoly(
    ploidy    = 4,
    pheno.file= pheno_out,
    geno.file = geno_out,
    format    = "numeric",
    n.traits  = length(keep),
    delim     = ","
  )
  
  # ---- LOCO kinship ----
  data_gp <- set.K(data_gp, LOCO = TRUE, n.core = max(1, parallel::detectCores()-1))
  
  # ---- Params ----
  N <- nrow(data_gp@pheno)
  params <- set.params(
    geno.freq  = 1 - 5/N,                    # light MAF-ish filter
    fixed      = if (is.null(fixed_formula)) NULL else fixed_formula,
    fixed.type = if (is.null(fixed_formula)) NULL else "factor"
  )
  
  # ---- Fit ----
  fit <- GWASpoly(
    data   = data_gp,
    models = models_use,
    traits = keep,
    params = params,
    n.core = max(1, parallel::detectCores()-1),
    quiet  = TRUE
  )
  
  fit_thr <- set.threshold(fit, method = "M.eff", level = 0.05)
  
  # ---- Save outputs ----
  for (tr in keep) {
    base <- file.path(out_dir, paste0("GWAS_", tag, "_", cohort_name, "_", cov_mode, "_", tr))
    write.GWASpoly(fit_thr, trait = tr, filename = paste0(base, "_scores.csv"), what = "scores")
    # Top QTL (coarse LD window; tweak for your genome)
    qtl <- get.QTL(fit_thr, traits = tr, models = "additive", bp.window = 5e6)
    fwrite(qtl, paste0(base, "_QTL.csv"))
    # Plots
    png(paste0(base, "_QQ.png"), width = 1200, height = 900, res = 150)
    print(qq.plot(fit_thr, trait = tr)); dev.off()
  }
  # Manhattan across all kept traits
  man_base <- file.path(out_dir, paste0("GWAS_", tag, "_", cohort_name, "_", cov_mode, "_Manhattan.png"))
  png(man_base, width = 1800, height = 1000, res = 150)
  print(manhattan.plot(fit_thr, traits = keep)); dev.off()
  
  invisible(TRUE)
}

# =============================================================
# 8) Run the full grid
# =============================================================
for (cohort_name in names(cohorts)) {
  dfC <- cohorts[[cohort_name]]
  message("\n=== Cohort: ", cohort_name, " | N=", nrow(dfC), " ===")
  
  for (cov_mode in cov_modes) {
    for (tag in names(trait_lists)) {
      trait_vec <- trait_lists[[tag]]
      run_one(df = dfC, trait_vec = trait_vec, cov_mode = cov_mode,
              cohort_name = cohort_name, tag = tag)
    }
  }
}

message("\nAll runs finished. Results in: ", normalizePath(out_dir))
# =============================================================
# End
# =============================================================

