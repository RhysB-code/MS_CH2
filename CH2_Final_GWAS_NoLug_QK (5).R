## =====================================================================================================================
## CH2_Final_GWAS_NoLug_QK.R   (Rhys Brock, MS thesis Ch. 2; 20 Nov 2026)  -- fresh rerun of this year's GWAS
##     source("C:/Users/RhysB/OneDrive/Desktop/CH2 GWAS/CH2_Final_GWAS_NoLug_QK.R")
## Needs in <work_dir>:  PHENO_MasterSheet_CURRENT.xlsx   and   Overlap_CaptureSeq_FlexSeq_SharedMarkers_GWASready_plusPF.csv
## RUNS (7):  BB10_pooled = field 10% budbreak, 2025-26 pooled, the 3 very early genotypes (A_2615T, A_2857T, A_2877T) set to DAY 45 (already done in GWAS_Phenotypes)
##            BB90_pooled = field 90% budbreak, 2025-26 pooled (day of year)
##            Exp1.1, Exp1.2, Exp2.1, Exp2.2 (each experiment alone) and LC_Both (ALL four experiments), low-chill BLUEs WITHOUT lug correction ("| no lug corr." columns)
## MODEL: GWASpoly ploidy 4, LOCO kinship, P3D, Q + K (Q = Cluster_2 of ONE K = 2 DAPC on the concordance-filtered markers, 60 PCs), M.eff Bonferroni alpha 0.05, MAF 0.05, geno.freq 0.98,
##        1-dom / 2-dom markers need >= 10 genotypes in the smaller class. K-only fits are OFF (set RUN_K_ONLY <- TRUE to add them).
## OUTPUT (every run is fitted ONCE with all four models; each folder gets its own models, so results are identical to separate fits):
##   MAIN set (all 6 runs; additive + 1-dom): Manhattan_<run>.png + QQ_<run>_by_chromosome.png go to  <work_dir>/FieldBudbreakGWAS/Final_Results  (field)  and  <work_dir>/LowChillGWAS/Final_Output  (low-chill)
##   ALTERNATIVE set (LOW-CHILL runs only; additive, general, 1-dom, 2-dom): its PNGs go to  <work_dir>/LowChillGWAS/Final_Output/5_Models
##   Workbooks / pptx / rds: <work_dir>/Final_GWAS_noLug (main) and <work_dir>/Final_GWAS_noLug_AllModels (alternative)
## ALSO (section 3b, runs before the GWAS fits): DAPC K=2 tables/plots + FF/PF histograms -> <work_dir>/DAPC_K=2 ; correlation matrices, gene density, marker density -> <work_dir>/Thesis_Figures
##   (gene density uses <work_dir>/GENO_Inputs_Seq_Files/Annotated Reference Genome.xlsx)
## FOLDER MAP (everything is written under <work_dir> = "C:/Users/RhysB/OneDrive/Desktop/CH2 GWAS"):
##   FieldBudbreakGWAS/Final_Results        Manhattan + QQ for field 10% and field 90%
##   LowChillGWAS/Final_Output              Manhattan + QQ for Exp 1.1, 1.2, 2.1, 2.2 and LC_Both (additive + 1-dom);  .../Final_Output/5_Models = the same plus general + 2-dom
##   DAPC_K=2                               DAPC tables/plots, % membership + coordinates, FF/PF histograms (field and low-chill)
##   Thesis_Figures                         correlation matrices (low-chill, field 10 vs 90, field vs low-chill), gene density, marker density
##   Final_GWAS_noLug / Final_GWAS_noLug_AllModels   results workbook, pptx, rds (main set / alternative set)
##   (to rename the low-chill results folder change LC_PNG_DIR in the SETTINGS block)
## When it finishes, paste the block between  COPY FROM HERE  and  END OF COPY  back to Claude.
## =====================================================================================================================
suppressMessages({ library(data.table); library(readxl); library(GWASpoly) })
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
HAVE_ADEGENET <- requireNamespace("adegenet", quietly = TRUE); if (HAVE_ADEGENET) suppressMessages(library(adegenet))

## ================= SETTINGS =================
if (!exists("work_dir")) work_dir <- "C:/Users/RhysB/OneDrive/Desktop/CH2 GWAS"
if (!exists("master_xlsx")) master_xlsx <- file.path(work_dir, "PHENO_MasterSheet_CURRENT.xlsx")   # to use another copy: master_xlsx <- "C:/.../PHENO_MasterSheet_CURRENT.xlsx" before source()
GENO_DIR <- file.path(work_dir, "GENO_Inputs_Seq_Files")   # genotype inputs live here
if (!exists("marker_file")) marker_file <- file.path(GENO_DIR, "Overlap_CaptureSeq_FlexSeq_SharedMarkers_GWASready_plusPF.csv")   # merged CaptureSeq + FlexSeq (+ PF) GWAS input
RUN_K_ONLY  <- FALSE      # TRUE = also fit every run with kinship only (doubles the run time)
CONC_CUTOFF <- 0.90; MIN_PAIRS <- 10; KEEP_UNTESTABLE <- FALSE; PAIR_R_MIN <- 0.83; PAIR_R_MAX <- 0.93
K <- 2; N_PCA <- 60; SEED <- 2026; MIN_MAF <- 0.05; MAX_GENO_FREQ <- 0.98; THRESH_LEVEL <- 0.05; MIN_CLASS_N <- 10
MODELS_ALL <- c("additive", "general", "1-dom", "2-dom")                  # fitted once
SETS <- list(main = c("additive", "1-dom-alt", "1-dom-ref"), alt = c("additive", "general", "1-dom-alt", "1-dom-ref", "2-dom-alt", "2-dom-ref"))   # what each folder reports
OUT_DIRS <- c(main = file.path(work_dir, "Final_GWAS_noLug"), alt = file.path(work_dir, "Final_GWAS_noLug_AllModels"))   # xlsx / pptx / rds
FIELD_PNG_DIR <- file.path(work_dir, "FieldBudbreakGWAS", "Final_Results"); LC_PNG_DIR <- file.path(work_dir, "LowChillGWAS", "Final_Output")
PNG_DIRS <- list(main_field = FIELD_PNG_DIR, main_lc = LC_PNG_DIR, alt_lc = file.path(LC_PNG_DIR, "5_Models"))
is_field <- function(r) grepl("^BB", r)
for (d_ in PNG_DIRS) dir.create(d_, recursive = TRUE, showWarnings = FALSE)
LOCUS_WINDOW <- 2.5e6; OVERLAP_WINDOW <- 1e5
KEY_MARKERS <- c("^Ra01:144240_", "^Ra01:625526_", "^Ra05:34864163_", "^Ra03:2190906_")
N_CORES <- max(1L, parallel::detectCores() - 1L)
RUNS <- data.table(Run = c("BB10_pooled", "BB90_pooled", "Exp1.1", "Exp1.2", "Exp2.1", "Exp2.2", "LC_Both"),
  Col = c("10% 2025-26 pooled, early-3 at day 45", "90% 2025-26 pooled", "BLUE Exp 1.1", "BLUE Exp 1.2", "BLUE Exp 2.1", "BLUE Exp 2.2", "BLUE both years (all 4 experiments)"),
  Label = c("Field 10% budbreak 2025-26 (day of year; early-3 = day 45)", "Field 90% budbreak 2025-26 (day of year)", "Low-chill Exp 1.1 (6 Dec 2024)", "Low-chill Exp 1.2 (18 Dec 2024)", "Low-chill Exp 2.1 (11 Dec 2025)",
            "Low-chill Exp 2.2 (5 Jan 2026)", "Low-chill ALL four experiments (experiment fixed)"), Group = "main")
## ============================================

REAL_CHROMS <- paste0("Ra0", 1:7); MODEL_ORDER <- c("additive", "general", "1-dom-alt", "1-dom-ref", "2-dom-alt", "2-dom-ref")
MCOL <- c(additive = "#0072B2", general = "#009E73", "1-dom-alt" = "#E6B800", "1-dom-ref" = "#E6B800", "2-dom-alt" = "#CC79A7", "2-dom-ref" = "#CC79A7")
MPCH <- c(additive = 15, general = 3, "1-dom-alt" = 16, "1-dom-ref" = 17, "2-dom-alt" = 18, "2-dom-ref" = 4)

FIG_RUNS <- RUNS[Group == "main"]$Run
find_file <- function(p) { if (file.exists(p)) return(p); h <- list.files(work_dir, paste0("^", gsub("([.])", "\\\\\\1", basename(p)), "$"), recursive = TRUE, full.names = TRUE); if (length(h)) h[1] else NA_character_ }
lambda_gc <- function(score) { s <- score[is.finite(score)]; if (length(s) < 20) return(NA_real_); median(qchisq(10^(-s), 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
skew <- function(v) { v <- v[is.finite(v)]; if (length(v) < 3 || sd(v) == 0) return(NA_real_); mean((v - mean(v))^3) / sd(v)^3 }
quiet <- function(expr) { invisible(capture.output(r <- suppressMessages(expr))); r }
safe <- function(what, expr) tryCatch(expr, error = function(e) { cat(sprintf("[diagnostic '%s' failed: %s]\n", what, conditionMessage(e))); NULL })
rnd <- function(m) { storage.mode(m) <- "numeric"; pmin(pmax(round(m), 0), 4) }
tmp_files <- character(0); on.exit(unlink(tmp_files), add = TRUE)

## ---- 1. phenotypes from the updated master ----
mx <- find_file(master_xlsx); if (is.na(mx)) stop("Master phenotype workbook not found: ", master_xlsx)
## read a temp COPY of the master so a copy that is open in Excel / locked by OneDrive sync still reads
mx_read <- tempfile(fileext = ".xlsx"); cp_ok <- isTRUE(suppressWarnings(file.copy(mx, mx_read, overwrite = TRUE))); if (cp_ok) tmp_files <- c(tmp_files, mx_read) else mx_read <- mx
cat(sprintf("Master sheet: %s | size %s bytes | copied to a temp file: %s\n", mx, format(file.size(mx), big.mark = ","), cp_ok))
if (is.na(file.size(mx)) || file.size(mx) < 50000) stop("The master sheet is empty or online-only (", file.size(mx), " bytes). In File Explorer right-click it > Always keep on this device, wait for the green check, and run again.")
if (!cp_ok) stop("R cannot open or copy the master sheet (another program has it locked, or OneDrive is not finished syncing it). Close Excel, or copy the file to a normal folder (e.g. C:/Users/RhysB/Documents) and run  master_xlsx <- \"C:/Users/RhysB/Documents/PHENO_MasterSheet_CURRENT.xlsx\"  before source().")
rd <- function(s) as.data.table(read_excel(mx_read, sheet = s, skip = 2, .name_repair = "minimal"))
PH <- rd("GWAS_Phenotypes"); PR <- rd("Marker_Pairs"); PH[, Selection := trimws(as.character(Selection))]
PH <- PH[Selection != "APF_630T"]; stopifnot(!anyDuplicated(PH$Selection))
ph <- data.table(Genotype = PH$Selection)
for (i in seq_len(nrow(RUNS))) { if (!RUNS$Col[i] %in% names(PH)) stop("Column not found in GWAS_Phenotypes: ", RUNS$Col[i]); ph[[RUNS$Run[i]]] <- suppressWarnings(as.numeric(PH[[RUNS$Col[i]]])) }
cat(sprintf("================ PHENOTYPES (master %s; low-chill BLUEs, no lug correction) ================\n", basename(mx)))
PHS <- rbindlist(lapply(RUNS$Run, function(r) { v <- ph[[r]]; data.table(Run = r, n = sum(is.finite(v)), mean = round(mean(v, na.rm = TRUE), 2), sd = round(sd(v, na.rm = TRUE), 2), min = round(min(v, na.rm = TRUE), 2), max = round(max(v, na.rm = TRUE), 2), skew = round(skew(v), 2)) })); print(PHS)

## ---- 2. markers: K = 2 filtered set ----
mf <- find_file(marker_file); if (is.na(mf)) stop("Marker file not found: ", marker_file)
G <- fread(mf, header = TRUE, check.names = FALSE); names(G)[1:3] <- c("Marker", "Chrom", "Position")
G <- G[Chrom %in% REAL_CHROMS]; G <- G[!duplicated(paste(Chrom, Position))]
col_use <- trimws(as.character(PH[["Seq name (marker column used in GWAS)"]]))
pl <- list(); pa_ <- trimws(as.character(PR[["Marker column A (used in GWAS)"]])); pb_ <- trimws(as.character(PR[["Marker column B (same genotype, other sample)"]])); pr_ <- suppressWarnings(as.numeric(PR[["Concordance r"]]))
for (i in seq_along(pa_)) if (!is.na(pa_[i]) && pa_[i] %in% names(G) && pb_[i] %in% names(G) && pa_[i] != pb_[i] && is.finite(pr_[i]) && pr_[i] >= PAIR_R_MIN && pr_[i] <= PAIR_R_MAX) pl[[length(pl) + 1]] <- data.table(a = pa_[i], b = pb_[i])
pairs <- unique(rbindlist(pl)); pairs <- pairs[!duplicated(paste(pmin(a, b), pmax(a, b)))]
if (nrow(pairs) < 5) stop("Fewer than 5 cross-platform duplicate pairs found in the Marker_Pairs sheet.")
ca <- pairs$a; cb <- pairs$b; PA <- rnd(as.matrix(G[, ..ca])); PB <- rnd(as.matrix(G[, ..cb])); called <- !is.na(PA) & !is.na(PB)
n_pairs <- rowSums(called); conc <- ifelse(n_pairs > 0, rowSums(PA == PB & called, na.rm = TRUE) / n_pairs, NA_real_)
keep_conc <- (n_pairs >= MIN_PAIRS & !is.na(conc) & conc >= CONC_CUTOFF) | (KEEP_UNTESTABLE & n_pairs < MIN_PAIRS)
cat(sprintf("\nMarkers in the shared-marker file (real chromosomes, one per position): %d | duplicate pairs used: %d | kept by the < %.0f%% concordance filter: %d\n", nrow(G), nrow(pairs), 100 * CONC_CUTOFF, sum(keep_conc)))
keepv <- tolower(as.character(PH[["Keep in GWAS (one entry per genotype)"]])) %in% c("true", "1", "yes"); drop_sel <- PH$Selection[!keepv]
for (i in which(keepv & col_use != PH$Selection & col_use %in% names(G))) set(G, j = PH$Selection[i], value = G[[col_use[i]]])   # genotype whose marker column has another name (e.g. Tupy_FLEX)
if ("A_2868TN_7" %in% PH$Selection || !("A_2868TN_2" %in% PH$Selection)) stop("A_2868TN_2 must be the A-2868TN genotype in the phenotype sheet (not _7).")
ph <- ph[!(Genotype %in% drop_sel)]; samples <- intersect(ph$Genotype, names(G)); ph <- ph[match(samples, Genotype)]
Mg <- rnd(as.matrix(G[keep_conc, ..samples])); Gd <- cbind(G[keep_conc, .(Marker, Chrom, Position)], as.data.table(Mg))
geno_out <- tempfile(fileext = ".csv"); tmp_files <- c(tmp_files, geno_out); fwrite(Gd, geno_out); cat("Genotypes in the DAPC:", length(samples), "| markers in the GWAS:", nrow(Gd), "\n")

## ---- 3. DAPC K = 2 on the filtered markers ----
frq <- rowMeans(Mg, na.rm = TRUE) / 4; miss <- rowMeans(is.na(Mg)); use <- which(is.finite(frq) & frq >= MIN_MAF & frq <= 1 - MIN_MAF & miss <= 0.10)
X <- t(Mg[use, , drop = FALSE]) / 4; rownames(X) <- samples; cm <- colMeans(X, na.rm = TRUE); nas <- which(is.na(X), arr.ind = TRUE); if (nrow(nas)) X[nas] <- cm[nas[, 2]]
X <- X[, apply(X, 2, sd) > 0, drop = FALSE]; N_DAPC_MARKERS <- ncol(X); cat(sprintf("Markers used in the DAPC (MAF 0.05, <= 10%% missing): %d  (Margaret: 6069)\n", N_DAPC_MARKERS))
DAPC_COUNTS <- rbindlist(lapply(c(0.02, 0.05, 0.10, 0.20, 0.50), function(mc) { u <- which(is.finite(frq) & frq >= MIN_MAF & frq <= 1 - MIN_MAF & miss <= mc)
  data.table(Max_missing = mc, DAPC_markers = sum(apply(Mg[u, , drop = FALSE], 1, function(v) isTRUE(sd(v, na.rm = TRUE) > 0)))) }))
pc <- prcomp(X, center = TRUE, scale. = FALSE); n_pca <- min(N_PCA, ncol(pc$x), nrow(X) - 1); PCs <- pc$x[, seq_len(n_pca), drop = FALSE]; set.seed(SEED)
if (HAVE_ADEGENET) { fc <- quiet(find.clusters(PCs, n.clust = K, n.pca = n_pca, scale = FALSE)); ord <- names(sort(table(fc$grp), decreasing = TRUE)); grp <- factor(match(as.character(fc$grp), ord), levels = seq_len(K))
  dp <- quiet(dapc(PCs, grp = grp, n.pca = n_pca, n.da = K - 1, center = TRUE, scale = FALSE)); P <- as.matrix(dp$posterior)[, as.character(seq_len(K)), drop = FALSE]
} else { km <- kmeans(PCs, centers = K, nstart = 25, iter.max = 1e5); ord <- order(table(km$cluster), decreasing = TRUE); grp <- factor(match(km$cluster, ord), levels = seq_len(K)); P <- predict(MASS::lda(PCs, grouping = grp))$posterior[, as.character(seq_len(K)), drop = FALSE] }
colnames(P) <- paste0("Cluster_", seq_len(K)); Qn <- "Cluster_2"; cat("DAPC K = 2 cluster sizes (1 = largest):", paste(table(grp), collapse = " / "), "\n")
ph_q <- cbind(ph, as.data.table(P)[, ..Qn])

## ---- 3b. DAPC outputs, histograms, correlation, gene / marker density (publication figures + tables) ----
## DAPC_K=2 folder: membership %, coordinates, cluster x habit tables, DAPC plots, FF (blue) vs PF (red) histograms of phenotypes and BLUEs for both studies.
## Thesis_Figures folder: correlation matrices (Godwin et al. 2025 / Alexander style), gene density and marker density on all chromosomes (Chizk et al. 2023 style).
## 2024 is left out: only 2025 and 2026 field data are used. Low-chill Exp 1.1 / 1.2 are the Dec 2024-start experiments that are GWAS runs; drop them from the figures by setting FIG_EXPS below.
DAPC_DIR <- file.path(work_dir, "DAPC_K=2"); FIG_DIR <- file.path(work_dir, "Thesis_Figures"); FIG_EXPS <- c("1.1", "1.2", "2.1", "2.2")
GENE_FILE <- file.path(GENO_DIR, "Annotated Reference Genome.xlsx"); HAB_COL <- c(FF = "#1F4FBF", PF = "#D62728"); CL_COL <- c("1" = "#009E73", "2" = "#CC79A7")
for (d_ in c(DAPC_DIR, FIG_DIR)) dir.create(d_, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2"); suppressMessages(library(ggplot2)); HAVE_PW <- requireNamespace("patchwork", quietly = TRUE)
th_pub <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey90", colour = "black"), strip.text = element_text(face = "bold"), legend.key = element_blank())
sv <- function(p, name, dir, w, h) { ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white"); ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h) }
tag_comb <- function(pl, ncol, nm, dir, w, h) { if (HAVE_PW) sv(patchwork::wrap_plots(pl, ncol = ncol) + patchwork::plot_annotation(tag_levels = "a"), nm, dir, w, h) }
HAB <- toupper(trimws(as.character(PH$Habit[match(samples, PH$Selection)]))); HAB[!HAB %in% c("FF", "PF")] <- NA; names(HAB) <- samples
LD1 <- if (exists("dp", inherits = FALSE)) as.numeric(dp$ind.coord[, 1]) else as.numeric(predict(MASS::lda(PCs, grouping = grp))$x[, 1]); if (mean(LD1[grp == 2]) < mean(LD1[grp == 1])) LD1 <- -LD1
pve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
MEM <- data.table(Genotype = samples, Habit = HAB, Assigned_cluster = as.integer(as.character(grp)), Pct_class_1 = round(100 * P[, "Cluster_1"], 1), Pct_class_2 = round(100 * P[, "Cluster_2"], 1),
                  PC1 = round(pc$x[, 1], 3), PC2 = round(pc$x[, 2], 3), LD1 = round(LD1, 3))
CxH <- MEM[!is.na(Habit), .N, by = .(Habit, Assigned_cluster)]; CxH[, Pct_of_habit := round(100 * N / sum(N), 1), by = Habit]; CxH[, Pct_of_cluster := round(100 * N / sum(N), 1), by = Assigned_cluster]
CS <- MEM[, .(n = .N, n_FF = sum(Habit == "FF", na.rm = TRUE), n_PF = sum(Habit == "PF", na.rm = TRUE), mean_pct_class_1 = round(mean(Pct_class_1), 1), mean_pct_class_2 = round(mean(Pct_class_2), 1), n_ambiguous_30_70 = sum(Pct_class_2 > 30 & Pct_class_2 < 70)), by = Assigned_cluster][order(Assigned_cluster)]
## --- phenotype tables (raw and BLUE), 2025-26 field and four low-chill experiments only ---
fld <- data.table(Genotype = PH$Selection, Habit = toupper(trimws(as.character(PH$Habit))), `Field 10% 2025 (DOY)` = suppressWarnings(as.numeric(PH[["10% 2025 (DOY)"]])), `Field 10% 2026 (DOY; early-3 = day 45)` = suppressWarnings(as.numeric(PH[["10% 2026, early-3 at day 45"]])),
  `Field 90% 2025 (DOY)` = suppressWarnings(as.numeric(PH[["90% 2025 (DOY)"]])), `Field 90% 2026 (DOY)` = suppressWarnings(as.numeric(PH[["90% 2026 (DOY)"]])))
PL <- safe("plot-level sheet", { x <- rd("LowChill_plot_level"); x[, Selection := trimws(as.character(Selection))]; x[, Exp := as.character(Experiment)]; x <- x[Exp %in% FIG_EXPS & Selection %in% PH$Selection & is.finite(suppressWarnings(as.numeric(`% budbreak`)))]
  data.table(Genotype = x$Selection, Habit = toupper(trimws(as.character(PH$Habit[match(x$Selection, PH$Selection)]))), Trait = paste("Low-chill Exp", x$Exp, "(% budbreak, plot)"), Value = as.numeric(x[["% budbreak"]])) })
HP <- rbindlist(list(melt(fld, id.vars = c("Genotype", "Habit"), variable.name = "Trait", value.name = "Value")[, Study := "Field earliness"], if (!is.null(PL)) PL[, Study := "Low-chill forcing"]), use.names = TRUE, fill = TRUE); HP[, Trait := as.character(Trait)]
bl <- data.table(Genotype = PH$Selection, Habit = toupper(trimws(as.character(PH$Habit))), `Field 10% BLUE 2025-26 (DOY)` = suppressWarnings(as.numeric(PH[["10% 2025-26 pooled, early-3 at day 45"]])), `Field 90% BLUE 2025-26 (DOY)` = suppressWarnings(as.numeric(PH[["90% 2025-26 pooled"]])))
for (e in FIG_EXPS) set(bl, j = paste0("Low-chill Exp ", e, " BLUE (%)"), value = suppressWarnings(as.numeric(PH[[paste0("BLUE Exp ", e)]]))); set(bl, j = "Low-chill ALL four experiments BLUE (%)", value = suppressWarnings(as.numeric(PH[["BLUE both years (all 4 experiments)"]])))
HB <- melt(bl, id.vars = c("Genotype", "Habit"), variable.name = "Trait", value.name = "Value")[, Trait := as.character(Trait)]; HB[, Study := ifelse(grepl("^Field", Trait), "Field earliness", "Low-chill forcing")]
HP <- HP[is.finite(Value) & Habit %in% c("FF", "PF")]; HB <- HB[is.finite(Value) & Habit %in% c("FF", "PF")]
HSUM <- rbind(HP[, .(Kind = "Phenotype (raw)", n = .N, mean = round(mean(Value), 2), sd = round(sd(Value), 2), min = round(min(Value), 2), max = round(max(Value), 2)), by = .(Study, Trait, Habit)],
              HB[, .(Kind = "BLUE (no lug corr.)", n = .N, mean = round(mean(Value), 2), sd = round(sd(Value), 2), min = round(min(Value), 2), max = round(max(Value), 2)), by = .(Study, Trait, Habit)])
writexl::write_xlsx(list(Membership_pct_coordinates = MEM, Cluster_by_habit = CxH, Cluster_summary = CS, Phenotype_BLUE_by_habit = HSUM), file.path(DAPC_DIR, "DAPC_K2_tables.xlsx"))
fwrite(MEM, file.path(DAPC_DIR, "DAPC_K2_membership_percent_and_coordinates.csv")); fwrite(CxH, file.path(DAPC_DIR, "DAPC_K2_cluster_by_habit.csv")); fwrite(HSUM, file.path(DAPC_DIR, "Phenotype_BLUE_summary_by_habit.csv"))
## --- DAPC plots ---
safe("DAPC plots", { md <- copy(MEM); md[, Cl := factor(Assigned_cluster)]; md[, Hab := factor(ifelse(is.na(Habit), "Unknown", Habit), levels = c("FF", "PF", "Unknown"))]
  hs <- scale_shape_manual(values = c(FF = 22, PF = 21, Unknown = 4)); hf <- scale_fill_manual(values = c(HAB_COL, Unknown = "grey60"))
  pa <- ggplot(md, aes(PC1, PC2, fill = Hab, shape = Hab)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") + geom_point(size = 2.6, colour = "black", stroke = 0.4, alpha = 0.85) + hs + hf +
    labs(x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("PC2 (%.1f%%)", pve[2]), fill = NULL, shape = NULL) + th_pub + theme(legend.position = c(0.92, 0.88), legend.background = element_blank())
  cen <- md[, .(LD1 = mean(LD1), PC1 = mean(PC1), PC2 = mean(PC2)), by = Cl]
  pb <- ggplot(md, aes(LD1, PC1, fill = Hab, shape = Hab)) + stat_ellipse(type = "norm", level = 0.95, linewidth = 1, inherit.aes = FALSE, data = md, mapping = aes(LD1, PC1, group = Cl, colour = Cl)) +
    geom_point(size = 2.6, colour = "black", stroke = 0.4, alpha = 0.85) + hs + hf + scale_colour_manual(values = CL_COL, guide = "none") + geom_label(data = cen, aes(LD1, PC1, label = Cl), inherit.aes = FALSE, fontface = "bold", size = 5, label.r = unit(0.2, "lines")) +
    labs(x = "LD1", y = sprintf("PC1 (%.1f%%)", pve[1]), fill = NULL, shape = NULL) + th_pub + theme(legend.position = "none")
  pc_ <- ggplot(md, aes(LD1, fill = Cl)) + geom_histogram(bins = 30, colour = "black", linewidth = 0.2, alpha = 0.85, position = "identity") + scale_fill_manual(values = CL_COL, name = "DAPC class") + labs(x = "LD1", y = "Genotypes") + th_pub + theme(legend.position = c(0.15, 0.85))
  sv(pa, "DAPC_K2_PCA_PC1_PC2", DAPC_DIR, 6.5, 5); sv(pb, "DAPC_K2_LD1_vs_PC1", DAPC_DIR, 6.5, 5); sv(pc_, "DAPC_K2_LD1_histogram", DAPC_DIR, 6.5, 4); tag_comb(list(pa, pb, pc_), 2, "DAPC_K2_Figure_combined", DAPC_DIR, 13, 9.5)
  ## percent membership of every genotype (structure-style bars), ordered within habit
  mb <- melt(md[, .(Genotype, Hab, Assigned_cluster, Pct_class_1, Pct_class_2)], id.vars = c("Genotype", "Hab", "Assigned_cluster"), variable.name = "Class", value.name = "Pct"); mb[, Class := factor(sub("Pct_class_", "", Class), levels = c("1", "2"))]
  ordg <- md[order(Hab, Pct_class_2), Genotype]; mb[, Genotype := factor(Genotype, levels = ordg)]
  pm <- ggplot(mb, aes(Genotype, Pct, fill = Class)) + geom_col(width = 1) + facet_grid(. ~ Hab, scales = "free_x", space = "free_x") + scale_fill_manual(values = CL_COL, name = "DAPC class") + scale_y_continuous(expand = c(0, 0)) + labs(x = "Genotype (ordered by % class 2 within habit)", y = "Membership (%)") +
    th_pub + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing = unit(0.15, "lines"))
  sv(pm, "DAPC_K2_percent_membership_by_genotype", DAPC_DIR, 11, 3.8)
  ph_ <- ggplot(CxH, aes(factor(Assigned_cluster), Pct_of_cluster, fill = Habit)) + geom_col(colour = "black", width = 0.6) + geom_text(aes(label = sprintf("%d (%.0f%%)", N, Pct_of_cluster)), position = position_stack(vjust = 0.5), colour = "white", fontface = "bold", size = 3.6) +
    scale_fill_manual(values = HAB_COL, name = NULL) + labs(x = "DAPC class", y = "Genotypes in class (%)") + th_pub
  sv(ph_, "DAPC_K2_habit_composition_of_classes", DAPC_DIR, 4.6, 4.4) })
## --- distribution histograms: FF blue / PF red, one figure per study x (phenotype, BLUE) ---
hist_fig <- function(d, ttl, xl) { d <- copy(d); d[, Trait := factor(Trait, levels = unique(Trait))]; mu <- d[, .(m = mean(Value)), by = .(Trait, Habit)]
  ggplot(d, aes(Value, fill = Habit)) + geom_histogram(bins = 25, position = "identity", alpha = 0.55, colour = "black", linewidth = 0.2) + geom_vline(data = mu, aes(xintercept = m, colour = Habit), linetype = 2, linewidth = 0.7, show.legend = FALSE) +
    facet_wrap(~ Trait, scales = "free", ncol = 2) + scale_fill_manual(values = HAB_COL, name = NULL) + scale_colour_manual(values = HAB_COL) + labs(x = xl, y = "Count") + th_pub + theme(legend.position = "top") }
safe("histograms", { for (st in c("Field earliness", "Low-chill forcing")) { tag <- if (st == "Field earliness") "Field" else "LowChill"
    a <- HP[Study == st]; b <- HB[Study == st]; nr <- function(d) max(1, ceiling(uniqueN(d$Trait) / 2))
    if (nrow(a)) sv(hist_fig(a, st, if (st == "Field earliness") "Day of year (field budbreak)" else "Budbreak (%)"), paste0("Histogram_", tag, "_phenotypes_FFblue_PFred"), DAPC_DIR, 9, 3 * nr(a) + 0.6)
    if (nrow(b)) sv(hist_fig(b, st, if (st == "Field earliness") "BLUE, day of year" else "BLUE, budbreak (%)"), paste0("Histogram_", tag, "_BLUEs_FFblue_PFred"), DAPC_DIR, 9, 3 * nr(b) + 0.6) } })
## --- correlation matrices (diagonal density, Pearson r above, scatter below; stars as in Godwin et al.) ---
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.1, ".", ""))))
corr_fig <- function(M, file_base) { M <- as.data.frame(M); nm <- names(M); n <- length(nm); out <- list()
  for (fmt in c("png", "pdf")) { f <- file.path(FIG_DIR, paste0(file_base, ".", fmt)); if (fmt == "png") png(f, width = 2400, height = 2200, res = 300) else pdf(f, width = 8, height = 7.3)
    par(mfrow = c(n, n), mar = c(0.25, 0.25, 0.25, 0.25), oma = c(3.2, 3.6, 3.2, 3.2), mgp = c(1.6, 0.5, 0), tcl = -0.25, cex.axis = 0.8)
    for (i in seq_len(n)) for (j in seq_len(n)) { xi <- M[[j]]; yi <- M[[i]]
      if (i == j) { v <- xi[is.finite(xi)]; dn <- density(v); plot(dn, main = "", xlab = "", ylab = "", axes = FALSE, bty = "o"); box()
      } else if (j > i) { ok <- is.finite(xi) & is.finite(yi); ct <- cor.test(xi[ok], yi[ok]); plot.new(); box(); text(0.5, 0.64, "Corr:", col = "grey40", cex = 1.1); text(0.5, 0.40, sprintf("%.3f%s", ct$estimate, stars(ct$p.value)), col = "grey40", cex = 1.3)
        text(0.5, 0.13, sprintf("n = %d", sum(ok)), col = "grey55", cex = 0.8)
        if (fmt == "png") out[[length(out) + 1]] <- data.table(Trait_1 = nm[i], Trait_2 = nm[j], n = sum(ok), r = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3))
      } else { ok <- is.finite(xi) & is.finite(yi); plot(xi[ok], yi[ok], pch = 16, cex = 0.55, col = rgb(0, 0, 0, 0.65), axes = FALSE, xlab = "", ylab = ""); box() }
      if (i == n && !(i == j && FALSE)) axis(1, labels = TRUE) else if (i == j) NULL
      if (j == 1 && i != 1) axis(2, labels = TRUE)
      usr <- par("usr"); if (i == 1) { hh <- diff(usr[3:4]) * 0.17; rect(usr[1], usr[4], usr[2], usr[4] + hh, col = "grey90", border = "black", xpd = NA); text(mean(usr[1:2]), usr[4] + hh / 2, nm[j], font = 2, cex = 0.62, xpd = NA) }
      if (j == n) { ww <- diff(usr[1:2]) * 0.17; rect(usr[2], usr[3], usr[2] + ww, usr[4], col = "grey90", border = "black", xpd = NA); text(usr[2] + ww / 2, mean(usr[3:4]), nm[i], srt = 270, font = 2, cex = 0.62, xpd = NA) } }
    invisible(dev.off()) }
  rbindlist(out) }
safe("correlation figures", { LCs <- c(paste("Exp", FIG_EXPS), "All four")
  M1 <- data.frame(sapply(FIG_EXPS, function(e) suppressWarnings(as.numeric(PH[[paste0("BLUE Exp ", e)]]))), `All four` = suppressWarnings(as.numeric(PH[["BLUE both years (all 4 experiments)"]])), check.names = FALSE); names(M1) <- LCs
  M2 <- data.frame(`Field 10% 2025` = fld[[3]], `Field 10% 2026` = fld[[4]], `Field 10% pooled` = suppressWarnings(as.numeric(PH[["10% 2025-26 pooled, early-3 at day 45"]])), `Low-chill all four BLUE` = suppressWarnings(as.numeric(PH[["BLUE both years (all 4 experiments)"]])), check.names = FALSE)
  M3 <- data.frame(M2[, 3, drop = FALSE], M1, check.names = FALSE); names(M3)[1] <- "Field 10% pooled"
  R1 <- corr_fig(M1, "Fig_Correlation_LowChill_BLUEs_noLug"); R2 <- corr_fig(M2, "Fig_Correlation_Field_vs_LowChill"); R3 <- corr_fig(M3, "Fig_Correlation_Field10_vs_each_LowChill_experiment")
  M4 <- data.frame(`Field 10% 2025` = fld[[3]], `Field 10% 2026` = fld[[4]], `Field 10% pooled` = suppressWarnings(as.numeric(PH[["10% 2025-26 pooled, early-3 at day 45"]])),
                   `Field 90% 2025` = fld[[5]], `Field 90% 2026` = fld[[6]], `Field 90% pooled` = suppressWarnings(as.numeric(PH[["90% 2025-26 pooled"]])), check.names = FALSE); R4 <- corr_fig(M4, "Fig_Correlation_Field_10_vs_90_budbreak")
  writexl::write_xlsx(list(LowChill_BLUEs = R1, Field_vs_LowChill = R2, Field10_vs_experiments = R3, Field_10_vs_90 = R4), file.path(FIG_DIR, "Correlation_tables.xlsx")) })
## --- gene density and marker density on all chromosomes (1 Mb bins, colour classes as in Chizk et al. / Alexander Supp. Fig. 1) ---
bin_counts <- function(pos, chrom, lens, bin = 1e6) rbindlist(lapply(REAL_CHROMS, function(ch) { nb <- ceiling(lens[[ch]] / bin); cnt <- tabulate(floor(pos[chrom == ch] / bin) + 1, nbins = nb); data.table(Chrom = ch, start = (seq_len(nb) - 1) * bin, end = pmin(seq_len(nb) * bin, lens[[ch]]), n = cnt) }))
dens_fig <- function(B, ttl, unit) { nz <- B$n[B$n > 0]; stp <- max(1, ceiling(as.numeric(quantile(nz, 0.99)) / 8 / 5) * 5); br <- c(0, stp * seq_len(8)); lab <- c("0", sprintf("%d-%d", c(1, br[2:7] + 1), br[2:8]), sprintf(">= %d", br[8] + 1))
  B[, cls := cut(n, c(-Inf, 0, br[2:8], Inf), labels = lab)]; pal <- setNames(c("#F0F0F0", "#3F8C0E", "#7FB000", "#C4DB00", "#FFFF00", "#FFC000", "#FF8000", "#FF4500", "#FF0000"), lab); B[, y := match(Chrom, rev(REAL_CHROMS))]
  ggplot(B) + geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = y - 0.32, ymax = y + 0.32, fill = cls)) + scale_fill_manual(values = pal, name = unit, drop = FALSE) + scale_y_continuous(breaks = seq_along(REAL_CHROMS), labels = rev(REAL_CHROMS), expand = c(0.02, 0.02)) +
    scale_x_continuous(position = "top", breaks = seq(0, 50, 5), labels = function(x) paste0(x, "Mb"), expand = c(0.005, 0)) + labs(x = NULL, y = NULL, title = ttl) + theme_minimal(base_size = 12) + theme(panel.grid = element_blank(), axis.line.x = element_line(colour = "black"), axis.ticks.x = element_line(colour = "black")) }
safe("gene and marker density", { gf <- find_file(GENE_FILE); mpos <- data.table(Chrom = as.character(Gd$Chrom), Position = as.numeric(Gd$Position)); allpos <- data.table(Chrom = as.character(G$Chrom), Position = as.numeric(G$Position)); lens <- setNames(as.list(rep(0, 7)), REAL_CHROMS)
  GN <- if (is.na(gf)) NULL else { g <- as.data.table(read_excel(gf, skip = 2, .name_repair = "minimal")); g <- g[, .(Chrom = as.character(Chromosome), s = as.numeric(`Start position`), e = as.numeric(`End position`))][Chrom %in% REAL_CHROMS & is.finite(s) & is.finite(e)]; g[, mid := (s + e) / 2]; g }
  if (is.null(GN)) cat("[gene annotation file not found:", GENE_FILE, "- gene density skipped]\n")
  for (ch in REAL_CHROMS) lens[[ch]] <- max(c(mpos[Chrom == ch]$Position, allpos[Chrom == ch]$Position, if (!is.null(GN)) GN[Chrom == ch]$e), na.rm = TRUE)
  BM <- bin_counts(mpos$Position, mpos$Chrom, lens); BA <- bin_counts(allpos$Position, allpos$Chrom, lens); BG <- if (!is.null(GN)) bin_counts(GN$mid, GN$Chrom, lens)
  pM <- dens_fig(copy(BM), sprintf("Distribution of %s GWAS markers across the seven chromosomes (SNPs per Mb)", format(nrow(Gd), big.mark = ",")), "Markers per Mb"); sv(pM, "Fig_Marker_density_all_chromosomes", FIG_DIR, 10, 4.6)
  sv(dens_fig(copy(BA), sprintf("Distribution of all %s shared markers before the concordance filter (per Mb)", format(nrow(allpos), big.mark = ",")), "Markers per Mb"), "Fig_Marker_density_all_shared_markers", FIG_DIR, 10, 4.6)
  if (!is.null(BG)) { pG <- dens_fig(copy(BG), sprintf("Gene density across the seven chromosomes (%s predicted genes per Mb)", format(nrow(GN), big.mark = ",")), "Genes per Mb"); sv(pG, "Fig_Gene_density_all_chromosomes", FIG_DIR, 10, 4.6)
    if (HAVE_PW) sv(patchwork::wrap_plots(list(pG, pM), ncol = 1) + patchwork::plot_annotation(tag_levels = "a"), "Fig_Gene_and_marker_density_combined", FIG_DIR, 10, 8.6)
    DN <- merge(BG[, .(Chrom, start, genes = n)], BM[, .(Chrom, start, markers = n)], by = c("Chrom", "start")); cat(sprintf("Genes vs markers per Mb: Pearson r = %.3f\n", cor(DN$genes, DN$markers))); fwrite(DN, file.path(FIG_DIR, "Gene_and_marker_counts_per_Mb.csv")) } })
cat("\nFigures and tables written to:\n  ", DAPC_DIR, "\n  ", FIG_DIR, "\n")

## ---- 4. GWAS per run: Q + K AND kinship only ----
fit_run <- function(rid, with_q) { mods <- if (is_field(rid)) c("additive", "1-dom") else MODELS_ALL   # general and 2-dom: low-chill runs only
  ok <- is.finite(ph_q[[rid]]); smp <- ph_q$Genotype[ok]; yv <- ph_q[[rid]][ok]
  pf <- tempfile(fileext = ".csv"); tmp_files <<- c(tmp_files, pf); pdt <- data.table(Genotype = smp); pdt[, (rid) := yv]; if (with_q) pdt <- cbind(pdt, ph_q[ok, ..Qn]); fwrite(pdt, pf, na = "NA")
  dat <- quiet(read.GWASpoly(ploidy = 4, pheno.file = pf, geno.file = geno_out, format = "numeric", n.traits = 1, delim = ",")); dat <- quiet(set.K(dat, LOCO = TRUE))
  prm <- if (with_q) set.params(MAF = MIN_MAF, geno.freq = MAX_GENO_FREQ, fixed = Qn, fixed.type = "numeric") else set.params(MAF = MIN_MAF, geno.freq = MAX_GENO_FREQ)
  fit <- quiet(GWASpoly(data = dat, models = mods, traits = rid, params = prm, n.core = N_CORES)); thr <- quiet(set.threshold(fit, method = "M.eff", level = THRESH_LEVEL))
  mp <- as.data.frame(fit@map); names(mp)[1:3] <- c("Marker", "Chrom", "Position"); mp$Chrom <- as.character(mp$Chrom); mp$Position <- as.numeric(mp$Position)
  th <- thr@threshold; th <- if (is.matrix(th)) th[1, ] else th; raw <- as.data.frame(fit@scores[[rid]])
  Mt <- Mg[match(mp$Marker, Gd$Marker), match(smp, samples), drop = FALSE]; nn <- rowSums(!is.na(Mt)); n0 <- rowSums(Mt == 0, na.rm = TRUE); n4 <- rowSums(Mt == 4, na.rm = TRUE)
  n1 <- rowSums(Mt <= 1, na.rm = TRUE); n2 <- rowSums(Mt <= 2, na.rm = TRUE)
  cls <- list("1-dom-alt" = pmin(n0, nn - n0), "1-dom-ref" = pmin(n4, nn - n4), "2-dom-alt" = pmin(n1, nn - n1), "2-dom-ref" = pmin(n2, nn - n2)); dl <- list(); sm <- list()
  for (mo in intersect(MODEL_ORDER, names(raw))) { v <- raw[[mo]]; if (mo %in% names(cls)) v[cls[[mo]] < MIN_CLASS_N] <- NA_real_; ib <- if (any(is.finite(v))) which.max(v) else NA_integer_
    dl[[mo]] <- data.table(Marker = mp$Marker, Chrom = mp$Chrom, Position = mp$Position, score = v, model = mo, thr = th[[mo]])
    sm[[mo]] <- data.table(Run = rid, Model = mo, n = length(smp), Threshold = round(th[[mo]], 2), N_above = sum(v >= th[[mo]], na.rm = TRUE), Top_score = round(v[ib], 2), Top_marker = mp$Marker[ib],
                           Chrom = mp$Chrom[ib], Position = mp$Position[ib], Lambda = round(lambda_gc(v), 3)) }
  list(df = rbindlist(dl)[is.finite(score)], summ = rbindlist(sm)) }
RES <- list(); SUMM <- list(); RESK <- list(); SUMMK <- list(); t0 <- Sys.time()
for (rid in RUNS$Run) { cat(sprintf("\n>>> %-11s n = %d ... Q+K (4 models) ", rid, sum(is.finite(ph_q[[rid]]))))
  o <- tryCatch(fit_run(rid, TRUE), error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }); if (!is.null(o)) { RES[[rid]] <- o$df; SUMM[[rid]] <- o$summ; cat(sprintf("done (%.1f min so far) ", as.numeric(difftime(Sys.time(), t0, units = "mins")))) }
  if (RUN_K_ONLY) { o2 <- tryCatch(fit_run(rid, FALSE), error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }); if (!is.null(o2)) { RESK[[rid]] <- o2$df; SUMMK[[rid]] <- o2$summ; cat("K only done") } }
  invisible(gc()) }
cat("\n"); if (!length(SUMM)) stop("No run completed.")
S_FULL <- rbindlist(SUMM); S_FULL <- merge(S_FULL, RUNS[, .(Run, Description = Label)], by = "Run"); setcolorder(S_FULL, c("Run", "Description"))
SK_FULL <- if (length(SUMMK)) rbindlist(SUMMK) else NULL; RES_FULL <- RES; RESK_FULL <- RESK

for (SET in names(SETS)) {
MS <- SETS[[SET]]; OUT_DIR <- OUT_DIRS[[SET]]; dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); unlink(list.files(OUT_DIR, "\\.(png|xlsx|pptx|rds)$", full.names = TRUE))
RUNS_SET <- if (SET == "alt") RUNS$Run[!is_field(RUNS$Run)] else RUNS$Run
RES <- lapply(RES_FULL[intersect(RUNS_SET, names(RES_FULL))], function(d) d[model %in% MS]); RESK <- lapply(RESK_FULL[intersect(RUNS_SET, names(RESK_FULL))], function(d) d[model %in% MS])
S <- S_FULL[Model %in% MS & Run %in% RUNS_SET][order(match(Run, RUNS$Run), match(Model, MODEL_ORDER))]; SK <- if (is.null(SK_FULL)) NULL else SK_FULL[Model %in% MS & Run %in% RUNS_SET][order(match(Run, RUNS$Run), match(Model, MODEL_ORDER))]
cat(sprintf("\n######## SET '%s': models %s -> %s\n", SET, paste(MS, collapse = ", "), OUT_DIR))
saveRDS(list(S = S, SK = SK, RES = RES, RESK = RESK, ph_q = ph_q, grp = grp, P = P, samples = samples, RUNS = RUNS), file.path(OUT_DIR, "CH2_Final_state.rds"))

## ---- 5. significant loci (Q + K) ----
LOC <- rbindlist(lapply(names(RES), function(rid) { d <- RES[[rid]][score >= thr][order(Chrom, Position)]; if (!nrow(d)) return(NULL)
  d[, locus := cumsum(c(TRUE, diff(Position) > LOCUS_WINDOW | Chrom[-1] != Chrom[-.N])), by = .(model)][, Run := rid]
  d[, {i <- which.max(score); .(Start_bp = min(Position), End_bp = max(Position), Markers = .N, Lead_marker = Marker[i], Lead_position = Position[i], Lead_score = round(score[i], 2))}, by = .(Run, model, Chrom, locus)][, locus := NULL] }), fill = TRUE)
if (!is.null(LOC) && nrow(LOC)) { setnames(LOC, "model", "Model"); LOC <- LOC[order(match(Run, RUNS$Run), Chrom, Start_bp)] } else LOC <- data.table(Note = "No marker above its model threshold in any run")

## ---- 5b. diagnostics ----
KVQ <- if (is.null(SK)) NULL else safe("Q+K vs K only", merge(S[, .(Run, Model, n, Threshold, QK_top = Top_score, QK_N_above = N_above, QK_lambda = Lambda)],
  SK[, .(Run, Model, K_top = Top_score, K_N_above = N_above, K_lambda = Lambda, K_top_marker = Top_marker)], by = c("Run", "Model"), all.x = TRUE)[order(match(Run, RUNS$Run), match(Model, MODEL_ORDER))])
TC <- safe("top score by chromosome", { tc <- function(L, fit) rbindlist(lapply(names(L), function(r) L[[r]][, .(Top = round(max(score), 1)), by = .(Chrom, Model = model)][, `:=`(Run = r, Fit = fit)]))
  W <- dcast(rbind(tc(RES, "Q+K"), if (length(RESK)) tc(RESK, "K only")), Run + Model + Fit ~ Chrom, value.var = "Top"); W[order(match(Run, RUNS$Run), match(Model, MODEL_ORDER), match(Fit, c("Q+K", "K only")))] })
PQ <- safe("phenotype by cluster", rbindlist(lapply(RUNS$Run, function(r) { ok <- is.finite(ph_q[[r]]); y <- ph_q[[r]][ok]; g <- as.integer(as.character(grp))[ok]
  data.table(Run = r, n_cl1 = sum(g == 1), n_cl2 = sum(g == 2), mean_cl1 = round(mean(y[g == 1]), 2), mean_cl2 = round(mean(y[g == 2]), 2), R2_phenotype_on_Q = round(summary(lm(y ~ ph_q$Cluster_2[ok]))$r.squared, 3)) })))
KM <- safe("key marker classes", rbindlist(lapply(KEY_MARKERS, function(p) { mi <- grep(p, Gd$Marker)[1]; if (is.na(mi)) { cat("[key marker not in the GWAS marker set:", p, "]\n"); return(NULL) }; dos <- Mg[mi, ]
  rbindlist(lapply(RUNS$Run, function(r) { y <- ph_q[[r]]; o <- list(Marker = Gd$Marker[mi], Run = r); for (d in 0:4) { ii <- which(dos == d & is.finite(y)); o[[paste0("n", d)]] <- length(ii); o[[paste0("mean", d)]] <- if (length(ii)) round(mean(y[ii]), 2) else NA_real_ }; as.data.table(o) })) }), fill = TRUE))
KC <- safe("key marker by cluster", rbindlist(lapply(KEY_MARKERS, function(p) { mi <- grep(p, Gd$Marker)[1]; if (is.na(mi)) return(NULL); dos <- Mg[mi, ]; g <- as.integer(as.character(grp))
  data.table(Marker = Gd$Marker[mi], mean_dosage_cl1 = round(mean(dos[g == 1], na.rm = TRUE), 2), mean_dosage_cl2 = round(mean(dos[g == 2], na.rm = TRUE), 2)) })))
KMEM <- safe("members of the dosage-4 class", { mi <- grep(KEY_MARKERS[1], Gd$Marker)[1]; if (is.na(mi)) NULL else { dos <- Mg[mi, ]; v <- intersect(c("LC_Y1", "LC_Both", "Exp1.1", "Exp1.2", "Exp2.1", "Exp2.2", "BB10_pooled", "BB90_pooled"), names(ph_q))
  d <- data.table(Marker = Gd$Marker[mi], Genotype = samples, DAPC_cluster = as.integer(as.character(grp)), Dosage = dos); d <- cbind(d, round(ph_q[, ..v], 2)); d <- d[Dosage == 4][order(-LC_Both)]; head(d, 60) } })
XT <- safe("field vs low-chill overlaps", { sg <- rbindlist(lapply(names(RES), function(r) RES[[r]][score >= thr][, Run := r])); fs <- sg[grepl("^BB", Run)]; ls <- sg[!grepl("^BB", Run)]
  if (!nrow(fs) || !nrow(ls)) NULL else { X <- merge(ls[, .(LC_Run = Run, LC_Model = model, Chrom, LC_marker = Marker, LC_pos = Position, LC_score = round(score, 2))],
      fs[, .(Field_Run = Run, Field_Model = model, Chrom, Field_marker = Marker, Field_pos = Position, Field_score = round(score, 2))], by = "Chrom", allow.cartesian = TRUE)[abs(LC_pos - Field_pos) <= OVERLAP_WINDOW]
    if (!nrow(X)) NULL else X[, .(LC_markers = uniqueN(LC_marker), LC_lead = LC_marker[which.max(LC_score)], LC_lead_score = max(LC_score), Field_lead = Field_marker[which.max(Field_score)], Field_lead_score = max(Field_score),
      Field_lead_model = Field_Model[which.max(Field_score)], Min_distance_bp = min(abs(LC_pos - Field_pos))), by = .(Chrom, LC_Run, LC_Model, Field_Run)][order(Chrom, LC_Run, LC_Model)] } })
EXC <- safe("experiment correlations", { v <- intersect(c("Exp1.1", "Exp1.2", "Exp2.1", "Exp2.2", "LC_Both", "LC_Y1"), names(ph)); m <- round(cor(as.matrix(ph[, ..v]), use = "pairwise.complete.obs"), 2); data.table(Run = rownames(m), as.data.table(m)) })

## ---- 6. Manhattan + QQ plots (style of the earlier K = 2 figures: shaded chromosomes, blue circles = additive, yellow triangles = 1-dom) ----
Rf <- RES[intersect(FIG_RUNS, names(RES))]; allp <- rbindlist(Rf); lens <- allp[, .(len = max(Position)), by = Chrom][match(REAL_CHROMS, Chrom)]
offs <- c(0, cumsum(lens$len))[seq_along(REAL_CHROMS)]; names(offs) <- REAL_CHROMS; mids <- offs + lens$len / 2
MC <- c(additive = "#4A9AD9", general = "#3BAA6E", "1-dom-alt" = "#F5B82E", "1-dom-ref" = "#F5B82E", "2-dom-alt" = "#C0504D", "2-dom-ref" = "#C0504D")
MP <- c(additive = 16, general = 15, "1-dom-alt" = 17, "1-dom-ref" = 17, "2-dom-alt" = 18, "2-dom-ref" = 18); MPC <- c(additive = 0.45, general = 0.45, "1-dom-alt" = 0.55, "1-dom-ref" = 0.55, "2-dom-alt" = 0.6, "2-dom-ref" = 0.6)
MAB <- c(additive = "add", general = "gen", "1-dom-alt" = "1d-alt", "1-dom-ref" = "1d-ref", "2-dom-alt" = "2d-alt", "2-dom-ref" = "2d-ref")
LEGLAB <- c(additive = "additive", general = "general", "1-dom" = "1-dom (alt and ref)", "2-dom" = "2-dom (alt and ref)"); fam <- function(m) sub("-(alt|ref)$", "", m)
short_lab <- c(BB10_pooled = "Field 10% budbreak 2025-26 pooled", BB90_pooled = "Field 90% budbreak 2025-26 pooled", Exp1.1 = "Low-chill Exp 1.1", Exp1.2 = "Low-chill Exp 1.2", Exp2.1 = "Low-chill Exp 2.1", Exp2.2 = "Low-chill Exp 2.2", LC_Both = "Low-chill ALL four experiments pooled")
png_dir <- function(rid) if (is_field(rid)) PNG_DIRS$main_field else if (SET == "alt") PNG_DIRS$alt_lc else PNG_DIRS$main_lc
for (rid in names(Rf)) { d <- copy(Rf[[rid]]); d[, x := Position + offs[Chrom]]; nn <- S[Run == rid]$n[1]; pd <- png_dir(rid); dir.create(pd, recursive = TRUE, showWarnings = FALSE)
  ttl <- sprintf("%s  |  n = %d  |  Q + K (DAPC K = 2)  |  DAPC K = 2, platform-discordant markers removed", if (rid %in% names(short_lab)) short_lab[[rid]] else RUNS[Run == rid]$Label, nn)
  ss <- S[Run == rid][match(intersect(MS, S[Run == rid]$Model), Model)]; sub_ <- paste0("markers above threshold: ", paste(sprintf("%s %d (>= %.2f)", ss$Model, ss$N_above, ss$Threshold), collapse = "  |  "))
  ymax <- max(d$score, na.rm = TRUE) * 1.07; ytk <- pretty(c(0, ymax)); ytk <- ytk[ytk <= ymax]
  png(file.path(pd, sprintf("Manhattan_%s.png", rid)), width = 2050, height = 1100, res = 150); par(mar = c(5, 5, 5.2, 1.2), bty = "n")
  plot(NA, xlim = c(0, max(offs + lens$len)), ylim = c(-ymax * 0.02, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  for (k in seq(2, length(REAL_CHROMS), by = 2)) rect(offs[k], par("usr")[3], offs[k] + lens$len[k], par("usr")[4], col = "#F0F0F0", border = NA)
  for (m in intersect(MODEL_ORDER, MS)) { dm <- d[model == m]; points(dm$x, dm$score, pch = MP[[m]], cex = MPC[[m]], col = MC[[m]]) }
  thm <- unique(d[, .(model, thr)]); for (i in seq_len(nrow(thm))) abline(h = thm$thr[i], lty = 2, lwd = 1, col = "grey35")
  axis(2, at = ytk, las = 1, cex.axis = 1.25, lwd = 0, lwd.ticks = 1); segments(par("usr")[1], par("usr")[3], par("usr")[1], par("usr")[4]); segments(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[3])
  axis(1, at = mids, labels = names(offs), tick = FALSE, cex.axis = 1.3, line = 0.6); mtext("Chromosome", 1, line = 3.6, cex = 1.3); mtext(expression(-log[10](italic(p))), 2, line = 3, cex = 1.4)
  mtext(ttl, 3, line = 2.7, adj = 0, font = 2, cex = 1.05); mtext(sub_, 3, line = 1.1, adj = 0, cex = 0.92, col = "grey30")
  fm <- unique(fam(intersect(MODEL_ORDER, MS))); fk <- ifelse(fm %in% c("1-dom", "2-dom"), paste0(fm, "-alt"), fm); legend("topright", bty = "n", cex = 1.3, legend = LEGLAB[fm], col = MC[fk], pch = MP[fk], pt.cex = 1.2); invisible(dev.off())
  png(file.path(pd, sprintf("QQ_%s_by_chromosome.png", rid)), width = 2050, height = 1090, res = 150); par(mfrow = c(2, 4), mar = c(4.6, 4.6, 2.2, 0.8), oma = c(0, 0, 2.6, 0), bty = "n")
  for (pn in c("All chromosomes", REAL_CHROMS)) { dd <- if (pn == "All chromosomes") d else d[Chrom == pn]; qq <- list(); mxx <- 0; mxy <- 0
    for (m in intersect(MODEL_ORDER, MS)) { v <- sort(dd[model == m]$score, decreasing = TRUE); if (!length(v)) next; e <- -log10(ppoints(length(v))); mxx <- max(mxx, e); mxy <- max(mxy, v); qq[[m]] <- list(e = e, v = v, l = lambda_gc(v)) }
    plot(NA, xlim = c(-0.05 * mxx, mxx * 1.03), ylim = c(-0.05 * mxy, mxy * 1.04), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = pn, cex.main = 1.4, font.main = 2, las = 1, cex.lab = 1.1, cex.axis = 1.05, xaxs = "i", yaxs = "i"); box(bty = "l", lwd = 1)
    abline(0, 1, col = "grey40", lty = 2); for (m in names(qq)) points(qq[[m]]$e, qq[[m]]$v, pch = MP[[m]], cex = MPC[[m]] * 0.9, col = MC[[m]])
    legend("topleft", bty = "n", cex = 0.78, y.intersp = 0.95, inset = c(0.02, 0.01), legend = sprintf("%-5s lambda %.2f", MAB[names(qq)], sapply(qq, function(q) q$l))) }
  mtext(ttl, 3, outer = TRUE, adj = 0, font = 2, cex = 0.95, line = 0.4); invisible(dev.off()) }

## ---- 7. workbook (ONE file) ----
PU <- ph[Genotype %in% samples, c("Genotype", RUNS$Run), with = FALSE]; for (r in RUNS$Run) set(PU, j = r, value = round(PU[[r]], 3)); setnames(PU, RUNS$Run, RUNS$Label)
PU <- PU[rowSums(!is.na(PU[, -1, with = FALSE])) > 0]
DM <- data.table(Genotype = samples, Assigned_cluster = as.integer(as.character(grp)), as.data.table(round(P, 4)))
NT <- data.table(Sheet = c("Summary", "Significant_loci", "QK_vs_Konly", "Top_by_chromosome", "Pheno_by_cluster", "Key_marker_classes", "Key_marker_by_cluster", "Key_marker_dosage4_members", "Field_vs_lowchill", "Exp_correlation", "Phenotypes_used", "DAPC_membership", "Method", "Exp 2.1 removals"),
  What = c("One row per run x model (Q + K): markers above the M.eff threshold (alpha 0.05), top score = -log10(p), position in bp, lambda.",
           sprintf("Markers above their threshold merged into loci (within %.1f Mb on a chromosome), Q + K; lead marker and score.", LOCUS_WINDOW / 1e6),
           "Same run fitted with Q + K and with kinship only: top score, markers above threshold, lambda.",
           "Highest score on each chromosome, Q + K and K only, per run and model.",
           "Phenotype mean in each DAPC cluster and R-squared of phenotype on Q: how much of each trait is confounded with population structure.",
           "Genotype-class sizes (n0..n4 = dosage 0..4) and mean phenotype per class at the key markers, per run.",
           "Mean allele dosage of the key markers in each DAPC cluster.", "Genotypes carrying dosage 4 at the first key marker (the 1-dom-ref class), with their phenotypes and DAPC cluster.",
           sprintf("Significant field and low-chill markers within %.0f kb of each other (same chromosome).", OVERLAP_WINDOW / 1e3),
           "Correlation among the four experiment BLUEs, the combined BLUE and Year 1.",
           "Trait values that went into GWASpoly, one row per genotype.", "K = 2 DAPC cluster and membership probabilities; Cluster_2 is the Q covariate in every Q + K run.",
           paste0("GWASpoly, ploidy 4, LOCO kinship, P3D, models ", paste(MS, collapse = ", "), ", Q + K, MAF 0.05, geno.freq 0.98, 1-dom/2-dom markers need >= 10 genotypes in the smaller class; K = 2 DAPC on ", N_DAPC_MARKERS, " markers (Margaret: 6069); ", nrow(Gd), " markers in the GWAS. Low-chill BLUEs WITHOUT lug correction."),
           "16 Exp 2.1 plots (APF_268, A_2615T, A_2956TN, A_2957TN, blocks A-D) are left out by the rule 'Year-1 high performer significantly worse in this Year-2 experiment'."))
sheets <- list(Summary = S, Significant_loci = LOC, QK_vs_Konly = KVQ, Top_by_chromosome = TC, Pheno_by_cluster = PQ, Key_marker_classes = KM, Key_marker_by_cluster = KC, Key_marker_dosage4_members = KMEM, Field_vs_lowchill = XT, Exp_correlation = EXC,
               Phenotypes_used = PU, DAPC_membership = DM, Notes = NT); sheets <- sheets[!sapply(sheets, is.null)]
writexl::write_xlsx(sheets, file.path(OUT_DIR, "CH2_Final_GWAS_Results.xlsx"))

## ---- 7b. one pptx with the thesis figures and the summary table ----
safe("pptx", { if (!requireNamespace("officer", quietly = TRUE)) try(install.packages("officer"), silent = TRUE)
  if (requireNamespace("officer", quietly = TRUE)) { library(officer); px <- read_pptx(); mst <- "Office Theme"
    px <- add_slide(px, layout = "Title Only", master = mst); px <- ph_with(px, "Summary of GWAS runs (K = 2 DAPC, Q + K)", ph_location_type("title"))
    st <- as.data.frame(S[Run %in% names(Rf), .(Run, Model, n, Threshold, N_above, Top_score, Chrom, Position, Lambda)]); px <- ph_with(px, st, location = ph_location(left = 0.3, top = 1.3, width = 9.4, height = 5.5))
    for (rid in names(Rf)) {
      px <- add_slide(px, layout = "Title Only", master = mst); px <- ph_with(px, sprintf("%s: Manhattan", RUNS[Run == rid]$Label), ph_location_type("title"))
      px <- ph_with(px, external_img(file.path(png_dir(rid), sprintf("Manhattan_%s.png", rid)), width = 9.2, height = 4.95), ph_location(left = 0.4, top = 1.5, width = 9.2, height = 4.95))
      px <- add_slide(px, layout = "Title Only", master = mst); px <- ph_with(px, sprintf("%s: QQ by chromosome (supplemental)", RUNS[Run == rid]$Label), ph_location_type("title"))
      px <- ph_with(px, external_img(file.path(png_dir(rid), sprintf("QQ_%s_by_chromosome.png", rid)), width = 9.2, height = 4.89), ph_location(left = 0.4, top = 1.5, width = 9.2, height = 4.88)) }
    print(px, target = file.path(OUT_DIR, "CH2_Final_GWAS_Figures.pptx")); cat("Wrote the pptx.\n") } else cat("officer package not available: pptx skipped (the PNGs are in the output folder).\n") })
cat("\nSaved to:", OUT_DIR, "\n")

## ---- 8. print ----
options(width = 250); pr <- function(title, x) { if (!is.null(x) && NROW(x)) { cat("\n--", title, "--\n"); print(as.data.table(x), nrows = 400) } }
cat(sprintf("\n\n################ COPY FROM HERE: set %s (%s), K = 2 DAPC, Q + K, no lug correction ################\n", SET, paste(MS, collapse = ", ")))
cat(sprintf("Markers in the GWAS: %d | markers in the DAPC: %d | genotypes: %d | DAPC K = 2 cluster sizes: %s | low-chill BLUEs: %s | adegenet: %s\n", nrow(Gd), N_DAPC_MARKERS, length(samples), paste(table(grp), collapse = " / "),
            "no lug correction", HAVE_ADEGENET))
pr("Phenotypes", PHS); pr("DAPC marker count by max missing-data cutoff (script uses 0.10; Margaret 6069)", DAPC_COUNTS)
pr("Q + K summary", S[, .(Run, Model, n, Threshold, N_above, Top_score, Chrom, Position, Top_marker, Lambda)]); pr("Significant loci (Q + K)", LOC)
pr("Q + K vs kinship only", KVQ); pr("Top score by chromosome", TC); pr("Phenotype by DAPC cluster", PQ)
pr("Key markers: class sizes (n0..n4 = dosage) and mean phenotype", KM); pr("Key markers: mean dosage by DAPC cluster", KC); pr("First key marker: genotypes in the dosage-4 class (the 1-dom-ref class), highest combined BLUE first", KMEM)
pr(sprintf("Field vs low-chill significant markers within %.0f kb", OVERLAP_WINDOW / 1e3), XT); pr("Correlation among experiments", EXC)
cat("################ END OF COPY ################\n")
}
