# ============================================================
# 05_univariate.R
# Linear mixed models for individual compounds identified by
# SIMPER as the most consistent contributors to swa dissimilarity.
#
# Focus compounds: nonanal and octanal
#   - Highest and second-highest consistency ratios across datasets
#   - Both are C6-C9 aliphatic aldehydes (lipoxygenase pathway)
#   - Directionally consistent: higher in endophyte-absent plants
#     in living tissue across GH and field contexts
#
# Model structure:
#   log(abundance + 0.01) ~ swa + covariates + (1 | pair)
#
#   - Log transformation: emission data are right-skewed with zeros
#   - Pseudocount 0.01 matches the CLR pipeline and prevents log(0)
#   - swa: fixed effect (factor, 0 = absent, 1 = present)
#   - year: fixed effect (only 2 levels; treated as fixed not random)
#   - pair/pair2: random intercept (accounts for within-pair correlation)
#   - collection.type: fixed effect in field living tissue models
#
# p-values via Satterthwaite approximation (lmerTest).
#
# INPUTS (from 01_clr_processing.R):
#   gh_raw_mat, gh_meta
#   field_live_raw_mat, field_live_meta
#   field_sens_raw_mat, field_sens_meta
# ============================================================

library(lme4)
library(lmerTest)

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "01_clr_processing_2.R"))


# ============================================================
# HELPER
# ============================================================

# Pull a single compound from a raw matrix, join metadata,
# and add a log-transformed column ready for modeling.
make_compound_df <- function(raw_mat, meta, compound) {
  if (!compound %in% colnames(raw_mat)) {
    stop("Compound '", compound, "' not found in matrix columns.")
  }
  meta |>
    rownames_to_column("row_id") |>
    mutate(
      abundance     = raw_mat[row_id, compound],
      log_abundance = log(abundance + 0.01),
      swa           = factor(swa, levels = c(0, 1),
                             labels = c("Absent", "Present"))
    )
}

# Fit lmer, print a clean summary, and return the model invisibly.
fit_and_print <- function(formula, data, label) {
  m <- lmer(formula, data = data, REML = TRUE)
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat(label, "\n")
  cat(strrep("=", 60), "\n")
  cat("Formula: "); print(formula)
  cat("\nFixed effects (Satterthwaite p-values):\n")
  print(coef(summary(m)), digits = 3)
  cat("\nRandom effects:\n")
  print(VarCorr(m))
  invisible(m)
}


# ============================================================
# GREENHOUSE
# ============================================================
# Model: log(abundance) ~ swa + year + (1 | pair)
# Year treated as fixed (2 levels only).
# Pair random intercept accounts for genetic/microsite variation
# between matched plant pairs.

for (cmpd in c("nonanal", "octanal")) {

  df <- make_compound_df(gh_raw_mat, gh_meta, cmpd)

  m <- fit_and_print(
    log_abundance ~ swa + year + (1 | pair),
    data  = df,
    label = paste("GH |", cmpd)
  )

  assign(paste0("m_gh_", cmpd), m)
}


# ============================================================
# FIELD — Living tissue (flowers + leaves)
# ============================================================
# Model: log(abundance) ~ swa + collection.type + year + (1 | pair2)
# collection.type accounts for tissue-type baseline differences.

for (cmpd in c("nonanal", "octanal")) {

  df <- make_compound_df(field_live_raw_mat, field_live_meta, cmpd)

  m <- fit_and_print(
    log_abundance ~ swa + collection.type + year + (1 | pair2),
    data  = df,
    label = paste("Field (living) |", cmpd)
  )

  assign(paste0("m_live_", cmpd), m)
}


# ============================================================
# FIELD — Senesced leaves (supplementary)
# ============================================================
# Model: log(abundance) ~ swa + year + (1 | pair2)
# collection.type omitted (single tissue type).

for (cmpd in c("nonanal", "octanal")) {

  df <- make_compound_df(field_sens_raw_mat, field_sens_meta, cmpd)

  m <- fit_and_print(
    log_abundance ~ swa + year + (1 | pair2),
    data  = df,
    label = paste("Field (senesced) |", cmpd)
  )

  assign(paste0("m_sens_", cmpd), m)
}


# ============================================================
# SUMMARY TABLE
# ============================================================
# Extract swa coefficient, SE, t, and p across all models.

extract_swa_row <- function(model, label) {
  s <- coef(summary(model))
  row <- s[grep("Present", rownames(s)), , drop = FALSE]
  data.frame(
    model     = label,
    estimate  = round(row[, "Estimate"], 3),
    se        = round(row[, "Std. Error"], 3),
    t         = round(row[, "t value"], 3),
    p         = round(row[, "Pr(>|t|)"], 4),
    row.names = NULL
  )
}

swa_summary <- bind_rows(
  extract_swa_row(m_gh_nonanal,      "GH | nonanal"),
  extract_swa_row(m_gh_octanal,      "GH | octanal"),
  extract_swa_row(m_live_nonanal,    "Field live | nonanal"),
  extract_swa_row(m_live_octanal,    "Field live | octanal"),
  extract_swa_row(m_sens_nonanal,    "Field senesced | nonanal"),
  extract_swa_row(m_sens_octanal,    "Field senesced | octanal")
)

cat("\n", strrep("=", 60), "\n", sep = "")
cat("SWA EFFECT SUMMARY (swa: Absent = reference)\n")
cat("Negative estimate = lower abundance when endophyte present\n")
cat(strrep("=", 60), "\n")
print(swa_summary, row.names = FALSE)
