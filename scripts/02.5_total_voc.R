# ============================================================
# 03_total_voc.R
# LMM tests for total volatile emission rate (ng/g/hr).
# Complements the compositional PERMANOVA in 02_permanova.R.
#
# INPUTS (from 01_clr_processing.R):
#   gh_wide, gh_meta
#   field_live_wide, field_live_meta
#
# OUTPUTS:
#   gh_voc_lmm, gh_voc_2022, gh_voc_2023
#   field_voc_lmm, field_voc_2022, field_voc_2023
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "02_permanova_2.R"))

library(lme4)
library(lmerTest)   # Satterthwaite df and p-values for lmer

gh_compounds <- colnames(gh_comp_mat)

# ============================================================
# GREENHOUSE
# ============================================================

gh_voc_df <- gh_wide |>
  mutate(
    row_id    = paste(plant.id, year, pair, sep = "~"),
    total_voc = rowSums(across(all_of(gh_compounds)))
  ) |>
  filter(row_id %in% rownames(gh_meta)) |>
  select(row_id, plant.id, year, swa, pair, total_voc) |>
  mutate(log_total = log(total_voc))

# Combined model: year as fixed effect, pair as random intercept
gh_voc_lmm <- lmer(
  log_total ~ swa + year + (1 | pair),
  data = gh_voc_df,
  REML = TRUE
)
cat("\n--- GH total VOC LMM (combined) ---\n")
print(summary(gh_voc_lmm))

# Year-stratified models
run_gh_voc_year <- function(year_val) {
  d <- gh_voc_df |> filter(year == year_val)
  lmer(log_total ~ swa + (1 | pair), data = d, REML = TRUE)
}

gh_voc_2022 <- run_gh_voc_year("Spring 2022")
gh_voc_2023 <- run_gh_voc_year("Spring 2023")

cat("\nGH total VOC LMM — Spring 2022:\n"); print(summary(gh_voc_2022))
cat("\nGH total VOC LMM — Spring 2023:\n"); print(summary(gh_voc_2023))


# ============================================================
# FIELD — living tissue (flowers + leaves)
# ============================================================

field_voc_df <- field_live_wide |>
  mutate(
    row_id    = paste(plant.id, year, collection.type, pair2, sep = "~"),
    total_voc = rowSums(across(all_of(field_compounds)))
  ) |>
  filter(row_id %in% rownames(field_live_meta)) |>
  select(row_id, plant.id, year, swa, collection.type, pair2, total_voc) |>
  mutate(log_total = log(total_voc))

# Combined model
field_voc_lmm <- lmer(
  log_total ~ swa + collection.type + year + (1 | pair2),
  data = field_voc_df,
  REML = TRUE
)
cat("\n--- Field total VOC LMM — living tissue (combined) ---\n")
print(summary(field_voc_lmm))

# Year-stratified models
run_field_voc_year <- function(year_val) {
  d <- field_voc_df |> filter(year == year_val)
  formula <- if (year_val == "2023") {
    log_total ~ swa + (1 | pair2)
  } else {
    log_total ~ swa + collection.type + (1 | pair2)
  }
  lmer(formula, data = d, REML = TRUE)
}

field_voc_2022 <- run_field_voc_year("2022")
field_voc_2023 <- run_field_voc_year("2023")

cat("\nField total VOC LMM — 2022:\n"); print(summary(field_voc_2022))
cat("\nField total VOC LMM — 2023:\n"); print(summary(field_voc_2023))


# ============================================================
# RESULTS SUMMARY
# ============================================================

extract_lmm_row <- function(model, term_str, label) {
  coef_tbl <- as.data.frame(coef(summary(model)))
  coef_tbl |>
    rownames_to_column("term") |>
    filter(grepl(term_str, term, fixed = TRUE)) |>
    transmute(
      dataset  = label,
      term     = term,
      estimate = round(Estimate, 3),
      se       = round(`Std. Error`, 3),
      df       = round(df, 1),
      t        = round(`t value`, 3),
      p        = round(`Pr(>|t|)`, 4)
    )
}

cat("\n========================================\n")
cat("TOTAL VOC LMM RESULTS SUMMARY\n")
cat("========================================\n")

bind_rows(
  extract_lmm_row(gh_voc_lmm,     "swa", "GH combined"),
  extract_lmm_row(gh_voc_2022,    "swa", "GH 2022"),
  extract_lmm_row(gh_voc_2023,   "swa", "GH 2023"),
  extract_lmm_row(field_voc_lmm, "swa", "Field combined"),
  extract_lmm_row(field_voc_2022,"swa", "Field 2022"),
  extract_lmm_row(field_voc_2023,"swa", "Field 2023")
) |> print(row.names = FALSE)