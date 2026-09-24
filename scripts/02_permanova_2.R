# ============================================================
# 02_permanova.R
# Betadisper + PERMANOVA for greenhouse and field datasets.
#
# Field analysis uses two datasets from 01_clr_processing.R:
#   field_live_* (flowers + leaves) -> main PERMANOVA
#   field_sens_* (senesced leaves)  -> supplementary PERMANOVA
#
# INPUTS (from 01_clr_processing.R):
#   gh_clr, gh_meta, gh_dist
#   field_live_clr, field_live_meta, field_live_dist
#   field_sens_clr, field_sens_meta, field_sens_dist
#
# OUTPUTS:
#   betadisp_gh_year_p, betadisp_gh_swa_p
#   betadisp_live_year_p, betadisp_live_type_p, betadisp_live_swa_p
#   betadisp_sens_year_p, betadisp_sens_swa_p
#   gh_permanova, gh_perm_2022, gh_perm_2023
#   field_permanova, field_perm_2022, field_perm_2023
#   field_sens_permanova
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "01_clr_processing_2.R"))

# Helper: extract betadisper p-value for the "Groups" term
get_betadisp_p <- \(bd) permutest(bd)$tab |>
  as.data.frame() |>
  rownames_to_column("term") |>
  filter(term == "Groups") |>
  pull(`Pr(>F)`)


# ============================================================
# GREENHOUSE — Betadisper
# ============================================================
# Betadisper tests whether within-group spread (dispersion) differs
# between groups. Significant dispersion differences can inflate
# PERMANOVA p-values — significant PERMANOVA terms should then be
# interpreted cautiously (may reflect variance, not just centroid shift).

bd_gh_year <- betadisper(gh_dist, gh_meta$year)
bd_gh_swa  <- betadisper(gh_dist, gh_meta$swa)

betadisp_gh_year_p <- get_betadisp_p(bd_gh_year)
betadisp_gh_swa_p  <- get_betadisp_p(bd_gh_swa)

cat("GH betadisper by year, p =", betadisp_gh_year_p, "\n")
cat("GH betadisper by swa,  p =", betadisp_gh_swa_p,  "\n")


# ============================================================
# GREENHOUSE — PERMANOVA
# ============================================================
# Tests whether VOC profiles differ by endophyte status (swa),
# after accounting for plant pair identity.
# by = "margin": marginal (Type III) test — term order does not affect results.
# strata = year: permutations restricted within years to control for year
#   as a blocking factor. Year is excluded from the formula because it is
#   fully collinear with strata and cannot be estimated separately (Df = 0).

set.seed(4721)
gh_permanova <- adonis2(
  gh_dist ~ swa + pair,
  data         = gh_meta,
  permutations = 999,
  strata       = gh_meta$year,
  by           = "margin"
)
cat("\n--- GH PERMANOVA (combined, strata = year) ---\n")
print(gh_permanova)

# Year-stratified models: confirm the swa effect is consistent across years.
run_gh_year_permanova <- function(year_val) {
  keep <- gh_meta$year == year_val
  d    <- vegdist(as.matrix(gh_clr)[keep, ], method = "euclidean")
  set.seed(4721)
  adonis2(d ~ swa + pair, data = gh_meta[keep, ],
          permutations = 999, by = "margin")
}

gh_perm_2022 <- run_gh_year_permanova("Spring 2022")
gh_perm_2023 <- run_gh_year_permanova("Spring 2023")

cat("\nGH PERMANOVA Spring 2022:\n"); print(gh_perm_2022)
cat("\nGH PERMANOVA Spring 2023:\n"); print(gh_perm_2023)


# ============================================================
# FIELD — Betadisper (living tissue: flowers + leaves)
# ============================================================

bd_live_year <- betadisper(field_live_dist, field_live_meta$year)
bd_live_type <- betadisper(field_live_dist, field_live_meta$collection.type)
bd_live_swa  <- betadisper(field_live_dist, field_live_meta$swa)

betadisp_live_year_p <- get_betadisp_p(bd_live_year)
betadisp_live_type_p <- get_betadisp_p(bd_live_type)
betadisp_live_swa_p  <- get_betadisp_p(bd_live_swa)

cat("\nField (living) betadisper by year,            p =", betadisp_live_year_p, "\n")
cat("Field (living) betadisper by collection type, p =", betadisp_live_type_p, "\n")
cat("Field (living) betadisper by swa,             p =", betadisp_live_swa_p,  "\n")


# ============================================================
# FIELD — PERMANOVA (living tissue: flowers + leaves)
# ============================================================
# collection.type (flowers vs. leaves) is included so the endophyte
# effect is tested after accounting for tissue type differences.
# Senesced leaves are excluded; see field_sens_permanova below.

set.seed(4721)
field_permanova <- adonis2(
  field_live_dist ~ swa + collection.type + pair2,
  data         = field_live_meta,
  permutations = 999,
  strata       = field_live_meta$year,
  by           = "margin"
)
cat("\n--- Field PERMANOVA — living tissue (combined, strata = year) ---\n")
print(field_permanova)

# Year-stratified models
run_field_year_permanova <- function(year_val) {
  keep <- field_live_meta$year == year_val
  d    <- vegdist(as.matrix(field_live_clr)[keep, ], method = "euclidean")
  formula <- if (year_val == "2023") {
    d ~ swa + pair2
  } else {
    d ~ swa + collection.type + pair2
  }
  set.seed(4721)
  adonis2(formula, data = field_live_meta[keep, ],
          permutations = 999, by = "margin")
}

field_perm_2022 <- run_field_year_permanova("2022")
field_perm_2023 <- run_field_year_permanova("2023")

cat("\nField PERMANOVA 2022:\n"); print(field_perm_2022)
cat("\nField PERMANOVA 2023:\n"); print(field_perm_2023)


# ============================================================
# FIELD — Betadisper + PERMANOVA (senesced leaves, supplementary)
# ============================================================
# Exploratory analysis to test whether endophyte status is
# associated with VOC profiles in senesced tissue. Reported
# separately from the main living-tissue PERMANOVA.

bd_sens_year <- betadisper(field_sens_dist, field_sens_meta$year)
bd_sens_swa  <- betadisper(field_sens_dist, field_sens_meta$swa)

betadisp_sens_year_p <- get_betadisp_p(bd_sens_year)
betadisp_sens_swa_p  <- get_betadisp_p(bd_sens_swa)

cat("\nField (senesced) betadisper by year, p =", betadisp_sens_year_p, "\n")
cat("Field (senesced) betadisper by swa,  p =", betadisp_sens_swa_p,  "\n")

set.seed(4721)
field_sens_permanova <- adonis2(
  field_sens_dist ~ swa + pair2,
  data         = field_sens_meta,
  permutations = 999,
  strata       = field_sens_meta$year,
  by           = "margin"
)
cat("\n--- Field PERMANOVA — senesced leaves (combined, strata = year) ---\n")
print(field_sens_permanova)


# ============================================================
# RESULTS SUMMARY
# ============================================================

extract_permanova_row <- function(aov_obj, row_term) {
  as.data.frame(aov_obj) |>
    rownames_to_column("term") |>
    filter(term == row_term) |>
    transmute(
      term,
      Df,
      F      = round(F, 3),
      R2_pct = round(R2 * 100, 2),
      p      = `Pr(>F)`
    )
}

gh_summary <- bind_rows(
  extract_permanova_row(gh_permanova, "swa"),
  extract_permanova_row(gh_permanova, "pair")
)

field_summary <- bind_rows(
  extract_permanova_row(field_permanova, "swa"),
  extract_permanova_row(field_permanova, "collection.type"),
  extract_permanova_row(field_permanova, "pair2")
)

field_sens_summary <- bind_rows(
  extract_permanova_row(field_sens_permanova, "swa"),
  extract_permanova_row(field_sens_permanova, "pair2")
)

cat("\n========================================\n")
cat("GREENHOUSE PERMANOVA (marginal, strata = year)\n")
cat("========================================\n")
print(gh_summary, row.names = FALSE)

cat("\n========================================\n")
cat("FIELD PERMANOVA — living tissue (marginal, strata = year)\n")
cat("========================================\n")
print(field_summary, row.names = FALSE)

cat("\n========================================\n")
cat("FIELD PERMANOVA — senesced leaves (marginal, strata = year)\n")
cat("========================================\n")
print(field_sens_summary, row.names = FALSE)

cat("\n--- GH betadisper (year / swa) ---\n")
cat("  year p =", betadisp_gh_year_p, "\n")
cat("  swa  p =", betadisp_gh_swa_p,  "\n")

cat("\n--- Field betadisper — living tissue (year / type / swa) ---\n")
cat("  year p =", betadisp_live_year_p, "\n")
cat("  type p =", betadisp_live_type_p, "\n")
cat("  swa  p =", betadisp_live_swa_p,  "\n")

cat("\n--- Field betadisper — senesced leaves (year / swa) ---\n")
cat("  year p =", betadisp_sens_year_p, "\n")
cat("  swa  p =", betadisp_sens_swa_p,  "\n")


