# ============================================================
# table1_permanova.R
# Table 1: PERMANOVA and betadisper results for GH and field
# VOC compositional analyses.
#
# Output is printed to console for copy-paste into Word/Excel.
#
# INPUTS (from 02_permanova.R):
#   gh_permanova, field_permanova, field_sens_permanova
#   betadisp_gh_year_p, betadisp_gh_swa_p
#   betadisp_live_year_p, betadisp_live_type_p, betadisp_live_swa_p
#   betadisp_sens_year_p, betadisp_sens_swa_p
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "02_permanova_2.R"))


# ============================================================
# HELPER
# ============================================================

# Format p-values: <0.001 for very small, exact to 3 decimal places otherwise.
fmt_p <- function(p) {
  ifelse(is.na(p), "—",
    ifelse(p < 0.001, "<0.001",
      formatC(p, digits = 3, format = "f")))
}

# Extract a named term row from an adonis2 object into a display-ready tibble.
extract_row <- function(aov_obj, term, term_label, dataset) {
  df_aov <- as.data.frame(aov_obj)
  if (!term %in% rownames(df_aov)) return(NULL)
  row <- df_aov[term, ]
  tibble(
    Dataset = dataset,
    Term    = term_label,
    Df      = row$Df,
    F       = round(row$F, 3),
    `R² (%)` = round(row$R2 * 100, 2),
    p       = fmt_p(row$`Pr(>F)`)
  )
}


# ============================================================
# BUILD PERMANOVA TABLE
# ============================================================

permanova_table <- bind_rows(

  # --- Greenhouse ---
  extract_row(gh_permanova, "swa",  "Endophyte (swa)",  "Greenhouse"),
  extract_row(gh_permanova, "pair", "Pair",              "Greenhouse"),

  # --- Field living tissue ---
  extract_row(field_permanova, "swa",             "Endophyte (swa)",  "Field — living tissue"),
  extract_row(field_permanova, "collection.type", "Tissue type",      "Field — living tissue"),
  extract_row(field_permanova, "pair2",           "Pair",             "Field — living tissue"),

  # --- Field senesced (supplementary) ---
  extract_row(field_sens_permanova, "swa",   "Endophyte (swa)", "Field — senesced leaves"),
  extract_row(field_sens_permanova, "pair2", "Pair",            "Field — senesced leaves")

)


# ============================================================
# BUILD BETADISPER TABLE
# ============================================================

betadisper_table <- tribble(
  ~Dataset,                    ~Group,          ~p,
  "Greenhouse",                "Year",          fmt_p(betadisp_gh_year_p),
  "Greenhouse",                "Endophyte",     fmt_p(betadisp_gh_swa_p),
  "Field — living tissue",     "Year",          fmt_p(betadisp_live_year_p),
  "Field — living tissue",     "Tissue type",   fmt_p(betadisp_live_type_p),
  "Field — living tissue",     "Endophyte",     fmt_p(betadisp_live_swa_p),
  "Field — senesced leaves",   "Year",          fmt_p(betadisp_sens_year_p),
  "Field — senesced leaves",   "Endophyte",     fmt_p(betadisp_sens_swa_p)
)


# ============================================================
# PRINT
# ============================================================

cat("\n")
cat("================================================================\n")
cat("TABLE 1A: PERMANOVA results (adonis2, marginal, strata = year)\n")
cat("Aitchison distance (CLR-transformed ng/g/hr)\n")
cat("================================================================\n")
print(as.data.frame(permanova_table), row.names = FALSE)

cat("\n")
cat("================================================================\n")
cat("TABLE 1B: Betadisper results (homogeneity of dispersion)\n")
cat("Non-significant p indicates PERMANOVA results reflect centroid\n")
cat("differences, not dispersion differences between groups.\n")
cat("================================================================\n")
print(as.data.frame(betadisper_table), row.names = FALSE)

cat("\n")
cat("Notes:\n")
cat("  - All PERMANOVAs: permutations = 999, by = 'margin', strata = year\n")
cat("  - GH: pair = matched plant pair identity (n = 30 pairs per year)\n")
cat("  - Field: pair2 = matched field pair identity\n")
cat("  - Field 2023 year-stratified model: collection.type omitted due to\n")
cat("    collinearity with pair2 (each pair sampled as single tissue type)\n")
cat("  - Senesced leaves analyzed separately (exploratory)\n")
