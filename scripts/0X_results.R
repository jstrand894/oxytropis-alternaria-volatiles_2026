# ============================================================
# 00_results.R
# Consolidated results output. Sources all analysis scripts
# and prints every key result in one place.
#
# Source order:
#   02_permanova.R  -> chains: 00_load_data -> 01_clr_processing -> 02
#   04_simper.R     -> re-chains 00 -> 01 (fast); adds simper objects
#   05_univariate.R -> re-chains 00 -> 01 (fast); adds lmer objects
#
# Note: 00_load_data and 01_clr_processing run 3x due to chaining.
# This is acceptable — they are fast. Permanova (slow) runs once.
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "02_permanova_2.R"))
source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "04_simper.R"))
source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "05_univariate.R"))


divider  <- function() cat(strrep("=", 64), "\n")
subdiv   <- function() cat(strrep("-", 64), "\n")
section  <- function(title) { cat("\n"); divider(); cat(title, "\n"); divider() }
subsect  <- function(title) { cat("\n"); subdiv();  cat(title, "\n"); subdiv() }

fmt_p <- function(p) {
  ifelse(is.na(p), "—",
    ifelse(p < 0.001, "<0.001",
      formatC(p, digits = 3, format = "f")))
}


# ============================================================
# SAMPLE SIZES
# ============================================================

section("SAMPLE SIZES")

cat("\nGreenhouse:\n")
print(gh_meta %>% count(year, swa))

cat("\nField — living tissue:\n")
print(field_live_meta %>% count(year, collection.type, swa))

cat("\nField — senesced leaves:\n")
print(field_sens_meta %>% count(year, swa))


# ============================================================
# BETADISPER
# ============================================================

section("BETADISPER (homogeneity of multivariate dispersion)")
cat("Non-significant p: PERMANOVA reflects centroid shift, not dispersion.\n")

subsect("Greenhouse")
cat("  Year:      p =", fmt_p(betadisp_gh_year_p), "\n")
cat("  Endophyte: p =", fmt_p(betadisp_gh_swa_p),  "\n")

subsect("Field — living tissue")
cat("  Year:        p =", fmt_p(betadisp_live_year_p), "\n")
cat("  Tissue type: p =", fmt_p(betadisp_live_type_p), "\n")
cat("  Endophyte:   p =", fmt_p(betadisp_live_swa_p),  "\n")

subsect("Field — senesced leaves")
cat("  Year:      p =", fmt_p(betadisp_sens_year_p), "\n")
cat("  Endophyte: p =", fmt_p(betadisp_sens_swa_p),  "\n")


# ============================================================
# PERMANOVA
# ============================================================

section("PERMANOVA (adonis2, marginal, strata = year, 999 permutations)")

subsect("Greenhouse — combined")
print(gh_permanova)

subsect("Greenhouse — Spring 2022")
print(gh_perm_2022)

subsect("Greenhouse — Spring 2023")
print(gh_perm_2023)

subsect("Field — living tissue (combined)")
print(field_permanova)

subsect("Field — living tissue 2022")
print(field_perm_2022)

subsect("Field — living tissue 2023 (collection.type omitted: collinear with pair2)")
print(field_perm_2023)

subsect("Field — senesced leaves (combined, supplementary)")
print(field_sens_permanova)


# ============================================================
# PERMANOVA SUMMARY TABLE
# ============================================================

section("PERMANOVA SUMMARY TABLE")

extract_row <- function(aov_obj, term, term_label, dataset) {
  df_aov <- as.data.frame(aov_obj)
  if (!term %in% rownames(df_aov)) return(NULL)
  row <- df_aov[term, ]
  tibble(
    Dataset  = dataset,
    Term     = term_label,
    Df       = row$Df,
    F        = round(row$F, 3),
    `R² (%)` = round(row$R2 * 100, 2),
    p        = fmt_p(row$`Pr(>F)`)
  )
}

permanova_table <- bind_rows(
  extract_row(gh_permanova,         "swa",             "Endophyte",    "Greenhouse"),
  extract_row(gh_permanova,         "pair",            "Pair",         "Greenhouse"),
  extract_row(field_permanova,      "swa",             "Endophyte",    "Field — living"),
  extract_row(field_permanova,      "collection.type", "Tissue type",  "Field — living"),
  extract_row(field_permanova,      "pair2",           "Pair",         "Field — living"),
  extract_row(field_sens_permanova, "swa",             "Endophyte",    "Field — senesced"),
  extract_row(field_sens_permanova, "pair2",           "Pair",         "Field — senesced")
)

print(as.data.frame(permanova_table), row.names = FALSE)


# ============================================================
# SIMPER
# ============================================================

section("SIMPER — Top 10 compounds by average dissimilarity")

extract_simper_print <- function(sim, n = 10) {
  summary(sim)[[1]] |>
    rownames_to_column("compound") |>
    arrange(desc(average)) |>
    slice_head(n = n) |>
    transmute(
      compound,
      avg_dissim  = round(average, 4),
      sd          = round(sd, 4),
      ratio       = round(ratio, 3),
      avg_absent  = round(ava, 3),
      avg_present = round(avb, 3),
      cumsum      = round(cumsum, 3)
    )
}

subsect("Greenhouse — combined")
print(extract_simper_print(simper_gh), row.names = FALSE)

subsect("Greenhouse — Spring 2022")
print(extract_simper_print(simper_gh_2022), row.names = FALSE)

subsect("Greenhouse — Spring 2023")
print(extract_simper_print(simper_gh_2023), row.names = FALSE)

subsect("Field — living tissue (combined)")
print(extract_simper_print(simper_live), row.names = FALSE)

subsect("Field — living tissue: flowers")
print(extract_simper_print(simper_live_flowers), row.names = FALSE)

subsect("Field — living tissue: leaves")
print(extract_simper_print(simper_live_leaves), row.names = FALSE)

subsect("Field — senesced leaves (supplementary)")
print(extract_simper_print(simper_sens), row.names = FALSE)

subsect("Cross-dataset: top 5 by consistency ratio")
print(ratio_summary, row.names = FALSE)


# ============================================================
# UNIVARIATE MODELS
# ============================================================

section("UNIVARIATE MIXED MODELS (lmer, Satterthwaite p-values)")
cat("Response: log(ng/g/hr + 0.01). Random effect: pair identity.\n")
cat("Focus compounds: nonanal and octanal (top SIMPER ratio compounds).\n")
cat("Negative estimate = lower abundance when endophyte present.\n\n")

print(swa_summary, row.names = FALSE)
