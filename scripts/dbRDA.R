# ============================================================
# 02b_dbrda.R
# dbRDA with Condition(pair) as alternative to adonis2 + pair term.
#
# Laura's suggestion (comments 55-56): formally partial out pair
# identity rather than estimating it as a model term, which avoids
# eating pair's df from the swa test and removes any ambiguity
# about whether pair is a strata or a model term.
#
# This script sources the same CLR outputs as 02_permanova.R
# and runs parallel models for GH and field (living tissue).
# Compare results to gh_permanova and field_permanova.
#
# OUTPUTS:
#   dbrda_gh, dbrda_field_live
#   anova_dbrda_gh, anova_dbrda_field_live
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "01_clr_processing_2.R"))


# ============================================================
# GREENHOUSE — dbRDA with Condition(pair)
# ============================================================
# Condition(pair) partials pair identity out of the ordination
# space before testing swa. Permutations restricted within year
# (same strata scheme as 02_permanova.R).

dbrda_gh <- capscale(
  gh_clr ~ swa + Condition(pair),
  data   = gh_meta,
  dist   = "euclidean"
)

set.seed(4721)
anova_dbrda_gh <- anova(
  dbrda_gh,
  permutations = 999,
  strata       = gh_meta$year,
  by           = "margin"
)

cat("--- GH dbRDA: swa after Condition(pair), strata = year ---\n")
print(anova_dbrda_gh)


# ============================================================
# FIELD — dbRDA with Condition(pair) — living tissue
# ============================================================

dbrda_field_live <- capscale(
  field_live_clr ~ swa + collection.type + Condition(pair2),
  data   = field_live_meta,
  dist   = "euclidean"
)

set.seed(4721)
anova_dbrda_field_live <- anova(
  dbrda_field_live,
  permutations = 999,
  strata       = field_live_meta$year,
  by           = "margin"
)

cat("\n--- Field dbRDA: swa + collection.type after Condition(pair2), strata = year ---\n")
print(anova_dbrda_field_live)

# ============================================================
# FIELD — dbRDA with Condition(pair) — senesced leaves
# ============================================================
dbrda_field_sens <- capscale(
  field_sens_clr ~ swa + Condition(pair2),
  data   = field_sens_meta,
  dist   = "euclidean"
)

set.seed(4721)
anova_dbrda_field_sens <- anova(
  dbrda_field_sens,
  permutations = 999,
  strata       = field_sens_meta$year,
  by           = "margin"
)

cat("\n--- Field dbRDA: swa after Condition(pair2), strata = year — senesced leaves ---\n")
print(anova_dbrda_field_sens)


# ============================================================
# COMPARISON: adonis2 swa result vs. dbRDA swa result
# ============================================================
# Expect similar F and p for swa. Df gain from removing pair
# as an estimated term shows up as larger Residual Df in dbRDA.

extract_term <- function(aov_obj, term) {
  as.data.frame(aov_obj) |>
    rownames_to_column("term") |>
    filter(term == !!term) |>
    transmute(term, Df, F = round(F, 3), p = `Pr(>F)`)
}

cat("\n========================================\n")
cat("GH: adonis2 vs. dbRDA — swa term\n")
cat("========================================\n")
bind_rows(
  extract_term(gh_permanova, "swa")        |> mutate(method = "adonis2"),
  extract_term(anova_dbrda_gh, "swa")      |> mutate(method = "dbRDA + Condition(pair)")
) |> print(row.names = FALSE)

cat("\n========================================\n")
cat("Field (living): adonis2 vs. dbRDA — swa term\n")
cat("========================================\n")
bind_rows(
  extract_term(field_permanova, "swa")          |> mutate(method = "adonis2"),
  extract_term(anova_dbrda_field_live, "swa")   |> mutate(method = "dbRDA + Condition(pair2)")
) |> print(row.names = FALSE)

cat("\n========================================\n")
cat("Field (senesced): adonis2 vs. dbRDA — swa term\n")
cat("========================================\n")
bind_rows(
  extract_term(field_sens_permanova, "swa")        |> mutate(method = "adonis2"),
  extract_term(anova_dbrda_field_sens, "swa")      |> mutate(method = "dbRDA + Condition(pair2)")
) |> print(row.names = FALSE)