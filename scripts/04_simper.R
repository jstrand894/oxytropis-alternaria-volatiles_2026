# ============================================================
# 04_simper.R
# SIMPER analysis identifying compounds driving VOC dissimilarity
# between endophyte-present (swa = 1) and endophyte-absent (swa = 0)
# plants, in greenhouse and field datasets.
#
# Note: SIMPER uses Bray-Curtis on raw (non-CLR) ng/g/hr values.
# CLR is not appropriate here because Bray-Curtis requires non-negative
# values. The raw matrices retain true zeros; no pseudocount is applied.
#
# Permutation testing is omitted: with paired designs and small n,
# SIMPER permutation p-values are unreliable. Results are reported
# as ranked compound contributions to average dissimilarity.
#
# INPUTS (from 01_clr_processing.R):
#   gh_raw_mat, gh_meta
#   field_live_raw_mat, field_live_meta
#   field_sens_raw_mat, field_sens_meta
#
# OUTPUTS:
#   simper_gh, simper_gh_2022, simper_gh_2023
#   simper_live, simper_live_flowers, simper_live_leaves
#   simper_sens
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "01_clr_processing_2.R"))


# ============================================================
# HELPER
# ============================================================

# Extract top n compounds from a simper() result.
# simper() returns one list element per group pair; with two swa levels
# there is always exactly one element.
#
# Columns:
#   avg_dissim  - average contribution to Bray-Curtis dissimilarity
#   sd          - standard deviation of contribution across sample pairs
#   ratio       - avg_dissim / sd; higher = more consistent contributor
#   avg_absent  - mean raw abundance in swa = 0 (endophyte absent)
#   avg_present - mean raw abundance in swa = 1 (endophyte present)
#   cumsum      - cumulative proportion of total dissimilarity

extract_simper <- function(sim, n = 10) {
  s <- summary(sim)[[1]]
  s |>
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

# Thin wrapper: subset a raw matrix and metadata, run simper on swa.
run_simper <- function(raw_mat, meta, keep = rep(TRUE, nrow(meta))) {
  simper(raw_mat[keep, ], meta$swa[keep])
}


# ============================================================
# GREENHOUSE SIMPER
# ============================================================

# Combined across years
simper_gh      <- run_simper(gh_raw_mat, gh_meta)
simper_gh_2022 <- run_simper(gh_raw_mat, gh_meta, gh_meta$year == "Spring 2022")
simper_gh_2023 <- run_simper(gh_raw_mat, gh_meta, gh_meta$year == "Spring 2023")

cat("============================================================\n")
cat("GH SIMPER — combined (top 10 compounds)\n")
cat("============================================================\n")
print(extract_simper(simper_gh), row.names = FALSE)

cat("\n------------------------------------------------------------\n")
cat("GH SIMPER — Spring 2022\n")
cat("------------------------------------------------------------\n")
print(extract_simper(simper_gh_2022), row.names = FALSE)

cat("\n------------------------------------------------------------\n")
cat("GH SIMPER — Spring 2023\n")
cat("------------------------------------------------------------\n")
print(extract_simper(simper_gh_2023), row.names = FALSE)


# ============================================================
# FIELD — Living tissue SIMPER (flowers + leaves)
# ============================================================

simper_live <- run_simper(field_live_raw_mat, field_live_meta)

simper_live_flowers <- run_simper(
  field_live_raw_mat, field_live_meta,
  field_live_meta$collection.type == "flowers"
)

simper_live_leaves <- run_simper(
  field_live_raw_mat, field_live_meta,
  field_live_meta$collection.type == "leaves"
)

cat("\n============================================================\n")
cat("Field (living) SIMPER — combined (top 10 compounds)\n")
cat("============================================================\n")
print(extract_simper(simper_live), row.names = FALSE)

cat("\n------------------------------------------------------------\n")
cat("Field (living) SIMPER — flowers only\n")
cat("------------------------------------------------------------\n")
print(extract_simper(simper_live_flowers), row.names = FALSE)

cat("\n------------------------------------------------------------\n")
cat("Field (living) SIMPER — leaves only\n")
cat("------------------------------------------------------------\n")
print(extract_simper(simper_live_leaves), row.names = FALSE)


# ============================================================
# FIELD — Senesced leaves SIMPER (supplementary)
# ============================================================

simper_sens <- run_simper(field_sens_raw_mat, field_sens_meta)

cat("\n============================================================\n")
cat("Field (senesced) SIMPER — combined (top 10 compounds)\n")
cat("============================================================\n")
print(extract_simper(simper_sens), row.names = FALSE)


# ============================================================
# CROSS-DATASET SUMMARY
# ============================================================
# Top compound by ratio (most consistent contributor) in each dataset.
# Useful for quickly scanning whether the same compounds appear
# across contexts or whether drivers are dataset-specific.

top_by_ratio <- function(sim, label, n = 5) {
  extract_simper(sim, n = nrow(summary(sim)[[1]])) |>
    arrange(desc(ratio)) |>
    slice_head(n = n) |>
    mutate(dataset = label) |>
    select(dataset, compound, avg_dissim, ratio, avg_absent, avg_present)
}

ratio_summary <- bind_rows(
  top_by_ratio(simper_gh,           "GH combined"),
  top_by_ratio(simper_gh_2022,      "GH 2022"),
  top_by_ratio(simper_gh_2023,      "GH 2023"),
  top_by_ratio(simper_live,         "Field live combined"),
  top_by_ratio(simper_live_flowers, "Field live flowers"),
  top_by_ratio(simper_live_leaves,  "Field live leaves"),
  top_by_ratio(simper_sens,         "Field senesced")
)

cat("\n============================================================\n")
cat("CROSS-DATASET: Top 5 compounds by consistency ratio\n")
cat("(ratio = avg_dissim / sd; higher = more consistent contributor)\n")
cat("============================================================\n")
print(ratio_summary, row.names = FALSE)







