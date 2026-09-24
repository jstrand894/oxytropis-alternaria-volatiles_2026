# ============================================================
# 03_pca.R
# PCA on CLR-transformed VOC data for GH and field datasets.
# Produces scored data frames and scaled loading arrows
# used by the figure scripts.
#
# INPUTS (from 01_clr_processing.R):
#   gh_clr, gh_meta
#   field_clr, field_meta   <- all three tissue types (flowers, leaves,
#                              senesced_leaves). PCA is run on all three
#                              for visualization. The main PERMANOVA uses
#                              living tissue only (field_live_*).
#
# OUTPUTS:
#   gh_pca, gh_pca_data, gh_pca_loadings, gh_pc1_var, gh_pc2_var
#   field_pca, field_pca_data, field_pca_loadings, field_pc1_var, field_pc2_var
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "01_clr_processing_2.R"))


# ============================================================
# GREENHOUSE PCA
# ============================================================

gh_pca    <- prcomp(as.matrix(gh_clr), scale. = TRUE)
gh_scores <- as.data.frame(gh_pca$x) %>%
  rownames_to_column("row_id")

gh_pca_data <- gh_meta %>%
  rownames_to_column("row_id") %>%
  left_join(gh_scores, by = "row_id") %>%
  mutate(
    swa  = factor(swa,  levels = c(0, 1), labels = c("Absent", "Present")),
    year = recode(year, "Spring 2022" = "2022", "Spring 2023" = "2023")
  )

# Variance explained by PC1 and PC2
gh_pca_imp <- summary(gh_pca)$importance |>
  as.data.frame() |>
  slice(2)  # Proportion of Variance row
gh_pc1_var <- round(gh_pca_imp$PC1 * 100, 1)
gh_pc2_var <- round(gh_pca_imp$PC2 * 100, 1)

# Top loading compounds (for interpretation)
gh_rot <- as.data.frame(gh_pca$rotation) |> rownames_to_column("compound")
gh_top_pc1 <- gh_rot |> slice_max(abs(PC1), n = 3) |> pull(compound)
gh_top_pc2 <- gh_rot |> slice_max(abs(PC2), n = 3) |> pull(compound)
cat("GH PC1 top compounds:", paste(gh_top_pc1, collapse = ", "), "\n")
cat("GH PC2 top compounds:", paste(gh_top_pc2, collapse = ", "), "\n")

# Loading arrows: scaled to the range of PC scores for biplot overlay.
# Only the 6 most influential compounds are retained for a clean biplot.
gh_pca_loadings <- as.data.frame(gh_pca$rotation[, 1:2]) %>%
  rownames_to_column("compound") %>%
  mutate(
    PC1 = PC1 * max(abs(gh_pca_data$PC1)) * 0.8,
    PC2 = PC2 * max(abs(gh_pca_data$PC2)) * 0.8
  ) %>%
  slice_max(order_by = sqrt(PC1^2 + PC2^2), n = 6)


# ============================================================
# FIELD PCA
# ============================================================
# Uses field_clr / field_meta which includes all three tissue types
# (flowers, leaves, senesced_leaves). All three are shown in the
# ordination so reviewers can see the full compositional spread.
# The main statistical test (PERMANOVA) uses living tissue only.

field_pca    <- prcomp(as.matrix(field_clr), scale. = TRUE)
field_scores <- as.data.frame(field_pca$x) %>%
  rownames_to_column("row_id")

field_pca_data <- field_meta %>%
  rownames_to_column("row_id") %>%
  left_join(field_scores, by = "row_id") %>%
  mutate(
    swa             = factor(swa, levels = c(0, 1), labels = c("Absent", "Present")),
    collection.type = recode(collection.type,
                             "flowers"         = "Flowers",
                             "leaves"          = "Leaves",
                             "senesced_leaves" = "Senesced Leaves"
    )
  )

field_pca_imp <- summary(field_pca)$importance |>
  as.data.frame() |>
  slice(2)  # Proportion of Variance row
field_pc1_var <- round(field_pca_imp$PC1 * 100, 1)
field_pc2_var <- round(field_pca_imp$PC2 * 100, 1)

field_rot <- as.data.frame(field_pca$rotation) |> rownames_to_column("compound")
field_top_pc1 <- field_rot |> slice_max(abs(PC1), n = 3) |> pull(compound)
field_top_pc2 <- field_rot |> slice_max(abs(PC2), n = 3) |> pull(compound)
cat("Field PC1 top compounds:", paste(field_top_pc1, collapse = ", "), "\n")
cat("Field PC2 top compounds:", paste(field_top_pc2, collapse = ", "), "\n")

field_pca_loadings <- as.data.frame(field_pca$rotation[, 1:2]) %>%
  rownames_to_column("compound") %>%
  mutate(
    PC1 = PC1 * max(abs(field_pca_data$PC1)) * 0.8,
    PC2 = PC2 * max(abs(field_pca_data$PC2)) * 0.8
  ) %>%
  slice_max(order_by = sqrt(PC1^2 + PC2^2), n = 6)
