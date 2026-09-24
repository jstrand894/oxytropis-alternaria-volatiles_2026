# ============================================================
# fig2_field_pca.R
# Figure 2: Field PCA biplot (2022 + 2023), faceted by tissue type
# Three panels: Flowers, Leaves, Senesced Leaves
#
# PCA is run on all three tissue types (field_clr / field_meta)
# so the ordination reflects the full compositional space.
# The main PERMANOVA uses living tissue only (see 02_permanova.R).
#
# INPUTS (from 03_pca.R):
#   field_pca_data, field_pc1_var, field_pc2_var
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "03_pca_2.R"))


# Set panel order: living tissue first, senesced last
field_pca_data <- field_pca_data %>%
  mutate(collection.type = factor(collection.type,
                                  levels = c("Flowers", "Leaves", "Senesced Leaves")))

FIELD_PCA_Plot <-
  ggplot(field_pca_data) +
  geom_point(aes(PC1, PC2, shape = year, fill = swa),
             size = 2.5, stroke = 0.4, alpha = 0.8, color = "black") +
  stat_ellipse(aes(PC1, PC2, color = swa, group = swa),
               level = 0.95) +
  scale_fill_manual(values  = c("Absent" = "#D55E00", "Present" = "#009E73")) +
  scale_color_manual(values = c("Absent" = "#D55E00", "Present" = "#009E73")) +
  scale_shape_manual(values = c("2022" = 21, "2023" = 24)) +
  facet_wrap(~ collection.type) +
  labs(
    x     = paste0("PC1 (", field_pc1_var, "%)"),
    y     = paste0("PC2 (", field_pc2_var, "%)"),
    fill  = "Endophyte",
    color = "Endophyte",
    shape = "Year"
  ) +
  guides(
    shape = guide_legend(override.aes = list(fill = "gray40")),
    fill  = guide_legend(override.aes = list(shape = 21))
  ) +
  theme_bw() +
  theme(
    panel.grid        = element_blank(),
    strip.background  = element_rect(fill = NA),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  )

FIELD_PCA_Plot

ggsave("figures/Fig2_Field_PCA.png", FIELD_PCA_Plot,
       width = 10, height = 4.5, dpi = 300)
