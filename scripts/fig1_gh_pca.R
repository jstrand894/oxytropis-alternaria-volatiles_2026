# ============================================================
# fig1_gh_pca.R
# Figure 1: Greenhouse PCA biplot (Spring 2022 + 2023)
#
# INPUTS (from 03_pca.R):
#   gh_pca_data, gh_pc1_var, gh_pc2_var
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "03_pca_2.R"))

GH_PCA_Plot <-
  ggplot(gh_pca_data) +
  geom_point(aes(PC1, PC2, shape = year, fill = swa),
             size = 2.5, stroke = 0.4, alpha = 0.8, color = "black") +
  stat_ellipse(aes(PC1, PC2, color = swa, group = swa),
               level = 0.95) +
  scale_fill_manual(values  = c("Absent" = "#D55E00", "Present" = "#009E73")) +
  scale_color_manual(values = c("Absent" = "#D55E00", "Present" = "#009E73")) +
  scale_shape_manual(values = c("2022" = 21, "2023" = 24)) +
  labs(
    x     = paste0("PC1 (", gh_pc1_var, "%)"),
    y     = paste0("PC2 (", gh_pc2_var, "%)"),
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
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  )

GH_PCA_Plot

ggsave("figures/Fig1_GH_PCA.png", GH_PCA_Plot,
       width = 6, height = 5, dpi = 300)
