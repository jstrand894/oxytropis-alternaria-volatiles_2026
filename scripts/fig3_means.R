# fig3_means.R

# ============================================================
# fig3_simper.R
# Figure 3: Paired mean VOC emission rates (E+ vs E-) for compounds
# flagged by SIMPER as top contributors to Bray-Curtis dissimilarity.
#
# Two panels: GH (combined) and Field living tissue (combined).
# Compounds ordered by SIMPER contribution rank (highest first).
# Points = raw mean ng/g/hr, error bars = 95% CI (t-based).
#
# Replaces the old contribution-magnitude bar chart. That version
# showed how much a compound moved SIMPER's dissimilarity metric,
# which looks like a real difference whether or not E+ and E- means
# actually overlap. This version plots the means directly so overlap
# is visible, consistent with the non-significant PERMANOVA/LMM
# results.
#
# INPUTS (from 04_simper.R, which sources 01_clr_processing_2.R):
#   simper_gh, simper_live
#   gh_raw_mat, gh_meta
#   field_live_raw_mat, field_live_meta
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "04_simper.R"))

# ============================================================
# COMPOUND NAME LOOKUP (unchanged from old version)
# ============================================================
name_lookup <- c(
  "Caryophyllene"      = "\u03b2-Caryophyllene",
  "Bpinene"            = "\u03b2-Pinene",
  "Bocimene"           = "\u03b2-Ocimene",
  "Dlimonene"          = "(R)-Limonene",
  "Methylsalicylate"   = "Methyl salicylate",
  "Hexenylisovalerate" = "(Z)-3-Hexenyl isovalerate",
  "Cisthreehexiso"     = "(Z)-3-Hexenyl isobutyrate",
  "Z3hex"              = "(Z)-3-Hexenyl acetate",
  "Sixmethyl"          = "6-Methyl-5-hepten-2-one",
  "Octanal"            = "Octanal",
  "Nonanal"            = "Nonanal",
  "Decanal"            = "Decanal",
  "Heptanal"           = "Heptanal",
  "Manoyl oxide"       = "Manoyl oxide"
)

# ============================================================
# STEP 1: pull top-N compound names (raw, un-prettified) from SIMPER,
# in rank order, per dataset.
# ============================================================

top_simper_compounds <- function(sim, n = 10) {
  summary(sim)[[1]] %>%
    rownames_to_column("compound") %>%
    arrange(desc(average)) %>%
    slice_head(n = n) %>%
    pull(compound)
}

gh_top_compounds   <- top_simper_compounds(simper_gh)
live_top_compounds <- top_simper_compounds(simper_live)

# ============================================================
# STEP 2: raw group means + 95% CI (t-based) per compound, per swa
# group, pulled from the same raw matrices SIMPER was run on.
# ============================================================

compound_group_means <- function(raw_mat, meta, compounds, label) {
  meta %>%
    rownames_to_column("row_id") %>%
    select(row_id, swa) %>%
    filter(!is.na(swa)) %>%
    tidyr::crossing(compound = compounds) %>%
    mutate(
      abundance = raw_mat[cbind(row_id, compound)],
      swa       = factor(swa, levels = c(0, 1), labels = c("E-", "E+"))
    ) %>%
    group_by(compound, swa) %>%
    summarise(
      n       = n(),
      mean    = mean(abundance, na.rm = TRUE),
      sd      = sd(abundance, na.rm = TRUE),
      se      = sd / sqrt(n),
      ci_half = ifelse(n > 1, qt(0.975, df = n - 1) * se, NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lower = mean - ci_half,
      ci_upper = mean + ci_half,
      dataset  = label
    )
}

gh_means   <- compound_group_means(gh_raw_mat, gh_meta, gh_top_compounds, "Greenhouse")
live_means <- compound_group_means(field_live_raw_mat, field_live_meta, live_top_compounds, "Field (living tissue)")

means_data <- bind_rows(gh_means, live_means) %>%
  mutate(
    compound = str_to_sentence(compound),
    compound = recode(compound, !!!name_lookup),
    dataset  = factor(dataset, levels = c("Greenhouse", "Field (living tissue)"))
  )

# ============================================================
# STEP 3: order compounds within each panel by SIMPER rank
# (decoupled levels per panel, same trick as the old script, so
# panels don't force a shared compound order).
# ============================================================

gh_levels <- str_to_sentence(gh_top_compounds) %>%
  recode(!!!name_lookup) %>%
  rev() %>%
  paste("Greenhouse", ., sep = "__")

live_levels <- str_to_sentence(live_top_compounds) %>%
  recode(!!!name_lookup) %>%
  rev() %>%
  paste("Field (living tissue)", ., sep = "__")

means_data <- means_data %>%
  mutate(
    compound_panel = paste(dataset, compound, sep = "__"),
    compound_panel = factor(compound_panel, levels = c(live_levels, gh_levels))
  )

# ============================================================
# PLOT
# ============================================================

SIMPER_Plot <-
  ggplot(means_data, aes(x = mean, y = compound_panel, color = swa)) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                position = position_dodge(width = 0.5),
                width = 0.25, linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  scale_color_manual(
    values = c("E-" = "#D55E00", "E+" = "#009E73"),
    name   = "Endophyte\nStatus"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_discrete(labels = \(x) sub(".*__", "", x)) +
  facet_wrap(~ dataset, scales = "free") +
  labs(
    x = expression("Mean emission rate (ng "*g^-1*" "*h^-1*") \u00b1 95% CI"),
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor    = element_blank(),
    strip.background    = element_rect(fill = NA),
    legend.background   = element_rect(color = "black", fill = "white", linewidth = 0.3),
    axis.text.y         = element_text(size = 9)
  )

SIMPER_Plot

ggsave("figures/Fig3_paired_means.png", SIMPER_Plot,
       width = 9, height = 5, dpi = 300)