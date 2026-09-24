library(dplyr)
library(ggplot2)

library(dplyr)

emm_df <- tibble(
  year      = rep(c(2023, 2024), each = 2),
  endophyte = rep(c("Not Present", "Present"), times = 2),
  
  # Panel A: Pollinator Visits
  pollinator_emm  = c(0.651, 0.629, 0.635, 0.705),
  
  # Panel B: Nectar % Sucrose
  sucrose_emm   = c(21, 20, 50, 49.5),
  sucrose_lower = c(-125, -130, -50, -75),
  sucrose_upper = c(125, 130, 200, 185),
  
  # Panel C: Nectar Volume (mm)
  volume_emm   = c(1.95, 1.93, 1.1, 1.0),
  volume_lower = c(-4, -4, -3, -4),
  volume_upper = c(4, 4, 7, 6)
) %>% 
  mutate(
    pollinator_lower = pollinator_emm - 0.05,
    pollinator_upper = pollinator_emm + 0.05,
    sucrose_lower = sucrose_emm - 150,
    sucrose_upper = sucrose_emm + 150,
    volume_lower = volume_emm - ifelse(year == 2024, 6.5, 5.5),
    volume_upper = volume_emm + ifelse(year == 2024, 6.5, 5.5)
  )
emm_df



emm_df$year      <- factor(emm_df$year)
emm_df$endophyte <- factor(emm_df$endophyte, levels = c("Not Present", "Present"))




emm_long <- emm_df %>%
  pivot_longer(
    cols = c(pollinator_emm, sucrose_emm, volume_emm,
             pollinator_lower, sucrose_lower, volume_lower,
             pollinator_upper, sucrose_upper, volume_upper),
    names_to = c("panel", ".value"),
    names_pattern = "(.+)_(emm|lower|upper)"
  ) %>%
  mutate(
    panel = factor(panel,
                   levels = c("pollinator", "sucrose", "volume"),
                   labels = c("(A) Pollinator Visits",
                              "(B) Nectar % Sucrose",
                              "(C) Nectar Volume (µL)")),
    year = factor(year),
    endophyte = factor(endophyte, levels = c("Not Present", "Present"))
  ) %>%
  mutate(endophyte = fct_recode(endophyte, "Absent" = "Not Present"))

# Plot
library(ggh4x)

nectar_plot <-
ggplot(emm_long, aes(x = year, y = emm, color = endophyte, shape = endophyte, 
                     fill = endophyte, group = endophyte)) +
  geom_point(size = 3, position = position_dodge(width = 0.3),
             fill = "white", stroke = 1.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.3)) +
  facet_wrap(~ panel, scales = "free_y", nrow = 1) +
  facetted_pos_scales(y = list(
    panel == "(A) Pollinator Visits" ~ scale_y_continuous(breaks = seq(0.55, 0.75, by = 0.05)),
    panel == "(B) Nectar % Sucrose"  ~ scale_y_continuous(breaks = seq(-200, 200, by = 100)),
    panel == "(C) Nectar Volume (µL)" ~ scale_y_continuous(breaks = seq(-8, 8, by = 4))
  )) +
  scale_color_manual(values = c("Absent" = "#D55E00", "Present" = "#009E73"),
                     name = "Endophyte") +
  scale_shape_manual(values = c("Absent" = 21, "Present" = 21),
                     name = "Endophyte") +
  scale_fill_manual(values = c("Absent" = "white", "Present" = "#009E73"),
                    name = "Endophyte") +
  labs(x = NULL, y = "Estimated Marginal Means") +
  theme_bw() +
  theme(
    strip.text = element_text(hjust = 0, face = "plain"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10)
  )

ggsave("figures/Fig5_nectar.png", nectar_plot,
       width = 9, height = 5, dpi = 300)
