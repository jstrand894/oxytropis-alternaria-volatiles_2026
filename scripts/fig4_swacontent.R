library(dplyr)
library(ggplot2)

df <- data.frame(
  flower_part = c("Leaf", "Calyx", "Petal", "Pistil", "Pollen", "Nectar"),
  swainsonine = c(0.2, 0.0807461, 0.02380279, 0.02353042, 0.02248456, 0.002),
  SE          = c(0.0000008, 0.010938, 0.00832, 0.009213, 0.00617, NA)
)

df$flower_part <- factor(df$flower_part, 
                         levels = c("Leaf", "Calyx", "Petal", 
                                    "Pistil", "Pollen", "Nectar"))

bar_colors <- colorRampPalette(c("#009E73", "#009E73"))(6)

swa_cont <-
df %>%
  ggplot(aes(x = flower_part, y = swainsonine, fill = flower_part)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = swainsonine - SE, ymax = swainsonine + SE), 
                width = 0.2) +
  scale_fill_manual(values = setNames(bar_colors, levels(df$flower_part))) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.position = "none"
  ) +
  labs(x = NULL,
       y = "Swainsonine content (mg/g dry weight)")
swa_cont

ggsave("figures/Fig4_swa_cont.png", swa_cont,
       width = 6, height = 5, dpi = 300)
