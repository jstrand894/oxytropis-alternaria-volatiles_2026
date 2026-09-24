# ============================================================
# 00_load_data.R
# Load and lightly clean raw VOC data from all datasets.
# This replaces knitr::knit() calls to the proc-docs Rmds so that
# objects are reliably created in the calling environment and paths
# are portable (all relative to the project root).
#
# All paths assume the working directory is the project root.
# In RStudio, open the .Rproj file to set this automatically.
#
# OUTPUTS:
#   joined.data   - Long-format GH VOC data, raw ng/g/hr (Spring 2022 + 2023)
#   field_comb    - Long-format field VOC data, raw ng/g/hr (2022 + 2023)
#   proc_2019     - CLR-processed 2019 field data (Megan; visualization only)
# ============================================================

library(tidyverse)
library(vegan)
library(compositions)
library(ggrepel)
library(indicspecies)
library(patchwork)
library(readxl)


# ============================================================
# GREENHOUSE — Spring 2022 & 2023
# (adapted from proc-docs/GH-data-proc.Rmd; paths corrected to relative)
# ============================================================

gh22spring <- read_excel("data-files/Spring 2022 GH Locoweed RESULTS.xlsx",
                         sheet = "RAW2")

gh23spring <- read_excel("data-files/Spring GH 23 Locoweed Volatiles.xlsx",
                         sheet = "RAW_PROC")

# Combine both years into one long-format data frame.
# Metadata (swa, pair) is joined from the raw sheets.
joined.data <-
  gh23spring %>%
  dplyr::select(plant.id, compound, ngghr) %>%
  group_by(plant.id, compound) %>%
  na.omit() %>%
  dplyr::summarise(ngghr = sum(ngghr), .groups = "drop") %>%
  distinct() %>%
  pivot_wider(names_from = "compound", values_from = "ngghr") %>%
  pivot_longer(-plant.id, names_to = "compound", values_to = "ngghr") %>%
  mutate(year = "Spring 2023") %>%
  left_join(
    gh23spring %>% dplyr::select(plant.id, swa, pair),
    by = "plant.id",
    relationship = "many-to-many"
  ) %>%
  distinct() %>%
  rbind(
    gh22spring %>%
      dplyr::select(plant.id, compound, ngghr) %>%
      group_by(plant.id, compound) %>%
      na.omit() %>%
      dplyr::summarise(ngghr = sum(ngghr), .groups = "drop") %>%
      distinct() %>%
      pivot_wider(names_from = "compound", values_from = "ngghr") %>%
      pivot_longer(-plant.id, names_to = "compound", values_to = "ngghr") %>%
      distinct() %>%
      mutate(year = "Spring 2022") %>%
      left_join(
        gh22spring %>% dplyr::select(plant.id, swa, pair),
        by = "plant.id",
        relationship = "many-to-many"
      )
  ) %>%
  mutate_at(vars(pair, swa), list(factor)) %>%
  distinct()


# ============================================================
# FIELD — 2022 & 2023
# (adapted from proc-docs/Field-data-proc.Rmd)
# ============================================================

field22 <- read_excel("data-files/2022 Field Locoweed.xlsx",
                      sheet = "RAW proc") %>%
  mutate(year = "2022")

field23 <- read_excel("data-files/2023 Field Locoweed.xlsx",
                      sheet = "proc") %>%
  mutate(year = "2023")

field_comb <-
  field22 %>%
  dplyr::select(sample, plant.id, compound, ngghr, year, swa, collection.type, pair2) %>%
  rbind(
    field23 %>%
      rename(collection.type = type) %>%
      dplyr::select(sample, plant.id, compound, ngghr, year, swa, collection.type, pair2)
  ) %>%
  mutate(collection.type = case_when(
    collection.type == "l" ~ "leaves",
    collection.type == "f" ~ "flowers",
    collection.type == "s" ~ "senesced_leaves",
    TRUE                   ~ collection.type
  )) %>%
  filter(!collection.type %in% c("0", NA, "NA"))


# ============================================================
# HISTORICAL FIELD — 2019 (Megan; supplemental visualization only)
# (adapted from proc-docs/Megan-data-proc.Rmd; paths corrected to relative)
# Not included in main PERMANOVA — different compound set, no pairing metadata.
# ============================================================

raw_2019 <- read_excel(
  "data-files/Locoweed Field Volatiles 2019 Summary.xlsx",
  sheet = "ng per g per hr",
  range = "A1:L41"
)

proc_2019 <-
  raw_2019 %>%
  unite("plant_id", plant, num, sep = "_") %>%
  dplyr::select(-c(biomass)) %>%
  pivot_longer(-c(collection, endophyte, plant_id)) %>%
  mutate(value = ifelse(value == 0, 0.01, value),
         value = value + 1) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  unite("unite", collection, endophyte, plant_id, sep = "_") %>%
  column_to_rownames("unite") %>%
  clr() %>%
  as.data.frame() %>%
  rownames_to_column("unite") %>%
  separate(unite,
           into = c("collection", "endophyte", "plant_id", "num"),
           sep  = "_") %>%
  pivot_longer(-c(plant_id, collection, endophyte, num),
               names_to  = "name",
               values_to = "value") %>%
  mutate(value = value + 2) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  mutate_at(vars(endophyte), list(factor)) %>%
  distinct() %>%
  na.omit()
