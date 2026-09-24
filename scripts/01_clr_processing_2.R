# ============================================================
# 01_clr_processing.R
# CLR transformation and Aitchison distance matrices.
# Applied once to raw ng/g/hr values — not to pre-transformed data.
#
# INPUTS (from 00_load_data.R):
#   joined.data, field_comb
#
# OUTPUTS:
#   gh_meta, gh_clr, gh_dist, gh_raw_mat         - GH data
#
#   field_meta, field_clr, field_dist,
#   field_raw_mat                                 - Field data (all 3 tissue types;
#                                                   used for PCA visualization only)
#
#   field_live_meta, field_live_clr,
#   field_live_dist, field_live_raw_mat           - Living tissue only (flowers + leaves;
#                                                   used for main PERMANOVA)
#
#   field_sens_meta, field_sens_clr,
#   field_sens_dist, field_sens_raw_mat           - Senesced leaves only
#                                                   (used for supplementary PERMANOVA)
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "00_load_data.R"))


# ============================================================
# GREENHOUSE
# ============================================================

# Compounds to exclude: internal standards, known contaminants,
# non-biogenic, synthetic/industrial, or straight-chain alkane contaminants.
gh_exclude <- c("tetradecane", "tetradec", "tetramethyl",
                "standard", "o+C5xe", "cyclohexene",
                "homosalate", "ethylhexylsalicylate", "diethylphthalate",
                "benzothiazole", "butylatedhydroxtoluene",
                "dodecylacrylate", "caprylyl acetate",
                "hexadecane", "tridecane", "pentadecane", "nanodecane",
                "p-trimethyl", "dimethylnona")

# Aggregate multiple measurements per plant x compound -> one value per combo.
gh_wide <- joined.data %>%
  filter(!compound %in% gh_exclude,
         !is.na(ngghr)) %>%
  group_by(plant.id, year, swa, pair, compound) %>%
  summarise(ngghr = sum(ngghr), .groups = "drop") %>%
  pivot_wider(names_from  = compound,
              values_from = ngghr,
              values_fill = 0) %>%
  arrange(plant.id, year)

# Separate metadata from compound abundances
gh_meta <- gh_wide %>%
  select(plant.id, year, swa, pair) %>%
  mutate(row_id = paste(plant.id, year, pair, sep = "~")) %>%
  column_to_rownames("row_id")

gh_comp_mat <- gh_wide %>%
  select(-plant.id, -year, -swa, -pair) %>%
  as.matrix()
rownames(gh_comp_mat) <- rownames(gh_meta)

# Replace true zeros with pseudocount before CLR.
# CLR requires strictly positive values; 0.01 matches the original pipeline.
gh_comp_mat[gh_comp_mat == 0] <- 0.01

# CLR transformation (applied once, to raw values).
# Euclidean distance on CLR data = Aitchison distance.
gh_clr <- clr(gh_comp_mat)

# Drop all-zero rows and samples missing swa/pair (e.g., plant 3A2 Spring 2023).
nonzero_gh <- rowSums(as.matrix(gh_clr)) != 0
gh_meta    <- gh_meta |> filter(nonzero_gh, !is.na(swa), !is.na(pair))
gh_clr     <- gh_clr[rownames(gh_meta), ]

# Aitchison distance matrix
gh_dist <- vegdist(as.matrix(gh_clr), method = "euclidean")

# Raw matrix aligned with filtered metadata (used in SIMPER -- Bray-Curtis
# requires non-negative values, so CLR cannot be used for SIMPER).
gh_raw_mat <- gh_wide |>
  mutate(row_id = paste(plant.id, year, pair, sep = "~")) |>
  filter(row_id %in% rownames(gh_meta)) |>
  arrange(match(row_id, rownames(gh_meta))) |>
  select(-plant.id, -year, -swa, -pair, -row_id) |>
  as.matrix()
rownames(gh_raw_mat) <- rownames(gh_meta)


# ============================================================
# FIELD — shared setup
# ============================================================

# Core compounds detected consistently across both years.
# Pre-selection focuses the analysis on biologically meaningful,
# reliably detected compounds; rare/trace compounds are excluded
# to reduce noise in the distance matrix.
field_compounds <- c("heptanal", "sixmethyl", "octanal", "z3hex",
                     "dlimonene", "nonanal", "decanal",
                     "methylsalicylate", "cisthreehexiso", "bocimene")

# All three tissue types retained here. Downstream objects are split into:
#   field_*      (all types)  -> PCA visualization
#   field_live_* (flowers + leaves only) -> main PERMANOVA
#   field_sens_* (senesced leaves only)  -> supplementary PERMANOVA
field_wide <- field_comb %>%
  filter(
    compound %in% field_compounds,
    !is.na(collection.type),
    !collection.type %in% c("0")
  ) %>%
  group_by(plant.id, year, swa, collection.type, pair2, compound) %>%
  summarise(ngghr = sum(ngghr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = compound,
              values_from = ngghr,
              values_fill = 0)


# ============================================================
# FIELD — All tissue types (PCA visualization)
# ============================================================

field_meta <- field_wide %>%
  select(plant.id, year, swa, collection.type, pair2) %>%
  mutate(row_id = paste(plant.id, year, collection.type, pair2, sep = "~")) %>%
  column_to_rownames("row_id")

field_comp_mat <- field_wide %>%
  select(all_of(field_compounds)) %>%
  as.matrix()
rownames(field_comp_mat) <- rownames(field_meta)

field_comp_mat[field_comp_mat == 0] <- 0.01

field_clr     <- clr(field_comp_mat)
nonzero_field <- rowSums(as.matrix(field_clr)) != 0
field_meta    <- field_meta |> filter(nonzero_field)
field_clr     <- field_clr[rownames(field_meta), ]

field_dist <- vegdist(as.matrix(field_clr), method = "euclidean")

field_raw_mat <- field_wide |>
  mutate(row_id = paste(plant.id, year, collection.type, pair2, sep = "~")) |>
  filter(row_id %in% rownames(field_meta)) |>
  arrange(match(row_id, rownames(field_meta))) |>
  select(all_of(field_compounds)) |>
  as.matrix()
rownames(field_raw_mat) <- rownames(field_meta)


# ============================================================
# FIELD — Living tissue only (flowers + leaves)
# Main PERMANOVA dataset. Senesced leaves are excluded because
# their VOC profiles are driven by degradation chemistry, not
# active plant-endophyte interactions, which would dilute the
# swa signal and conflate tissue senescence with treatment effects.
# ============================================================

field_live_wide <- field_wide %>%
  filter(collection.type != "senesced_leaves")

field_live_meta <- field_live_wide %>%
  select(plant.id, year, swa, collection.type, pair2) %>%
  mutate(row_id = paste(plant.id, year, collection.type, pair2, sep = "~")) %>%
  column_to_rownames("row_id")

field_live_comp_mat <- field_live_wide %>%
  select(all_of(field_compounds)) %>%
  as.matrix()
rownames(field_live_comp_mat) <- rownames(field_live_meta)

field_live_comp_mat[field_live_comp_mat == 0] <- 0.01

field_live_clr  <- clr(field_live_comp_mat)
nonzero_live    <- rowSums(as.matrix(field_live_clr)) != 0
field_live_meta <- field_live_meta |> filter(nonzero_live)
field_live_clr  <- field_live_clr[rownames(field_live_meta), ]

field_live_dist <- vegdist(as.matrix(field_live_clr), method = "euclidean")

field_live_raw_mat <- field_live_wide |>
  mutate(row_id = paste(plant.id, year, collection.type, pair2, sep = "~")) |>
  filter(row_id %in% rownames(field_live_meta)) |>
  arrange(match(row_id, rownames(field_live_meta))) |>
  select(all_of(field_compounds)) |>
  as.matrix()
rownames(field_live_raw_mat) <- rownames(field_live_meta)


# ============================================================
# FIELD — Senesced leaves only (supplementary PERMANOVA)
# Exploratory: tests whether endophyte status is associated with
# VOC profiles in senesced tissue, reported separately.
# ============================================================

field_sens_wide <- field_wide %>%
  filter(collection.type == "senesced_leaves")

field_sens_meta <- field_sens_wide %>%
  select(plant.id, year, swa, collection.type, pair2) %>%
  mutate(row_id = paste(plant.id, year, collection.type, pair2, sep = "~")) %>%
  column_to_rownames("row_id")

field_sens_comp_mat <- field_sens_wide %>%
  select(all_of(field_compounds)) %>%
  as.matrix()
rownames(field_sens_comp_mat) <- rownames(field_sens_meta)

field_sens_comp_mat[field_sens_comp_mat == 0] <- 0.01

field_sens_clr  <- clr(field_sens_comp_mat)
nonzero_sens    <- rowSums(as.matrix(field_sens_clr)) != 0
field_sens_meta <- field_sens_meta |> filter(nonzero_sens)
field_sens_clr  <- field_sens_clr[rownames(field_sens_meta), ]

field_sens_dist <- vegdist(as.matrix(field_sens_clr), method = "euclidean")

field_sens_raw_mat <- field_sens_wide |>
  mutate(row_id = paste(plant.id, year, collection.type, pair2, sep = "~")) |>
  filter(row_id %in% rownames(field_sens_meta)) |>
  arrange(match(row_id, rownames(field_sens_meta))) |>
  select(all_of(field_compounds)) |>
  as.matrix()
rownames(field_sens_raw_mat) <- rownames(field_sens_meta)

