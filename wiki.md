# Locoweed Volatiles 2026 — Project Wiki

## Table of Contents

1. [Research Overview](#research-overview)
2. [Datasets](#datasets)
3. [Repository Structure](#repository-structure)
4. [Data Processing Pipeline](#data-processing-pipeline)
5. [Statistical Approach](#statistical-approach)
6. [Key Analytical Decisions](#key-analytical-decisions)
7. [Known Data Issues](#known-data-issues)

---

## Research Overview

This project asks whether *Oxytropis sericea* (white locoweed) plants harboring the endophytic fungus *Alternaria oxytropis* (formerly *Undifilum oxytropis*) differ in their volatile organic compound (VOC) profiles from plants without the endophyte. Endophyte status is confirmed by seed wash assay (SWA), with SWA-positive plants designated E+ and SWA-negative plants designated E-.

Plants were spatially paired (one E+, one E-) to control for microsite and genetic variation. Analyses are conducted separately for greenhouse and field collections, with year treated as a blocking factor throughout.

---

## Datasets

| Dataset | Years | Context | Tissue | File(s) |
|---|---|---|---|---|
| Greenhouse (Jackson) | Spring 2022, Spring 2023 | Whole-plant VOC collection; plants paired by seed family / planting block | Whole plant | `Spring 2022 GH Locoweed RESULTS.xlsx`, `Spring GH 23 Locoweed Volatiles.xlsx` |
| Field (Jackson) | 2022, 2023 | Co-occurring paired plants at a natural population | Flowers and leaves separately | `2022 Field Locoweed.xlsx`, `2023 Field Locoweed.xlsx` |
| Historical field (Megan) | 2018, 2019 | Supplemental visualization only; not included in main PERMANOVA | -- | `Locoweed Field Volatiles 2018.xlsx`, `Locoweed Field Volatiles 2019 Summary.xlsx` |

The historical dataset is excluded from primary analyses because it uses a different compound set and lacks pairing metadata, making it non-comparable to the 2022--2023 collections.

---

## Repository Structure

```
locoweed-volatiles-2026/
├── archive/
│   └── Locoweed All.Rmd            # Original analysis script (deprecated; see Known Issues)
├── data-files/
│   ├── 2022 Field Locoweed.xlsx
│   ├── 2023 Field Locoweed.xlsx
│   ├── Locoweed Field Volatiles 2018.xlsx
│   ├── Locoweed Field Volatiles 2019 Summary.xlsx
│   ├── Spring 2022 GH Locoweed RESULTS.xlsx
│   └── Spring GH 23 Locoweed Volatiles.xlsx
├── figures/
│   ├── Fig1_GH_PCA.png
│   ├── Fig2_Field_PCA.png
│   └── Fig3_Swainsonine.png
├── proc-docs/
│   ├── GH-data-proc.Rmd            # Greenhouse data import and cleaning
│   ├── GH-data-proc.md             # Knitted output
│   ├── Field-data-proc.Rmd         # Field data import and cleaning
│   └── Megan-data-proc.Rmd         # 2019 historical data processing
├── locoweed_volatiles_analysis.R   # Main analysis script
├── locoweed-volatiles-2026.Rproj
├── results_summary.txt
└── README.md
```

---

## Data Processing Pipeline

The main analysis script (`locoweed_volatiles_analysis.R`) calls the three processing documents via `knitr::knit()` before any analysis runs:

```r
knitr::knit("proc-docs/GH-data-proc.Rmd",    envir = globalenv())
knitr::knit("proc-docs/Megan-data-proc.Rmd", envir = globalenv())
knitr::knit("proc-docs/Field-data-proc.Rmd", envir = globalenv())
```

These produce three key objects used downstream:

| Object | Source | Description |
|---|---|---|
| `joined.data` | GH-data-proc.Rmd | Long-format greenhouse VOC data, raw ng/g/hr, both years |
| `field_comb` | Field-data-proc.Rmd | Long-format field VOC data, raw ng/g/hr, both years |
| `proc_2019` | Megan-data-proc.Rmd | CLR-processed 2019 field data; visualization only |

> **Note:** `GH-data-proc.Rmd` also produces `gh.joined.data.clr` and `field_comb_clr_final`, but these are not used in the main analysis. CLR is re-derived correctly from raw values in `locoweed_volatiles_analysis.R`.

### Compound Exclusions

**Greenhouse:** The following compound categories are excluded prior to analysis via string matching against compound names.

| Category | Examples |
|---|---|
| Synthetic UV filters / plasticizers | homosalate, ethylhexylsalicylate, diethylphthalate |
| Industrial / rubber additives | benzothiazole |
| SPME fiber contaminants | butylatedhydroxytoluene (BHT) |
| Synthetic acrylates | dodecylacrylate |
| Straight-chain alkane contaminants | hexadecane, tridecane, pentadecane, nanodecane |
| Anthropogenic aromatics | p-trimethyl (trimethylbenzene) |
| Not established as biogenic in *Oxytropis* | caprylyl acetate |

**Field:** Analysis is restricted to 10 compounds detected consistently across both years:

```
heptanal, sixmethyl, octanal, z3hex, dlimonene,
nonanal, decanal, methylsalicylate, cisthreehexiso, bocimene
```

---

## Statistical Approach

### 1. Compositional Data Transformation (CLR)

Raw VOC emission rates (ng/g/hr) are transformed using centered log-ratio (CLR) transformation:

> CLR(x)_i = log( x_i / geometric_mean(x) )

CLR is applied once, to raw values, after replacing true zeros with a pseudocount of 0.01. This pseudocount is consistent with the original processing pipeline and is required because CLR is undefined for zero values.

Euclidean distance on CLR-transformed data equals Aitchison distance, which is the appropriate dissimilarity metric for compositional data such as VOC profiles.

### 2. PERMANOVA

Multivariate VOC differences are tested using `adonis2` from the vegan package.

**Greenhouse model:**
```r
adonis2(gh_dist ~ swa + pair, data = gh_meta,
        permutations = 999, strata = gh_meta$year, by = "margin")
```

**Field model:**
```r
adonis2(field_dist ~ swa + collection.type + pair2, data = field_meta,
        permutations = 999, strata = field_meta$year, by = "margin")
```

- `by = "margin"`: marginal (Type III) testing -- each term is tested after all others, so results do not depend on term order
- `strata = year`: permutations restricted within years to account for year as a blocking factor; year is not included as a formula term because it is fully collinear with strata (Df = 0)
- Year-stratified PERMANOVAs (2022 and 2023 separately) are run in addition to the combined model to confirm consistency of the endophyte effect across years

### 3. Betadisper

Before interpreting PERMANOVA results, within-group dispersion is tested using `betadisper` + `permutest`. If groups differ in spread rather than centroid location, PERMANOVA p-values can be inflated.

Betadisper is run separately for each grouping variable (year, swa, collection type). A significant result (p < 0.05) is noted as a potential confounder and interpreted accordingly.

### 4. PCA

Principal components analysis is run on CLR-transformed data using `prcomp(scale. = TRUE)`. Biplots show individual plant scores colored by year, with endophyte status encoded by point shape, and loading arrows for the six most influential compounds (selected by Euclidean length in PC1-PC2 space).

Variance explained by PC1 and PC2 is extracted from the correct per-dataset PCA object -- greenhouse labels use `summary(gh_pca)`, field labels use `summary(field_pca)`.

### 5. SIMPER

SIMPER analysis (`simper`, vegan) is run on raw ng/g/hr values using Bray-Curtis dissimilarity, because CLR-transformed values can be negative and are not appropriate for Bray-Curtis. This identifies which compounds most contribute to dissimilarity between E+ and E- groups.

A parallel CLR-based compound ranking is also computed: each compound's contribution to squared Aitchison distance between group centroids is expressed as a percentage of the total squared centroid distance. This provides an Aitchison-consistent analog to SIMPER for ranking compounds by their contribution to compositional separation.

### 6. IndVal

Indicator value analysis (`multipatt`, indicspecies) identifies compounds significantly associated with E+ or E- status. `IndVal.g` is used with 999 permutations. This is useful for identifying candidate compounds for follow-up chemical ecology work.

---

## Key Analytical Decisions

### CLR applied once, to raw values

Earlier versions of the analysis (archived in `Locoweed All.Rmd`) applied transformations sequentially:

1. `Field-data-proc.Rmd` applied `log2()` before CLR, producing `CLR(log2(x))`
2. The main analysis then squared values and re-applied CLR, resulting in a triple transformation with no valid compositional interpretation

The current script applies CLR exactly once, to raw ng/g/hr values, before any distance matrix is computed. The CLR outputs from the processing Rmds (`gh.joined.data.clr`, `field_comb_clr_final`) are deliberately not used.

### `by = "margin"` used for both datasets

The archived script used `by = "terms"` (sequential, Type I) for the greenhouse PERMANOVA and `by = "margin"` for the field PERMANOVA. Sequential testing conflates term order with effect estimation -- the endophyte effect would be estimated before accounting for pair, which is not the intended test. Marginal testing is used consistently here for both datasets.

### `strata = year` vs. `strata = pair`

For a paired design, restricting permutations within pairs (`strata = pair`) would most directly preserve the pairing structure. However, `strata = year` is used here as a reasonable conservative choice given sample sizes. This is a known limitation and may be revisited if sample sizes permit within-pair permutation.

### Pseudocount value

Zeros are replaced with 0.01 prior to CLR. This value is carried forward from the original processing pipeline for consistency. The choice of pseudocount can influence CLR results, particularly for compounds with many zeros, and should be noted in any methods section.

### `set.seed(4721)`

A consistent seed is set before all permutation-based analyses (PERMANOVA, betadisper, IndVal) to ensure reproducibility.

---

## Known Data Issues

The following issues were identified during audit of `Locoweed All.Rmd` and require resolution before publication. The main analysis script corrects the statistical issues; the data integrity issues below still require manual verification.

### Action Required

#### 1. Undocumented correction factors in field bar plot data

The archived script applies unexplained numeric correction factors to 2022 field data in the bar plot pipeline:

| Sample | Correction applied |
|---|---|
| 2022 leaves | `mean_ngghr / 3`, `se / 3` |
| 2022 flowers | `mean_ngghr / 22`, `se / 22` |
| 2022 bocimene (flowers) | `mean_ngghr * 19`, `se * 5` |

No comments explain what these correct for (e.g., unit conversions, collection volume differences, instrument calibration). These factors must be traced back to a documented source and applied at the data processing stage rather than silently in a plot pipeline before publication.

#### 2. *d*-limonene relabeled as linalool in field plot

The archived field plot pipeline contains:

```r
mutate(compound = if_else(compound == "dlimonene", "linalool", compound))
```

These are chemically distinct compounds. This relabeling is either a legitimate compound identification correction (e.g., a misidentified peak in the original GC-MS run) or an error. The compound identity must be confirmed against the original chromatograms before publication.

#### 3. Hardcoded summary values in field data

Several rows are manually inserted into the field dataset via `bind_rows` with hardcoded mean and SE values (e.g., dlimonene 2023 E-, mean = 6.21, se = 1.5). No derivation or source is given for these values. They should be replaced with programmatic calculations from the raw data.

### Noted but Lower Priority

#### 4. Absolute paths in `GH-data-proc.Rmd`

`GH-data-proc.Rmd` reads from absolute paths pointing to a OneDrive directory outside the project root. `Field-data-proc.Rmd` correctly uses relative paths. The absolute paths in `GH-data-proc.Rmd` must be updated to point to `data-files/` within the project before this script will run on any other machine or after directory reorganization.

#### 5. Incomplete chunk in archived script

The "Greenhouse Table" chunk in `Locoweed All.Rmd` contains only `joined.data %>%` with no further code. This is a remnant of an incomplete draft and does not affect the current analysis.

#### 6. Inconsistent variable naming in archived script

The archived script mixes dot-separated (`permanova.output.field`) and underscore-separated (`permanova_output_gh`) naming conventions. The current script uses underscores throughout for consistency.
