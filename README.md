# Closest Cousin Cell Line Finder

Identifies the cancer cell line whose drug-response profile most closely
matches a given patient's drug-response profile, based on binarized
(sensitive vs. resistant) IC50 calls compared using K-Nearest Neighbors
with Hamming distance.

## Motivation

Preclinical cell line models are often selected based on tissue type or
genomic similarity to a patient's tumor. This project instead selects a
"cousin" cell line based on **functional drug-response similarity** —
i.e., a cell line that responds to the same panel of drugs in the same
sensitive/resistant pattern as the patient — which may better reflect a
patient's real phenotypic drug behavior for downstream experimental or
translational use.

## Method

1. **Data cleaning** — cell line names and drug names are standardized
   (removing formatting inconsistencies such as spacing, hyphens, and
   salt-form suffixes like "...HCl") so the same entity isn't
   inadvertently split into multiple labels.
2. **Outlier removal** — IC50 outliers are removed independently for each
   drug (1.5 × IQR rule), since different drugs have very different IC50
   ranges and a single global threshold would be inappropriate.
3. **Binarization** — for each drug, the median LN_IC50 across the cell
   line cohort is used as a sensitivity threshold: cell lines (and later,
   the patient) at or below the median are labeled "sensitive" (1),
   otherwise "resistant" (0). The patient's response is binarized using
   the **same cell-line-derived medians**, so both are on a consistent
   reference scale.
4. **Matrix construction** — a cell line × drug binary response matrix is
   built. Drugs tested in fewer than 50% of cell lines are dropped;
   remaining missing values are filled with the per-drug median response.
5. **Matching** — the patient's binary response vector is compared against
   all cell lines using **Hamming distance** (the fraction of drugs where
   two binary vectors disagree), via `sklearn.neighbors.NearestNeighbors`.
   The closest cell lines are returned, ranked by similarity.

## Output

- A ranked table of the top N closest cell lines, with distance and
  similarity percentage.
- A list of the specific drugs (if any) where the patient and their single
  closest "cousin" cell line disagree on sensitive/resistant classification.
- An overall percentage similarity in response across all shared drugs.

## Input Data

This repository does not include raw data files.

**Cell line drug response file** (Excel), expected columns:
| Column | Description |
|---|---|
| `CELL_LINE_NAME` | Name of the cell line |
| `DRUG_NAME` | Name of the drug tested |
| `LN_IC50` | Natural-log-transformed IC50 value |

Data of this format is commonly sourced from
[GDSC (Genomics of Drug Sensitivity in Cancer)](https://www.cancerrxgene.org).

**Patient drug response file** (Excel), expected columns:
| Column | Description |
|---|---|
| `DRUG` | Name of the drug |
| `Response` | Patient's IC50 (or IC50-equivalent) value for that drug |

Place both files under a local `data/` directory (not tracked in this
repo), or update the file paths at the top of `closest_cousin.py`.

## Requirements

```
pandas
numpy
scikit-learn
openpyxl
```

## Usage

```bash
python closest_cousin.py
```

## Notes / Limitations

- Binarization at the cohort median is a simplification of a continuous
  drug-response measurement; borderline IC50 values near the median may
  flip sensitivity calls with small changes in the underlying data.
- Hamming distance weighs all drugs equally; it does not account for
  differing clinical importance or confidence across drugs.
- Results depend on the overlap between drugs tested on the patient and
  drugs available in the cell line panel (`common_drugs`); a small overlap
  will reduce the reliability of the match.
