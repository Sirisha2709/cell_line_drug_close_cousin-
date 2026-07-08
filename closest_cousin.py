"""
Closest Cousin Cell Line Finder
--------------------------------
Identifies the cell line whose drug-response pattern most closely matches
a given patient's drug-response pattern, using binarized (sensitive/resistant)
IC50 calls and K-Nearest Neighbors with Hamming distance.

Motivation: preclinical cell line models are often chosen based on tissue
type or genomic similarity alone. This approach instead selects a "cousin"
cell line based on its functional drug-response profile, which may better
reflect real phenotypic behavior for a specific patient.

Expected input files:

1. Cell line drug response file (Excel), with columns:
   - CELL_LINE_NAME : name of the cell line (e.g. from GDSC)
   - DRUG_NAME      : name of the drug tested
   - LN_IC50        : natural-log-transformed IC50 value

2. Patient drug response file (Excel), with columns:
   - DRUG           : name of the drug tested on/for the patient
   - Response       : the patient's IC50 (or IC50-equivalent) value for that drug

Data source note: cell line IC50 data of this kind is commonly sourced from
GDSC (Genomics of Drug Sensitivity in Cancer, https://www.cancerrxgene.org).
Raw data files are not included in this repository -- see README for details.
"""

import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors


# ---------------------------------------------------------------------------
# Config -- update these paths to point at your local data files
# ---------------------------------------------------------------------------
CELLLINE_FILE = "data/cellline_ic50.xlsx"
PATIENT_FILE = "data/patient_ic50.xlsx"
N_NEIGHBORS = 10


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def standardised_drug_name(name):
    """
    Normalize drug names so the same compound isn't treated as different
    drugs due to formatting differences (case, spacing, hyphens, or salt
    forms such as '...HCl').
    """
    return str(name).lower().strip().replace("-", "").replace(" ", "").split("hcl")[0].strip()


def remove_outliers_iqr(df, column="LN_IC50"):
    """
    Remove outliers from each drug's IC50 distribution independently, using
    the standard 1.5*IQR rule. Done per-drug (not globally) since different
    drugs have very different IC50 ranges.
    """
    cleaned = []
    for drug, group in df.groupby("DRUG"):
        q1 = group[column].quantile(0.25)
        q3 = group[column].quantile(0.75)
        iqr = q3 - q1
        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr
        cleaned.append(group[(group[column] >= lower) & (group[column] <= upper)])
    return pd.concat(cleaned)


# ---------------------------------------------------------------------------
# 1. Load and clean cell line drug response data
# ---------------------------------------------------------------------------
df_raw = pd.read_excel(CELLLINE_FILE)

# Strip spaces/hyphens from cell line names for consistent matching
df_raw["CELL_LINE_NAME"] = df_raw["CELL_LINE_NAME"].str.replace(r"[-\s]", "", regex=True)

df = df_raw.rename(columns={"CELL_LINE_NAME": "CELL_LINE", "DRUG_NAME": "DRUG", "LN_IC50": "LN_IC50"})
df["DRUG"] = df["DRUG"].apply(standardised_drug_name)

# Remove per-drug outliers before computing reference medians
df_clean = remove_outliers_iqr(df, column="LN_IC50")

# Median LN_IC50 per drug across all (cleaned) cell lines -- used as the
# sensitivity/resistance threshold for both cell lines and the patient
drug_medians = df_clean.groupby("DRUG")["LN_IC50"].median()

# Binarize: 1 = sensitive (IC50 at or below the drug's median), 0 = resistant
df_clean["RESPONSE"] = df_clean.apply(
    lambda row: 1 if row["LN_IC50"] <= drug_medians[row["DRUG"]] else 0, axis=1
)

# ---------------------------------------------------------------------------
# 2. Build the cell line x drug binary response matrix
# ---------------------------------------------------------------------------
drug_cellline_matrix = df_clean.pivot_table(
    index="CELL_LINE",
    columns="DRUG",
    values="RESPONSE",
    aggfunc="first",
)

# Keep only drugs tested in at least 50% of cell lines
drug_cellline_matrix = drug_cellline_matrix.dropna(
    axis=1,
    thresh=int(0.5 * drug_cellline_matrix.shape[0]),
)

# Fill any remaining gaps with the per-drug median response
drug_cellline_matrix = drug_cellline_matrix.fillna(drug_cellline_matrix.median())

print(f"Cell lines: {drug_cellline_matrix.shape[0]}")
print(f"Drugs: {drug_cellline_matrix.shape[1]}")

# ---------------------------------------------------------------------------
# 3. Load and binarize patient drug response data
# ---------------------------------------------------------------------------
df_pat = pd.read_excel(PATIENT_FILE)
df_pat["DRUG"] = df_pat["DRUG"].apply(standardised_drug_name)
df_pat = df_pat.rename(columns={"Response": "LN_IC50"})

# Binarize the patient's response using the SAME per-drug medians computed
# from the cell line cohort, so patient and cell line calls are on the same
# reference scale. Drugs not present in the cell line cohort get NaN.
df_pat["RESPONSE"] = df_pat.apply(
    lambda row: (1 if row["LN_IC50"] <= drug_medians[row["DRUG"]] else 0)
    if row["DRUG"] in drug_medians.index
    else np.nan,
    axis=1,
)
patient = df_pat.set_index("DRUG")["RESPONSE"]

# ---------------------------------------------------------------------------
# 4. Restrict to drugs available for both patient and cell lines
# ---------------------------------------------------------------------------
common_drugs = drug_cellline_matrix.columns.intersection(patient.index)
drug_cellline_matrix_subset = drug_cellline_matrix[common_drugs]
patient_vector = patient[common_drugs].values

# ---------------------------------------------------------------------------
# 5. Find nearest cell lines by drug-response pattern (Hamming distance)
# ---------------------------------------------------------------------------
# Hamming distance is used because the response vectors are binary
# (sensitive/resistant); it simply counts the fraction of drugs where the
# patient and a given cell line disagree.
knn = NearestNeighbors(n_neighbors=N_NEIGHBORS, metric="hamming")
knn.fit(drug_cellline_matrix_subset.values)

dists, idxs = knn.kneighbors([patient_vector])

# ---------------------------------------------------------------------------
# 6. Report results
# ---------------------------------------------------------------------------
results = pd.DataFrame({
    "Cell_Line": [drug_cellline_matrix_subset.index[i] for i in idxs[0]],
    "Distance": dists[0].round(3),
    "Similarity": ((1 - dists[0]) * 100).round(1),
})

closest = results.iloc[0]["Cell_Line"]

comparison = pd.DataFrame({
    "Patient": patient[common_drugs],
    "Closest_cousin": drug_cellline_matrix.loc[closest, common_drugs],
})
comparison["Match"] = comparison["Patient"] == comparison["Closest_cousin"]

print("\n==== Drugs where patient and closest cousin disagree ====")
print(comparison[comparison["Match"] == False])

print("\nPercentage Similarity in Response:", comparison["Match"].sum() / len(comparison) * 100)

print("\n==== Closest Cousins (top {}) ====".format(N_NEIGHBORS))
print(results.to_string(index=False))
