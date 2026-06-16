import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from umap import UMAP
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import multipletests
import plotly.express as px
from pathlib import Path
from typing import Literal
import argparse as ap


"""
This script is for analyzing the result of AMRFinderPlus and creating visualizations.

Here are functions for:

1. Merging AMRFinderPlus result with metadata
2. Plotting frequency of AMR classes, subclasses, and genes
3. Plotting AMR burden (number of unique AMR genes per isolate)
4. Plotting MDR levels (number of unique AMR classes per isolate)
5. Plotting MDR levels per collection
6. Plotting comparison of AMR classes between studies
7. PCA and UMAP of AMR presence/absence profiles
8. Heatmap of presence/absence of AMR genes across isolates
9. 3D UMAP of AMR profiles colored by metadata
10. Fold-change plots comparing two collections
"""

ANALYSIS_DIR = Path("1_data/3_results/amrfinder_analysis")
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)


# =====================
# MERGING AMRFINDERPLUS RESULTS
# =====================
def merge_amrfinder_results_metadata(report :str, metadata :str, output_dir :str) -> pd.DataFrame:
    """
    Merges AMRFinderPlus report with metadata and saves the merged DataFrame to a CSV file.
    Parameters:
    - report: Path to the AMRFinderPlus result file (CSV format).
    - metadata: Path to the metadata file (Excel format).
    - output: Path to the output merged CSV file.
    """
    out_path = output_dir / "amr_metadata_merged.csv"


    # Load the AMRFinderPlus report and metadata
    amr_df = pd.read_csv(report, sep="\t")
    amr_df["Sanger_ID"] = amr_df["Name"].str.split("_contigs").str[0]
    meta_df = pd.read_excel(metadata)

    # Merge
    merged_df = pd.merge(amr_df, meta_df, on='Sanger_ID', how='left')

    # Save the merged DataFrame to a CSV file
    merged_df.to_csv(out_path, index=False)
    print(f"Merged AMRFinderPlus results with metadata saved to {out_path}")

    return merged_df


# =====================
# PLOTTING
# =====================

def plot_frequency(df: pd.DataFrame, column: str, filename: str, title: str, output_dir: str):
    """
    Generic frequency plot for any categorical AMR column.
    Example usage:
        plot_frequency(df, "Class", "class_frequency.png", "AMR Class Frequency")
        plot_frequency(df, "Subclass", "subclass_frequency.png", "AMR Subclass Frequency")
    """

    # number of isolates
    total_isolates = df["Sanger_ID"].nunique()
    freq = (
        df.groupby("Sanger_ID")[column]
        .unique()
        .explode()
        .value_counts()
        .sort_values(ascending=False))
    percent = (freq / total_isolates * 100).round(1)


    height = max(6, len(freq) * 0.35)
    plt.figure(figsize=(12, height))
    ax = sns.barplot(x=freq.values, y=freq.index, color="#1447E6")

    # labels: "n (x%)"
    labels = [f"{n} ({p}%)" for n, p in zip(freq, percent)]
    # add labels at the end of bars
    for i, (value, label) in enumerate(zip(freq.values, labels)):
        ax.text(value + 0.5, i, label, va='center', fontsize=9)

    plt.title(title)
    plt.xlabel("Count", fontweight="bold")
    plt.ylabel(column, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_dir / "plots" / filename, dpi=300)
    plt.close()



def plot_amr_burden(df: pd.DataFrame, output_dir: str):
    """
    Plot distribution of AMR burden (number of unique AMR genes per isolate).
    """
    burden = df.groupby("Sanger_ID")["Element symbol"].nunique()

    plt.figure(figsize=(8, 5))
    sns.histplot(burden, bins=20, kde=True, color="purple")
    plt.title("AMR burden per isolate")
    plt.xlabel("Number of AMR genes", fontweight="bold")
    plt.ylabel("Number of isolates", fontweight="bold")
    plt.tight_layout()
    plt.grid(alpha=0.3)
    plt.savefig(output_dir / "plots" / "amr_burden.png", dpi=300)
    plt.close()


def plot_mdr_levels(df: pd.DataFrame, output_dir: str):
    """
    Plot distribution of MDR levels:
    MDR3 = isolates with 3 AMR classes
    MDR4 = isolates with 4 AMR classes
    etc.
    """

    counts = (
        df.groupby("MDR_level")["Sanger_ID"]
           .nunique()
           .reset_index(name="count")
    )

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x=counts["MDR_level"], y=counts["count"], color="#1447E6")

    plt.title("Distribution of MDR levels (unique AMR classes per isolate)")
    plt.xlabel("Number of AMR classes (MDR level)", fontweight="bold")
    plt.ylabel("Number of isolates", fontweight="bold")
    plt.tight_layout()
    plt.grid(alpha=0.3)
    plt.savefig(output_dir / "plots" / "mdr_levels.png", dpi=300)
    plt.close()

def plot_mdr_per_collection(df: pd.DataFrame, output_dir: str):
    """
    Plot distribution of MDR levels per collection.
    """

    # Count how many genomes have a given MDR level in each collection
    counts = (
        df.groupby(["Collection_name", "MDR_level"])["Sanger_ID"]
           .nunique()
           .reset_index(name="count")
    )

    plt.figure(figsize=(12, 6))
    sns.barplot(
        data=counts,
        x="MDR_level",
        y="count",
        hue="Collection_name",
        palette="tab20"
    )
    plt.title("MDR level distribution per collection")
    plt.xlabel("MDR level", fontweight="bold")
    plt.ylabel("Number of genomes", fontweight="bold")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_dir / "plots" / "mdr_per_collection.png", dpi=300)
    plt.close()


def plot_study_comparison(df: pd.DataFrame, output_dir: str):
    """
    Plot how many isolates in each study have each AMR class.
    Counts unique classes per isolate, not number of genes.
    """

    # Unique classes per isolate and study
    unique_classes = (
        df.groupby(["Sanger_ID", "Collection_name"])["Class"]
        .unique()
        .explode()
        .reset_index()
    )

    # Count how many isolates have each class in each study
    study_class = (
        unique_classes.groupby(["Collection_name", "Class"])["Sanger_ID"]
        .nunique()
        .reset_index(name="Isolate_count")
    )

    plt.figure(figsize=(14, 7))
    sns.barplot(
        data=study_class,
        x="Collection_name",
        y="Isolate_count",
        hue="Class"
    )

    plt.title("Number of isolates with each AMR class per study")
    plt.xlabel("Study", fontweight="bold")
    plt.ylabel("Number of isolates", fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_dir / "plots" / "study_class_comparison.png", dpi=300)
    plt.close()


def plot_pca_amr(df: pd.DataFrame, color_by: str, output_dir: str):
    """
    PCA of AMR presence/absence profiles with optional coloring.
    """

    # presence/absence matrix
    pivot = (
        df.assign(present=1)
          .pivot_table(
              index="Sanger_ID",
              columns="Element symbol",
              values="present",
              aggfunc=lambda x:1,
              fill_value=0
          )
    )

    # PCA
    pca = PCA(n_components=2)
    coords = pca.fit_transform(pivot)
    pca_df = pd.DataFrame(coords, columns=["PC1", "PC2"], index=pivot.index)

    # coloring points
    if color_by is not None:
        # add metadata to PCA dataframe
        meta = df.drop_duplicates("Sanger_ID").set_index("Sanger_ID")
        pca_df[color_by] = meta[color_by]

    plt.figure(figsize=(8, 6))

    if color_by is None:
        sns.scatterplot(data=pca_df, x="PC1", y="PC2")
    else:
        sns.scatterplot(
            data=pca_df,
            x="PC1",
            y="PC2",
            hue=color_by,
            palette="tab20",
            s=60
        )

    plt.title(f"PCA of AMR profiles (colored by {color_by})" if color_by else "PCA of AMR profiles")
    plt.tight_layout()
    plt.savefig(output_dir / "plots" / f"pca_amr_{color_by if color_by else 'plain'}.png", dpi=300)
    plt.close()


def plot_heatmap_presence_absence(df: pd.DataFrame, output_dir: str):
    """
    Heatmap of presence/absence of AMR genes across isolates.
    Rows = isolates, columns = genes, color = presence/absence.
    """
    pivot = df.pivot_table(
        index="Sanger_ID",
        columns="Element symbol",
        values="Start",
        aggfunc=lambda x: 1,   # 1 when present, 0 when absent
        fill_value=0
    )

    plt.figure(figsize=(30, 50))
    sns.heatmap(pivot, cmap="viridis", cbar=False)
    plt.title("Presence/absence heatmap of AMR genes")
    plt.tight_layout()
    plt.savefig(output_dir / "plots" / "heatmap_presence_absence.png", dpi=300)
    plt.close()


def plot_umap_amr(df: pd.DataFrame, color_by: str, output_dir: str):
    """
    UMAP of AMR presence/absence profiles with optional coloring.
    Parameters:
    - df: DataFrame containing AMR data and metadata
    - color_by: column name in metadata to color points by (e.g. "Collection_name", "MDR_level", "ST")
    """

    # presence/absence matrix
    pivot = (
        df.assign(present=1)
          .pivot_table(
              index="Sanger_ID",
              columns="Element symbol",
              values="present",
              aggfunc="max",
              fill_value=0
          )
    )

    # UMAP
    umap = UMAP(
        n_neighbors=15,
        min_dist=0.1,
        metric="hamming",
        random_state=42
    )
    coords = umap.fit_transform(pivot)

    umap_df = pd.DataFrame(coords, columns=["UMAP1", "UMAP2"], index=pivot.index)

    # Color points
    if color_by is not None:
        meta = df.drop_duplicates("Sanger_ID").set_index("Sanger_ID")
        umap_df[color_by] = meta[color_by]

    plt.figure(figsize=(8, 6))

    if color_by is None:
        sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", s=60)
    else:
        sns.scatterplot(
            data=umap_df,
            x="UMAP1",
            y="UMAP2",
            hue=color_by,
            palette="tab20",
            s=60
        )

    plt.title(f"UMAP of AMR profiles (colored by {color_by})" if color_by else "UMAP of AMR profiles")
    plt.tight_layout()
    plt.savefig(output_dir / "plots" / f"umap_2d_amr_{color_by if color_by else 'plain'}.png", dpi=300)
    plt.close()

    return umap_df


def plot_umap_3d(df: pd.DataFrame, color_by: str, output_dir: str):
    """
    3D UMAP of AMR presence/absence profiles colored by a metadata column.
    Parameters:
    - df: DataFrame containing AMR data and metadata
    - color_by: column name in metadata to color points by (e.g. "Collection_name")
    """
    # pivot presence/absence
    pivot = (
        df.assign(present=1)
          .pivot_table(
              index="Sanger_ID",
              columns="Element symbol",
              values="present",
              aggfunc="max",
              fill_value=0
          )
    )

    # UMAP 3D
    reducer = UMAP(
        n_components=3,
        random_state=42,
        metric="hamming"
    )
    embedding = reducer.fit_transform(pivot)

    # df with results
    umap_df = pd.DataFrame({
        "Sanger_ID": pivot.index,
        "UMAP1": embedding[:, 0],
        "UMAP2": embedding[:, 1],
        "UMAP3": embedding[:, 2]
    })

    # add metadata for coloring
    umap_df = umap_df.merge(
        df.drop_duplicates("Sanger_ID")[["Sanger_ID", color_by]],
        on="Sanger_ID",
        how="left"
    )

    # interactive 3D plot with plotly
    fig = px.scatter_3d(
        umap_df,
        x="UMAP1",
        y="UMAP2",
        z="UMAP3",
        color=color_by,
        hover_name="Sanger_ID",
        title=f"3D UMAP colored by {color_by}"
    )

    fig.write_html(output_dir / "plots" / f"umap_3d_{color_by}.html")
    print(f"Saved 3D UMAP to umap_3d_{color_by}.html")


def plot_foldchange(df, collection_a, collection_b, output_dir: str, col="Class"):
    """
    Plot log2 fold-change of proportion of isolates with each AMR class between two collections.
    Parameters:
    - df: DataFrame containing AMR data and metadata
    - collection_a: name of the first collection (e.g. "HICF") - new
    - collection_b: name of the second collection (e.g. "Prospective_study") - old/reference
    - col: column to compare (e.g. "Class", "Subclass")
    """

    # Count how many isolates have each class in each collection
    counts = (
        df.groupby(["Collection_name", col])["Name"]
        .nunique()
        .reset_index(name="n")
    )

    # Count total number of genomes in each collection
    totals = (
        df.groupby("Collection_name")["Name"]
        .nunique()
        .reset_index(name="total_genomes")
    )
    counts = counts.merge(totals, on="Collection_name")
    # Proportions
    counts["prop"] = counts["n"] / counts["total_genomes"]
    pivot = counts.pivot(index=col, columns="Collection_name", values="prop").fillna(0)

    # pseudo-count to avoid division by zero, should be small compared to the smallest non-zero proportion
    eps = 0.01

    # log2 fold-change
    pivot["log2FC"] = np.log2((pivot[collection_a] + eps) / (pivot[collection_b] + eps))

    # sorting
    pivot = pivot.sort_values("log2FC", ascending=True)

    # plotting
    plt.figure(figsize=(10, 8))
    plt.barh(pivot.index, pivot["log2FC"], color="#0066FF")
    plt.margins(x=0.2)
    plt.xlabel(f"Log2FC ({collection_a} / {collection_b})", fontweight="bold")
    plt.ylabel(f"{col}", fontweight="bold")

    # etiqutes with log2FC values
    for y, x in zip(pivot.index, pivot["log2FC"]):
        plt.text(x + 0.05, y, f"{x:.2f}", va="center")

    plt.tight_layout()
    plt.grid(alpha=0.3)
    plt.savefig(output_dir / "plots" / f"foldchange_{col}.png", dpi=300)
    plt.close()

    return pivot

def test_amr_change(df, collection_a, collection_b, direction, output_dir: str, col="Class") -> pd.DataFrame:
    """
    Test if the proportion of isolates with each AMR class is significantly different between two collections using z-test for proportions.
    Parameters:
    - df: DataFrame containing AMR data and metadata
    - collection_a: name of the first collection (e.g. "HICF") - new
    - collection_b: name of the second collection (e.g. "Prospective_study") - old
    - col: column to compare (e.g. "Class", "Subclass")
    - direction: "smaller" means p1 < p2, "larger" means p1 > p2, "two-sided" means p1 != p2
    """

    # presence/absence per genome
    pa = (
        df.assign(present=1)
          .pivot_table(
              index="Sanger_ID",
              columns=col,
              values="present",
              aggfunc="max",
              fill_value=0
          )
    )

    # add metadata for collections
    meta = df.drop_duplicates("Sanger_ID")[["Sanger_ID", "Collection_name"]]
    pa = pa.merge(meta, on="Sanger_ID")

    results = []

    for class_name in pa.columns:
        if class_name in ["Sanger_ID", "Collection_name"]:
            continue

        # liczby dla 2012
        k1 = pa.loc[pa["Collection_name"] == collection_a, class_name].sum()
        n1 = (pa["Collection_name"] == collection_a).sum()

        # liczby dla 2019
        k2 = pa.loc[pa["Collection_name"] == collection_b, class_name].sum()
        n2 = (pa["Collection_name"] == collection_b).sum()

        # test proporcji
        stat, pval = proportions_ztest([k1, k2], [n1, n2], alternative=direction)

        results.append({
            "Class": class_name,
            f"{collection_a}_n": n1,
            f"{collection_a}_k": k1,
            f"{collection_a}_prop": k1 / n1 if n1 > 0 else None,
            f"{collection_b}_n": n2,
            f"{collection_b}_k": k2,
            f"{collection_b}_prop": k2 / n2 if n2 > 0 else None,
            "z_stat": stat,
            "p_value": pval
        })

    results_df = pd.DataFrame(results)
    mask = results_df["p_value"].notna()
    results_df.loc[mask, "p_adj"] = multipletests(results_df.loc[mask, "p_value"], method="fdr_bh")[1]

    results_df = results_df.sort_values("p_adj")
    results_df.to_csv(output_dir / "tests" / f"amr_change_test_{col}_{collection_a}_{collection_b}_{direction}.csv", index=False)

    return results_df


def main():
    parser = ap.ArgumentParser(description="Analysis of AMRFinderPlus results")
    parser.add_argument("-r", required=True, help="Path to the AMRFinderPlus result file (CSV format)")
    parser.add_argument("-m", required=True, help="Path to the metadata file (Excel format)")
    parser.add_argument("-o", required=True, help="Directory to save analysis results and plots")

    args = parser.parse_args()

    output_dir = Path(args.o)
    output_dir.mkdir(parents=True, exist_ok=True)

    merged_df = merge_amrfinder_results_metadata(args.r, args.m, output_dir)
    merged_df["MDR_level"] = merged_df.groupby("Sanger_ID")["Class"].transform(lambda x: len(set(x)))
    burden = (
        merged_df.groupby("Sanger_ID")["Element symbol"]
        .nunique()
        .rename("amr_burden")
        )
    merged_df = merged_df.merge(burden, on="Sanger_ID")
    merged_df.to_csv(output_dir / "amr_metadata_mdr_burden.csv", index=False)

    output_dir.joinpath("plots").mkdir(exist_ok=True)
    output_dir.joinpath("tests").mkdir(exist_ok=True)


    plot_frequency(merged_df, "Class", "class_frequency.png", "AMR Class Frequency", output_dir=output_dir)
    plot_frequency(merged_df, "Subclass", "subclass_frequency.png", "AMR Subclass Frequency", output_dir=output_dir)
    plot_frequency(merged_df, "Element symbol", "gene_frequency.png", "AMR Gene Frequency", output_dir=output_dir)

    plot_amr_burden(merged_df, output_dir=output_dir)
    plot_mdr_levels(merged_df, output_dir=output_dir)
    plot_mdr_per_collection(merged_df, output_dir=output_dir)

    plot_study_comparison(merged_df, output_dir=output_dir)

    plot_pca_amr(merged_df, color_by="Collection_name", output_dir=output_dir)
    plot_pca_amr(merged_df, color_by="MDR_level", output_dir=output_dir)


    plot_umap_amr(merged_df, color_by="Collection_name", output_dir=output_dir)
    plot_umap_amr(merged_df, color_by="MDR_level", output_dir=output_dir)
    plot_umap_amr(merged_df, color_by="amr_burden", output_dir=output_dir)
    plot_umap_amr(merged_df, color_by="Collection_geographical_range", output_dir=output_dir)
    umap_df = plot_umap_amr(merged_df, color_by="ST", output_dir=output_dir)
    # Checking if UMAP clusters correspond to STs
    #upper_cluster = umap_df[umap_df["UMAP2"] > umap_df["UMAP2"].quantile(0.75)]
    #print(upper_cluster["ST"].value_counts())

    
    plot_heatmap_presence_absence(merged_df, output_dir=output_dir)
    

    plot_umap_3d(merged_df, color_by="Collection_name", output_dir=output_dir)
    plot_umap_3d(merged_df, color_by="MDR_level", output_dir=output_dir)
    plot_umap_3d(merged_df, color_by="Collection_geographical_range", output_dir=output_dir)
    plot_umap_3d(merged_df, color_by="ST", output_dir=output_dir)
    plot_umap_3d(merged_df, color_by="amr_burden", output_dir=output_dir)


    plot_foldchange(merged_df, collection_a="HICF", collection_b="Prospective_study", output_dir=output_dir, col="Class",)

    test_amr_change(merged_df, collection_a="HICF", collection_b="Prospective_study", direction="two-sided", output_dir=output_dir)
    test_amr_change(merged_df, collection_a="HICF", collection_b="Prospective_study", direction="smaller", output_dir=output_dir)
    test_amr_change(merged_df, collection_a="HICF", collection_b="Prospective_study", direction="larger", output_dir=output_dir)

    print(f"DONE! All analyses completed and plots saved to {output_dir}")



if __name__ == "__main__":
    main()