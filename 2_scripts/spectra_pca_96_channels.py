import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import argparse as ap
import os

"""
Perform PCA on 96 channels spectra for all genomes
"""


def filter_supertable(df: pd.DataFrame, mut_num: int) -> tuple:
    """Remove isolates with too few mutations."""
    
    # 96 channels
    mutation_cols = [col for col in df.columns if "[" in col and "]" in col]
    
    df["total_mut"] = df[mutation_cols].sum(axis=1)

    df_filt = df[df["total_mut"] >= mut_num].copy()
    return df_filt

def pca_96_channels(df_filt: pd.DataFrame) -> pd.DataFrame:

    channels_96 = [col for col in df_filt.columns if "[" in col and "]" in col]
    X = df_filt[channels_96].fillna(0)

    pca = PCA(n_components=5)
    pcs = pca.fit_transform(X)

    df_filt["PC1"] = pcs[:,0]
    df_filt["PC2"] = pcs[:,1]
    df_filt["PC3"] = pcs[:,2]
    df_filt["PC4"] = pcs[:,3]
    df_filt["PC5"] = pcs[:,4]

    explained = pca.explained_variance_ratio_

    df_exp = pd.DataFrame({
        "PC": ["PC1","PC2","PC3","PC4","PC5"],
        "explained_variance": explained})


    return df_filt, df_exp

def plot_pca_colored(df_pca: pd.DataFrame, output_dir: str):
    """Plot PCA (PC1 vs PC2) colored by Collection_name."""
    plt.figure(figsize=(7,6))

    # unique cohorts
    cohorts = df_pca["Collection_name"].unique()

    for cohort in cohorts:
        subset = df_pca[df_pca["Collection_name"] == cohort]
        plt.scatter(
            subset["PC1"], 
            subset["PC2"], 
            label=cohort, 
            s=20, 
            alpha=0.7
        )

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA on 96-channel spectra (colored by Collection_name)")
    plt.legend(markerscale=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/pca_colored.png", dpi=300)
    plt.close()


def main():
    parser = ap.ArgumentParser(description="Perform PCA on 96 channels.")
    parser.add_argument("--supertable", required=True, help="Path to supertable .csv file.")
    parser.add_argument("--mut_num", type=int, required=True, help="Minimal number of mutations")
    parser.add_argument("--output_dir", required=True, help="Where to save output")

    args = parser.parse_args()
    # create output directory if missing
    os.makedirs(args.output_dir, exist_ok=True)

    # load table
    df = pd.read_csv(args.supertable)

    # filter
    df_filt = filter_supertable(df, args.mut_num)

    # PCA
    df_pca, df_exp = pca_96_channels(df_filt)

    plot_pca_colored(df_pca, args.output_dir)


    # save PCA scatter
    # plt.figure(figsize=(6,5))
    # plt.scatter(df_pca["PC1"], df_pca["PC2"])
    # plt.xlabel("PC1")
    # plt.ylabel("PC2")
    # plt.title("PCA on 96-channel mutational spectra")
    # plt.tight_layout()
    # plt.savefig(f"{args.output_dir}/pca_96_channels.png", dpi=300)
    # plt.close()

    # save table with PCs
    df_pca.to_csv(f"{args.output_dir}/supertable_with_pca.csv", index=False)
    df_exp.to_csv(f"{args.output_dir}/pca_explained_variance.csv", index=False)


if __name__ == "__main__":
    main()


