import pandas as pd
import statsmodels.formula.api as smf
import argparse as ap
import os
import matplotlib.pyplot as plt
import seaborn as sns




def fit_signature_models(df: pd.DataFrame, signatures: list):
    """
    Fit linear models for mutational signatures:
    C>A, C>T, T>G
    """

    results = {}

    for sig in signatures:
        # Model 1: signature ~ AMR_burden
        m1 = smf.ols(f"{sig} ~ amr_burden", data=df).fit()
        results[f"{sig}_model1"] = m1.summary().as_text()

        # Model 2: signature ~ AMR_burden + Collection_name
        m2 = smf.ols(f"{sig} ~ amr_burden + C(Collection_name)", data=df).fit()
        results[f"{sig}_model2"] = m2.summary().as_text()

        # Model 3: signature ~ AMR_burden * Collection_name
        m3 = smf.ols(f"{sig} ~ amr_burden * C(Collection_name)", data=df).fit()
        results[f"{sig}_model3"] = m3.summary().as_text()

    return results

def plot_signature_vs_amr(df, sig, output_dir):
    """Scatterplot signatures vs AMR_burden with regression line."""
    plt.figure(figsize=(6,5))
    sns.regplot(
        data=df,
        x="amr_burden",
        y=sig,
        scatter_kws={"alpha":0.5, "s":20},
        line_kws={"color":"red"}
    )
    plt.title(f"{sig} vs AMR_burden")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/{sig}_vs_amr.png", dpi=300)
    plt.close()

def plot_signature_boxplot(df, sig, output_dir):
    """Boxplot signatures per cohort."""
    plt.figure(figsize=(10,5))
    sns.boxplot(
        data=df,
        x="Collection_name",
        y=sig
    )
    plt.xticks(rotation=45, ha="right")
    plt.title(f"{sig} distribution across cohorts")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/{sig}_boxplot.png", dpi=300)
    plt.close()



def main():
    parser = ap.ArgumentParser(description="Fit linear models for mutational signatures.")
    parser.add_argument("--table6", required=True, help="CSV with 6-channel spectra + metadata")
    parser.add_argument("--mut_num", type=int, required=True, help="Minimal number of mutations")
    parser.add_argument("--signatures", required=True, help="Comma-separated list of signatures to analyze")
    parser.add_argument("--output_dir", required=True, help="Where to save model summaries")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(f"{args.output_dir}/models", exist_ok=True)
    os.makedirs(f"{args.output_dir}/plots", exist_ok=True)

    SIG_LIST = args.signatures.split(",")

    df = pd.read_csv(args.table6)
    df = df.rename(columns={
        "C>A": "C_A",
        "C>T": "C_T",
        "T>G": "T_G",
        "C>G": "C_G",
        "T>A": "T_A",
        "T>C": "T_C"
    })

    df["total_mut"] = df[["C_A", "C_T", "T_G", "C_G", "T_A", "T_C"]].sum(axis=1)
    df = df[df["total_mut"] >= args.mut_num]

    results = fit_signature_models(df, SIG_LIST)
    # Plotting
    for sig in SIG_LIST:
        plot_signature_vs_amr(df, sig, args.output_dir)
        plot_signature_boxplot(df, sig, args.output_dir)

    out_path = f"{args.output_dir}/models/signature_models.txt"
    with open(out_path, "w") as f:
        for name, summary in results.items():
            f.write(f"=== {name} ===\n")
            f.write(summary)
            f.write("\n\n")

    print(f"Saved signature models to: {out_path}")


if __name__ == "__main__":
    main()
