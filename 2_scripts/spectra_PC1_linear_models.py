import pandas as pd
import statsmodels.formula.api as smf
import argparse as ap
import os

"""
Fit models
"""

def fit_pc_models(df: pd.DataFrame):
    """
    Fit linear models for PC1.
    Returns dict: model_name -> summary_text
    """

    results = {}

    # Model 1: PC1 ~ AMR_burden
    m1 = smf.ols("PC1 ~ amr_burden", data=df).fit()
    results["PC1_model1"] = m1.summary().as_text()

    # Model 2: PC1 ~ AMR_burden + Study
    m2 = smf.ols("PC1 ~ amr_burden + C(Collection_name)", data=df).fit()
    results["PC1_model2"] = m2.summary().as_text()

    # Model 3: PC1 ~ AMR_burden * Collection_name (interaction)
    m3 = smf.ols("PC1 ~ amr_burden * C(Collection_name)", data=df).fit()
    results["PC1_model3"] = m3.summary().as_text()

    return results

def main():
    parser = ap.ArgumentParser(description="Fit linear models for PC1.")
    parser.add_argument("--pca_table", required=True, help="CSV with PC1/PC2 and metadata")
    parser.add_argument("--output_dir", required=True, help="Where to save model summaries")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.pca_table)

    # Fit models
    results = fit_pc_models(df)

    # Save results
    out_path = f"{args.output_dir}/PC1_linear_models.txt"
    with open(out_path, "w") as f:
        for name, summary in results.items():
            f.write(f"=== {name} ===\n")
            f.write(summary)
            f.write("\n\n")

    print(f"Saved model summaries to: {out_path}")


if __name__ == "__main__":
    main()