import pandas as pd
from pathlib import Path
import argparse as ap

"""
For generating files indicating to which group(label) each leaf belongs
(for calculating splitted spectra)
"""

def generate_labels(amr_metadata_file, metadata_col, output_dir):
    """
    Create labels file: Branch, group
    """
    df = pd.read_csv(amr_metadata_file)
    labels_df = df[["Name", metadata_col]].copy()
    labels_df["Name"] = labels_df["Name"].astype(str) + ".fasta"
    #labels_df.columns = ["Branch", "Label"]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    out_path = output_dir / f"labels_{metadata_col}.csv"
    labels_df.to_csv(out_path, index=False)

    return out_path

def main():
    parser = ap.ArgumentParser(description="Generate labels for MutTui analysis")
    parser.add_argument("-i", "--input", required=True, help="Path to AMR metadata CSV file")
    parser.add_argument("-c", "--column", required=True, help="Metadata column to use for labels")
    parser.add_argument("-o", "--output", required=True, help="Output directory for labels CSV")
    args = parser.parse_args()

    label_file = generate_labels(args.input, args.column, args.output)
    print(f"Labels generated and saved to {label_file}")

if __name__ == "__main__":
    main()