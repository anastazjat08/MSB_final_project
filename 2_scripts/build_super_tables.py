import pandas as pd
from pathlib import Path
import argparse as ap

"""
Module for building tabels of spectra for all labels (genomes) for downstream analysis
"""


def collapse_96_to_6(df):
    six = {"C>A": 0, "C>G": 0, "C>T": 0, "T>A": 0, "T>C": 0, "T>G": 0}

    for col in df.columns:
        if "[" in col and "]" in col:
            sub = col.split("[")[1].split("]")[0]
            six[sub] += df[col].iloc[0]

    return pd.DataFrame([six])


def build_master_table(
    spectra_dir,
    amr_metadata_file,
    output_dir,
    name_col_meta="Sanger_ID",
    class_col="Class",
    build_profiles=False,
    collapse_to_six=False
):

    spectra_dir = Path(spectra_dir)
    type_split = spectra_dir.name

    spectra_files = [
        f for f in spectra_dir.glob("mutational_spectrum_label_*.csv")
        if not f.name.endswith(f"{type_split}.csv")]
    spectra_list = []

    for f in spectra_files:
        df = pd.read_csv(f)
        df = df.set_index("Substitution").T

        label = f.stem.replace("mutational_spectrum_label_", "")

        if collapse_to_six:
            df6 = collapse_96_to_6(df)
            df6["Label_ID"] = label
            spectra_list.append(df6)
        else:
            df["Label_ID"] = label
            spectra_list.append(df)

    spectra = pd.concat(spectra_list, ignore_index=True)

    # OUTPUT DIR for this split type
    out_dir = Path(output_dir) / type_split
    out_dir.mkdir(parents=True, exist_ok=True)

    # If no profiles → save and exit
    if not build_profiles:
        suffix = "6" if collapse_to_six else "96"
        output_path = out_dir / f"{type_split}_{suffix}.csv"
        spectra.to_csv(output_path, index=False)
        return spectra

    # Load metadata
    meta = pd.read_csv(amr_metadata_file)
    meta[name_col_meta] = meta[name_col_meta].astype(str)

    # Build Profile = set(Class) per genome
    profiles = (
        meta.groupby(name_col_meta)[class_col]
            .apply(lambda x: set(x.dropna().astype(str)))
            .reset_index()
            .rename(columns={class_col: "Profile"})
    )

    # Merge ONLY Label_ID with Profile
    master = spectra.merge(
        profiles,
        left_on="Label_ID",
        right_on="Sanger_ID",
        how="left"
    ).drop(columns=[name_col_meta])

    # Add metadata for future fitting models
    meta_cols = ["Sanger_ID", "MDR_level", "amr_burden", "ST", "Collection_name"]
    meta_sub = meta[meta_cols].copy()

    master = master.merge(
        meta_sub,
        left_on="Label_ID",
        right_on="Sanger_ID",
        how="left"
    ).drop(columns=["Sanger_ID"])

    suffix = "6" if collapse_to_six else "96"
    output_path = out_dir / f"{type_split}_{suffix}.csv"
    master.to_csv(output_path, index=False)

    return master


def main():
    parser = ap.ArgumentParser(description="Build master tables of mutational spectra.")
    parser.add_argument("--spectra_dir", required=True, help="Folder with MutTui spectra")
    parser.add_argument("--amr_metadata", required=True, help="Metadata CSV")
    parser.add_argument("--output_dir", required=True, help="Where to save output")
    parser.add_argument("--build_profiles", action="store_true", help="Add AMR profiles")
    parser.add_argument("--collapse_to_six", action="store_true", help="Collapse 96->6 channels")

    args = parser.parse_args()

    build_master_table(
        spectra_dir=args.spectra_dir,
        amr_metadata_file=args.amr_metadata,
        output_dir=args.output_dir,
        build_profiles=args.build_profiles,
        collapse_to_six=args.collapse_to_six
    )


if __name__ == "__main__":
    main()
