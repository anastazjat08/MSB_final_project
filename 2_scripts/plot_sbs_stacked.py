import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap

"""
Stacked barplot to summarise SBS for metadata category (e.g., Collection_name, ST)
"""

def plot_stacked_bars(df: pd.DataFrame, outfile: str):

    # set Label_ID as index
    df = df.set_index("Label_ID")

    # compute proportions
    df_prop = df.div(df.sum(axis=1), axis=0)

    # plot
    plt.figure(figsize=(12, 6))
    df_prop.plot(kind="bar", stacked=True, edgecolor="black")

    plt.ylabel("Proportion of SBS")
    plt.xlabel("Group")
    plt.title("SBS composition across groups")
    plt.xticks(rotation=45, ha="right")
    plt.legend(title="Mutation type", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    # save to file
    plt.savefig(outfile, dpi=300)
    plt.close()


def main():
    parser = ap.ArgumentParser(description="Plot stacked SBS barplot")
    parser.add_argument("--input", required=True, help="Input CSV with SBS counts")
    parser.add_argument("--output", required=True, help="Output PNG file")

    args = parser.parse_args()

    # read input
    df = pd.read_csv(args.input)

    # plot
    plot_stacked_bars(df, args.output)


if __name__ == "__main__":
    main()
