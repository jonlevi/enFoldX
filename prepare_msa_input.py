import pandas as pd
import json
import argparse
import os


AF3_BASE_DICT = {
    "name": "",
    "sequences": [{"protein": {"id": "X", "sequence": ""}}],
    "modelSeeds": [1],
    "dialect": "alphafold3",
    "version": 1,
}


def main(args):

    try:
        seq_df = pd.read_csv(args.sequences_file)
    except FileNotFoundError as e:
        print(f"Error: The file {args.sequences_file} was not found.")
    except Exception as e:
        print(f"Error reading CSV file: {e}")

    for name, col in [
        ("alpha", args.alpha_col),
        ("beta", args.beta_col),
        ("mhc", args.mhc_col),
        ("peptide", args.peptide_col),
    ]:
        assert (
            col in seq_df.columns
        ), f"Missing {name} sequence column which was to {col}. Please check it is set properly and try again."
        assert (
            seq_df[col].notna().all() and (seq_df[col] != "").all()
        ), f"{col} contains missing values"

    # Extract Unique Sequences
    unique_alphas = seq_df[args.alpha_col].unique()
    unique_betas = seq_df[args.beta_col].unique()
    unique_mhcs = seq_df[args.mhc_col].unique()

    if not os.path.exists(args.output_dir):
        print(f"{args.output_dir} does not exist... Creating new directory")
        os.makedirs(args.output_dir)

    MSA_input_path = args.output_dir

    if not os.path.exists(MSA_input_path):
        os.mkdir(MSA_input_path)

    chain_mapper_path = os.path.join(args.output_dir, "chain_ids_to_sequences.txt")

    with open(chain_mapper_path, "w") as f:
        f.write("ID\tSEQ")
        f.write("\n")
        for i, a in enumerate(unique_alphas):
            f.write(f"alpha_{i}\t{a}\n")
        for i, b in enumerate(unique_betas):
            f.write(f"beta_{i}\t{b}\n")
        for i, b in enumerate(unique_mhcs):
            f.write(f"mhc_{i}\t{b}\n")

    id_mapper = pd.read_table(chain_mapper_path)

    for _, row in id_mapper.iterrows():
        dd = AF3_BASE_DICT.copy()
        dd["name"] = row["ID"]
        dd["sequences"][0]["protein"]["sequence"] = row["SEQ"]
        json_path = os.path.join(MSA_input_path, f"{row['ID']}.json")
        with open(json_path, "w") as json_file:
            json.dump(dd, json_file, indent=4)

    print(f"MSA input JSON files successfully written to {MSA_input_path}")
    print(f"Chain ID map written to {chain_mapper_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-s",
        "--sequences-file",
        type=str,
        required=True,
        help="Path to input file with TCR-pMHC sequences",
    )

    parser.add_argument(
        "-a",
        "--alpha-col",
        type=str,
        required=False,
        default="TRA_aa",
        help="Name of column containing TCRa chain sequence in <--sequences-file>",
    )

    parser.add_argument(
        "-b",
        "--beta-col",
        type=str,
        required=False,
        default="TRB_aa",
        help="Name of column containing TCRb chain sequence in <--sequences-file>",
    )

    parser.add_argument(
        "-m",
        "--mhc-col",
        type=str,
        required=False,
        default="M_aa",
        help="Name of column containing MHC chain sequence in <--sequences-file>",
    )

    parser.add_argument(
        "-p",
        "--peptide-col",
        type=str,
        required=False,
        default="peptide",
        help="Name of column containing peptide sequence in <--sequences-file>",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Directory to place output JSON files",
    )

    args = parser.parse_args()

    main(args)
