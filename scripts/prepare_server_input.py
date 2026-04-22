import pandas as pd
import json
import argparse
import os
from pathlib import Path
import re
import tqdm
import copy


def get_next_json_path(output_dir: str) -> str:
    output_dir = Path(output_dir)
    base_name = "server_input"
    pattern = re.compile(rf"^{base_name}(?:_(\d+))?\.json$")

    used = set()

    for f in output_dir.iterdir():
        if f.is_file():
            m = pattern.match(f.name)
            if m:
                idx = 0 if m.group(1) is None else int(m.group(1))
                used.add(idx)

    n = 0
    while n in used:
        n += 1

    filename = f"{base_name}.json" if n == 0 else f"{base_name}_{n}.json"
    return str(output_dir / filename)


JSON_MAX_LENGTH = 100

AF3_BASE_DICT = {
    "name": "AGTRNIDYL",
    "dialect": "alphafoldserver",
    "version": 1,
    "dialect": "alphafold3",
    "version": 1,
}


def make_chain(chain_id, sequence):
    return {
        "proteinChain": {
            "id": chain_id,
            "sequence": sequence,
            "count": 1,
        }
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

    if not os.path.exists(args.output_dir):
        print(f"{args.output_dir} does not exist... Creating new directory")
        os.makedirs(args.output_dir)

    if seq_df.shape[0] > JSON_MAX_LENGTH:
        raise Exception(
            f"Current input CSV is larger than {JSON_MAX_LENGTH} rows. AF3 Server JSON can only be run in batches of {JSON_MAX_LENGTH} at a time. Please update the CSV and run again."
        )

    full_json = []

    for i, (idx, row) in tqdm.tqdm(enumerate(seq_df.iterrows())):
        dd = AF3_BASE_DICT.copy()
        dd["name"] = f"index_{idx}"
        dd["modelSeeds"] = [args.af3_seed]
        dd["sequences"] = [
            make_chain("A", row[args.alpha_col]),
            make_chain("B", row[args.beta_col]),
            make_chain("P", row[args.peptide_col]),
            make_chain("M", row[args.mhc_col]),
        ]

        full_json.append(dd)

    json_path = get_next_json_path(args.output_dir)
    with open(json_path, "w") as json_file:
        json.dump(full_json, json_file, indent=4)

    print(
        f"Server input JSON files successfully written to {json_path}. You can upload that file on https://alphafoldserver.com/"
    )


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
        "--af3-seed",
        type=int,
        required=False,
        default=1,
        help="AF3 Seed for Folding. Note that with the server, you can only input one seed at a time",
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
