import json
import pandas as pd
import tqdm
import os
import argparse


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
        ("complex_id", args.id_col),
    ]:
        assert (
            col in seq_df.columns
        ), f"Missing {name} column which was to {col}. Please check it is set properly and try again."
        assert (
            seq_df[col].notna().all() and (seq_df[col] != "").all()
        ), f"{col} contains missing values"

    assert seq_df[
        args.id_col
    ].is_unique, f"ID Column ({args.id_col}) entries are not unique!"

    try:
        chain_map_df = pd.read_table(args.chain_id_map)
    except FileNotFoundError as e:
        print(f"Error: The file {args.chain_id_map} was not found.")
    except Exception as e:
        print(f"Error reading CSV file: {e}")

    def map_seq_to_id(seq):
        my_row = chain_map_df[chain_map_df["SEQ"] == seq]
        assert (
            my_row.shape[0] == 1
        ), f"duplicate or missing entries found in chain_map for seq {seq}"
        return my_row["ID"].values[0]

    seq_df["alpha_id"] = seq_df[args.alpha_col].apply(map_seq_to_id)
    seq_df["beta_id"] = seq_df[args.beta_col].apply(map_seq_to_id)
    seq_df["mhc_id"] = seq_df[args.mhc_col].apply(map_seq_to_id)

    assert os.path.exists(args.MSA_output_dir), f"{args.MSA_output_dir} not found"

    if not os.path.exists(args.output_dir):
        print(f"{args.output_dir} does not exist... Creating new directory")
        os.makedirs(args.output_dir)

    fold_input_path = args.output_dir

    def get_msa(id):
        folder_path = os.path.join(args.MSA_output_dir, id)
        json_path = os.path.join(folder_path, f"{id}_data.json")
        with open(json_path, "r") as jf:
            msa_data = json.load(jf)["sequences"][0]["protein"]
            return (
                msa_data["unpairedMsa"],
                msa_data["templates"],
                msa_data["sequence"],
            )

    chain_base_dict = {
        "id": "",
        "sequence": "",
        "modifications": [],
        "pairedMsa": "",
        "unpairedMsa": "",
        "templates": None,
    }
    try:
        model_seeds = [int(s) for s in args.af3_seeds.split(",")]
    except Exception as e:
        print(
            "Could not process the input for af3-seeds correctly. Please set it to a comma-separated-string of integers like 1,2,3"
        )
        raise e

    base_dict = {
        "name": "",
        "sequences": None,
        "modelSeeds": model_seeds,
        "dialect": "alphafold3",
        "version": 2,
    }

    # Write the files
    for i, (idx, row) in tqdm.tqdm(enumerate(seq_df.iterrows())):

        a_id = row["alpha_id"]
        b_id = row["beta_id"]

        mhc_id = row["mhc_id"]
        mhc_seq = row[args.mhc_col]

        peptide_seq = row[args.peptide_col]

        this_af_input = base_dict.copy()
        this_af_input["sequences"] = []
        this_af_input["name"] = row[args.id_col]

        j_path = os.path.join(fold_input_path, f"{this_af_input['name']}.json")

        # alpha chain
        alpha_chain = chain_base_dict.copy()
        alpha_chain["id"] = "A"
        alpha_chain["sequence"] = row[args.alpha_col]
        a_msa, a_templates, a_seq = get_msa(a_id)

        assert a_seq == row[args.alpha_col]

        alpha_chain["unpairedMsa"] = a_msa
        alpha_chain["templates"] = a_templates

        this_af_input["sequences"].append({"protein": alpha_chain})

        # beta chain
        beta_chain = chain_base_dict.copy()
        beta_chain["id"] = "B"
        beta_chain["sequence"] = row[args.beta_col]
        b_msa, b_templates, b_seq = get_msa(b_id)

        assert b_seq == row[args.beta_col]

        beta_chain["unpairedMsa"] = b_msa
        beta_chain["templates"] = b_templates

        this_af_input["sequences"].append({"protein": beta_chain})

        # peptide (no MSA or templates)
        # use random peptide from netMHC calls
        peptide_chain = chain_base_dict.copy()
        peptide_chain["id"] = "P"
        peptide_chain["sequence"] = peptide_seq
        peptide_chain["templates"] = []
        peptide_chain["unpairedMsa"] = ""

        this_af_input["sequences"].append({"protein": peptide_chain})

        # MHC chain
        mhc_chain = chain_base_dict.copy()
        mhc_chain["id"] = "M"
        mhc_chain["sequence"] = mhc_seq
        m_msa, m_templates, m_seq = get_msa(mhc_id)

        assert mhc_seq == m_seq

        mhc_chain["unpairedMsa"] = m_msa
        mhc_chain["templates"] = m_templates

        this_af_input["sequences"].append({"protein": mhc_chain})

        with open(j_path, "w") as newf:
            json.dump(this_af_input, newf)


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
        "-i",
        "--id-col",
        type=str,
        required=False,
        default="complex_id",
        help="Name of column containing unique complex ID in <--sequences-file>",
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
        "--chain-id-map",
        type=str,
        required=True,
        help="File that maps chain IDs to sequences from MSA run. Should be an output from the <prepare_msa_input> step",
    )

    parser.add_argument(
        "--MSA-output-dir",
        type=str,
        required=True,
        help="Directory with MSA outputs. \
            This directory should contain one subdirectory per ID in the chain-id-map,\
            each containing a single file named <ID>_data.json",
    )

    parser.add_argument(
        "--af3-seeds",
        type=str,
        required=False,
        help="comma-separated string of seeds to input to AF3; default is 10 seeds",
        default="1,13,17,21,42,133,177,213,315,1001",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Directory to place output JSON files to pass to AlphaFold3",
    )

    args = parser.parse_args()

    main(args)
