import json
import pandas as pd
import tqdm
import os
import argparse
import numpy as np
from enum import Enum


class enFoldXModel(Enum):
    VDJdbFebruaryDecoyRidge = 1
    VDJdbFebruaryPermutedRidge = 2
    VDJdbAllDecoyRidge = 3
    VDJdbAllPermutedRidge = 4
    # TODO: finish list of all models included in the repo
    # TODO: do we separate training data from model type or keep them all as one list?


def load_model(model_type):
    return None


def load_scaler(model_type):
    return None


def main(args):

    try:
        features_df = pd.read_csv(args.features_file)
    except FileNotFoundError as e:
        print(f"Error: The file {args.features_file} was not found.")
    except Exception as e:
        print(f"Error reading CSV file: {e}")

    # TODO: assert features file is formated correctly
    assert features_df.columns == []

    # TODO: get user input to determine what model and scaler to use
    model_choice = args.model
    assert model_choice in [x for x in enFoldXModel]

    # TODO: run scaling - currently in pseudocode
    scaler = load_scaler(model_choice)
    scaled_data = scaler.scale(features_df)

    # TODO: run model prediction - currently in pseudocode
    model = load_model(model_choice)
    model_predictions = model.predict(scaled_data)

    if not os.path.exists(args.output_dir):
        print(f"{args.output_dir} does not exist... Creating new directory")
        os.makedirs(args.output_dir)

    print(f"Writing results to {os.path.join(args.output_dir,args.output_filename)}...")

    # TODO: write results
    model_predictions.to_csv(os.path.join(args.output_dir, args.output_filename))

    print(f"Output features written to {args.output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-f",
        "--features-file",
        type=str,
        required=True,
        help="Path to input file with enFoldX features",
    )

    parser.add_argument(
        "-m",
        "--model",
        type=str,
        required=False,
        default="XXX",  # TODO: set default model choice
        help="Choose which pre-trained enFoldX model to use",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Directory to place output CSV files with features",
    )

    parser.add_argument(
        "-of",
        "--output-filename",
        type=str,
        required=False,
        default="enFoldX_predictions.csv",
        help="Filename for output CSV files with features",
    )

    args = parser.parse_args()

    main(args)
