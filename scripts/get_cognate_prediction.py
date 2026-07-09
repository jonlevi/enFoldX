import pandas as pd
import os
import argparse
import numpy as np
from enum import Enum
import pickle


class enFoldXModel(str, Enum):
    human = "models/enFoldX_human_vFebSept_DecoyPerm.pkl"
    human_decoy = "models/enFoldX_human_vFebSept_Decoy.pkl"
    mouse = "models/enFoldX_mouse_vFeb_Decoy.pkl"
    human_1seed = "models/enFoldX_human_vFebSept_DecoyPerm_1seed.pkl"


def load_model_bundle(model_type: str):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_dir = os.path.dirname(script_dir)
    model_path = os.path.join(repo_dir, enFoldXModel[model_type].value)

    with open(model_path, "rb") as f:
        bundle = pickle.load(f)

    return bundle["scaler"], bundle["model"], bundle["features"]


def predict_from_model(data_predict, featurecols, scaler, model):
    X_test = data_predict[featurecols].copy()
    X_test_scaled = scaler.transform(X_test)
    # y_pred = model.predict(X_test_scaled)
    y_proba = model.predict_proba(X_test_scaled)[:, 1]
    return y_proba


def main(args):

    try:
        features_df = pd.read_csv(args.features_file)
    except FileNotFoundError as e:
        print(f"Error: The file {args.features_file} was not found.")
    except Exception as e:
        print(f"Error reading CSV file: {e}")

    model_choice = args.model
    assert model_choice in [x.name for x in enFoldXModel]

    scaler, model, features = load_model_bundle(model_choice)

    missing_features = set(features) - set(features_df.columns)
    if missing_features:
        raise ValueError(f"Missing required features: {sorted(missing_features)}")

    y_proba = predict_from_model(
        data_predict=features_df, featurecols=features, scaler=scaler, model=model,
    )

    non_feature_cols = [c for c in features_df.columns if c not in features]
    if non_feature_cols:
        model_predictions = features_df[non_feature_cols].copy()
    else:
        model_predictions = pd.DataFrame({"index": features_df.index})

    model_predictions["enFoldX_score"] = y_proba

    os.makedirs(args.output_dir, exist_ok=True)
    output_path = os.path.join(args.output_dir, args.output_filename)
    model_predictions.to_csv(output_path, index=False)

    print(f"Model predictions written to {output_path}")


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
        choices=[x.name for x in enFoldXModel],
        default="human",
        help="Pre-trained enFoldX model to use.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Directory to save the output CSV file with enFoldX scores.",
    )

    parser.add_argument(
        "-of",
        "--output-filename",
        type=str,
        default="enFoldX_predictions.csv",
        help="Filename for output CSV file with enFoldX scores.",
    )

    args = parser.parse_args()
    main(args)
