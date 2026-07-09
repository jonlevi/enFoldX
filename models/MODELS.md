# Models

This document summarizes the machine learning models included in this repository.

## Model Files

Each model is distributed as a serialized Python pickle (`.pkl`) file containing all objects required for inference. Every pickle stores a dictionary with the following entries:

- **`scaler`** – A fitted feature scaler (e.g., `sklearn.preprocessing.StandardScaler`) used to transform the input features before prediction.
- **`model`** – A fitted machine learning model (e.g., `sklearn.linear_model.LogisticRegression`) used to generate predictions.
- **`features`** – A list of feature names expected by the model. Input data should contain these features before scaling and prediction.

---

## Model 1: enFoldX_human_vFebSept_DecoyPerm

**Description:**  
This is a L2-regularized logistic regression enFoldX model trained on VDJdb Human data (09/2025 and 02/2025 releases, records with VDJdb score 1+), with Decoy and Permuted negatives, as described in the manuscript. This is our recommended model make cognate/non-cognate predictions on all human TCR:pMHC complexes. This model was trained with features from an 10-seed ensemble of structures.

---

## Model 2: enFoldX_human_vFebSept_Decoy

**Description:**  
This is a L2-regularized logistic regression enFoldX model trained on VDJdb Human data (09/2025 and 02/2025 releases, records with VDJdb score 1+), with Decoy negatives, as described in the manuscript. This is our recommended model for inference on mutational scan data. This model was trained with features from an 10-seed ensemble of structures.


---

## Model 3: enFoldX_mouse_vFeb_Decoy

**Description:**  
This is a L2-regularized logistic regression enFoldX model trained on VDJdb Mouse data (02/2025 releases, records with VDJdb score 1+), with Decoy negatives, as described in the manuscript. Please use it for all mouse TCR:pMHC data. This model was trained with features from an 10-seed ensemble of structures.


---

## Model 4: enFoldX_human_vFebSept_DecoyPerm_1seed

**Description:**  
This is a L2-regularized logistic regression enFoldX model trained on VDJdb Human data (09/2025 and 02/2025 releases, records with VDJdb score 1+), with Decoy and Permuted negatives, with features extracted from AF3 samples for just 1 seed. It does slightly worse than the 10-seed models (above) when used on test sets with 10-seed-based features, but if you only ran enFoldX pipeline with 1 seed (for example, using AlphaFold server) and so only have features for 1 seed (ensembled across 5 samples), please use this model.
