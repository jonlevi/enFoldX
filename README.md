# enFoldX: Leveraging ensembles of in silico structures for complex binding prediction


enFoldX (**En**semble of **Fold**ed Comple**X**es) enables the binding prediction of large datasets of input sequences by extracting distributions of high dimensional confidence features from ensembles of predicted structures.

This repo provides code for a pipeline to run enFoldX. This implementation of enFoldX uses AlphaFold3 and runs predictions for complex binding for TCR-peptide-MHC complexes (for class I MHC).


## enFoldX

![Project Banner](https://github.com/jonlevi/af3_tcr_pipeline/blob/main/enfoldx_diagram.png) 


## Terms of Use
By using enFoldX, you are agreeing to the terms set in the enFoldX [Terms of Use](TERMS_OF_USE.md), and the terms of use set by AlphaFold3 (see below).

## Issues / Contact
If you have any issues running enFoldX or any questions about the code or the method, please file an issue on the repository or [reach out directly via email](mailto:levinej4@mskcc.org).

## Installation
Clone the repository:

```bash
git clone https://github.com/jonlevi/enFoldX.git
cd enFoldX
```
Alternatively, you can download the zip file directly using the "download" option on the repository options. Either clone or download should complete in <5 minutes.

The python requirements to run these scripts are fairly minimal and you can run
```bash
conda env create -f environment.yml
conda activate enfoldx_env
```
to create a `enfoldx_env` virtual environment, and then activate it. The environment includes the basic toolkit (numpy, pandas, etc.) as well as a few tools that we interact with (biopython etc.).

## Getting Started
You can run enFoldX with AF3 locally installed, or by using the [AF3 server](https://alphafoldserver.com/). We include more detailed instructions for both of these options [below](#run-alphafold3-predictions).

**We highly recommend using local installation of AF3 for enFoldX, as it allows for the incorporation of multiple seeds at the same time**. Currently, AlphaFold Server runs just one seed for each job. If you want to sample multiple seeds, you need to run multiple identical job submissions with different input seeds, and then collect the results yourself, but this is rather hard to scale up. You also cannot separate the MSA steps from the diffusion steps on the AF3 server, which also limits your ability to re-use MSAs, which is critical for efficient prediction of the same TCR against multple potential peptides.

We provide step-by-step tutorials for how to use enFoldX for either the installed AF3 as for the server AF3:

To run enFoldX, you need to essentially run 3 main steps:
1) [Prepare sequence inputs](#prepare-sequence-inputs)
2) [Run Structure Predictions and Extract Features](#run-alphafold3-predictions)
3) [Predict binding](#predict-binding)

The first step is detailed below. The second step for running AF3 can be found in the companion tutorial files, depending on if you are using the local installation or the server.

## Step 1: Prepare sequence inputs

### TCRs
In order to run this pipeline, you will need a file that contains one row per TCR-pMHC complex that you wish to predict. If you already have full length TCR sequences, you can skip this step. For instructions on how to go from `V`, `J`, `CDR3` calls to full length TCR sequences, have a look at our [tutorial](https://github.com/jonlevi/enFoldX/blob/main/docs/format_tcr_sequences.md). (Note that full length includes leader/constant/framework sequences, not just the variable regions).

### MHC
You will also need the full length sequences for any MHC/HLA chains you want to model. For MHC sequence information, you can look up the allele in [Uniprot](https://www.uniprot.org/uniprotkb) or IPD-IMGT/HLA at https://www.ebi.ac.uk/ipd/imgt/hla/alleles/. For convenience, we include the [sequences of many common HLA alleles](https://github.com/jonlevi/enFoldX/blob/main/MHC_sequences/) together with our repo.

### Output
In order to continue with the next steps, you need a CSV that has columns containig the full-length amino acids sequences for the TCRa, TCRb, MHC, and peptide sequence per row. It also must contain the sequences of the CDR3a and the CDR3b, which must be substrings of the TCRa and TCRb respectively. Your CSV should look something like this:

| TRA_aa   | TRA_CDR3         | TRB_aa   | TRB_CDR3         | M_aa     | peptide   |
|----------|------------------|----------|------------------|----------|-----------|
| MDSSPG... | CALGDPPNTGKLTF  | MGSRL... | CASTSGVGQDTQYF   | MVPCTL... | TVYGFCLL |
| MLILS...  | CAMRSSGGSNAKLTF | MGAMA... | CASSGGANTGQLYF   | MVPCTL... | ASNENMETM |

etc.

(Note: the column names above are the default names for all of the scripts in the enFoldX code, but you can override with custom column names by passing in the apporpirate flags to each script with
`[-a <ALPHA_COL_NAME>] [-b <BETA_COL_NAME>] [-m <MHC_COL_NAME>] [-p <PEPTIDE_COL_NAME>] [-cdr3a <CDR3A_COL_NAME>] [-cdr3b <CDR3B_COL_NAME>]`)

## Step 2: Run AlphaFold3 Predictions

There are currently 2 ways to run AlphaFold3 predictions:
- [Tutorial for local installation of AF3](docs/TUTORIAL.md)
- [Tutorial for AF3 server](docs/TUTORIAL_SERVER_BASED.md)

After following the steps in either of those tutorials, you should have a CSV called `ensemble_features.csv` that can be used to predict binding below with our pre-trained models. The output also contains `all_structures_features.csv`, `avg_features.csv`, and `best_features.csv` which contains the full ensemble, just the average values over the ensemble, or the highest-ranked structure from the ensemble, respectively. For the pre-trained enFoldXs models, you will use `ensemble_features.csv`. Once you have that CSV, you can continue to the next step here.

## Step 3: Predict Binding
TODO: @olga to add pre-trained models here, as well as any other advice about scaling etc.
We include an example notebook at `analysis_notebooks/enFoldX_ML_notebook.ipynb` for how we used the features to predict specificity. The enFoldX_ML Jupyter notebook provides an example of training a Ridge regression model on the cross-reactivity dataset utilizing features prepared by the enFoldX AF3 pipeline. It includes preparation of the feature ensemble across all predicted structures for the same TCR/epitope pair and 5-fold cross-validation. The notebook needs fairly standard Python packages to be installed: numpy, pandas, sklearn, matplotlib, and seaborn.In terms of the pipeline, this notebook starts with the feature set that is created at the very end of the last step (Step 5). This notebook uses data from the cross-reactivity data that can be found as part of [a previous paper](https://www.nature.com/articles/s41586-022-04735-9).

Of course, all of these performance evaluations rely on the availability of labeled data. For the use of this tool in studying TCR-pMHCs for which the true label of binder/non-binder is not known, we recommend training on a relevant, related labeled data set, and then using the resultant model to predict specificity for the novel data. For example, you can run enFoldX to create features for all of the human data in VDJdb, then train a Logistic model on that labelled data, and then use the model to predict specificity for a new human dataset.


## Other Useful Things 
###  Access to labeled mutational scan data used in our paper

If you are looking for the labeled mutational scan datasets used for validation our paper, please visit the ```manuscript/data/mutational_scan_data``` folder. The other datasets used in the paper can be accessed directly from their original publications, as described in the methods section.

### RMSD Calculations
You can also run the rmsd calculation that we use for the specific scenario of comparing two different TCR:pMHC complexes to each other. We use this in the first section of our paper when comparing structural diversity across an ensemble of predictions. Run the script passing in the paths for structure 1 and 2, like this example:
```bash
python ./scripts/rmsd_calculation.py -s1 examples/af3_fold_outputs/index_0/seed-1_sample-0/model.cif -s2 examples/af3_fold_outputs/index_0/seed-1_sample-1/model.cif
```

