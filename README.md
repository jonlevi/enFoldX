# enFoldX: Leveraging ensembles of predicted structures for complex binding prediction


enFoldX (**En**semble of **Fold**ed Comple**X**es) enables the binding prediction of large datasets of input sequences by extracting distributions of high dimensional confidence features from ensembles of predicted structures.

This repo provides code for a pipeline to run enFoldX. This implementation of enFoldX uses AlphaFold3 and runs predictions for complex binding for TCR-peptide-MHC complexes (for class I MHC). 

![Project Banner](https://github.com/jonlevi/af3_tcr_pipeline/blob/main/enfoldx_diagram.png) 


## Terms of Use
By using enFoldX, you are agreeing to the terms set in the enFoldX [Terms of Use](TERMS_OF_USE.md), and the terms of use set by AlphaFold3 (see below).

## Installing this repository
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
You can run EnFoldX with AF3 locally installed, or AF3 server. 

**We highly recommend using local installation of AF3 for EnFoldX, as it allows for the incorporation of multiple seeds at a time**. Currently, AlphaFold Server runs just one seed for each job. If you want to sample multiple seeds, you could theoretically run multiple identical job submissions with different input seeds, and then collect the results yourself, but this is rather hard to scale up.

Nevertheless, we provide step-by-step tutorials for how to use EnFoldX for either the installed AF3 as for the server AF3:

To run EnFoldX, you need to essentially run 3 main steps:
1) [Prepare sequence inputs](#prepare-sequence-inputs)
2) [Run Structure Predictions and Extract Features](#run-alphafold3-predictions)
3) [Predict binding](#predict-binding)

The first step is detailed below. The rest of them can be found in the companion tutorial files, depending on if you are using the local installation or the server.

## Prepare sequence inputs

### TCRs
In order to run this pipeline, you will need a file that contains one row per TCR-pMHC complex that you wish to predict. If you already have full length TCR sequences, you can skip this step. For instructions on how to go from `V`, `J`, `CDR3` calls to full length TCR sequences, have a look at our [tutorial](https://github.com/jonlevi/enFoldX/blob/main/docs/format_tcr_sequences.md). (Note that full length includes leader/constant/framework sequences, not just the variable regions).

### MHC
You will also need the full length sequences for any MHC/HLA chains you want to model. For MHC sequence information, you can look up the allele in [Uniprot](https://www.uniprot.org/uniprotkb) or IPD-IMGT/HLA at https://www.ebi.ac.uk/ipd/imgt/hla/alleles/. For convenience, we include the [sequences of many common HLA alleles](https://github.com/jonlevi/enFoldX/blob/main/MHC_sequences/) together with our repo.

### Output
In order to continue with the next steps, you need a CSV that has a alpha, beta, MHC, and peptide sequence per row. It should look something like this:

| TRA_aa | TRB_aa | M_aa | peptide |
|--------|--------|------|---------|
| MDSSPG... | MGSRL... | MVPCTL... | TVYGFCLL |
| MLILS...  | MGAMA... | MVPCTL... | ASNENMETM |
etc.


## Run AlphaFold3 Predictions

There are currently 2 ways to run AlphaFold3 predictions:
a) [local installation of AF3](docs/TUTORIAL.md)
b) [AF3 server](docs/TUTORIAL_SERVER_BASED.md)

**We highly recommend using local installation of AF3 for EnFoldX, as it allows for the incorporation of multiple seeds at a time**. Currently, AlphaFold Server runs just one seed for each job. If you want to sample multiple seeds, you could theoretically run multiple identical job submissions with different input seeds, and then collect the results yourself, but this is rather hard to scale up.

Nevertheless, in those tutorials, we break down how to use EnFoldX for either choice. After following the steps in those tutorials, you should have a CSV called `ensemble_features.csv` that can be used to predict binding below. 

## Predict Binding
TODO: @olga to add pre-trained models here, as well as any other advice about scaling etc.
We include an example notebook at `analysis_notebooks/enFoldX_ML_notebook.ipynb` for how we used the features to predict specificity. The enFoldX_ML Jupyter notebook provides an example of training a Ridge regression model on the cross-reactivity dataset utilizing features prepared by the enFoldX AF3 pipeline. It includes preparation of the feature ensemble across all predicted structures for the same TCR/epitope pair and 5-fold cross-validation. The notebook needs fairly standard Python packages to be installed: numpy, pandas, sklearn, matplotlib, and seaborn.In terms of the pipeline, this notebook starts with the feature set that is created at the very end of the last step (Step 5). This notebook uses data from the cross-reactivity data that can be found as part of [a previous paper](https://www.nature.com/articles/s41586-022-04735-9).

Of course, all of these performance evaluations rely on the availability of labeled data. For the use of this tool in studying TCR-pMHCs for which the true label of binder/non-binder is not known, we recommend training on a relevant, related labeled data set, and then using the resultant model to predict specificity for the novel data. For example, you can run enFoldX to create features for all of the human data in VDJdb, then train a Logistic model on that labelled data, and then use the model to predict specificity for a new human dataset.


## Other Misc. Code
### RMSD Calculations
You can also run the rmsd calculation that we use for the specific scenario of comparing two different TCR:pMHC complexes to each other. We use this in the first section of our paper when comparing structural diversity across an ensemble of predictions. Run the script passing in the paths for structure 1 and 2, like this example:
```bash
python rmsd_calculation.py -s1 examples/af3_fold_outputs/index_0/seed-1_sample-0/model.cif -s2 examples/af3_fold_outputs/index_0/seed-1_sample-1/model.cif
```

