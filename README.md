# enFoldX: Leveraging ensembles of predicted structures for complex binding prediction


enFoldX (**En**semble of **Fold**ed Comple**x**es) enables the evaluation of large datasets of predicted structures simultaneously as it extracts distributions of high dimensional confidence features, including customized biophysical features. 

This repo provides code for a pipeline to run enFoldX. The current implementation of enFoldX uses AlphaFold3 and runs predictions for complex binding betweens TCRs and peptide-MHC complexes. 

![Project Banner](https://github.com/jonlevi/af3_tcr_pipeline/blob/main/enfoldx_diagram.png) 

## 📚 Table of Contents
- [Terms of use](#terms-of-use)
- [Installation](#installation)
- [Usage](#usage)

## Terms of Use
By using enFoldX, you are agreeing to the terms set in the enFoldX [Terms of Use](TERMS_OF_USE.md), and the terms of use set by AlphaFold3 (see below).

## Getting Started
To run EnFoldX, you need to essentially run 3 main steps:
1) Prepare sequence inputs
2) Run AlphaFold3 (see below)
3) Extract Features and Predict Binding from AlphaFold3 Results

Each of these steps is detailed below:

##  Prepare sequence inputs

### TCRs
In order to run this pipeline, you will need a CSV that contains one row for every TCR-pMHC complex you wish to predict. If you already have full length TCR sequences, you can skip this step. Usually, you will only have `V`,`J`, and `CDR3` calls from sequencing. We prefer to use the full length TCR sequences that come from the output of the useful tool [stitchr](https://jamieheather.github.io/stitchr/index.html), which reconstructs full-length TCR sequences from minimal input data. We will be using the command-line utility `thimble`, a wrapper included with `stitchr`, for batch processing of TCR alpha and beta chains. Before you begin, ensure you have installed `stitchr` and the necessary data (for the genome of interest) by following the instructions [here] (https://jamieheather.github.io/stitchr/installation.html#installation).

To start, you should have a CSV that has paired chain TCR information for each TCR you are interested in. These should include the 6 columns: 
`[TRAV,	TRAJ,	TRA_CDR3,	TRBV,	TRBJ,	TRB_CDR3]`. Thimble expects the input file to have a specific format (see [here](https://github.com/JamieHeather/stitchr/blob/main/templates/input_template_TRA-TRB.tsv) for a template).

To run thimble, you can pass in the file with the TCR calls and the appropriate species. For example:
`thimble -i my_human_tcrs.tsv -o stitchr_out_human_tcrs.tsv -r ab -s HUMAN`. The output will include columns containing the full AA sequences of the TCRa and TCRb chains using 1-letter amino acid codes

### MHC
## Get MHC/HLA Sequence Information
You will also need the full length sequences for any MHC chains you want to model. For MHC sequence information, you can look up the allele in [Uniprot](https://www.uniprot.org/uniprotkb) or IPD-IMGT/HLA at https://www.ebi.ac.uk/ipd/imgt/hla/alleles/. For convenience, we include the sequences of a few common HLA alleles in this repo at <>



## Installation
### Installing AlphaFold3
If you wish to use enFoldX with AlphaFold3 predictions, you must install AlphaFold3 following the instructions at https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md 

After the installation you should have:
- path to the AF3 parameters
- path to downloaded databases
- path to the docker/singularity container for running alphafold

Please keep in mind that you must additional comply with the [TERMS OF USE](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md) set by Deepmind in order to use AlphaFold3. Please also keep in mind that in order to run AlphaFold3 efficiently, you may need access to specialized hardware (see [https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md](https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md)). 

### Installing this repository
Clone the repository:

```bash
git clone https://github.com/jonlevi/enFoldX.git
cd enFoldX
```
Alternatively, you can download the zip file directly using the "download" option on the repository options. Either clone or download should complete in <5 minutes.

### python reqs
The python requirements to run these scripts are minimal so instead of forcing you to install a specific list of reqiurements, just make sure you run these scripts in any environment that has:
- numpy
- pandas
- tqdm

### Configure AlphaFold Paths:
Open `af3_config.sh` and replace the stubbed out paths with the correct paths that contain the AlphaFold3 installation from above:
```
ALPHAFOLD_DIR="path/to/alphafold3"
DATABASE_DIR="path/to/public_databases"
WEIGHTS_DIR="path/to/weights"
CONTAINER_PATH="path/to/container"
```

## Usage
The example code here is for running complex binding predictions for a TCR:pmHC complex, although enFoldX can be adapted to any protein complex with some trivial changes.

There are a few sequential steps to take in going from a TCR-pMHC sequence --> set of features:
1) Format JSONs for running AlphaFold3 MSA for each unique sequence
2) Run MSA on Input Sequences
3) Collect MSA results and Format JSONs per TCR-pMHC for running AlphaFold3 Folding
4) Run AlphaFold3 Folding on Input Sequences
5) Extract Features from Predicted Structures and Confidence Metadata

The key idea of this pipeline is to set things up so that we can parallelize as much as possible. The MSA is run on each unique sequence independently and thus can be parallelized across each one. The Folding is run on each TCR-pMHC complex independently, and thus can be parallelized as well. This README goes through a full example of going from TCR:pMHC sequences --> enFoldX features. All of the data for the examples can be found in `examples/`

### Step 0: 

### Step 1: Format JSONs for running AlphaFold3 MSA for each unique sequence
```bash
usage: prepare_msa_input.py [-h] -s SEQUENCES_FILE [-a ALPHA_COL] [-b BETA_COL] [-m MHC_COL] [-p PEPTIDE_COL] -o OUTPUT_DIR
options:
  -h, --help            show this help message and exit
  -s SEQUENCES_FILE, --sequences-file SEQUENCES_FILE
                        Path to input file with TCR-pMHC sequences
  -a ALPHA_COL, --alpha-col ALPHA_COL
                        Name of column containing TCRa chain sequence in <--sequences-file> (default: A_seq)
  -b BETA_COL, --beta-col BETA_COL
                        Name of column containing TCRb chain sequence in <--sequences-file> (default: B_seq)
  -m MHC_COL, --mhc-col MHC_COL
                        Name of column containing MHC chain sequence in <--sequences-file> (default: M_seq)
  -p PEPTIDE_COL, --peptide-col PEPTIDE_COL
                        Name of column containing peptide sequence in <--sequences-file> (default: P_seq)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to place output JSON files
```
#### Tutorial Example:
```python prepare_msa_input.py -s examples/example_input_tcr_pmhcs.csv -o examples/af3_msa_inputs```

This script will write out one JSON per unique TCRa, TCRb, and MHC sequence in the data, as well as a metadata file mapping the sequences to their chain IDs/JSON paths. 

You can see what a successful run looks like by looking at the output for the tutorial in `examples/af3_msa_inputs`


### Step 2: Run MSA on Input Sequences
You will run the AlphaFold3 container on `--no-run-inference` mode once per MSA JSON that was created. It is beneficial to run this in parallel. 
An example for what this might look like using [slurm](https://slurm.schedmd.com/documentation.html) can be found for the tutorial data in `slurm_af3_msa_example.sh`. If you use this submission script, make sure to adjust the partition name, the array size, and the input/output paths to match your preferences/requirements. 

You can see what a successful run looks like by looking at the output for the tutorial above in `examples/af3_msa_outputs` (although these JSON files are very large and can be difficult to view using regular IDEs or the github file browser). There will be one JSON per unique chain (except for peptides which don't get MSA), and a metadata chain mapping file to map the MSA output to the original sequences. (You will need this chain_id_map for the next step)

### Step 3: Collect MSA results and Format JSONs per TCR-pMHC for running AlphaFold3 Folding
```bash
usage: prepare_fold_input.py [-h] -s SEQUENCES_FILE [-a ALPHA_COL] [-b BETA_COL] [-m MHC_COL] [-p PEPTIDE_COL] --chain-id-map CHAIN_ID_MAP --MSA-output-dir MSA_OUTPUT_DIR [--af3-seeds AF3_SEEDS]
                             -o OUTPUT_DIR

options:
  -h, --help            show this help message and exit
  -s SEQUENCES_FILE, --sequences-file SEQUENCES_FILE
                        Path to input file with TCR-pMHC sequences
  -a ALPHA_COL, --alpha-col ALPHA_COL
                        Name of column containing TCRa chain sequence in <--sequences-file> (default: A_seq)
  -b BETA_COL, --beta-col BETA_COL
                        Name of column containing TCRb chain sequence in <--sequences-file> (default: B_seq)
  -m MHC_COL, --mhc-col MHC_COL
                        Name of column containing MHC chain sequence in <--sequences-file> (default: M_seq)
  -p PEPTIDE_COL, --peptide-col PEPTIDE_COL
                        Name of column containing peptide sequence in <--sequences-file> (default: P_seq)
  --chain-id-map CHAIN_ID_MAP
                        File that maps chain IDs to sequences from MSA run. Should be an output from the <prepare_msa_input> step
  --MSA-output-dir MSA_OUTPUT_DIR
                        Directory with MSA outputs. This directory should contain one subdirectory per ID in the chain-id-map, each containing a single file named <ID>_data.json
  --af3-seeds AF3_SEEDS
                        comma-separated string of seeds to input to AF3 (default: 1,2,5,10)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to place output JSON files to pass to AlphaFold3
```
You still need the original sequences file from Step 0 and you will also need:

1) the path that the chain id map was written to in Step 1 and 
2) the directory that AlphaFold3 wrote the MSA output to. (Note this is the output of step 2, not the output of step 1)

You can choose how many [random seeds](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md#random-seeds) to use in the AF3 inference. By default, we use 5 seeds per complex. You can add more or remove them by changing the `--af3-seeds` flag. Note that each additional seeds adds about an additional 1-2 minutes at minimum per complex to the runtime, but your predictive accuracy may improve with more seeds. See the our paper for more!

#### Example:
Following our previous tutorial data:
```
python prepare_fold_input.py -s examples/example_input_tcr_pmhcs.csv -o examples/af3_fold_inputs --chain-id-map examples/af3_msa_inputs/chain_ids_to_sequences.txt --MSA-output-dir examples/af3_msa_outputs
```
You can see what a successful run looks like by looking at the output for the tutorial in `examples/af3_fold_inputs` (although these JSON files are very large and can be difficult to view using regular IDEs or the github file browser).

### Step 4: Run AlphaFold3 Folding on Input Sequences
You will run the AlphaFold3 container on `--norun_data_pipeline` mode once per TCR-pMHC input JSON that was created. This requires a GPU to run, see the AlphaFold3 documentation for more detaiils. It is beneficial to run this in parallel. 
An example for what this might look like using [slurm](https://slurm.schedmd.com/documentation.html) can be found for the tutorial data in `slurm_af3_fold_example.sh`. There will be one directory per row in the original sequences CSV file. Each directory will contain 5 X Nseeds subdirectories for each prediction, and will contain the 3D structure .cif file as well as the confidence metadata JSONs. The output will also contain a ranking scores CSV file that ranks the outputs structures, and also a copy of the results for the "best" ranked structure is saved at the top level.

You can see what a successful run looks like by looking at the output for the example above in `examples/af3_fold_outputs`. For more information on the output of AlphaFold3, see [their docs](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md). 

### Step 5: Extract Features from Predicted Structures and Confidence Metadata
```bash
usage: extract_features.py [-h] -s SEQUENCES_FILE [-a ALPHA_COL] [-b BETA_COL] [-m MHC_COL] [-p PEPTIDE_COL] --af-output-dir AF_OUTPUT_DIR [--af3-seeds AF3_SEEDS] -o OUTPUT_DIR

options:
  -h, --help            show this help message and exit
  -s SEQUENCES_FILE, --sequences-file SEQUENCES_FILE
                        Path to input file with TCR-pMHC sequences
  -a ALPHA_COL, --alpha-col ALPHA_COL
                        Name of column containing TCRa chain sequence in <--sequences-file> (default: A_seq)
  -b BETA_COL, --beta-col BETA_COL
                        Name of column containing TCRb chain sequence in <--sequences-file> (default: B_seq)
  -m MHC_COL, --mhc-col MHC_COL
                        Name of column containing MHC chain sequence in <--sequences-file> (default: M_seq)
  -p PEPTIDE_COL, --peptide-col PEPTIDE_COL
                        Name of column containing peptide sequence in <--sequences-file> (default: P_seq)
  --af-output-dir AF_OUTPUT_DIR
                        Directory with AlphaFold3 outputs. This directory should contain one subdirectory per row in the sequences-file
  --af3-seeds AF3_SEEDS
                        comma-separated string of seeds that were input to AF3 (default: 1,2,5,10)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to place output CSV files with features
```
You still need the original sequences file from Step 0 and you will also need:

1) the directory that AlphaFold3 wrote the output to  (this is the output of step 4)
2) the choice of random seeds that you used in step 3 (if you did not use the default), as that determines the file system structure for the outputs

This script will output 3 CSV files in the output directory:
1) all_structures_features.csv - this contains a row for every predicted structure. So for example, if you had 5 seeds, then you will have 25 rows per TCR-pMHC input (5 seeds x 5 samples/per seed)
2) avg_features.csv -  this contains a row for every TCR-pMHC input, where each feature is the mean of that feature across all structures predicted for that TCR-pMHC
3) best_features.csv - this contains a row for every TCR-pMHC input, where each feature is taken from the structure that was ranked the highest by AF3 for that TCR-pMHC
4) ensemble_features.csv -  this contains a row for every TCR-pMHC input, where each feature is saved once with the mean of that feature across all structures predicted for that TCR-pMHC, and once with the standard deviation. This is the "structure ensemble" from the paper.

#### Example
Following our previous example data:
```
python extract_features.py -s examples/example_input_tcr_pmhcs.csv --af-output-dir examples/af3_fold_outputs -o examples/enfoldx_extracted_features
```
You can see the example output in `examples/enfoldx_extracted_features`

## Binding Classification from enFoldX Features
We include an example notebook at `analysis_notebooks/enFoldX_ML_notebook.ipynb` for how we used the features to predict specificity. The enFoldX_ML Jupyter notebook provides an example of training a Ridge regression model on the cross-reactivity dataset utilizing features prepared by the enFoldX AF3 pipeline. It includes preparation of the feature ensemble across all predicted structures for the same TCR/epitope pair and 5-fold cross-validation. The notebook needs fairly standard Python packages to be installed: numpy, pandas, sklearn, matplotlib, and seaborn.In terms of the pipeline, this notebook starts with the feature set that is created at the very end of the last step (Step 5). This notebook uses data from the cross-reactivity data that can be found as part of [a previous paper](https://www.nature.com/articles/s41586-022-04735-9).

Of course, all of these performance evaluations rely on the availability of labeled data. For the use of this tool in studying TCR-pMHCs for which the true label of binder/non-binder is not known, we recommend training on a relevant, related labeled data set, and then using the resultant model to predict specificity for the novel data. For example, you can run enFoldX to create features for all of the human data in VDJdb, then train a Logistic model on that labelled data, and then use the model to predict specificity for a new human dataset.


## Other Misc. Code
### RMSD Calculations
You can also run the rmsd calculation that we use for the specific scenario of comparing two different TCR:pMHC complexes to each other. We use this in the first section of our paper when comparing mutant peptides to wild-type peptides for the same TCR and MHC. You will first need to `pip install biopython`, and then open `rmsd_calculation.py`, change the paths of s1 and s2 to the two structures you want to compare, and run `python rmsd_calculation.py`. 
### Feature Clustering
You can see how we did our feature analysis by looking at the example notebook at `analysis_notebooks/feature_clustering.ipynb`, where we employ our feature importance scoring and feature clustering for the mouse VDJdb dataset.

