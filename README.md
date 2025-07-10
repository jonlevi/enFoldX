# enFoldX: Leveraging ensembles of predicted structures for complex binding prediction


enFoldX (**En**semble of **Fold**ed Comple**x**es)enables the evaluation of large datasets of predicted structures simultaneously as it extracts distributions of high dimensional confidence features, including customized biophysical features. We found these features varied significantly between TCR:pMHC pairs of differing T cell reactivity, implying that binders and non-binders are separable in feature space. We then use the resultant ensemble of features to create a binding model, analogous to an energy function in statistical physics. The enFoldX approach generalizes in a dataset agnostic way - for example, models trained on human data can classify binder TCR:pMHC pairs in mouse data. Moreover, we found, for the first time, that our structure-guided models trained on publicly available data can generalize to discriminate epitopes that differ by only one mutation away from self, the resolution needed for predicting mutation-derived cancer neoantigens.

You can see our paper here <>:

This repo provides code for a pipeline to run enFoldX. The current implementation of enFoldX uses AlphaFold3 and runs predictions for complex binding betweens TCRs and peptide-MHC complexes. 

![Project Banner](https://github.com/jonlevi/af3_tcr_pipeline/blob/main/enfoldx_diagram.pdf) 

## 📚 Table of Contents
- [Terms of use](#terms-of-use)
- [Installation](#installation)
- [Usage](#usage)

## Terms of Use
By using enFoldX, you are agreeing to the terms set in the enFoldX [Terms of Use](TERMS_OF_USE.md).

## Installation
### Installing AlphaFold3
If you wish to use enFoldX with AlphaFold3 predictions, you must install AlphaFold3 following the instructions at https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md 

After the installation you should have:
- path to the AF3 parameters
- path to downloaded databases
- path to the docker/singularity container for running alphafold

Please keep in mind that you must additional comply with the [TERMS OF USE](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md) set by Deepmind in order to use AlphaFold3.

### Installing this repository
Clone the repository:

```bash
git clone https://github.com/jonlevi/enFoldX.git
cd enFoldX
```

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

### Step 0: Format Sequences Input
In order to run this pipeline, you will need a CSV that contains one row for every TCR-pMHC complex you wish to predict. Each row should have four columns containing the sequences of the TCRa, TCRb, MHC, and Peptide chains using 1-letter amino acid codes. By default, the pipeline assumes that the columns are named ["A_seq","B_seq","M_seq","P_seq"] respectively, although you can pass in custom column names using the optional flags to each script. If you are starting from VDJ+CDR3 calls, we recommend you use [stitchr](https://jamieheather.github.io/stitchr/index.html) to get full length TCRa and TCRb sequences. For MHC sequence information, you can look up the allele in [Uniprot](https://www.uniprot.org/uniprotkb) or IPD-IMGT/HLA at https://www.ebi.ac.uk/ipd/imgt/hla/alleles/.


Our tutorial example file can be found in `examples/example_input_tcr_pmhcs.csv`

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

### From these features, you can then train models to predict binding. We do not include code for that in this repository, but you can see our paper for more on train/test splits and on publicly available data. 

## Other Misc. Code
You can also run the rmsd calculation that we use for the specific scenario of comparing two different TCR:pMHC complexes to each other. We use this in the first section of our paper when comparing mutant peptides to wild-type peptides for the same TCR and MHC. You will first need to `pip install biopython`, and then open `rmsd_calculation.py`, change the paths of s1 and s2 to the two structures you want to compare, and run `python rmsd_calculation.py`. 

