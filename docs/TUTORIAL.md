# Tutorial for enFoldX with locally installed AlphaFold3

## Locally Installing AlphaFold3
To run AlphaFold3 predictions on your own compute system, you first must install AlphaFold3 following the instructions at https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md. As is outlined in detail there, you will need to request access to the AF3 parameters, and install about 1TB of data along with a docker/singularity container. 

After the installation you should have the following handy:
- path to the AF3 parameters
- path to downloaded databases
- path to the docker/singularity container for running alphafold

Please keep in mind that you must additionally comply with the [TERMS OF USE](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md) set by Deepmind in order to use AlphaFold3. Please also keep in mind that in order to run AlphaFold3, you may need access to specialized hardware (see [https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md](https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md)). 


### Configure AlphaFold Paths:
Open `af3_config.sh` and replace the stubbed out paths with the correct paths that contain the AlphaFold3 installation from above:
```
ALPHAFOLD_DIR="path/to/alphafold3"
DATABASE_DIR="path/to/public_databases"
WEIGHTS_DIR="path/to/weights"
CONTAINER_PATH="path/to/container"
```

### Prediction Pipeline

There are a few sequential steps in the pipeline to go from a TCR-pMHC sequence --> predicted structure ensemble:
1) Format JSONs for running AlphaFold3 MSA for each unique sequence
2) Run MSA on Input Sequences (CPU)
3) Collect MSA results and Format JSONs per TCR-pMHC for running AlphaFold3 Folding
4) Run AlphaFold3 Container on Input JSONs (GPU)

The key idea of this pipeline is to set things up so that we can parallelize as much as possible. The MSA is run on each unique sequence independently and thus can be parallelized across each one. By splitting up the MSA and folding steps, we can re-use TCR MSAs when possible, as often a user will want to fold a TCR with many potential antigens. The folding is run on each TCR-pMHC complex independently, and thus can be parallelized as well if you have access to multiple compute nodes. This section goes through a full example of going from TCR:pMHC sequences --> enFoldX features. All of the data for the examples can be found in `examples/`

#### Pipeline Step 1: Format JSONs for running AlphaFold3 MSA for each unique sequence
```bash
usage: prepare_msa_input.py [-h] -s SEQUENCES_FILE [-a ALPHA_COL] [-b BETA_COL] [-m MHC_COL] [-p PEPTIDE_COL] -o
                            OUTPUT_DIR

options:
  -h, --help            show this help message and exit
  -s SEQUENCES_FILE, --sequences-file SEQUENCES_FILE
                        Path to input file with TCR-pMHC sequences (default: None)
  -a ALPHA_COL, --alpha-col ALPHA_COL
                        Name of column containing TCRa chain sequence in <--sequences-file> (default: TRA_aa)
  -b BETA_COL, --beta-col BETA_COL
                        Name of column containing TCRb chain sequence in <--sequences-file> (default: TRB_aa)
  -m MHC_COL, --mhc-col MHC_COL
                        Name of column containing MHC chain sequence in <--sequences-file> (default: M_aa)
  -p PEPTIDE_COL, --peptide-col PEPTIDE_COL
                        Name of column containing peptide sequence in <--sequences-file> (default: peptide)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to place output JSON files (default: None)
```
*Example*:
```python ./scripts/prepare_msa_input.py -s examples/MEL_data/MEL_enfoldx_input.csv -o examples/af3_msa_inputs```

This script will write out one JSON per unique TCRa, TCRb, and MHC sequence in the data, as well as a metadata file mapping the sequences to their chain IDs/JSON paths. 

You can see what a successful run looks like by looking at the output for the tutorial in `examples/af3_msa_inputs`


#### Pipeline Step 2: Run MSA on Input Sequences
You will run the AlphaFold3 container on `--no-run-inference` mode once per MSA JSON that was created. It is beneficial to run this in parallel across your dataset. Note that each job runs on a CPU node, however you can allocate multiple threads for AF3 to speed up jackhmmer. For example, in slurm, to add 8 threads you could use `#SBATCH --cpus-per-task=8`.

An example for what a submission script might look like using [slurm](https://slurm.schedmd.com/documentation.html) can be found for the tutorial data in `slurm_af3_msa_example.sh`. If you use this submission script, make sure to adjust the partition name, the array size, the compute requests, and the input/output paths to match your preferences/requirements. Note, the array size will be larger than the size of your input CSV, since you should be running one job per unique TCRa, TCRb, and MHC. You can check what the array size for slurm should be by looking at the the number of input jsons created in Step 1 (e.g. `ls examples/af3_msa_inputs/*.json | wc -l`). For our MEL TCRs, the array should be 5, 2 for each TCR and 1 for the HLA. You would run the script by invoking ``` sbatch ./examples/slurm_af3_msa_example.sh ```.

You can see what a successful MSA run looks like by looking at the output for the tutorial above in `examples/af3_msa_outputs` (although these JSON files are very large and can be difficult to view using regular IDEs or the github file browser). There will be one JSON per unique chain (except for peptides which don't get MSA), and a metadata chain mapping file to map the MSA output to the original sequences. (You will need this chain_id_map for the next step). 

Troubleshooting tip: One common error we have seen is that if your TCR sequences have a "*" character in them, which is sometimes used for a stop codon, then AlphaFold will fail. Make sure to remove those TCRs, as AlphaFold cannot model them.

#### Pipeline Step 3: Collect MSA results and Format JSONs per TCR-pMHC for running AlphaFold3 Folding
```bash
usage: prepare_fold_input.py [-h] -s SEQUENCES_FILE [-i ID_COL] [-a ALPHA_COL] [-b BETA_COL] [-m MHC_COL] [-p PEPTIDE_COL]
                             --chain-id-map CHAIN_ID_MAP --MSA-output-dir MSA_OUTPUT_DIR [--af3-seeds AF3_SEEDS] -o
                             OUTPUT_DIR

options:
  -h, --help            show this help message and exit
  -s SEQUENCES_FILE, --sequences-file SEQUENCES_FILE
                        Path to input file with TCR-pMHC sequences (default: None)
  -i ID_COL, --id-col ID_COL
                        Name of column containing unique complex ID in <--sequences-file> (default: complex_id)
  -a ALPHA_COL, --alpha-col ALPHA_COL
                        Name of column containing TCRa chain sequence in <--sequences-file> (default: TRA_aa)
  -b BETA_COL, --beta-col BETA_COL
                        Name of column containing TCRb chain sequence in <--sequences-file> (default: TRB_aa)
  -m MHC_COL, --mhc-col MHC_COL
                        Name of column containing MHC chain sequence in <--sequences-file> (default: M_aa)
  -p PEPTIDE_COL, --peptide-col PEPTIDE_COL
                        Name of column containing peptide sequence in <--sequences-file> (default: peptide)
  --chain-id-map CHAIN_ID_MAP
                        File that maps chain IDs to sequences from MSA run. Should be an output from the
                        <prepare_msa_input> step (default: None)
  --MSA-output-dir MSA_OUTPUT_DIR
                        Directory with MSA outputs. This directory should contain one subdirectory per ID in the chain-id-
                        map, each containing a single file named <ID>_data.json (default: None)
  --af3-seeds AF3_SEEDS
                        comma-separated string of seeds to input to AF3; default is 10 seeds (default:
                        1,13,17,21,42,133,177,213,315,1001)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to place output JSON files to pass to AlphaFold3 (default: None)
```
You still need the original sequences file from Step 0 and you will also need:

1) the path that the chain id map was written to in Step 1 and 
2) the directory that AlphaFold3 wrote the MSA output to. (Note this is the output of step 2, not the output of step 1)

You can choose how many [random seeds](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md#random-seeds) to use in the AF3 inference. By default, this pipeline uses 5 seeds per complex. You can add more or remove them by changing the `--af3-seeds` flag. Note that each additional seeds adds about an additional 1-2 minutes at minimum per complex to the runtime, but your predictive accuracy may improve with more seeds. See the our paper for more!

*Example*:
Following our previous tutorial data:
```
python ./scripts/prepare_fold_input.py -s examples/MEL_data/MEL_enfoldx_input.csv -o examples/af3_fold_inputs --chain-id-map examples/af3_msa_inputs/chain_ids_to_sequences.txt --MSA-output-dir examples/af3_msa_outputs
```
You can see what a successful run looks like by looking at the output for the tutorial in `examples/af3_fold_inputs` (although these JSON files are very large and can be difficult to view using regular IDEs or the github file browser). 

#### Pipeline Step 4: Run AlphaFold3 Folding on Input Sequences
You will run the AlphaFold3 container on `--norun_data_pipeline` mode once per TCR-pMHC input JSON that was created. This requires a GPU to run, see the AlphaFold3 documentation for more details. It is beneficial to run this in parallel. 
An example for what this might look like using [slurm](https://slurm.schedmd.com/documentation.html) can be found for the tutorial data in `slurm_af3_fold_example.sh`. There will be one directory per row in the original sequences CSV file. Each directory will contain 5 X Nseeds subdirectories for each prediction, and will contain the 3D structure .cif file as well as the confidence metadata JSONs. The output will also contain a ranking scores CSV file that ranks the outputs structures, and also a copy of the results for the "best" ranked structure is saved at the top level. You can invoke the script by running ``` sbatch ./examples/slurm_af3_fold_example.sh ```.

You can see what a successful run looks like by looking at the output for the example above in `examples/af3_fold_outputs`. For more information on the output of AlphaFold3, see [their docs](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md). 

## Extract features from AlphaFold3 results

#### Pipeline Step 5: Run AlphaFold3 Feature Extraction on Output Structures
```bash
usage: extract_features.py [-h] -s SEQUENCES_FILE --af-output-dir AF_OUTPUT_DIR [--af3-seeds AF3_SEEDS] -o OUTPUT_DIR
                           [-i ID_COL] [-a ALPHA_COL] [-b BETA_COL] [-cdr3a CDR3ALPHA_COL] [-cdr3b CDR3BETA_COL]
                           [-m MHC_COL] [-p PEPTIDE_COL]

options:
  -h, --help            show this help message and exit
  -s SEQUENCES_FILE, --sequences-file SEQUENCES_FILE
                        Path to input file with TCR-pMHC sequences (default: None)
  --af-output-dir AF_OUTPUT_DIR
                        Directory with AlphaFold3 outputs. This directory should contain one subdirectory per row in the
                        sequences-file (default: None)
  --af3-seeds AF3_SEEDS
                        Comma-separated string of seeds that were input to AF3 (default:
                        1,13,17,21,42,133,177,213,315,1001)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to place output CSV files with features (default: None)
  -i ID_COL, --id-col ID_COL
                        Name of column containing unique complex ID in <--sequences-file> (default: complex_id)
  -a ALPHA_COL, --alpha-col ALPHA_COL
                        Name of column containing TCRa chain sequence in <--sequences-file> (default: TRA_aa)
  -b BETA_COL, --beta-col BETA_COL
                        Name of column containing TCRb chain sequence in <--sequences-file> (default: TRB_aa)
  -cdr3a CDR3ALPHA_COL, --cdr3alpha-col CDR3ALPHA_COL
                        Name of column containing CDR3a chain sequence in <--sequences-file> (default: TRA_CDR3)
  -cdr3b CDR3BETA_COL, --cdr3beta-col CDR3BETA_COL
                        Name of column containing CDR3b chain sequence in <--sequences-file> (default: TRB_CDR3)
  -m MHC_COL, --mhc-col MHC_COL
                        Name of column containing MHC chain sequence in <--sequences-file> (default: M_aa)
  -p PEPTIDE_COL, --peptide-col PEPTIDE_COL
                        Name of column containing peptide sequence in <--sequences-file> (default: peptide)
```
You still need the original sequences file from Step 0 and you will also need:

1) the directory that AlphaFold3 wrote the output to  (this is the output of step 4)
2) the choice of random seeds that you used in step 3 (if you did not use the default), as that determines the file system structure for the outputs

This script will output 4 CSV files in the output directory:
1) all_structures_features.csv - this contains a row for every predicted structure. So for example, if you had 5 seeds, then you will have 25 rows per TCR-pMHC input (5 seeds x 5 samples/per seed)
2) avg_features.csv -  this contains a row for every TCR-pMHC input, where each feature is the mean of that feature across all structures predicted for that TCR-pMHC
3) best_features.csv - this contains a row for every TCR-pMHC input, where each feature is taken from the structure that was ranked the highest by AF3 for that TCR-pMHC
4) ensemble_features.csv -  this contains a row for every TCR-pMHC input, where each feature is saved once with the mean of that feature across all structures predicted for that TCR-pMHC, and once with the standard deviation. This is the "structure ensemble" from the paper.

#### Example
Following our previous example data:
```
python ./scripts/extract_features.py -s examples/MEL_data/MEL_enfoldx_input.csv --af-output-dir examples/af3_fold_outputs -o examples/enfoldx_extracted_features
```
You can see the example output in `examples/enfoldx_extracted_features`
