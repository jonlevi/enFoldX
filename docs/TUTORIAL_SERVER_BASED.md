### Server Based Predictions
Go to [https://alphafoldserver.com/](https://alphafoldserver.com/) and sign-in with a gmail account and agree to any necessary terms. If you want, you can manually add the sequences using the GUI, but we provide an easier way to generate multiple jobs using a single input JSON, using a handy python script included in this repo. Please note, you must process your data in batches of 100 or fewer rows at a time, since the server only accepts JSONs with 100 jobs max at a time. Also note that although you upload a JSON with many jobs, the current quota for AF3 server will only allow ~20 jobs per person per day.
To generate the JSON that can be uploaded to the AF3 server, you can run the handy script:
```bash
usage: prepare_server_input.py [-h] -s SEQUENCES_FILE [-a ALPHA_COL] [-b BETA_COL] [-m MHC_COL] [-p PEPTIDE_COL] [--af3-seed AF3_SEED] -o OUTPUT_DIR

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
  --af3-seed AF3_SEED   AF3 Seed for Folding. Note that with the server, you can only input one seed at a time (default: 1)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to place output JSON files (default: None)
```
#### Example
Following our previous example data:
```
python prepare_server_input.py -s examples/example_input_tcr_pmhcs.csv -o examples/af3_server_inputs
```
You can see the example output in `examples/af3_server_inputs`.

You can now use the "Upload JSON" button in the top right of alphafoldserver.com and upload the JSON that was output from the previous script. This will pre-populate one job per row of your input CSV into your draft jobs on the bottom of the screen. You must submit each job individually by clicking on the job row, and then clicking the blue button that says "continue and preview job". Once the jobs complete you can download the outputs and continue with the steps below.

### Extract features from AlphaFold3 Server Results
Once the jobs finish, you can select as many jobs as you ran and click the download icon. You should get a large zip file with a name that looks something like "fold_\<date\>.zip". This zip can be used directly by the script to extract results:
```bash
usage: extract_features_server.py [-h] -s SEQUENCES_FILE [-a ALPHA_COL] [-b BETA_COL] [-m MHC_COL] [-p PEPTIDE_COL] -z ZIP_FILE -o OUTPUT_DIR

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
  -z ZIP_FILE, --zip-file ZIP_FILE
                        Zip with AlphaFold3 Server Outputs. This directory should contain one subdirectory per job downloaded (default: None)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to place output CSV files with features (default: None)
```

#### Example
Following the example above, I selected 2/4 jobs to download and ran:
```
python extract_features_server.py -o ./examples/enfoldx_extracted_features_server -z ./examples/af3_server_outputs/folds_2026_04_21_20_08.zip -s ./examples/example_input_tcr_pmhcs.csv
```
You can see the example download in `examples/af3_server_outputs`, and the example extracted features in `examples/enfoldx_extracted_features_server`
