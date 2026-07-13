# Tutorial for enFoldX with AlphaFold3 server predictions

All of the data for the examples can be found either in `examples/` or the [downloaded tutorial data](https://github.com/jonlevi/enFoldX/releases/download/v0.1.0/enfoldx_tutorial_data.zip).

**WARNING: At the moment, pre-trained enFoldX models cannot be used with server predictions**

## Server Based Predictions
Firs, go to [https://alphafoldserver.com/](https://alphafoldserver.com/) and sign-in with a gmail account and agree to any necessary terms.

If you want, you can manually add the sequences using the GUI and run them directly, but for effiency, we provide an easier way to input multiple jobs using a single input JSON, using a handy python script included in this repo. Please note, you must process your data in batches of 100 or fewer rows at a time, since the server only accepts JSONs with 100 jobs max at a time. Also note that although you upload a JSON with many jobs, the current quota for AF3 server will only allow ~20 jobs per person per day, so it may take time for the jobs to complete.
To generate the JSON that can be uploaded to the AF3 server, you can run the handy script:
```bash
usage: prepare_server_input.py [-h] -s SEQUENCES_FILE [-a ALPHA_COL] [-b BETA_COL] [-m MHC_COL] [-p PEPTIDE_COL]
                               [-i ID_COL] [--af3-seed AF3_SEED] -o OUTPUT_DIR

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
  -i ID_COL, --id-col ID_COL
                        Name of column containing unique complex ID in <--sequences-file> (default: complex_id)
  --af3-seed AF3_SEED   AF3 Seed for Folding. Note that with the server, you can only input one seed at a time (default:
                        1)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to place output JSON file (default: None)
```
### Example
Following our previous example data:
```
python ./scripts/prepare_server_input.py -s examples/MEL_data/MEL_enfoldx_input.csv -o examples/af3_server_inputs
```
You can see the output in `examples/af3_server_inputs` or the example tutorial data in `enfoldx_tutorial_data/af3_server_inputs`.

You can now use the "Upload JSON" button in the top right of alphafoldserver.com and upload the JSON that was output from the previous script. This will pre-populate one job per row of your input CSV into your draft jobs on the bottom of the screen. You must submit each job individually by clicking on the job row, and then clicking the blue button that says "continue and preview job". Once the jobs complete you can download the outputs and continue with the steps below.

### Extract features from AlphaFold3 Server Results
Once the jobs finish, you can select as many jobs as you ran and click the download icon. You should get a large zip file with a name that looks something like "fold_\<date\>.zip". Following the example above, I selected all 18 MEL TCR jobs to download in one zip. If you want to follow along with the tutorial data you can download it from [this link](https://github.com/jonlevi/enFoldX/releases/download/v0.1.0/folds_2026_07_09_17_57.zip) and move it to `examples/af3_server_outputs`.

This zip can be used directly by the script to extract results:
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

```
python ./scripts/extract_features_server.py  -s examples/MEL_data/MEL_enfoldx_input.csv -o ./examples/enfoldx_extracted_features_server -z ./examples/af3_server_outputs/folds_2026_07_09_17_57.zip
```
The example extracted features were then written to `examples/enfoldx_extracted_features_server`

#### A note on multiple seeds with AF3 Server
A key result from our manuscript is that enFoldX works better with multiple random seeds. Unfortunately, the current server-based workflow only lets you run a single seed at a time. If you are using the GUI version, you can "Specify a seed" using the input box. If you are using the script you can specify a seed using ```--af3-seed```. If you want to create a multi-seed ensemble, you should run all of the steps above multiple times for the same input, changing the seed each time. Then, instead of using the ```ensemble_features.csv```, you should use the ```all_structures_features.csv```. You can then create your own ```ensemble_features.csv``` by running this snippet:

```
all_results_df = pd.read_csv('all_structures_features.csv',index_col=0)
ensemble = all_results_df.groupby(by=["original_index"], dropna=False).agg(
        ["mean", "std"]
    )
ensemble.columns = [f"{col}_{stat}" if stat=='std' else col for col, stat in ensemble.columns]
columns_to_drop = [
        col
        for col in ensemble.columns
        if col.startswith(("af3_seed", "af3_sample", "af3_ranking"))
]
ensemble = ensemble.drop(columns=columns_to_drop).reset_index()
ensemble.to_csv("ensemble_features.csv")
```
