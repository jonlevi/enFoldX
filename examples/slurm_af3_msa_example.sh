#!/bin/bash
#
#SBATCH --job-name=af3_msa
#SBATCH --partition=componc_cpu
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --array=1-5

# ^ set array size to number of input JSONs

# load config
source ./af3_config.sh

# set input path to directory of JSONs
json_dir="examples/af3_msa_inputs"

# set output path to desired directory
OUTDIR="examples/af3_msa_outputs"
mkdir -p "$OUTDIR"

# j is the array index for slurm. Grab the jth JSON filename
j=$SLURM_ARRAY_TASK_ID
filepath=$(ls -1 ${json_dir}/*.json | sed -n "${j}p")
filename=$(basename "$filepath")

# invoke docker or singularity depending on how the container is built/installed
# mount the paths for AF3 (saved in config file) and inputs and outputs
# running with --norun_inference does MSA only
singularity --debug exec --nv \
    --bind "$ALPHAFOLD_DIR:/root/alphafold3" \
    --bind "$DATABASE_DIR:/root/public_databases" \
    --bind "$WEIGHTS_DIR:/root/models" \
    --bind "$json_dir:/root/af_input" \
    --bind "$OUTDIR:/root/af_output" \
    "$CONTAINER_PATH" \
    python /root/alphafold3/run_alphafold.py \
        --json_path=/root/af_input/$filename \
        --model_dir=/root/models \
        --db_dir=/root/public_databases \
        --output_dir=/root/af_output \
        --norun_inference