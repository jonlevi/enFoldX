#!/bin/bash
#
#SBATCH --job-name=af3_fold
#SBATCH --partition=<partition>
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-gpu=24G
#SBATCH --array=<1-N>

# set path to directory of JSONs
json_dir="path/to/input_jsons"

j=$SLURM_ARRAY_TASK_ID
filename=$(ls -1 "${json_dir}" | sed -n "${j}p")

# docker or singularity depending on how the container is built/installed
singularity --debug exec --nv \
    --bind /path/to/alphafold3:/root/alphafold3 \
    --bind $json_dir:/root/af_input \
    --bind /path/to/output_directory:/root/af_output \
    --bind /path/to/public_databases:/root/public_databases \
    --bind //path/to/weights:/root/models \
    /path/to/alphafold3-container \
    python /root/alphafold3/run_alphafold.py \
   --json_path=/root/af_input/$filename \
   --model_dir=/root/models \
   --db_dir=/root/public_databases \
   --output_dir=/root/af_output \
   --norun_data_pipeline

