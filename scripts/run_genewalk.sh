#!/bin/sh
#SBATCH -p short
#SBATCH -J genewalk
#SBATCH -t 0-12:00
#SBATCH -c 4
#SBATCH --mem=8000

genes="$1"
filename=$(basename "$genes")
project="${filename%.*}"

set -eux

genewalk --project "$project" --base_folder $(pwd) --id_type ensembl_id --nproc 4 --genes "$genes"
