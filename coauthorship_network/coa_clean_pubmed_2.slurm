#!/bin/bash
#SBATCH --job-name=coA_clean_2
#SBATCH --output=clean_loc_bash/clean_loc_2.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --time=6-00:00:00
#SBATCH --qos=long
#SBATCH --mail-user=EUS16@pitt.edu
#SBATCH --mail-type=END,FAIL

module purge
module load python/3.7.0

python clean_pubmed_data.py -i 2 > clean_loc_out/clean_loc_out_2.txt
