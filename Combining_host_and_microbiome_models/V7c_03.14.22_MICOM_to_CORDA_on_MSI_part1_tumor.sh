#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 12:00:00
#SBATCH --array=1-1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

cd Merging_Burns_MICOM_and_CORDA/
export PATH=/home/blekhman/adamo010/miniconda3/envs/gurobi_micom/bin:$PATH
export GRB_LICENSE_FILE=/panfs/roc/msisoft/gurobi/license/gurobi.lic
file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Burns2015_tumor_pickle_file_sublist1.txt)
python 03.14.22_MICOM_to_CORDA_on_MSI_part1_tumor_V4a.py $file
