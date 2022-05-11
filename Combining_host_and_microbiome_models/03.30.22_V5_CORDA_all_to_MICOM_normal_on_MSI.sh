#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 12:00:00
#SBATCH --array=1-44%8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

cd Merging_Burns_MICOM_and_CORDA/
export PATH=/home/blekhman/adamo010/miniconda3/envs/gurobi_micom/bin:$PATH
export GRB_LICENSE_FILE=/panfs/roc/msisoft/gurobi/license/gurobi9_2021.2.lic
file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Burns2015_normal_pickle_file_list.txt)
python 03.30.22_CORDA_all_to_MICOM_normal_on_MSI.py $file
