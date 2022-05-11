#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH --array=1-2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

cd Combining_Niccolai_models/
source activate gurobi_micom
export GRB_LICENSE_FILE=/panfs/roc/msisoft/gurobi/license/gurobi9_2021.2.lic
file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Niccolai2020_OTU_abundances_subset1.txt)
python 10.22.21_V01_building_Niccolai2020_micom_communities_MSI.py $file
