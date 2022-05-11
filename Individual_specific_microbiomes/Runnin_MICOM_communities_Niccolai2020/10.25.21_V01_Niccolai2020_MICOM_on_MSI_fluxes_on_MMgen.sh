#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 12:00:00
#SBATCH --array=1-80%8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

cd Combining_Niccolai_models/
source activate gurobi_micom
export GRB_LICENSE_FILE=/panfs/roc/msisoft/gurobi/license/gurobi9_2021.2.lic
file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Niccolai2020_pickle_file_list.txt)
python 10.25.21_V01_Niccolai2020_MICOM_on_MSI_fluxes_on_MMgen.py $file


