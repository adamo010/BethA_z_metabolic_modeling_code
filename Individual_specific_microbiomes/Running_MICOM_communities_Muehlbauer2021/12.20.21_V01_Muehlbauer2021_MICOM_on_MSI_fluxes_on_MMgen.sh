#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 12:00:00
#SBATCH --array=1-12%6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu

cd Combining_Muehlbauer_models/
source activate gurobi_micom
export GRB_LICENSE_FILE=/panfs/roc/msisoft/gurobi/license/gurobi9_2021.2.lic
file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Muehlbauer2021_pickle_file_list.txt)
python 12.20.21_V01_Muehlbauer2021_MICOM_on_MSI_fluxes_on_MMgen_TDpointone.py $file


