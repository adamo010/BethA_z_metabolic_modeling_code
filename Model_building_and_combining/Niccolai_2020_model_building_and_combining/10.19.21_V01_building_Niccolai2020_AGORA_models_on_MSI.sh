#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=32gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

cd Combining_Niccolai_models/
source activate gurobi_micom
export GRB_LICENSE_FILE=/panfs/roc/msisoft/gurobi/license/gurobi9_2021.2.lic
python 10.19.21_V01_building_Niccolai2020_AGORA_models_on_MSI.py Niccolai2020_AGORA_models_that_need_OTU_builds.pickle 

