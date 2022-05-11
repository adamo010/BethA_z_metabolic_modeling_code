#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=32gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

cd Combining_Muehlbauer_models/
source activate gurobi_micom
export GRB_LICENSE_FILE=/panfs/roc/msisoft/gurobi/license/gurobi9_2021.2.lic
python 12.17.21_V01_building_Muehbauer2021_AGORA_models_on_MSI.py Muehlbauer2021_AGORA_models_that_need_OTU_builds_subset2.pickle 



