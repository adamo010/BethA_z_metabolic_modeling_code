#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 20:00:00
#SBATCH --array=1-1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

cd Combining_Niccolai_models/
source activate gurobi_micom
export GRB_LICENSE_FILE=/panfs/roc/msisoft/gurobi/license/gurobi9_2021.2.lic
#file= T0_B_PR1_CRC_36_community.pickle
python 11.02.21_V01a_Niccolai2020_MICOM_on_MSI_fluxes_on_MMgen_TDpointtwofive.py T0_B_PR1_CRC_36_community.pickle


