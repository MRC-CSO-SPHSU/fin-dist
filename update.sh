#!/bin/bash 

# Utilities
# This bash executable will update the fst dataset and upload it to the cluster

set -euo pipefail

cd "$(dirname "$0")"

HPC_USER=dd198b
HPC_HOST=mars-login.ice.gla.ac.uk
HPC_DIR=/users/dd198b/Documents/GitHub/

R -e 'source(here::here("R", "import_cleaning.R")); import_data(force = TRUE) |> clean_data() |> preproc_data()'
scp -r data/cache/raw_data.fst $HPC_USER@$HPC_HOST:$HPC_DIR/financial_distress/data/cache/raw_data.fst
