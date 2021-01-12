# Locate LOH in WGS data

## Install nextflow

curl -s https://get.nextflow.io | bash

# Make sample plan for each sample
python bin/write_files.py -d /data/kdi_prod/project_result/948/01.00/Analysis/D050 -o sample_plan.csv


conda env export --from-history -n name_of_your_env -f environment.yml
