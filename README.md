# Locate LOH in WGS data

## Install conda

## Create nextflow conda environment
```{bash}
conda create -n nextflow nextflow pandas -c bioconda
conda activate nextflow
```

# Make sample plan for each sample
```{bash}
python bin/write_files.py -d /data/kdi_prod/project_result/948/01.00/Analysis/D050 -o sample_plan.csv
```

```{bash}
conda env export --from-history -n name_of_your_env -f environment.yml
```
