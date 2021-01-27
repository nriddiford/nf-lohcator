#!/usr/bin/sh

echo "Submitting nf-lohcator for sample plan: $1"
out_dir=$(pwd)
user=$(whoami)

mkdir -p $out_dir/log

qsub -v VAR1=$out_dir,VAR2=$1 -o $out_dir/log/${1%%.*}.log -j oe -N ${user}_submit_nf-lohcator $out_dir/bin/submit_nextflow.pbs
