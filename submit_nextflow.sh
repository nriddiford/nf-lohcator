#!/usr/bin/sh

echo "Submitting nf-lohcator for sample plan: $1"
out_dir=$(pwd)
user=$(whoami)

mkdir -p $out_dir/log

ts=$(date +%d_%m_%Y_%H%M%S)

logfile=${1%%.*}_${ts}.log

echo $logfile

qsub -v VAR1=$out_dir,VAR2=$1 -o $out_dir/log/$logfile -j oe -N ${user}_submit_nf-lohcator $out_dir/bin/submit_nextflow.pbs
