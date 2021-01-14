#!/bin/sh
set -euo pipefail

b_green="\033[1;32m"
b_purple="\033[1;35m"

c_reset="\033[0m"

mkdir -p log


install_miniconda() {
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm -rf ~/miniconda3/miniconda.sh
    ~/miniconda3/bin/conda init bash
    ~/miniconda3/bin/conda init zsh
}


install_nextflow(){
    conda install --yes -c bioconda nextflow
}


check_installaton(){
    program=$1
    install_command=$2
    if ! command -v $program &> /dev/null
    then
        echo -e "  ${b_purple}→${c_reset} $program is not installed. Installing now..."
        echo "$install_command"
        $install_command
        exit
    else
        echo -e "  ${b_green}✔${c_reset} $program already installed"
    fi
}


check_installaton "conda" "install_miniconda > log/miniconda_install.log 2>&1"


check_installaton "nextflow" "install_nextflow > log/nextflow_install.log 2>&1"
