#!/bin/sh
set -euo pipefail

mkdir -p log


check_installaton(){
    program=$1
    install_command=$2
    if ! command -v $program &> /dev/null
    then
        echo "$program is not installed. Installing now..."
        echo "$install_command"
        exit
    else
        echo "$program already installed"
    fi
}


check_installaton "conda4" "install_miniconda > log/miniconda_install.log 2>&1"


install_miniconda () {
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm -rf ~/miniconda3/miniconda.sh
    ~/miniconda3/bin/conda init bash
    ~/miniconda3/bin/conda init zsh
}
