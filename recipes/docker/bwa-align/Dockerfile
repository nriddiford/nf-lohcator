FROM nfcore/base:1.9
LABEL authors="Nick Riddiford" \
      description="Docker image containing all software requirements for the nf-bwa-mem pipeline"

# Install the conda environment
COPY ../../conda/bwa.yml /
RUN conda env create -f ../../conda/bwa.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-bwa-mem-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-bwa-mem-1.0dev > nf-bwa-mem-1.0dev.yml
