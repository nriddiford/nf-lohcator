FROM nfcore/base:1.9
LABEL authors="Nick Riddiford" \
      description="Docker image containing all software requirements for the bwa-align processes"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

RUN echo "source activate $(head -1 environment.yml | cut -d' ' -f2)" > ~/.bashrc

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/$(head -1 environment.yml | cut -d' ' -f2)/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name $(head -1 environment.yml | cut -d' ' -f2) > $(head -1 environment.yml | cut -d' ' -f2).yml
