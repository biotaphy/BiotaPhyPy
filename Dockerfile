FROM conda/miniconda3:latest as biotaphypy

LABEL maintainer="Biotaphy <github.com/biotaphy>"
LABEL description="A Docker image for Biotaphy Tools"

# RUN apt-get update && apt-get install -y .....

RUN conda update -n base -c conda-forge conda && \
    conda install -y -c conda-forge gdal libspatialindex rtree

ENV PROJ_LIB=/usr/local/share/proj/

# RUN conda update -n base -c conda-forge conda && \
#     conda env create -f environment.yml
# SHELL ["conda", "run", "-n", "lmpy", "/bin/bash", "-c"]

# Copy in our library and install
WORKDIR /biotaphy
COPY . .
RUN pip install .
