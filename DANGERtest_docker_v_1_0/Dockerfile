############################################################
# Dockerfile to build DANGER test v1.0
############################################################
FROM mambaorg/micromamba:1.3.1
COPY env_base.yaml /tmp/env_base.yaml

# Create conda environment

RUN micromamba install -y -n base -f /tmp/env_base.yaml && \
    micromamba clean --all --yes

RUN micromamba create -y -n get_exp_permutation_data_env && \
    micromamba install -y -n get_exp_permutation_data_env -c anaconda python=3.8 pandas=1.5.2 && \
    micromamba clean --all --yes

RUN micromamba create -y -n calcDANGERindex_env && \
    micromamba install -y -n calcDANGERindex_env -c anaconda pandas=1.5.2 scipy=1.10.0 && \
    micromamba install -y -n calcDANGERindex_env -c conda-forge matplotlib=3.6.3 && \
    micromamba clean --all --yes

RUN micromamba create -y -n matplotlib_venn_env && \
    micromamba install -y -n matplotlib_venn_env -c conda-forge matplotlib-venn=0.11.5 && \
    micromamba clean --all --yes

RUN micromamba create -y -n topGO && \
   micromamba install -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-topgo=2.50.0 && \
   micromamba install -y -n topGO -c anaconda -c conda-forge -c bioconda bioconductor-rgraphviz=2.42.0 && \
   micromamba clean --all --yes

# File Author / Maintainer
MAINTAINER Kazuki Nakamae

# Copy DANGER test scripts
WORKDIR /tmp
COPY get_exp_permutation_data.py .
COPY get_multi_exp_permutation_data.sh .
COPY get_mm_permutation_data.py .
COPY get_multi_mm_permutation_data.sh .
COPY annotateGOv2.R .
COPY calcDANGERindex_v2.py .
COPY drawVennDiagram.py .
COPY dindex_test.py .
COPY makeDANGERenrichmentTablev2.R .
COPY 15ggsearch.sh .
COPY 15parseggsearch_linux.sh .
COPY 15parseggsearch.pl .
COPY 15mkannotbl.pl .
USER root
RUN chmod 755 get_multi_exp_permutation_data.sh
RUN chmod 755 get_multi_mm_permutation_data.sh
RUN chmod 755 15ggsearch.sh
RUN chmod 755 15parseggsearch_linux.sh
RUN chmod 755 15parseggsearch.pl
RUN chmod 755 15mkannotbl.pl
COPY dangertest_v1.sh .
RUN chmod 755 dangertest_v1.sh
COPY dangertest_v2.sh .
RUN chmod 755 dangertest_v2.sh

