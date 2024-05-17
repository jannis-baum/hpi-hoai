FROM ubuntu:latest

RUN apt-get update

# ------------------------------------------------------------------------------
# install required tools -------------------------------------------------------
ARG TOOLS_DIR=/tools
WORKDIR ${TOOLS_DIR}
# sra-tools
ADD https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.0/sratoolkit.3.1.0-ubuntu64.tar.gz sratoolkit.tar.gz
RUN tar -vxzf sratoolkit.tar.gz
ENV PATH="${PATH}:${TOOLS_DIR}/sratoolkit.3.1.0-ubuntu64/bin"
