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

# samtools
RUN apt install -y bzip2 build-essential libncurses5-dev libncursesw5-dev zlib1g-dev libbz2-dev liblzma-dev
ADD https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 samtools-1.20.tar.bz2
RUN tar -vxf samtools-1.20.tar.bz2
WORKDIR ${TOOLS_DIR}/samtools-1.20
RUN ./configure --prefix=${TOOLS_DIR}/samtools-build
RUN make
RUN make install
WORKDIR ${TOOLS_DIR}
ENV PATH="${PATH}:${TOOLS_DIR}/samtools-build/bin"

# hisat 2
RUN apt install -y unzip
ADD https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download hisat2-2.2.1.zip
RUN unzip hisat2-2.2.1.zip
ENV PATH="${PATH}:${TOOLS_DIR}/hisat2-2.2.1"

# stringtie
ADD https://github.com/gpertea/stringtie/releases/download/v2.2.3/stringtie-2.2.3.Linux_x86_64.tar.gz stringtie-2.2.3.Linux_x86_64.tar.gz
RUN tar -vxzf stringtie-2.2.3.Linux_x86_64.tar.gz
ENV PATH="${PATH}:${TOOLS_DIR}/stringtie-2.2.3.Linux_x86_64"
