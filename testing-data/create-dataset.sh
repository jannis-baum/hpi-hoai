#!/bin/sh

if [ "$#" -ne 1 -o -z "$1" ]; then
    echo "usage: create-dataset <SRA accession>"
    exit 1
fi

N_CORES=`getconf _NPROCESSORS_ONLN`
# samtools wants to know number of "additional threads", i.e. one less
N_ADD_THREADS=$(($N_CORES - 1))

SCRIPT_PATH=`dirname "$0"`; SCRIPT_PATH=`eval "cd \"$SCRIPT_PATH\" && pwd"`
REF_DIR="$SCRIPT_PATH/data/reference"

REF_GTF="$REF_DIR/Homo_sapiens.GRCh38.104.gtf"
if ! test -f "$REF_GTF"; then
    echo "Reference GTF not found. Please download to $REF_DIR as follows"
    echo '```'
    echo "wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz"
    echo "gunzip Homo_sapiens.GRCh38.104.gtf.gz"
    echo '```'
    echo "exiting"
    exit 1
fi

REF_HT2="$REF_DIR/genome_index_GRCh38"
if ! test -f "$REF_HT2.1.ht2"; then
    echo "HISAT2 index not found. Please download and build in $REF_DIR as follows"
    echo '```'
    echo "wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    echo "gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    echo "hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa genome_index_GRCh38"
    echo '```'
    echo "exiting"
    exit 1
fi

_bold() {
    printf "\033[1m$1\033[0m\n"
}
_run() {
    _bold "$1 ..."
    if ! eval "$2"; then
        _bold "'$1' exited with non-zero code, stopping"
        exit 1
    fi
    _bold "Success"
}

_run "Downloading SRA" \
    "prefetch $1; cd $1"

_run "Splitting files" \
    "fastq-dump --split-files $1.sra"

_run "Aligning data to reference genome" \
    "hisat2 --threads $N_CORES -x $REF_HT2 -1 $1_1.fastq -2 $1_2.fastq -S output.sam"

_run "Converting to BAM" \
    "samtools view --threads $N_ADD_THREADS -bS output.sam > output.bam"

_run "Sorting" \
    "samtools sort --threads $N_ADD_THREADS output.bam -o output_sorted.bam"

_run "Indexing" \
    "samtools index --threads $N_ADD_THREADS output_sorted.bam"

_run "Assempling transcripts and quantifying splice events" \
    "stringtie -p $N_CORES output_sorted.bam -G $REF_GTF -o output.gtf -A gene_abundances.tsv"
