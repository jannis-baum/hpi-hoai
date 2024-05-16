# Creating testing datasets

This file describes the process and requirements to create datasets with ground
truth RNA splicing distributions.

The idea here is to

1. Use RNA sequencing data from a real experiment that quantifies read counts,
   i.e. how much of the given RNA was read
2. Align the sequencing results to a reference genome to identify what part of
   the genome the RNA was transcribed from
3. Build a splicing distribution

## Installing required tools

- [sra-tools/fastq-dump](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
- [hisat2](http://daehwankimlab.github.io/hisat2/download/). Here we also have
  to build an index:
  ```sh
  wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

  hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa genome_index_GRCh38
  ```
- [samtools](http://www.htslib.org/download/), which we have to build ourselves
  according to the instructions on the download page
- [stringtie](https://github.com/gpertea/stringtie/releases). Here we will also
  need a reference GTF file:
  ```sh
  wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
  gunzip Homo_sapiens.GRCh38.104.gtf.gz
  ```

## Finding suitable datasets

The National Center for Biotechnology (NCBI) has a Sequence Read Archive (SRA)
where we can search for RNA-Seq data with paired layout. In the [search
portal](https://www.ncbi.nlm.nih.gov/sra) we query

```plain
((RNA-Seq[Strategy]) AND PAIRED[Layout]) AND "Homo sapiens"[orgn]
```

And select a suitable RNA sequencing experiment.

## Creating splice site distribution data

- pick an RNA-Seq dataset, `SRRXXX`
- load data as fastq
  ```sh
  prefetch SRRXXX
  fastq-dump --split-files SRRXXX/SRRXXX.sra
  ```
- align data to reference genome
  ```sh
  hisat2 -x genome_index_GRCh38 -1 SRRXXX_1.fastq -2 SRRXXX_2.fastq -S output.sam
  ```
- convert SAM to BAM, sort, and index
  ```sh
  samtools view -bS output.sam > output.bam
  samtools sort output.bam -o output_sorted.bam
  samtools index output_sorted.bam
  ```
- assemble transcripts and quantify splicing events
  ```sh
  stringtie output_sorted.bam -G Homo_sapiens.GRCh38.104.gtf -o output.gtf -A gene_abundances.tsv
  ```
