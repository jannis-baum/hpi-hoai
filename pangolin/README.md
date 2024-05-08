# Pangolin evaluation

## Installation

- install requirements.txt in a virtual env with Python version specified in
  `.python-version`
  - you may need to `pip install setuptools==57` separately to make it work
- run `pip install .` in the `Pangolin` repo directory

## Running the example

- download the data
  ```sh
  mkdir data; cd data
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
  wget https://www.dropbox.com/sh/6zo0aegoalvgd9f/AAA9Q90Pi1UqSzX99R_NM803a/gencode.v38lift37.annotation.db
  ```
- run pangolin
  ```
  pangolin Pangolin/examples/brca.vcf data/GRCh37.primary_assembly.genome.fa.gz data/gencode.v38lift37.annotation.db brca.pangolin
  ```
