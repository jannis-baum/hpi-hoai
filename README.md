# Evaluating Predictions of Splice Site Usage Probabilities

This repository is part of my work for the course [Hands-on Artificial
Intelligence for Digital
Health](https://hpi.de/digital-health-cluster/teaching/summer-term-2024/hands-on-artificial-intelligence-for-digital-health.html)
at HPI.
It contains the implementation of a pipeline to compute per-base splice site
usage probabilities from genetic variants with
[SpliceAI](https://github.com/Illumina/SpliceAI) and
[Pangolin](https://github.com/tkzeng/Pangolin), and to evaluate and compare them
to RNA-seq results based on visualizations and a suitable error metric.
In the context of this work, I executed the pipeline at the example of 100
breast cancer patients from the [TCGA-BRCA
study](https://www.cancerimagingarchive.net/collection/tcga-brca/) and the BRCA1
gene.

## Repository overview

- [`evaluation.ipynb`](evaluation.ipynb) is the main Notebook and produces the
  core of what is discussed in the report of my work, as well as all figures.
- [`tissue-selection.ipynb`](tissue-selection.ipynb) represents an appendix to
  the report for selecting a suitable tissue-specific Pangolin model for a given
  use case, in this case breast tissue and the BRCA1 gene.
- [`filter_files.py`](filter_files.py) and [`get_manifest.py`](get_manifest.py)
  are used to build a data repository from the [GDC Data
  Portal](https://portal.gdc.cancer.gov). More details can be found
  [below](#data).
- The [`script/`](script/) directory holds code supporting the Notebooks, e.g. a
  generalized interface to SpliceAI, Pangolin and baseline predictors in
  [`predictions.py`](script/predictions.py), reconstruction of genomes from
  variants in [`sequences.py`](script/sequences.py), and data loading and
  caching in [`data_loading.py`](script/data_loading.py)

## Setup

To use this repository to reproduce my results or investigate new data, follow
the steps below.

### Python

Create a new Python 3.10.14 environment. Then

```sh
pip install pip==23 setuptools==57
pip install -r requirements.txt
pyensembl install --release 111 --species homo_sapiens
```

> **Note:** The specific `pip` and `setuptools` versions are necessary to fix a
> [deprecation problem](https://github.com/jamescasbon/PyVCF/issues/334) of
> [PyVCF](https://github.com/jamescasbon/PyVCF), a dependency of Pangolin.

### Data

To set up the data required to run the pipeline in
[`evaluation.ipynb`](evaluation.ipynb) for your study, first create a `data/`
directory at the root of the repository.

```sh
mkdir data/
```

Next, download a reference genome sequence suitable for your study. In this
example, we will use the [reference sequence the GDC
pipeline](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
uses. You can download it from
[here](https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834).
Place it in the `data/` directory and name it `ch38.fa`

Now we need to get matched sequencing data in VCF format and RNA sequencing
results in the form of STAR2 splice junction counts for each individual
considered in the study. The scripts [`filter_files.py`](filter_files.py) and
[`get_manifest.py`](get_manifest.py) are made to support retrieving this data
from the [GDC Data Portal](https://portal.gdc.cancer.gov) by following the steps
below:

1. Go to the [Cohort
   Builder](https://portal.gdc.cancer.gov/analysis_page?app=CohortBuilder) and
   select the data source(s) you want to use. For this example I selected
   *TCGA-BRCA1* under *Project*, and filtered for cases where *Raw Simple
   Somatic Mutations* from *WGS* are available under *Available Data*. The
   filtering here is not necessary but significantly speeds up the next step.
2. In the website's navigation bar, click *Repository* and then *JSON* at the
   top-left of the table to download a list of the selected files. Save the JSON
   at `data/files.json`.

   > **Note:** The maximal download speed of this list seems to drop
   > significantly with the number of selected files so this may take a while.

3. Click *Add All Files to Cart* at the top-right of the table. Navigate to
   *Cart* at the top-right of the website, then open the drop-down menu
   *Download Cart* at the top-left, and select *Manifest*. Save the downloaded
   file at `data/manifest.txt`.
4. Run `filter_files.py` to filter the downloaded list for relevant files and
   create a data index for your study:

   ```sh
   python filter_files.py data/files.json --out data/info.csv
   ```

5. Use `get_manifest.py` to filter the manifest file for the relevant files:

   ```sh
   python get_manifest.py data/manifest.txt data/info.csv --out data/manifest.txt
   ```

6. Use the GDC Data Transfer Tool to download all files specified in the
   `data/manifest.txt` file. You can save them anywhere in the `data/`
   directory; the subdirectory structure is not relevant.

Now you have all data to run the implemented pipeline for your study!
