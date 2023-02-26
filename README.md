# PlasmidCompiler
This project simplifies the process of building plasmids and packaging them up into a single annoated genbank file that can be easily used on benchling or another workbench.  Although this task is commonly perfomred by hand in a browser or desktop application, performing senstive tasks like this with python allows us to be eaxct with our changes and to verify the changes.

The project aims to be a swiss army kinfe of common pipelines used in gene editing. The current pipeline consists of three main steps: (1) downloading gene data from NCBI, (2) discovering mutations within the genome using GATK4, and (3) predicting off-targets for a set of gRNA candidates using DeepCRISPR's DCModelOntar model.

## Prerequisites
- Python 3.7+
- Conda
- GATK4
- tenserflow 1.3.2 (todo: tf2 update)

## Installation
Clone the repository

`git clone https://github.com/SweetWaterAI/PlasmidCompiler`

Install required packages using the command below

`pip install -r requirements.txt`

Add GATK4 to the path

`export PATH=$PATH:/path/to/gatk4`


### Example Usage

download
This function downloads gene data for a given gene name from NCBI and saves all fields as annotations in a genbank file. If the gene data has already been downloaded, the function exits without re-downloading the data.

Usage:

`$ python main.py download <gene_name> [--email <email>]`

#### Arguments:

- gene_name (required): Name of the gene to download.
- email (optional): Email address to be used for NCBI queries.
### discover
This function uses GATK4 to find all mutations of a given gene within a genome. The function requires a reference genome file, the name of the gene of interest, and a path to the output VCF file.

Usage:
`$ python main.py discover <reference_file> <gene_name> <output_file>`

####Arguments:

- reference_file (required): Path to the reference genome file.
- gene_name (required): Name of the gene of interest.
- output_file (required): Path to the output VCF file.


### predict
This function predicts off-targets for a set of gRNA candidates and target sequences using DeepCRISPR's DCModelOntar model. The function requires the path to the gRNA candidates file, the path to the target sequences file, the path to the full genome sequence file, and the path to the DCModelOntar model file.

Usage:

`$ python main.py predict <gsRNA_candidates> <targets> <genome> [--output <output>] --model_path <model_path> [--pam <pam>]`

#### Arguments:

- gsRNA_candidates (required): Path to a file containing gRNA sequences.
- targets (required): Path to a file containing target sequences.
- genome (required): Path to the full genome sequence.
- output (optional): Path to the output file. Default is "ontar_results.tsv".
model_path (required): Path to the DCModelOntar model file.
- pam (optional): PAM sequence used by the gRNA. Default is "NGG".
### Reproduceablity / Flexablity
Commands can be saved within a .sh file and re-run to re-build a specific plasmid when needed. A team should be able to build any number of therapies or tests out of the same directory.  The code will take in fasta or genebank files.
