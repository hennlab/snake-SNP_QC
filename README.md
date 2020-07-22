# snake-SNP_QC
Henn Lab SNP Array QC pipeline from GenomeStudio Output
```https://sites.google.com/site/thehennlab/snp_array_qc```

# How to run pipeline:

### 1. Set up

Conda environment:
- plink
- python 2.7
- pandas 

Activating conda environment
```bash
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate SNP-QC
```

You must export path every time you run the pipeline
```bash
export PATH=/share/hennlab/progs/ZCall/GenomeStudio/:$PATH
```

### 2. Inputs to pipeline

The zcall portion of the pipeline requires two versions of the genome studio report output file: One containing the X and Y columns, and one without them. Please name these two files according to the following syntax:
- {dataset}.txt (file without xy)
- {dataset}.xy.txt (file containing xy)
It also requires a file from genome studio containing the rare SNPs, please name as follows:
- {dataset}.rarevariants
