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
It also requires a ped and map file from genome studio containing the common variants, please name as follows:
- {dataset}.commonvar.ped, {dataset}.commonvar.map

Folder set up:
Please create a separate folder for the snakefile to run in. The folder structure should be as follows:
  - Subdirectory GenomeStudio should contain required output files from GenomeStudio in the correct format
    GenomeStudio/{dataset}.txt
    GenomeStudio/{dataset}.xy.txt
    GenomeStudio/{dataset}.rarevariants
- Subdirectory scripts should contain all scripts needed for pipeline
    scripts/


Running the pipeline
config:
- paths to following files
- name of dataset (needs to match naming scheme of {dataset}.txt )
- strand information file
    /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NC_H3Africa/StrandReport/H3Africa_2017_20021485_A3_StrandReport_F.txt
### Scripts which need to be sped up
scripts/update_rsID_bim.py
