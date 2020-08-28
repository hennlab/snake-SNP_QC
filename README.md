# Henn Lab SNP-QC pipeline Steps

```bash https://sites.google.com/site/thehennlab/snp_array_qc```


### 1. Genome Studio (Austin)

Inputs to snakefile pipeline from zCall (must be named with the following syntax) Please name these two files according to the following syntax:
- {dataset}.txt (file without xy)
- {dataset}.xy.txt (file containing xy)
- {dataset}.rarevariants
- {dataset}.commonvar.ped, {dataset}.commonvar.map


### 2. Set up snakemake pipeline
After Genome Studio is finished, the zCall and other steps are run through the automated snakemake portion of the pipeline. To set up the snakemake, follow the following steps:

#### 2.1 Activate the conda environment

The conda environment must contain
- plink
- python 2.7
- pandas

```bash
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate SNP-QC
```

#### 2.2 Set up working directory

Please create a working directory for the snakefile to run in. The folder structure should be as follows:
```
|__ Snakefile
|__ config.yaml
|__ GenomeStudio
    |__ {dataset}.txt
    |__ {dataset}.xy.txt
    |__ {dataset}.rarevariants
    |__ {dataset}.commonvar.ped
    |__ {dataset}.commonvar.map
|__ scripts
    |__het_filter.R
    |__AB_to_top_allele_v2.R
    |__convertReportToTPED.py
    |__excess_het.R
    |__find_cg_at.py
    |__find_duplicates_bim_v3.R
    |__match_against_1000g.py
    |__update_rsID_bim_v2.R
    |__calibrateZ.py
    |__findBetas.r
    |__zCall.py
    |__findMeanSD.py
    |__dropSamplesFromReport.py
    |__findThresholds.py

```
The last 6 scripts are from the zcall program. They can be symbolically linked from this location on augrabies `/share/hennlab/progs/ZCall/GenomeStudio` into your scripts folder, or downloaded from the [zCall github page](https://github.com/jigold/zCall)

All other scripts can be found in this github.

The config file requires four arguments:
- dataset: the name of the dataset, must match the input genome studio scripts {dataset}.txt, {dataset}.xy.txt, {dataset}.rarevariants
- dbsnp: path to the dbsnp file
- legend: path to the illumina strand legend file

Below is an example format for the config file:

```yaml
dataset: cdb_80

dbsnp: /share/hennlab/reference/dbSNP/snp144_snpOnly_noContig.bed

legend: /share/hennlab/reference/1000g_legend_forQC/combined_autosome_X_XY_1000GP_Phase3_GRCh37_SNPonly.legend
```

### 3. Running the pipeline

```bash  
/share/hennlab/progs/miniconda3/bin/snakemake --config config.yaml -p -j 10
```

#### 3.1 Detailed description of each step in the snakefile portion of the pipeline



### 4. Manual steps still required after running snakefile
