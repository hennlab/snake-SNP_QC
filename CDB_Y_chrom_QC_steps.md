# Steps used to extract Y chromosome data and QC for CDB snp-array data

### 1. Download and open Genome Studio
### 2. Load New Project
-
Choose the "Genotyping" option to open the Project Wizard
- Name the project "CDB_Y_chrom"
- Choose the "Use sample sheet to load sample intensities" option
- Input the path to three files provided by the genotyping facility
    - Sample sheet (.csv)
    ```
    /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NCTB_H3Africa/UCDavis_Henn_H3Africa_2020.02_SampleSheet.csv
    ```
    - Directory for all your sample iDAT files
    ```
    /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NCTB_H3Africa/iDATs
    ```
    - Directory for the Manifest (.bpm)
    ```
    /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NCTB_H3Africa/Manifest
    ```

- include path to cluster file
```
/share/hennlab/data/snp-array/SAfrica_IlluminaArrays/UCDavis_Henn_H3Africa_2020.02_Final.egt
```
- Use default Gen Call Threshold of 0.15, click "Finish" to load project into GenomeStudio.

### 3. Sample QC, SNP QC, and filtering inside GenomeStudio

- Once samples are loaded, updated SNP statistics
- Update sample statistics by clicking the calculator icon in the Samples Table (bottom left)
- Add gender data into sample table using File > Import Phenotype Information from file > select file with gender info from desktop. This file must have two columns: the first listing the index of the sample in the sample sheet, and the second listing the gender of the sample at that index. The syntax must be "Male", "Female", or "Unknown".  
  - For some reason, GenomeStudio wasn't reading the sample sheet in correctly if the gender column was anything other than "Unknown" when loading the project, but luckily there is a way to add that information into the sample table after the loading in the project.
- Filter the sample table (bottom left) by Males only using the "filter rows" button (Looks like an upside triangle with a line underneath it)
- Filter the SNP table (upper right, second tab) by Y chromosome only using the filter rows button there.
- From the same Samples Table, Plot 'p10 GC' against 'Call Rate' and exclude all samples that do not lie in the major group. We chose to filter by a call rate threshold of 0.9
- Using the filter rows button, apply the following filters:
  - Rep errors <= 2
  - cluster sep >0.02  
  - AB freq <= 0.05

### 4. Write to plink file
- The requires you to install a plug-in available at: http://support.illumina.com/array/array_software/genomestudio/downloads.html
- Using the Report Wizard (Analysis > Reports > Report Wizard... from the top menu), export these to a PLINK formatted report
  - Change ForwardStrand param to true
  - Remove all non-visible SNPs, zeroed SNPs and excluded individuals

### 5. Plink QC outside of GenomeStudio
- Convert to bed/bim/fam format
```
cd /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NCTB_H3Africa/Y_chrom/
  plink --file GenomeStudio_Output/CDB_Y_chrom --make-bed --out CDB_Y_chrom
```
- Remove Duplicates:
  ```
  plink --bfile CDB_Y_chrom --list-duplicate-vars suppress-first --out CDB_Y_chrom
  wc -l CDB_Y_chrom.dupvar # Only 3 duplicates found
  plink --bfile CDB_Y_chrom --exclude CDB_Y_chrom.dupvar --make-bed --out CDB_Y_chrom.noDups
  ```
- Update SNP ids
```
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate SNP-QC
Rscript update_rsID_bim_v2.R --bim CDB_Y_chrom.noDups.bim --rsid /share/hennlab/reference/dbSNP/snp144_snpOnly_noContig.bed --out CDB_Y_chrom.noDups.rsids.bim # script from snake-SNP_QC pipeline on GitHub
cp {params.input}.bed {params.output}.bed
cp {params.input}.fam {params.output}.fam
```

### 6. Snappy
