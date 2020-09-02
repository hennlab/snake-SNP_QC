# Henn Lab SNP-QC pipeline Steps

```https://sites.google.com/site/thehennlab/snp_array_qc```


### 1. GenomeStudio 2.0
This is a free Windows-only GUI program from Illumina that you will have to create an account to download. 
    *Input:* 
            -.idat files for every sample
            -.bpm file for the project, 
            -a sample sheet for the project, 
            -and possibly a .egt (cluster) file, depending on the array used.
    *Output:* 
            -a plot of Call Rate vs. p10 GC, 
            -a plink file of common variants, 
            -a list of SNPs with poor cluster separation or excessive replicate errors, 
            -a list of rare SNPs that need to be replaced using zCall, 
            -and the Full Data Table for input to zCall

#### 1.1 Make new project
- choose the "Genotyping" option to open the Project Wizard
- name the project and put it wherever you want on your system
- choose the "Use sample sheet to load sample intensities" option
- input the path to three files provided by the genotyping facility
    - Sample sheet (.csv)
    - Directory for all your iDAT files 
    - Directory for the Manifest (.bpm)
- if you were provided with a cluster file (.egt) you can include the path to that
- Under the "Project Creation Actions" section:
    - if your system does not have much RAM, select 'Cluster SNPs' option
    - Use default Gen Call Threshold of 0.15
- click "Finish" and your project will be loaded into GenomeStudio. This could take several minutes to an hour depending on the number of samples and the speed of your system.


#### 1.2 Exclude outlier samples
        - Update sample statistics by clicking the calculator icon in the Samples Table (bottom left)
 - Plot 'p10 GC' against 'Call Rate' and exclude all samples that do not lie in the major group
 Ex. choose a Call Rate threshold < 0.9
        - Save an image of the plot for future reference
        - Update SNP statistics for all SNPs when prompted

        
#### 1.3 Apply some QC filters
*These are some standard QC metrics we use for genotype data. Depending on the quality of your data, you can adjust things accordingly*
 - Call freq < 0.85 (will generally exclude all of chr Y)
 - Rep errors > 2
 - Cluster separation values < 0.02
 - Heterozygote rate >= 0.80
        The line in the filter window should look similar to this: [ !("Call Freq" < 0.85 ) AND !("Rep Errors" > 2 ) AND !("<all columns>.Cluster Sep" < 0.02 ) AND !("AB Freq >= 0.80 ) ]

      
#### 1.4 Write out common variants as a PLINK file - these are done!
        - Add a minor allele frequency filter to capture the common variants
        One liner: [ !("Call Freq" < 0.85 ) AND !("Rep Errors" > 2 ) AND !("<all columns>.Cluster Sep" < 0.02 ) AND !("AB Freq >= 0.80 ) AND ("Minor Freq" >= 0.05 ) ]
            - Using the Report Wizard, export these to a PLINK formatted report
            - Plug-in available at: http://support.illumina.com/array/array_software/genomestudio/downloads.html
	    - Change ForwardStrand param to true, remove all non-visible SNPs and excluded individuals


#### 1.5 Identify SNPs with poor cluster separation or excessive replicate errors
        - Clear the previous filter, and apply the Rep Errors > 2 and Cluster Sep < 0.02 filters only, both branching from an OR statement
        One liner: [ [ ("Rep Errors" > 2 ) OR ("<all.columns>.Cluster Rep" < 0.02 ) ] ]
        - Display only the 'Name' sub column, and export the resulting table to a file to be used later.


#### 1.6 Identify rare SNPs to be replaced by zCall:
        - Clear the previous filter, and apply the Minor allele frequency < 0.05 filter
        One liner: [ ("Minor Freq" < 0.05 ) ]
        - This will display a list of SNPs whose calls should be replaced by the results of zCall. Export the resulting table to a file to be used later.


#### 1.7 Export data for zCall
        - Clear all filters
 - Click 'Full Data Table' tab (to the left of the 'SNP Table' tab)
 - Click 'Column Chooser' icon
 - Display 'Name', 'Chr', 'Position', and all sample columns, and 'GType', 'X', 'Y' subcolumns
 - Click OK and click 'Export displayed data to file' icon. Do not include individuals who were outliers on the Call Rate vs. p10 GC plot.

 
#### 1.8 Save project, can exit





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

The config file requires three arguments:
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
/share/hennlab/progs/miniconda3/bin/snakemake --configfile config.yaml -p -j 10
```

#### 3.1 Detailed description of each step in the snakefile portion of the pipeline

- rule **convert_tped**: convert genome studio output (without xy) to tped/tfam format

- rule **remove_het**: identify heterozygous samples with plink. plink needs the genome studio format without the xy columns
- rule **dropSamples**: identify samples with excess heterozygosity and remove them from the report. The script ```dropSamplesFromReport.py``` needs the genome studio report format with the xy columns
- rule **findMeanSD**: Calculate means of homozygote clusters
- rule **findBetas**: Calculate beta and p-values from the mean and SD using weighted linear regression. If this doesn't work, make sure you have removed individuals with excessive No Calls. Use the 'dropSamplesFromReport.py' script to remove them, then try again
- rule **findThresholds**: Find the thresholds t_x and t_y for each SNP using findThresholds.py, running multiple times with different values of z (3-10)
- rule **calibrateZ**: Run the calibration script to calculate threshold stats
- rule **global_concordance**: Recall SNPs using the best threshold - pick the one with highest global concordance
- rule **extract_rare**: Extract the rare variants from the zCalled plink files - Extract rare SNPs from .tped file and apply the genotype filter
- rule **more_QC**: Apply additional QC filters: Remove SNPs that have suspiciously high heterozygosity, poor cluster separation, or excessive replicate errors, append the output to the list of SNPs with poor cluster separation or too many replicate errors, delete the top line, and remove these from the plink file,
- rule **AB_convert_rare**: Convert Illumina A and B allele to Illumina Forward allele.
    - Note: There are still some loci assigned “NA” as allele code, because their SNP names are not found in the strand files. It’s technically possible to research them again using the original chr and pos, yet probably not worth it as :1) they could be indels at the same locus with a SNP, updating alleles according to chr:pos is not guaranteed to retrieve the right genotypes; 2) they could simply be duplicated loci assigned with a different SNP name. I suggest removing them.
        - Note: The current version of this script assumes Illumina ForwardAllele strand reports and not TopBottom strand reports.
- rule **make_commonvar**: convert common variant file from ped/map to bed/bim/fam.
- rule **remove_dups**: Using plink remove duplicates from common variants file and rare variants file
- rule **merge_variants**: Merge common and filtered rare variant files
- rule **update_snp_ids**: Update SNP names to dbsnp
        - Note: Often a chr:position will have more than one rsid. We allow the code to remove these duplicates randomly, so there is no guarantee that a particular rsid for any chr/pos will be included. If you want specific rsids to be included you will need to modify this behavior.
- rule **find_dups_bim**: Some loci have the same SNP names and are duplicated, this step removes them
    - Note: Some duplicates removed in this step will be different versions of the same variant (eg. alleles are AT and CG or 0A and TA in the duplicates). This program removes duplicates without regard to what allele they are. , which could result in missing data for some individuals at the affected snps. This could be an issue in smaller datasets or when you are interested in one of the affected snps.
- rule **align_strand**: Align strand to 1000 genome reference. T
    - The script has three outputs:
        1) Indel_[outfile].txt: a bim file of indels
        2) NonMatching_[outfile].txt: a bim file containing loci not found among 1000 genome (phase 3) variants (matched by positions), or those with different coding alleles than 1000 genome (tri-allelic, for example). To merge smoothly with other platforms, these loci are to be removed.
        3) FlipStrand_[outfile].txt: a bim file containing loci to flip. .bim format is directly recognizable as to follow --exclude or --keep in plink. No need to cut the second column out.
- rule **flip**: Flip strands
- rule **remove_non_matching**: remove non-matching loci
- rule **remove_cg_at**: Remove A/T, G/C loci using plink


### 4. Manual steps still required after running snakefile


- Update individuals' IDs from Illumina plate # to IDs in ethnographic data
Also check duplicated IDs.

- Gender check based on X heterozygosity; heterozygosity check for contamination
plink
Remove individuals you think that are contaminated, swapped etc.

- Geno (and MAF filter)
Note: sometimes fam file might have NAs in the parent column. This will set the corresponding individual to a non-founder, and the individuals be ignored automatically by plink in --freq or --maf commands. Either manually correct them, or add --nonfounders when running plink.
