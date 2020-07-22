"""
Title: snake-SNP_QC
Description: Henn lab pipeline for QC of SNP array data
Authors:
"""

# -------
# SET-UP
# -------

DATASET = config['dataset']

Z = [i for i in range(3,11)]

def get_Z_num(wildcards):
    files = expand("{dataset}_outputThreshold_{z}.txt", dataset=DATASET, z=Z)
    return files

# -------
# ZCALL
# -------

rule convert_tped:
  input:
    "{dataset}.txt"
  output:
    "{dataset}.tped",
    "{dataset}.tfam"
  params:
    prefix = DATASET
  shell:
    """
    convertReportToTPED.py -R {input} -O {params.prefix}
    """

rule remove_het:
  input:
    "{dataset}.tped",
    "{dataset}.tfam"
  output:
    "{dataset}.het"
  params:
    prefix = DATASET
  shell:
    """
    plink --tfile {params} --het --out {params}
    """

rule dropSamples:
  input:
    het = "{dataset}.het",
    xy_report = "{dataset}.xy.txt"
  output:
    drop = "{dataset}.SamplesToDrop.txt",
    new_report = "{dataset}.xy.drop.txt"
  shell:
    """
    Rscript scripts/excess_het.R {input.het} {output.drop}
    python dropSamplesFromReport.py {input.xy_report} {output.drop} > {output.new_report}
    """

rule findMeanSD:
  input:
    "{dataset}.xy.drop.txt"
  output:
    "{dataset}.meanSD.txt"
  shell:
    """
    python scripts/findMeanSD.py -R {input} > {output}
    # python scripts/findMeanSD.py -R zCall_correctorder.xy.drop.txt > data_for_zCall_meanSD.txt
    """

rule findBetas:
  input:
    "{dataset}.meanSD.txt"
  output:
    "{dataset}.betas.txt"
  shell:
    """
    Rscript scripts/findBetas.r {input} {output} 1
    # Rscript scripts/findBetas.r data_for_zCall_meanSD.txt data_for_zCall.betas.txt 1
    """

rule findThresholds:
  input:
    betas = "{dataset}.betas.txt",
    report = "{dataset}.xy.drop.txt",
  output:
    "{dataset}.thresholds.{Z}.txt"
  params: DATASET
  shell:
    """
    for Z in `seq 3 10` ; do python scripts/findThresholds.py -B {input.betas} -R {input.report} -Z $Z -I 0.2 > {output} ; done
    """

# for Z in `seq 3 10` ; do python scripts/findThresholds.py -B data_for_zCall.betas.txt -R zCall_correctorder.xy.drop.txt -Z $Z -I 0.2 > zCall.thresholds.${Z}.txt ; done

rule calibrateZ:
  input:
    report = "{dataset}.xy.drop.txt",
    thresh = "{dataset}.thresholds.{Z}.txt"
  output:
    "{dataset}.concordance.stats{Z}.txt",
  shell:
    """
    python scripts/calibrateZ.py -R {input.report} -T {input.thresh} > {output}
    """
# for Z in `seq 3 10`; do python scripts/calibrateZ.py -R zCall_correctorder.xy.drop.txt -T zCall.thresholds.${Z}.txt > zCall.xy.concordance.stats.${Z}.txt ; done

# this rule needs to go into a shell script ?
rule global_concordance:
  input:
    report = "{dataset}.xy.drop.txt",
    conc_z = "{dataset}.concordance.stats{Z}.txt"
    thresh = "{dataset}.thresholds.{Z}.txt"
  output:
    multiext("{dataset}.zCalled", ".tped", ".tfam")
  params:
    stats = "{dataset}.xy.concordance.stats",
    out = "{dataset}.zCalled"
  shell:
    """
    BESTZ = grep -oP '(?<=Global Concordance: ).*' *concordance.stats* | sort -sn -t ":" -k 2,2 | tail -n 1 | awk -F ":" '{{print$1}}' | sed -e s/[^0-9]//g
    python scripts/zCall.py -R {input.report} -T {params}.${BESTZ}.txt -O {params.out}
    """
# python scripts/zCall.py -R zCall_correctorder.xy.drop.txt -T zCall.thresholds.6.txt -O Austin_zCalled

rule extract_rare:
  input:
    t_in = multiext("{dataset}.zCalled", ".tped", ".tfam")
    rare = "{dataset}.rarevariants"
  output:
    multiext("{dataset}.zCalled.rare-filt", ".ped", ".map")
  params:
    input = "{dataset}.zCalled"
    output = "{dataset}.zCalled.rare-filt"
  shell:
    """
    plink --tfile {params.input} --extract {input.rare} --geno 0.85 --not-chr 0 --update-alleles --recode --out {params.output}
    """
# plink --tfile Austin_zCalled --extract /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NC_H3Africa/CDB/Output/rare_variants.txt --geno 0.85 --not-chr 0 --update-alleles --recode --out zCalled.rare-filt


rule more_QC:
  input:
    multiext("{dataset}.zCalled.rare-filt", ".ped", ".map")
  output:
    plink = "{dataset}.zCalled.rare-filt.plink.hwe",
    badsnp = multiext("{dataset}.zCalled.rare-filt.badsnp-filt", ".bed", ".fam", ".bim")
  params:
    infile = "{dataset}.zCalled.rare-filt",
    hwe = "{dataset}.zCalled.rare-filt.plink",
    badsnp = "{dataset}.zCalled.rare-filt.badsnp-filt"
  shell:
    """
    plink --file {params.infile} --hardy --out {params.hwe}
    Rscript scripts/het_filter.R {output.plink} over 0.8
    plink --file {params.infile} --exclude het_filter_out.txt --make-bed --out {params.badsnp}
    """

# plink --file zCalled.rare-filt --hardy --out Zcalled.rare-filt.plink.hwe
# Rscript scripts/het_filter.R plink.hwe over 0.8
# plink --file zCalled.rare-filt --exclude het_filter_out.txt --make-bed --out zCalled.rare-filt.badsnp-filt

rule AB_convert_rare:
# Rscript scripts/AB_to_top_allele_v2.R --bim zCalled.rare-filt.badsnp-filt.bim --strand /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NC_H3Africa/StrandReport/H3Africa_2017_20021485_A3_StrandReport_F.txt --out zCalled.rare-filt.badsnp-filt.conv.bim

rule merge_variants:
  input:
    com = multiext("{dataset}.commonvar", ".ped", ".map"),
    rare = multiext("{dataset}.zCalled.rare-filt.badsnp-filt", ".bed", ".fam", ".bim")
  output:
    multiext("{dataset}.commonvar", ".bed", ".bim", ".fam"),
    multiext("{dataset}.rareCommonMerged", ".bed", ".bim", ".fam")
  params:
    com = "{dataset}.commonvar",
    rare = "{dataset}.zCalled.rare-filt.badsnp-filt",
    merge = "{dataset}.rareCommonMerged"
  shell:
    """
    plink --file {params.com} --make-bed --out {params.com}
    plink --bfile {params.com} --bmerge {params.rare} --make-bed --out {params.merge}
    """

 # plink --file cdb_80_commonvar --make-bed --out cdb_80_commonvar
#  plink --bfile cdb_80_commonvar --bmerge zCalled.rare-filt.badsnp-filt --make-bed --out cdb_80.rareCommonMerged

# merging with austins new bim file
# plink --bfile cdb_80_commonvar --bmerge zCalled.rare-filt.badsnp-filt.conv --make-bed --out cdb_80.rareCommonMerged


# this step would be skipped
rule convert_strand:
  input:
    multiext("{dataset}.rareCommonMerged", ".bed", ".bim", ".fam"),
  output:
    "{dataset}.rareCommonMerged.strand.bim"
  shell:
    """
    python scripts/AB_to_top_allele.py --bim cdb_80.rareCommonMerged.bim --strand /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NC_H3Africa/StrandReport/H3Africa_2017_20021485_A3_StrandReport_F.txt --pos y --out cdb_80.rareCommonMerged.strand.bim
    """
# python scripts/AB_to_top_allele_.py --bim cdb_80.rareCommonMerged.bim --strand /share/hennlab/data/snp-array/SAfrica_IlluminaArrays/CDB_NC_H3Africa/StrandReport/H3Africa_2017_20021485_A3_StrandReport_F.txt --pos y --out cdb_80.rareCommonMerged.strand.bim


# /share/hennlab/reference/dbSNP/snp144_snpOnly_noContig.bed

rule update_snp_ids:
  input:
    "{dataset}.rareCommonMerged.bim"
  output:
    "{dataset}.rareCommonMerged.dbsnp.bim"

# python scripts/update_rsID_bim.py --bim cdb_80.rareCommonMerged.bim --bed /share/hennlab/reference/dbSNP/snp144_snpOnly_noContig.bed  --format T--codechr T --out cdb_80.rareCommonMerged.dbsnp.bim


rule find_dups_bim:
  input:
  output:
  shell:

# Rscript scripts/find_duplicates_bim_v2.R cdb_80.rareCommonMerged.dbsnp.bim
# cut -f2 NonDup_cdb_80.rareCommonMerged.dbsnp.bim > NonDup_SNPS.txt
# cp cdb_80.rareCommonMerged.bed cdb_80.rareCommonMerged.dbsnp.bed
# cp cdb_80.rareCommonMerged.fam cdb_80.rareCommonMerged.dbsnp.fam
# plink --bfile cdb_80.rareCommonMerged.dbsnp --extract NonDup_SNPS.txt --make-bed --out cdb_80.rareCommonMerged.dbsnp.noDup

rule align_strand:
  input
  output
  shell:

# python match_against_1000g.py --bim cdb_80.rareCommonMerged.dbsnp.noDup --legend /share/hennlab/reference/1000g_legend_forQC/combined_autosome_X_XY_1000GP_Phase3_GRCh37_SNPonly.legend --out MidFixOfOutput
