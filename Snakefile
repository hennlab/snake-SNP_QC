"""
Title: snake-SNP_QC
Description: Henn lab pipeline for QC of SNP array data
Authors:
"""

# -------
# SET-UP
# -------

DATASET = config['dataset']
STRAND = config['strand']
DBSNP = config['dbsnp']
LEGEND = config['legend']
Z = [i for i in range(3,11)]

def get_Z_num(wildcards):
    files = expand("{dataset}_outputThreshold_{z}.txt", dataset=DATASET, z=Z)
    return files

# -------
# ZCALL
# -------

rule convert_tped:
  input:
    "GenomeStudio/{dataset}.txt"
  output:
    "zCall/{dataset}.tped",
    "zCall/{dataset}.tfam"
  params:
    "zCall/"+DATASET
  shell:
    """
    scripts/convertReportToTPED.py -R {input} -O {params.prefix}
    """

rule remove_het:
  input:
    "zCall/{dataset}.tped",
    "zCall/{dataset}.tfam"
  output:
    "zCall/{dataset}.het"
  params:
    "zCall/"+DATASET
  shell:
    """
    plink --tfile {params} --het --out {params}
    """

rule dropSamples:
  input:
    het = "zCall/{dataset}.het",
    xy_report = "GenomeStudio/{dataset}.xy.txt"
  output:
    drop = "zCall/{dataset}.SamplesToDrop.txt",
    new_report = "zCall/{dataset}.xy.drop.txt"
  shell:
    """
    Rscript scripts/excess_het.R {input.het} {output.drop}
    python scripts/dropSamplesFromReport.py {input.xy_report} {output.drop} > {output.new_report}
    """

rule findMeanSD:
  input:
    "zCall/{dataset}.xy.drop.txt"
  output:
    "zCall/{dataset}.meanSD.txt"
  shell:
    """
    python scripts/findMeanSD.py -R {input} > {output}
    """

rule findBetas:
  input:
    "zCall/{dataset}.meanSD.txt"
  output:
    "zCall/{dataset}.betas.txt"
  shell:
    """
    Rscript scripts/findBetas.r {input} {output} 1
    """

rule findThresholds:
  input:
    betas = "zCall/{dataset}.betas.txt",
    report = "zCall/{dataset}.xy.drop.txt",
  output:
    "zCall/{dataset}.thresholds.{Z}.txt"
  shell:
    """
    python scripts/findThresholds.py -B {input.betas} -R {input.report} -Z {Z} -I 0.2 > {output}
    """

rule calibrateZ:
  input:
    report = "zCall/{dataset}.xy.drop.txt",
    thresh = "zCall/{dataset}.thresholds.{Z}.txt"
  output:
    "zCall/{dataset}.concordance.stats{Z}.txt",
  shell:
    """
    python scripts/calibrateZ.py -R {input.report} -T {input.thresh} > {output}
    """

# this rule needs to go into a shell script ?
rule global_concordance:
  input:
    report = "zCall/{dataset}.xy.drop.txt",
    conc_z = "zCall/{dataset}.concordance.stats{Z}.txt"
    thresh = "zCall/{dataset}.thresholds.{Z}.txt"
  output:
    multiext("zCall/{dataset}.zCalled", ".tped", ".tfam")
  params:
    stats = "zCall/{dataset}.xy.concordance.stats",
    out = "zCall/{dataset}.zCalled"
  shell:
    """
    BESTZ = grep -oP '(?<=Global Concordance: ).*' *concordance.stats* | sort -sn -t ":" -k 2,2 | tail -n 1 | awk -F ":" '{{print$1}}' | sed -e s/[^0-9]//g
    python scripts/zCall.py -R {input.report} -T {params}.${BESTZ}.txt -O {params.out}
    """
# python scripts/zCall.py -R zCall_correctorder.xy.drop.txt -T zCall.thresholds.6.txt -O Austin_zCalled

rule extract_rare:
  input:
    t_in = multiext("zCall/{dataset}.zCalled", ".tped", ".tfam")
    rare = "GenomeStudio/{dataset}.rarevariants"
  output:
    multiext("merging/{dataset}.zCalled.rare-filt", ".ped", ".map")
  params:
    input = "zCall/{dataset}.zCalled"
    output = "merging/{dataset}.zCalled.rare-filt"
  shell:
    """
    plink --tfile {params.input} --extract {input.rare} --geno 0.85 --not-chr 0 --update-alleles --recode --out {params.output}
    """

rule more_QC:
  input:
    multiext("merging/{dataset}.zCalled.rare-filt", ".ped", ".map")
  output:
    plink = "merging/{dataset}.zCalled.rare-filt.plink.hwe",
    badsnp = multiext("merging/{dataset}.zCalled.rare-filt.badsnp-filt", ".bed", ".fam", ".bim")
  params:
    infile = "merging/{dataset}.zCalled.rare-filt",
    hwe = "merging/{dataset}.zCalled.rare-filt.plink",
    badsnp = "merging/{dataset}.zCalled.rare-filt.badsnp-filt"
  shell:
    """
    plink --file {params.infile} --hardy --out {params.hwe}
    Rscript scripts/het_filter.R {output.plink} over 0.8
    plink --file {params.infile} --exclude het_filter_out.txt --make-bed --out {params.badsnp}
    """

rule AB_convert_rare:
  input:
    bim = "merging/{dataset}.zCalled.rare-filt.badsnp-filt.bim",
    bed = "merging/{dataset}.zCalled.rare-filt.badsnp-filt.bed",
    fam = "merging/{dataset}.zCalled.rare-filt.badsnp-filt.fam"
  output:
    bim = "merging/{dataset}.zCalled.rare-filt.badsnp-filt.conv.bim",
    bed = "merging/{dataset}.zCalled.rare-filt.badsnp-filt.conv.bed",
    fam = "merging/{dataset}.zCalled.rare-filt.badsnp-filt.conv.fam"
  shell:
    """
    Rscript scripts/AB_to_top_allele_v2.R --bim {input.bim} --strand {STRAND} --out {output}
    cp {input.bed} {output.bed}
    cp {input.fam} {output.fam}
    """

rule make_commonvar:
  input:
    com = multiext("GenomeStudio/{dataset}.commonvar", ".ped", ".map"),
  output:
    multiext("merging/{dataset}.commonvar", ".bed", ".bim", ".fam"),
  params:
    "merging/{dataset}.commonvar"
  shell:
    """
    plink --file {params} --make-bed --out {params}
    awk '{{$1 = $2; print}}' {params}.fam > test ; mv test {params}.fam
    """

rule remove_dups:
  input:
    com = multiext("merging/{dataset}.commonvar", ".bed", ".bim", ".fam"),
    rare = multiext("merging/{dataset}.zCalled.rare-filt.badsnp-filt.conv.bim")
  output:
    multiext("merging/{dataset}.commonvar.noDups", ".bed", ".bim", ".fam"),
    multiext("merging/{dataset}.zCalled.rare-filt.badsnp-filt.conv.noDups", ".bed", ".bim", ".fam")
  params:
    com_in = "merging/{dataset}.commonvar",
    rare_in = "merging/{dataset}.zCalled.rare-filt.badsnp-filt.conv",
    com_out = "merging/{dataset}.commonvar.noDups",
    rare_out = "merging/{dataset}.zCalled.rare-filt.badsnp-filt.conv.noDups"
  shell:
    """
    plink --bfile {params.com_in} --list-duplicate-vars suppress-first --out merging/common
    plink --bfile {params.rare_in} --list-duplicate-vars suppress-first --out merging/rare
    plink --bfile {params.com_in} --exclude merging/common.dupvar --make-bed --out {params.com_out}
    plink --bfile {params.rare_in} --exclude merging/rare.dupvar --make-bed --out {params.rare_out}
    """


rule merge_variants:
  input:
    multiext("merging/{dataset}.commonvar.noDups", ".bed", ".bim", ".fam"),
    multiext("merging/{dataset}.zCalled.rare-filt.badsnp-filt.conv.noDups", ".bed", ".bim", ".fam")
  output:
    multiext("merging/{dataset}.rareCommonMerged", ".bed", ".bim", ".fam")
  params:
    com = "merging/{dataset}.commonvar.noDups",
    rare = "merging/{dataset}.zCalled.rare-filt.badsnp-filt.conv.noDups",
    merge = "merging/{dataset}.rareCommonMerged"
  shell:
    """
    plink --bfile {params.com} --bmerge {params.rare} --merge-mode 2 --make-bed --out {params.merge}
    """

rule update_snp_ids:
  input:
    multiext("merging/{dataset}.rareCommonMerged", ".bed", ".bim", ".fam")
  output:
    multiext("QC/{dataset}.rareCommonMerged.dbsnp", ".bed", ".bim", ".fam")
  params:
    "merging/{dataset}.rareCommonMerged"
  shell:
    """
    python scripts/update_rsID_bim.py --bim {params}.bim --bed {DBSNP} --format T--codechr T --out {params}.dbsnp.bim
    cp {params}.bed {params}.dbsnp.bed
    cp {params}.fam {params}.dbsnp.fam
    """

rule find_dups_bim:
  input:
    multiext("QC/{dataset}.rareCommonMerged.dbsnp", ".bed", ".bim", ".fam")
  output:
    multiext("QC/{dataset}.rareCommonMerged.dbsnp.noDup", ".bed", ".bim", ".fam"),
    NonDup = "NonDup_{dataset}.rareCommonMerged.dbsnp.bim",
    txt = "QC/NonDup_SNPS.txt"
  params:
    input = "QC/{dataset}.rareCommonMerged.dbsnp",
    out = "QC/{dataset}.rareCommonMerged.dbsnp.noDup"
  shell:
    """
    Rscript scripts/find_duplicates_bim_v2.R {params}.bim
    cut -f2 {output.NonDup} > {output.txt}
    plink --bfile {params.input} --extract {output.txt} --make-bed --out {params.out}
    """

rule align_strand:
  input:
    "QC/{dataset}.rareCommonMerged.dbsnp.noDup.bim"
  output:
    "QC/Indel_{dataset}.txt",
    "QC/NonMatching_{dataset}.txt",
    "QC/FlipStrand_{dataset}.txt"
  shell:
    """
    python scripts/match_against_1000g.py --bim {input} --legend {LEGEND} --out {DATASET}
    """

rule flip:
  input:
    flip = "QC/FlipStrand_{dataset}.txt",
    multiext("QC/{dataset}.rareCommonMerged.dbsnp.noDup", ".bim", ".bed", ".fam")
  output:
    snps = "QC/SNPs_to_flip.txt",
    multiext("QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip", ".bim", ".bed", ".fam")
  params:
    input = "QC/{dataset}.rareCommonMerged.dbsnp.noDup"
    out = "QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip"
  shell:
    """
    cut -f2 {input.flip} > {output.snps}
    plink --bfile {params.input} --flip {output.snps} --make-bed --out {params.output}
    mv {params.input}.bed {params.out}.bed
    mv {params.input}.fam {params.out}.fam
    """

rule remove_non_matching:
  input:
    bfile = multiext("QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip", ".bim", ".bed", ".fam"),
    match = "QC/NonMatching_{dataset}.txt"
  output:
    multiext("QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip.match", ".bim", ".bed", ".fam")
  params:
    input = "QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip",
    out = "QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip.match"
  shell:
    """
    plink --bfile {params.in} --exclude {input.match} --make-bed --out {params.out}
    """

rule remove_cg_at:
  input:
    multiext("QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip.match", ".bim", ".bed", ".fam")
  output:
    multiext("QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip.match.noCGAT", ".bim", ".bed", ".fam"),
    txt = "QC/cg_at_loci.txt"
  params:
    input = "QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip.match",
    out = "QC/{dataset}.rareCommonMerged.dbsnp.noDup.flip.match.noCGAT"
  shell:
    """
    python scripts/find_cg_at.py {params.input}.bim > {output.txt}
    plink --bfile {params.input} --exclude {output.txt} --make-bed --out {params.out}
    """
