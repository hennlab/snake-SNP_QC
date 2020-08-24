# This script accepts the output of the plink function --hardy (.hwe file)
# It returns a list of SNPs that meet the user-specified criterion
# The arguments are: Input (the .hwe file), 'over' or 'under' to specify the direction of the filter, and a numeric threshold
# Optional fourth argument to specify the output file
# Ex.Rscript het_filter.R plink.hwe over 0.8
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("At least three arguments must be supplied", call.=FALSE)
} else if (length(args)>4) {
  stop("Too many arguments", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  args[4] = "het_filter_out.txt"
}
HWEfile <- read.table(args[1], header=T, stringsAsFactors=F)
if (args[2]=="over") {
  SNPnames <- HWEfile[which(HWEfile$O.HET. > args[3]), "SNP"]
} else if (args[2]=="under") {
  SNPnames <- HWEfile[which(HWEfile$O.HET. < args[3]), "SNP"]
}
write.table(SNPnames, args[4], quote=F, row.names=F, col.names=F)
