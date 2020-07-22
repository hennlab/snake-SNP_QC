# AB_to_top_allele_v2.R
# Author: Austin Reynolds
# An Rscript to replace the original python version by Meng Lin
#
# This script takes a bimfile of rare alleles called in "AB" format from zCall
# and coverts them to ACGT format as listed in an Illumina Strand Report.
#
# v2 speeds up the conversion process using the data.table package and eliminating
# the ifelse logic loop of the previous python script
#
# We are currently working with a ForwardAllele strand report for the H3Africa Array
# which does not contain top or bottom strand information, so alignment to something like
# 1000G is necessary to ensure your data is compatible with other things


# Define functions
pkgTest <- function(x){
  #A function to check if a package is installed, 
  #optionally install it, and load it for use
  if (!require(x,character.only = TRUE, quietly = T)){
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE, quietly = T)) stop("Package not found")
  }
}

# Call required packages
pkgTest("data.table")
pkgTest("optparse")

# Define input options
option_list <- list(
  make_option(c("-b", "--bim"), action="store", default="",
              help="Input path to bim file output from zCall"),
  make_option(c("-s", "--strand"), action="store", default="",
              help="Input path to Illumina Strand file"),
  make_option(c("-v", "--verbose"), action="store_false",default=TRUE, 
              help="Turn off verbose output"),
  make_option(c("-o", "--out"), action="store",
              help = "Output file path")
)
# Add user input options to list
opt<-parse_args(OptionParser(option_list=option_list))

# Read in bimfile
bimfile <- fread(opt$bim)
# Assign column names for bimfile
colnames(bimfile)<- c("chr","rsid","dist","pos","Allele1","Allele2")
#add an order column just in case the combining and replacement screw this up
#otherwise the bimfile won't match the bedfile and cause downstream errors
bimfile$ordercolumn<-seq(1,nrow(bimfile)) 

# Read in strandfile
strandfile <- fread(opt$strand)

#get the overlap between the bimfile and strandfile for easy allele replacement
#this is done based on SNP names because in our case these are shared and chr/pos notation is not between files
#but this won't always be true. This is the most likely source of future bugs...
overlap<-strandfile[bimfile, on = c(SNP_Name="rsid")]

#print number of overlapping sites
if (opt$verbose==TRUE) {
  print(paste(nrow(overlap),"/",nrow(bimfile)," sites overlap", 
              sep = ""))
}

# Make logical subsets of bim alleles 1 and 2 that correspond to either "A" or "B"
# there is probably a more elegant way to do this
A1A <- overlap$Allele1 %in% "A"
A1B <- overlap$Allele1 %in% "B"
A2A <- overlap$Allele2 %in% "A"
A2B <- overlap$Allele2 %in% "B"

# Sequentially replace "AB" allele notation with corresponding "ACGT" notation from the strand file
# Since "A" is present in both notations we have to be careful not to overwrite and introduce errors
overlap[A1A,"Allele1"] <- overlap[A1A, "Forward_Allele1"]
overlap[A1B,"Allele1"] <- overlap[A1B, "Forward_Allele2"]
overlap[A2A,"Allele2"] <- overlap[A2A, "Forward_Allele1"]
overlap[A2B,"Allele2"] <- overlap[A2B, "Forward_Allele2"]

# Remove all the extra columns to get back to bim format
newbimfile<-overlap[order(ordercolumn),c("chr","SNP_Name","dist","pos","Allele1","Allele2")]

# Write converted bimfile to out
write.table(newbimfile, opt$out, quote=F, row.names=F, col.names=F, sep="\t")
