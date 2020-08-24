# find_duplicates_bim_v3.R
# Author: Austin Reynolds
# 12 Aug 2020
#------------
#This script takes a bimfile as input and identifies duplicates based on chromosome and position information.
#The program then outputs a file in a plink-compatible format that can be used with the 
#--extract command to remove duplicate snps from your dataset.
#------------
#Note to user. Many of the duplicates being removed, at least in the case of the H3Africa array
#that this program was validated on, are different versions of the same variant 
#(eg. alleles are AT and CG or 0A and TA in the duplicates). 
#These are a very small number of cases in our validation dataset (~1500/2.1m), so this program 
#removes duplicates without regard to what allele they are. 
#This could result in missing data for some individuals at the affected snps, which
#could be an issue in smaller datasets or when you are interested in one of the affected snps.
#------------

# Define functions
pkgTest <- function(x){
  if (!require(x,character.only = TRUE, quietly = T)){
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE, quietly = T)) stop("Package not found")
  }
}

# Call packages
pkgTest("optparse")
pkgTest("data.table")
pkgTest("dplyr")

# User-defined options
option_list <- list(
  make_option(c("-b", "--bim"), action="store", default="",
              help="Input bim files"),
  make_option(c("-s", "--silent"), action="store_false",default=TRUE, 
              help="Turn off verbose output"),
  make_option(c("-o", "--out"), action="store",
              help = "Output filename")
)
opt<-parse_args(OptionParser(option_list=option_list))

# Read in the bimfile
bimfile <- fread(opt$bim)
# Assign colnames
colnames(bimfile)<-c("chr","rsid","cm","pos","ref","alt")

# Remove duplicate snps from the dataset
newbimfile <- bimfile %>% distinct(chr,pos, .keep_all= TRUE)

# Unless --silent option is chosen, report the number of duplicates
if (opt$silent==TRUE) {
  numdups <- nrow(bimfile) - nrow(newbimfile)
  cat(paste(numdups,"duplicates to remove.\n"))
}

# Print chr/snps to keep in the bim file
write.table(newbimfile$rsid, opt$out, quote=F, row.names=F, col.names=F)



#to identify duplicate columns
#dupe = bimfile[,c('chr','pos')] # select columns to check duplicates
#duplicates <- bimfile[duplicated(dupe) | duplicated(dupe, fromLast=TRUE),]
