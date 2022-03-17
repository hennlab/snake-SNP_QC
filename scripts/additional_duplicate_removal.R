# additional_duplicate_removal.R
# Author: Austin Reynolds
# 10 Mar 2022
#------------
#This script takes a bimfile with chr/pos duplicates and removes them.
#This was made because plink will remove chr/pos duplicates with the same ref/alt alleles
#but not those with different ref/alt alleles
#------------
#Note to user. We allow the code to remove these duplicates, so there is no guarantee 
#that a particular rsid for any chr/pos will be included. if you want specific rsids 
#to be included you will need to modify this.
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

# Assign column names to bed and bimfiles
colnames(bimfile)<-c("chr","rsid","cm","pos","ref","alt")

# get list of variants with chr0
zeros <- bimfile[which(bimfile$chr==0),]

#get list of duplicate variants
dups <- bimfile[duplicated(bimfile, by = c("chr", "pos"))]

#combine the snps to remove
to_remove <- rbind(zeros, dups)

# Unless --silent option is chosen, report the number of duplicates
if (opt$silent==TRUE) {
  cat(paste(nrow(dups),"duplicates identified for removal.\n"))
}

# Write to remove snps to file
write.table(to_remove$rsid, opt$out, quote=F, row.names=F, col.names=F,sep = "\t")
