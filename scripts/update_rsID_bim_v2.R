# update_rsID_bim_v2.R
# Author: Austin Reynolds
# 13 Aug 2020
#------------
#This script takes a bimfile with incorrect or no rsids and a bedfile with desired rsids as input
#and outputs a new bimfile with correct rsids
#------------
#Note to user. Often a position will have more than one rsid. We allow the code
#to remove these duplicates, so there is no guarantee that a particular rsid for any chr/pos will be included.
#if you want specific rsids to be included you will need to modify this.
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
  make_option(c("-r", "--rsid"), action="store", default="",
              help="Input rsids in bed format"),
  make_option(c("-s", "--silent"), action="store_false",default=TRUE, 
              help="Turn off verbose output"),
  make_option(c("-o", "--out"), action="store",
              help = "Output filename")
)
opt<-parse_args(OptionParser(option_list=option_list))

# Read in the bed and bimfile
rsids <- fread(opt$rsid)
bimfile <- fread(opt$bim)

# Assign column names to bed and bimfiles
colnames(bimfile)<-c("chr","bimid","cm","pos","ref","alt")
colnames(rsids)<-c("chr","V2","pos","rsid")

#plink uses numeric chromosome designations but
#dbsnp often uses chrN designations so this changes those to numeric
rsids$chr<-gsub("chr","",rsids$chr)
rsids$chr<-gsub("XY", 25, rsids$chr)
rsids$chr<-gsub("X", 23, rsids$chr)
rsids$chr<-gsub("Y", 24, rsids$chr)
rsids$chr<-gsub("MT", 26, rsids$chr)
rsids$chr<-as.numeric(rsids$chr)

# Add chromosome and position as the keys for these data tables
rsids<-data.table(rsids,key = c("chr","pos"))
bimfile<-data.table(bimfile,key = c("chr","pos"))

# Find matching and nonmatching rows between the bed and bimfiles
matches<-merge(rsids, bimfile, all=FALSE)
nonmatches<-bimfile[!rsids]

# Remove extra columns for final bimfile
matches<-matches[,c("chr","rsid","cm","pos","ref","alt")]

# Remove duplicate rsids from matches
matches<-unique(matches, by = c("chr","pos"))

# Rename nonmatches columns for final bimfile
colnames(nonmatches)<-c("chr","rsid","cm","pos","ref","alt")

# Combine matches and nonmatches
newbimfile<-rbind(matches,nonmatches)

# Sort by chr/pos
newbimfile<-data.table(newbimfile,key = c("chr","pos"))

# Unless --silent option is chosen, report the number of duplicates
if (opt$silent==TRUE) {
  cat(paste(nrow(matches),"out of",nrow(bimfile),"rsids corrected.\n"))
}

# Write new bimfile
write.table(newbimfile, opt$out, quote=F, row.names=F, col.names=F)


