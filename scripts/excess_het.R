# usage 
# Rscript excess_het.R infile outfile 

#read in plink output
args <- commandArgs(trailingOnly = TRUE)

x<-read.table(args[1] ,header = T)
output = args[2]

#Calculate mean heterozygosity as (N–O)/N,
#where N is the number of non-missing genotypes
#and O is the observed number of homozygous genotypes for a given individual)

x$meanhet<-((x$N.NM.-x$O.HOM.)/x$N.NM.)

#get mean and 3sd upper and lower limits for the dataset
meanhet <- mean(x$meanhet)
upper.limit <- mean(x$meanhet)+(sd(x$meanhet)*3)
lower.limit <- mean(x$meanhet)-(sd(x$meanhet)*3)

#find individuals >3sd away from mean
dropSamples = x$IID[which(x$meanhet>upper.limit | x$meanhet<lower.limit)]

# write out 
write.table(dropSamples, file = output, sep = "\t", row.names = F, col.names = F, quote = FALSE)


