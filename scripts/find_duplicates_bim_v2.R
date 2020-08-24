library(data.table)

args=commandArgs(TRUE)
bim.file <- args[1]

bim <- fread(bim.file)
index<-which(duplicated(bim$V2))
if(length(index!=0)){
bim$V2 <- as.character(bim$V2)
bim$V2[index] <- paste(bim$V2[index],"_DUP",sep="")
non_dup <- bim[-index,]
write.table(bim, bim.file, col.names=F,row.names=F,quote=F,sep='\t')
write.table(non_dup, paste("NonDup_",bim.file,sep=""), col.names=F,row.names=F,quote=F,sep='\t')
}else{
    print("No duplicates detected!")
}

