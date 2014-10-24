###This is the R script that is inserted into the bash script in order to run the final annotations
###This is just a test script

options(stringsAsFactors=F)

#specific data
dat <- read.delim("TEST.tsv")
clinvar <- read.delim("TEST.tsv.tempann.hg19_multianno.txt")

#annotation table
annot <- read.delim("Final_Gene_Annotation_Table.txt")

#prepare annotations
clinmat <- match(paste(dat[,1], dat[,2], dat[,4], dat[,5]), paste(clinvar[,1], clinvar[,2], clinvar[,4], clinvar[,5]))
clinvar <- clinvar[clinmat,6]
clinvar[is.na(clinvar)] <- "."

RVIS.specific <- rep(".", nrow(dat))
whi <- which(dat[,"VariantClass"]=="synonymousSNV")
RVIS.specific[whi] <- annot[whi,"TOLERANCE_SYNONYMOUS"]
whi <- which(dat[,"VariantClass"]=="nonsynonymousSNV")
RVIS.specific[whi] <- annot[whi,"TOLERANCE_MISSENSE"]
whi <- which(dat[,"VariantClass"]%in%c("stopgain", "stoploss"))
RVIS.specific[whi] <- annot[whi,"TOLERANCE_NONSENSE"]
whi <- which(grepl("splicing", dat[,"VariantFunction"]))
RVIS.specific[whi] <- annot[whi,"TOLERANCE_SPLICE"]
whi <- which(grepl("^frameshift", dat[,"VariantClass"]))
RVIS.specific[whi] <- annot[whi,"TOLERANCE_FRAMESHIFT"]

#find columns to split dat and annot at
GT1col <- grep("PredictionSummary", colnames(dat))+1
INFOcol <- grep("INFO", colnames(dat))-1
col1 <- grep("TOLERANCE_ALL_DALY", colnames(annot))
col2 <- grep("TOLERANCE_FRAMESHIFT", colnames(annot))+1

#add leading space to genotypes to deal with Excel date format issue
for(i in GT1col:(INFOcol-1)){
  dat[,i] <- paste(" ", dat[,i], sep="")
}

annmat <- match(dat[,"Gene"], annot[,"GENE"])

out <- cbind(dat[,1:INFOcol], clinvar, annot[annmat,1:col1], RVIS.specific, annot[annmat,col2:ncol(annot)])
write.table(out, "TESTOUTTSV", sep="\t", col.names=T, row.names=F, quote=F)










    
  
