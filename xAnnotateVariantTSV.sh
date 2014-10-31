#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N AnnTSV

# The purpose of this script to annotate a tab separated files containing lists of target variants. The file firsts gets ClinVar annotations from Annovar. The R script adds annotations from ClinVar at the variants level, and then annotations from a custom annotation table at the gene level.
# The input can either be a single tsv file or a list of tsv files. In the latter case the script should be called as an array job (qsub -t 1-X, where X is the number of lines in the file), each array job will run on a single tsv.

#get arguments
while getopts i:o: opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        o) OutFil="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

#check for array and if appropriate get the file name from the provided list
if [[ ! -z $SGE_TASK_ID ]] && [[ $SGE_TASK_ID != "undefined" ]];then
    InpFil=`head -n $SGE_TASK_ID $InpFil | tail -n 1`
fi

#if not output name has been specified derive it from the input file name
if [[ -z $OutFil ]]; then 
    OutFil=${InpFil/tsv/annot.tsv}
else
    OutFil=$OutFil.annot.tsv
fi

# location of annovar database adn the custom annotation table
AnnDir=/ifs/scratch/c2b2/af_lab/ads2202/bin/annovar
AnnTab=/ifs/scratch/c2b2/ys_lab/yshen/WENDY/VariantsAnnotationPipeline/Final_Gene_Annotation_Table.txt

echo $InpFil
echo $OutFil

#check input has variants to annotate
LinNum=`cat $InpFil | wc -l`
if [[ $LinNum -lt 2 ]];then
    echo "There are no variants to annotate"
    exit
fi


TempInp=$InpFil.tempann #name for temporary files

#cut the first 8 column of the tsv to simulate vcf format
tail -n +2 $InpFil | cut -f 1-8 > $TempInp

#convert to annovar format - annovar adjusts the way indels are represented
CMD1="convert2annovar.pl -includeinfo -format vcf4 $TempInp > $TempInp.2; mv -f $TempInp.2 $TempInp"

#annotate the new file
CMD2="table_annovar.pl $TempInp $AnnDir/humandb/ -buildver hg19 -remove -protocol clinvar_20140929 -operation f -nastring . -otherinfo -outfile $TempInp"

#Run commands
echo "Start Annovar at `date`"
echo $CMD1
eval $CMD1
echo $CMD2
eval $CMD2


#annotate using R
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)
#get data
dat.all <- read.delim("$InpFil")
dat <- dat.all[,1:(grep("FILTER", colnames(dat.all)))]
  #add leading space to genotypes to deal with Excel date format issue
GT1col <- grep("PredictionSummary", colnames(dat))+1
for(i in GT1col:ncol(dat)){
  dat[,i] <- paste(" ", dat[,i], sep="")
}
#get annotation table and extract relevant lines
annot <- read.delim("$AnnTab")
annmat <- match(dat[,"Gene"], annot[,"GENE"])
annot <- annot[annmat,]
if(any(is.na(annot[,1]))) { annot[is.na(annot[,1]),] <- "." } #replace missing data (NA) with "."
#get clinvar and match by variatns (CHR_POS_REF_ALT)
CLINVAR <- read.table("$TempInp.hg19_multianno.txt", skip=1)
  #need to match to specific variant, use the original variants that are at the end of the annovar output
clinmat <- match(paste(dat[,1], dat[,2], dat[,4], dat[,5]), paste(CLINVAR[,7], CLINVAR[,8], CLINVAR[,10], CLINVAR[,11]))
CLINVAR <- CLINVAR[clinmat,6]
if(any(is.na(CLINVAR))){CLINVAR[is.na(CLINVAR)] <- "."} #replace missing data (NA) with "."
#Create TOLERANCE.specific column with variant types
TOLERANCE.specific <- rep(".", nrow(dat))
whi <- which(dat[,"VariantClass"]=="synonymousSNV")
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_SYNONYMOUS"]}
whi <- which(dat[,"VariantClass"]=="nonsynonymousSNV")
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_MISSENSE"]}
whi <- which(dat[,"VariantClass"]%in%c("stopgain", "stoploss"))
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_NONSENSE"]}
whi <- which(grepl("splicing", dat[,"VariantFunction"]))
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_SPLICE"]}
whi <- which(grepl("^frameshift", dat[,"VariantClass"]))
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_FRAMESHIFT"]}
#find columns to split annot at 
col1 <- grep("TOLERANCE_ALL_DALY", colnames(annot))
col2 <- grep("TOLERANCE_FRAMESHIFT", colnames(annot))+1
out <- cbind(dat, CLINVAR, annot[,2:col1], TOLERANCE.specific, annot[,col2:ncol(annot)])
splicord <- grepl("splicing", out[,"VariantFunction"]) #sort splicing to the nonsense to the top
nonseord <- grepl("stop", out[,"VariantClass"]) # then stopgain/stoploss
indelfsord <- grepl("^frameshift", out[,"VariantClass"]) # then frameshift indels
indelnford <- grepl("nonframeshift", out[,"VariantClass"]) # then non-frameshift indels
predord <- match(out[,"PredictionSummary"], c("Med", "High")) #then high likelihood deleterious missense (sort order is "decreasing", so "M" would come above "H")
for(i in 1:ncol(out)){
  cnt <- which(nchar(out[,i])>32760) 
  if(length(cnt)>0){
    out[cnt,i] <- paste(substr(out[cnt,i], 1, 32740), ":::---TRUNCATED_FOR_EXCEL")
  }
}
ord <- order(predord, splicord, nonseord, indelfsord, indelnford, out[,"CADDscore"], out[,"GERP.."], decreasing=T) #final sort by CADD and then GERP
write.table(out[ord,], "$OutFil", sep="\t", col.names=T, row.names=F, quote=F)
RSCRIPT

if [[ $? -gt 0 ]]; then
    echo "Annotation Failed"
    exit
fi


###This is specific for the REGENERON PROJECT as the samples names were extended by Regeneron, this just removes the additional bits
# e.g our sample DD68 would be COL_CHUNG_DD68_Diabetes_123456789, this reverts it to DD68
cat $OutFil | awk '{ if ( $1 == "Chromosome" ){
    gsub ( /COL.CHUNG_/, "");
    gsub ( /_PAH_[0-9]*/, "");
    gsub ( /_FPPH_ADDED_[0-9]*/, "");
    gsub ( /_OBESITY_[0-9]*/, "");
    gsub ( /_RARE_OR_MENDELIAN_DISEASE_[0-9]*/, "");
    gsub ( /_WES_STUDY__RARE_OR_MENDELIAN__[0-9]*/, "");
    gsub ( /_DIABETES_[0-9]*/, "");
    gsub ( /_CHD_[0-9]*/, "");
    gsub ( /GERP../, "GERP++");
    gsub ( /prediction/, "");
    gsub ( /\.1/, "");
    };
    print }' > $OutFil.temp
mv $OutFil.temp $OutFil

# remove temporary files
rm -f $InpFil.tempann*
echo "Annotation complete"
