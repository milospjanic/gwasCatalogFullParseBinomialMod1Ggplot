#!/bin/bash

SECONDS=0

#download gwas catalog and create bed file with chr;position;position+1;proxy_gene;phenotype
echo "
#download gwas catalog 

wget https://www.dropbox.com/s/43d9vh10k44jqbu/full.mod.crossmap.plus.snpsnpinteract

#awk -F\"\t\" '{if (\$12!=\"\") print \$12\"\t\"\$13\"\t\"\$15\"\t\"\$8}' gwascatalog.txt > tmp
#awk -F\"\t\" '{print \$1\"\t\"\$2\"\t\"\$2+1\"\t\"\$3\"\t\"\$4}' tmp > GwasCatalog.bed
#rm tmp" > GwasCatalog2Bed.sh
#chmod 775 GwasCatalog2Bed.sh
#./GwasCatalog2Bed.sh

cp full.mod.crossmap.plus.snpsnpinteract GwasCatalog.bed

#tail -n +2 GwasCatalog.bed > tmp
#mv tmp GwasCatalog.bed


#adding CardiogramPlusC4D to GWASCatalog

wget https://stanfordmedicine.box.com/shared/static/pqxkuzwgv8bhl8ne05a3ohlmzgwbir28.bed
mv pqxkuzwgv8bhl8ne05a3ohlmzgwbir28.bed CARDIOGRAMplusC4DleadSNPs.bed

awk '{print $1,$2,$3,$4,"CardiogramPlusC4D"}' FS='\t' OFS='\t' CARDIOGRAMplusC4DleadSNPs.bed  > CARDIOGRAMC4Dplusnovel.txt.tmp
sed -i 's/^chr//g' CARDIOGRAMC4Dplusnovel.txt.tmp
cat CARDIOGRAMC4Dplusnovel.txt.tmp >> GwasCatalog.bed
rm CARDIOGRAMC4Dplusnovel.txt.tmp
sed -i 's/\//-/g' GwasCatalog.bed 

awk -F'\t' '
    FNR == 1 { header = $0;next }
    !seen[$5]++ {}
    {
        print > ("GWASCatalogPhenotype_"$5".txt");
    }
' GwasCatalog.bed

echo "Gwas Catalog number of SNP-phenotype associations:" > report.txt
wc -l GwasCatalog.bed >> report.txt

echo "Gwas Catalog number of SNP-phenotype associations:" 
wc -l GwasCatalog.bed 

echo "Gwas Catalog number of SNP-phenotype associations per category:" >> report.txt
find -name "GWASCatalogPhenotype_*" -type f | rename 's/ /_/g'

for file in `find . -name "GWASCatalogPhenotype_*.txt"`
do
  filename=$(basename $file .txt)
  echo Phenotype: ${filename##GWASCatalogPhenotype_} >> report.txt
wc -l < "$file" >> report.txt
done
echo "***Finished parsing Gwas Catalog SNP-phenotype associations per category***"


for file in `find . -name "GWASCatalogPhenotype_*.txt"`
do
cut -f1-3 $file | sort -k1,1V -k2,2n | uniq > $file.cut.sort.uniq
done

#substitute 23 24 with X Y, add chr
echo ***Converting to GWAS Phenotype specific bed files*** >>report.txt
echo ***Converting to GWAS Phenotype specific bed files*** 
for file in `find . -name "GWASCatalogPhenotype_*.txt.cut.sort.uniq"`
do
sed -i 's/^23/X/g' $file
sed -i 's/^24/Y/g' $file
sed  's/^/chr/g' $file >  $file.chrXY
rm $file
done

echo "***Finished creating Gwas Catalog phenotype specific bed files***"

#overlap with input
echo ***Overlapping Phenotype SNPs with input bed***
echo ***Overlapping Phenotype SNPs with input bed***>>report.txt
for file in `find . -name "GWASCatalogPhenotype_*.txt.cut.sort.uniq.chrXY"`
do
bedtools intersect -a $file -b $1 > $file.overlap
filename=$(basename $file .txt.cut.sort.uniq.chrXY)
echo Overlap: ${filename##GWASCatalogPhenotype_} >>report.txt
wc -l < $file.overlap >>report.txt
done
echo ***Finished overlapping GWAS Catalog Phenotype specific bed with input bed***

echo ***Downloading hg19 chromosome sizes***

#get the size of hg19
wget https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
hg19=$(cat hg19.chrom.sizes | awk '{ sum+=($2)} END {print sum}')
echo Human Genome size version hg19: $hg19

#bedtools merge on input bed
sort -k1,1V -k2,2n $1 > tmp
bedtools merge -n -i tmp > $1
wc -l $1
rm tmp

echo ***Selecting FDR 10% overlaps***
#select top overlaps
find -name "*sort.uniq.chrXY.overlap" -type f -exec wc -l {} + | grep -v total$ | sort -rn  | head -n5

find -name "*sort.uniq.chrXY.overlap" -type f -exec wc -l {} + | grep -v total$ | sort -rn  | head -n100 > top
awk '{print $2}' FS=' ' OFS='\t' top > top2
tr '\n' ' ' < top2 >top3
for i in $(cat top3); do cp $i $i.select;done

sed 's/.overlap//g' top3 > top4
for i in $(cat top4); do cp $i $i.select;done

echo ***Start calculating coverage***
#calculate coverage
cov=$(cat $1 | awk '{ sum+=($3-$2)} END {print sum}')
echo Coverage of BED file $cov >>report.txt
fra=$(cat $1 | awk '{ sum+=($3-$2)} END {print sum/"'"$hg19"'"}')
echo Fraction of hg19 $fra >>report.txt

echo ***Overlapping Phenotype SNPs with input bed to calculate coverage***
#make overlaps with input bed ranges
echo ***Overlapping Phenotype SNPs with input bed to calculate coverage***>>report.txt
for file in `find . -name "GWASCatalogPhenotype_*.txt.cut.sort.uniq.chrXY.select"`
do	
bedtools intersect -wb -a $file -b $1 > $file.overlap.input.int
cut -f4-6 $file.overlap.input.int > $file.overlap.input.int.cut
done

echo "***Creating R script***"

#create  R script
touch script.R
echo "#!/usr/bin/Rscript" > script.R
echo "existingDF <- as.data.frame(matrix(seq(4),nrow=1,ncol=4))" >>script.R
#calculate coverage and fraction per category, load into R script
i=1
echo GWAS Catalog Phenotype     Total SNPs      Overlap Fold change     Fraction of hg19        Peak coverage > output.table.txt

for file in `find . -name "GWASCatalogPhenotype_*.txt.cut.sort.uniq.chrXY.select.overlap.input.int.cut"`
do
let i=i+1
#removing input spaces and '
var2a=$(echo $(basename $file .txt.cut.sort.uniq.chrXY.select.overlap.input.int.cut) | tr ' ' "_")
var2=${var2a//[^[a-zA-Z0-9_]]/}

var3="$(cat "$file" | awk '{ sum+=($3-$2)} END {print sum}')"

fra=$(cat "$file" | awk '{ sum+=($3-$2)} END {print sum/"'"$hg19"'"}')

#calculating fold change
overlap="$(wc -l $(basename $file .txt.cut.sort.uniq.chrXY.select.overlap.input.int.cut).txt.cut.sort.uniq.chrXY.overlap.select | cut -f1 -d ' ')"
total="$(wc -l $(basename $file .txt.cut.sort.uniq.chrXY.select.overlap.input.int.cut).txt.cut.sort.uniq.chrXY.select | cut -f1 -d ' ')"
fold="$(awk 'BEGIN {print (100*"'"$overlap"'"/"'"$total"'")}')"
filename=$(basename $file .txt.cut.sort.uniq.chrXY.select.overlap.input.int.cut)

echo GWAS Catalog Phenotype: ${filename##GWASCatalogPhenotype_} Total SNPs: $total Overlap: $overlap Fold change: $fold Fraction of hg19 $fra Peak coverage: $var3
echo ${filename##GWASCatalogPhenotype_}	$total	$overlap	$fold	$fra	$var3 >> output.table.txt 

echo  "x<-dbinom ($(wc -l $(basename $file .txt.cut.sort.uniq.chrXY.select.overlap.input.int.cut).txt.cut.sort.uniq.chrXY.overlap.select | cut -f1 -d ' '), $(wc -l $(basename $file .txt.cut.sort.uniq.chrXY.select.overlap.input.int.cut).txt.cut.sort.uniq.chrXY.select | cut -f1 -d ' '), "$fra")" >> script.R
echo "y<-"$fold"">> script.R
echo "name<-\""$var2"\"">>script.R
echo "s<-$(wc -l $(basename $file .txt.cut.sort.uniq.chrXY.select.overlap.input.int.cut).txt.cut.sort.uniq.chrXY.select | cut -f1 -d ' ')">>script.R
echo "z<-c(name,x,y,s)">>script.R
echo "existingDF <- rbind(existingDF,z)">>script.R
done


#finishing R script

echo "existingDF<-existingDF[-1,]">>script.R
echo "existingDF[,2:4]<-sapply(existingDF[,2:4], as.numeric)">>script.R #you can test if numeric with sapply(existingDF, mode)
echo "existingDF<-transform(existingDF, V2=-log(V2))" >> script.R

echo "data<-existingDF" >>script.R
echo "rownames(data) <- data\$V1" >>script.R
echo "data<-data[,2:4]">>script.R

#adding categories
echo "k<-dim (data)" >>script.R
echo "rep <-rep(\"Other\", k[1])">>script.R
echo "data\$V5 <- rep">>script.R

echo "data=transform(data, V5=ifelse(grepl("\"Schizophrenia\"", rownames(data), ignore.case=T)==TRUE, \"Brain\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Bipolardisorder\"", rownames(data), ignore.case=T)==TRUE, \"Brain\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Multiplesclerosis\"", rownames(data), ignore.case=T)==TRUE, \"Brain\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Parkinsonsdisease\"", rownames(data), ignore.case=T)==TRUE, \"Brain\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Alzheimersdisease\"", rownames(data), ignore.case=T)==TRUE, \"Brain\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Rheumatoidarthritis\"", rownames(data), ignore.case=T)==TRUE, \"Chronic Inflammatory\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Ulcerativecolitis\"", rownames(data), ignore.case=T)==TRUE, \"Chronic Inflammatory\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Crohnsdisease\"", rownames(data), ignore.case=T)==TRUE, \"Chronic Inflammatory\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Lupus\"", rownames(data), ignore.case=T)==TRUE, \"Chronic Inflammatory\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"CardiogramplusC4D\"", rownames(data), ignore.case=T)==TRUE, \"Cardiovascular\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Myocardialinfarction\"", rownames(data), ignore.case=T)==TRUE, \"Cardiovascular\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Coronaryarterycalcification\"", rownames(data), ignore.case=T)==TRUE, \"Cardiovascular\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Artery\"", rownames(data), ignore.case=T)==TRUE, \"Cardiovascular\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Coronary\"", rownames(data), ignore.case=T)==TRUE, \"Cardiovascular\", V5))">>script.R
echo "data=transform(data, V5=ifelse(grepl("\"Cancer\"", rownames(data), ignore.case=T)==TRUE, \"Cancer\", V5))">>script.R

#removing 0s
echo "data[, 1:3] <- sapply(data[,1:3], as.numeric)">>script.R
echo "row_sub = apply(data[,1:3], 1, function(y) all(y != 0))">>script.R
echo "data<-data[row_sub,]">>script.R
echo "data[, 1:3] <- sapply(data[,1:3], as.numeric)">>script.R
echo "colnames(data)<-c(\"LogP\", \"FC\", \"Phenotype SNPs\", \"Category\")">>script.R

#short row names
echo "rownames(data) = gsub(\"GWASCatalogPhenotype_\",\"\",rownames(data))">>script.R

#select top
echo "data <- data[order(data\$LogP,decreasing = TRUE), ]">>script.R
echo "data <- data[1:25, ]">>script.R
echo "print(\"***Final table - top 25 categories - for plot - output.pdf***\")" >> script.R
echo "data">>script.R

#making ggplot2 graph
echo "library(ggplot2)" >> script.R
echo "library(wesanderson)">>script.R
echo "library(directlabels)">>script.R
echo "ymax<-max(data\$LogP,na.rm = TRUE)">>script.R
echo "ymin<-min(data\$LogP,na.rm = TRUE)">>script.R
echo "xmax<-max(data\$FC,na.rm = TRUE)">>script.R
echo "xmin<-min(data\$FC,na.rm = TRUE)">>script.R
echo Finishing R script 
echo Writing ggplot2 code for plot in output.pdf
echo "p<- ggplot(data, aes(x=data\$FC, y=data\$LogP,label=row.names(data))) + geom_point(shape=19, alpha=1/8, color=\"red\", aes(size=data\$\"Phenotype SNPs\"), max_size=max(data\$\"Phenotype SNPs\")) + xlab(\"Fold change\") + ylab(\"-log P-value\") + ggtitle (\"GWAS SNPs enrichment - binomial test\") + geom_dl(aes(label=row.names(data)), method=list(\"last.bumpup\"), col=\"blue\", alpha=1/2)+ylim(ymin, ymax) +xlim(xmin-15, xmax)" >>script.R
echo "pdf(\"output.pdf\",width=10, height=8)">>script.R
echo "print(p+ geom_dl(aes(colour = data\$\"Category\"), method=list(\"last.bumpup\")) + scale_colour_hue(name=\"Category\") + labs(size=\"Phenotype SNPs\", color=\"Category\") + scale_size(range = c(0,50)) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face=\"bold\")) + scale_color_manual(values = c(wes_palette(\"Cavalcanti\"), wes_palette(\"Royal1\"), wes_palette(\"GrandBudapest\"), wes_palette(\"Royal2\"), wes_palette(\"Darjeeling\"), wes_palette(\"Zissou\")))) ">> script.R
echo "dev.off()">>script.R


chmod 775 script.R
./script.R
rm script.R
rm GWASCatalogPhenotype*
rm top*
rm hg19.chrom.sizes
rm gwascatalog.txt

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
