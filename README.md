# SALL2_conserved_network_analysis
Analysis of SALL2 isoform specific ChIP-seq data
```
###############################################
###############################################
### Isoform Expression, related to Figure 1 ###
###############################################
###############################################

mkdir Fig1 && cd Fig1

##### Import RNA-seq quantification  #####

wget -O SALL2_normal_tumor.isoforms https://usegalaxy.org/datasets/bbd44e69cb8906b53c31efa2219fa05e/display?to_ext=tabular

##### R enviroment #####
R 

#install.packages("dplyr")
library(dplyr)
#install.packages("tidyverse")
library(tidyverse)
require(ggplot2)
require(gridExtra)
data1<-read.table("SALL2_normal_tumor.isoforms", header = TRUE, sep = "\t", row.names=1)
data2<-as.matrix(data1)
dim(data2)

breaksList = seq(0, 30, by = 0.01)

# install.packages("pheatmap")
library(pheatmap)
library("RColorBrewer")
pdf("SALL2_Isoform_Expression.pdf", width=10, height=10)
plot2<-pheatmap(data2, main="SALL2 Isoform Expression", cluster_rows=TRUE, cluster_cols=FALSE, show_colnames=TRUE, fontsize=10.5, border_color=NA, cellwidth=30, cellheight=15, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)  
dev.off()
```
```
#################################################################
#################################################################
#### ChIPseeker Analysis: GBM vs ENCODE, related to Figure 2 ####
#################################################################
#################################################################

mkdir Fig1 && cd Fig1

##### Import ChIP-seq quantification and GO terms #####

wget -O ENCODE-SALL2-eGFP-MAPPED_COORDINATES_hg19.bed https://usegalaxy.org/datasets/bbd44e69cb8906b5f5ca0fb6b2260df1/display?to_ext=bed
wget -O GSM1306364_MGG8TPC.SALL2.bed https://usegalaxy.org/datasets/bbd44e69cb8906b59973e08722beb4ba/display?to_ext=bed
wget -O GO_short_E1A_HEK293.tabular https://usegalaxy.org/datasets/bbd44e69cb8906b5f7c0b7bd67f405bd/display?to_ext=tabular
wget -O GO_MGG8TPC.tabular https://usegalaxy.org/datasets/bbd44e69cb8906b554a7c1428ef1fb5a/display?to_ext=tabular

##### R enviroment #####
R  

##### ChIPseeker Analysis #####
# install.packages("tidyverse")
library(tidyverse)
# install.packages("dplyr")
library(dplyr)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("ChIPseeker")

library("ChIPseeker")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)

short_E1A_HEK293 <- read.table("ENCODE-SALL2-eGFP-MAPPED_COORDINATES_hg19.bed", header=F)
colnames(short_E1A_HEK293)[1] <- "chr"
colnames(short_E1A_HEK293)[2] <- "start"
colnames(short_E1A_HEK293)[3] <- "end"
head(short_E1A_HEK293)
short_E1A_HEK293 <- makeGRangesFromDataFrame(short_E1A_HEK293)
peakAnno <- annotatePeak(short_E1A_HEK293, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(peakAnno, file="short_E1A_HEK293_peaks.txt", sep="\t")

MGG8TPC <- read.table("GSM1306364_MGG8TPC.SALL2.bed", header=F)
colnames(MGG8TPC)[1] <- "chr"
colnames(MGG8TPC)[2] <- "start"
colnames(MGG8TPC)[3] <- "end"
head(MGG8TPC)
MGG8TPC <- makeGRangesFromDataFrame(MGG8TPC)
peakAnno <- annotatePeak(MGG8TPC, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(peakAnno, file="MGG8TPC_peaks.txt", sep="\t")


files <- GRangesList("short_E1A_HEK293" = short_E1A_HEK293, "MGG8TPC" = MGG8TPC)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

pdf("AnnoBar.pdf", width=9, height=3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("DistToTSS.pdf", width=9, height=3)
plotDistToTSS(peakAnnoList)
dev.off()


#### Plotting GO terms ####

require(ggplot2)
require(gridExtra)

data1<-read.table("GO_short_E1A_HEK293.tabular", header = TRUE, sep = "\t", row.names=1)
dim(data1)
data2<-read.table("GO_MGG8TPC.tabular", header = TRUE, sep = "\t", row.names=1)
dim(data2)

### Plotting Heatmap GO: ENCODE

subset <- c("negative_log10_of_adjusted_p_value")
data1.1<-data1[subset]
names(data1.1)[names(data1.1) == "negative_log10_of_adjusted_p_value"] <- "-log10(p_value)"

breaksList = seq(0, 15, by = 0.001)

library(pheatmap)
library("RColorBrewer")
pdf("GO_ENCODE.pdf", width=6, height=3)
plot1<-pheatmap(data1.1, main="short_E1A HEK293", cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=TRUE, fontsize=10.5, border_color=NA, cellwidth=25, cellheight=40, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList) 
dev.off()
###

### Plotting Heatmap GO: MGG8TPC

subset <- c("negative_log10_of_adjusted_p_value")
data2.1<-data2[subset]
names(data2.1)[names(data2.1) == "negative_log10_of_adjusted_p_value"] <- "-log10(p_value)"

breaksList = seq(0, 15, by = 0.001)

library(pheatmap)
library("RColorBrewer")
pdf("GO_MGG8TPC.pdf", width=9, height=9)
plot2<-pheatmap(data2.1, main="MGG8TPC", cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=TRUE, fontsize=10.5, border_color=NA, cellwidth=25, cellheight=13, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList) 
dev.off()
###
```
```
##########################################################################
##########################################################################
#### ChIPseeker Analysis: WT vs E1A-KO in HEK293, related to Figure 4 ####
##########################################################################
##########################################################################

mkdir Fig4 && cd Fig4

##### Import ChIP-seq quantification and GO terms #####

wget -O MACS2_callpeak_WT_p_0.005.bed https://usegalaxy.org/datasets/bbd44e69cb8906b52eb38042e78ca3e9/display?to_ext=interval
wget -O MACS2_callpeak_E1_p_0.005.bed https://usegalaxy.org/datasets/bbd44e69cb8906b51728a5e4733670ee/display?to_ext=interval
wget -O GO_WT.tabular https://usegalaxy.org/datasets/bbd44e69cb8906b5ab65f56fc668ddc9/display?to_ext=tabular
wget -O GO_E1.tabular https://usegalaxy.org/datasets/bbd44e69cb8906b526a771b6bf68b31e/display?to_ext=tabular

##### R enviroment #####
R

##### ChIPseeker Analysis #####
# install.packages("tidyverse")
library(tidyverse)
# install.packages("dplyr")
library(dplyr)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("ChIPseeker")

library("ChIPseeker")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)

WT <- read.table("MACS2_callpeak_WT_p_0.005.bed", header=F)
colnames(WT)[1] <- "chr"
colnames(WT)[2] <- "start"
colnames(WT)[3] <- "end"
head(WT)
WT <- makeGRangesFromDataFrame(WT)
peakAnno <- annotatePeak(WT, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(peakAnno, file="WT_peaks.txt", sep="\t")

E1 <- read.table("MACS2_callpeak_E1_p_0.005.bed", header=F)
colnames(E1)[1] <- "chr"
colnames(E1)[2] <- "start"
colnames(E1)[3] <- "end"
head(E1)
E1 <- makeGRangesFromDataFrame(E1)
peakAnno <- annotatePeak(E1, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(peakAnno, file="E1_peaks.txt", sep="\t")

files <- GRangesList("WT" = WT, "E1" = E1)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

pdf("AnnoBar.pdf", width=9, height=3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("DistToTSS.pdf", width=9, height=3)
plotDistToTSS(peakAnnoList)
dev.off()


#### Plotting GO terms ####

require(ggplot2)
require(gridExtra)

data1<-read.table("GO_WT.tabular", header = TRUE, sep = "\t", row.names=1)
dim(data1)
data2<-read.table("GO_E1.tabular", header = TRUE, sep = "\t", row.names=1)
dim(data2)

### Plotting Heatmap GO: WT

subset <- c("negative_log10_of_adjusted_p_value")
data1.1<-data1[subset]
names(data1.1)[names(data1.1) == "negative_log10_of_adjusted_p_value"] <- "-log10(p_value)"

breaksList = seq(0, 15, by = 0.01)

library(pheatmap)
library("RColorBrewer")
pdf("GO_WT_HEK293.pdf", width=9, height=9)
plot1<-pheatmap(data1.1, cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=TRUE, fontsize=10.5, border_color=NA, cellwidth=25, cellheight=13, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
dev.off()
###

### Plotting Heatmap GO: E1

subset <- c("negative_log10_of_adjusted_p_value")
data2.1<-data2[subset]
names(data2.1)[names(data2.1) == "negative_log10_of_adjusted_p_value"] <- "-log10(p_value)"

breaksList = seq(12, 35, by = 0.01)

library(pheatmap)
library("RColorBrewer")
pdf("GO_E1_HEK293.pdf", width=9, height=9)
plot2<-pheatmap(data2.1, cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=TRUE, fontsize=10.5, border_color=NA, cellwidth=25, cellheight=13, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList) 
dev.off()
###
```
```
#################################################
#################################################
### ChIP-seq correlation, related to Figure 4 ###
#################################################
#################################################

#### Using fastq-dump to download and recover fastq reads from SRA files. 

prefetch --max-size 800G -O ./ SRR11184885  # WT-HEK293 
prefetch --max-size 800G -O ./ SRR11184886  # E1-HEK293
prefetch --max-size 800G -O ./ SRR11184887  # KO-HEK293
prefetch --max-size 800G -O ./ SRR577512    # H3K4me3 rep1 ENCODE
prefetch --max-size 800G -O ./ SRR577513    # H3K4me3 rep2 ENCODE
prefetch --max-size 800G -O ./ SRR9604322   # H3K27me3 PRJNA551302
prefetch --max-size 800G -O ./ DRR014667    # H3K27ac PRJDB2619
wget https://www.encodeproject.org/files/ENCFF002ABA/@@download/ENCFF002ABA.fastq.gz # H3K27ac ENCODE rep1
wget https://www.encodeproject.org/files/ENCFF002ABC/@@download/ENCFF002ABC.fastq.gz # H3K27ac ENCODE rep2	

#### Downloading hg19 build 

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
chmod 755 twoBitToFa
./twoBitToFa hg19.2bit hg19.fa
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
chmod 755 genePredToGtf
gunzip refGene.txt.gz
cut -f 2- refGene.txt | ./genePredToGtf file stdin -source=hg19_Ref hg19.gtf

#### Converting SRA files to fastq.gz

SRA= ls -1 *.sra
for SRA in *.sra; do fastq-dump --gzip ${SRA}
done

#### Trimming downloaded Illumina datasets with fastp.
a= ls -1 *.fastq.gz
for a in *.fastq.gz; do fastp -w 16 -i ${a} -o ${a}.fastp
gzip ${a}.fastp
done


#### Aligning illumina datasets againts reference with minimap, using n threads."   screen -r -167871.pts-0.X39932CORE
samtools faidx hg19.fa
b= ls -1 *.fastq.gz.fastp.gz
for b in *.fastq.gz.fastp.gz; do minimap2 -ax sr hg19.fa ${b} > ${b}.sam -t 55
samtools sort ${b}.sam > ${b}.sam.sorted.bam -@ 55
rm ${b}.sam
rm ${b}
done

#### Renaming files in bash
for filename in *.bam; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fastq.gz.fastp.gz.sam.sorted//g')";  done

#### Indexing BAM files
f= ls -1 *.bam
for f in *.bam; do samtools index ${f}; done

#### BigWig files         # https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
bamCoverage -b SRR11184885.bam -o WT-HEK293.bw -p 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220
bamCoverage -b SRR11184886.bam -o E1-HEK293.bw -p 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220
bamCoverage -b SRR11184887.bam -o KO-HEK293.bw -p 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220
bamCoverage -b SRR577512.bam -o H3K4me3_rep1_ENCODE.bw -p 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220
bamCoverage -b SRR577513.bam -o H3K4me3_rep2_ENCODE.bw -p 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220
bamCoverage -b SRR9604322.bam -o H3K27me3.bw -p 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220
bamCoverage -b DRR014667.bam -o H3K27ac.bw -p 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220
bamCoverage -b ENCFF002ABA.bam -o H3K27ac_rep1_ENCODE.bw -p 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220
bamCoverage -b ENCFF002ABC.bam -o H3K27ac_rep2_ENCODE.bw -p 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220


#### Plot correlations
wget -O MACS2_callpeak_WT_p_0.005.bed https://usegalaxy.org/datasets/bbd44e69cb8906b52eb38042e78ca3e9/display?to_ext=interval
wget -O MACS2_callpeak_E1_p_0.005.bed https://usegalaxy.org/datasets/bbd44e69cb8906b51728a5e4733670ee/display?to_ext=interval

multiBigwigSummary BED-file --BED MACS2_callpeak_WT_p_0.005.bed --bwfiles WT-HEK293.bw E1-HEK293.bw H3K4me3_rep1_ENCODE.bw H3K4me3_rep2_ENCODE.bw H3K27me3.bw H3K27ac.bw H3K27ac_rep1_ENCODE.bw H3K27ac_rep2_ENCODE.bw --labels WT-HEK293 E1-HEK293 H3K4me3_rep1_ENCODE H3K4me3_rep2_ENCODE H3K27me3 H3K27ac H3K27ac_rep1_ENCODE H3K27ac_rep2_ENCODE --numberOfProcessors 55 --binSize 50 -o WT_BigWig.npz

multiBigwigSummary BED-file --BED MACS2_callpeak_E1_p_0.005.bed --bwfiles WT-HEK293.bw E1-HEK293.bw H3K4me3_rep1_ENCODE.bw H3K4me3_rep2_ENCODE.bw H3K27me3.bw H3K27ac.bw H3K27ac_rep1_ENCODE.bw H3K27ac_rep2_ENCODE.bw --labels WT-HEK293 E1-HEK293 H3K4me3_rep1_ENCODE H3K4me3_rep2_ENCODE H3K27me3 H3K27ac H3K27ac_rep1_ENCODE H3K27ac_rep2_ENCODE --numberOfProcessors 55 --binSize 50 -o E1_BigWig.npz

# Pearson_BigWig_WT
plotCorrelation -in WT_BigWig.npz --corMethod pearson --skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o WT_PearsonCorr_readCounts_BigWig.pdf --outFileCorMatrix WT_PearsonCorr_readCounts_BigWig.tab --removeOutliers
# Spearman_BigWig_E1
plotCorrelation -in WT_BigWig.npz --corMethod spearman --skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o WT_SpearmanCorr_readCounts_BigWig.pdf --outFileCorMatrix WT_SpearmanCorr_readCounts_BigWig.tab --removeOutliers

# Pearson_BigWig_E1
plotCorrelation -in E1_BigWig.npz --corMethod pearson --skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o E1_PearsonCorr_readCounts_BigWig.pdf --outFileCorMatrix E1_PearsonCorr_readCounts_BigWig.tab --removeOutliers
# Spearman_BigWig_E1
plotCorrelation -in E1_BigWig.npz --corMethod spearman --skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o E1_SpearmanCorr_readCounts_BigWig.pdf --outFileCorMatrix E1_SpearmanCorr_readCounts_BigWig.tab --removeOutliers
```
