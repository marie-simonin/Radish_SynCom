# Radish_SynCom
Study performed at INRAE - IRHS : Script to reproduce bioinformatics and figures of the article

# Title: Transmission of synthetic seed bacterial communities to radish seedlings: impact on microbiota assembly and plant phenotype

# Bioinformatic analysis on raw reads - Denoising with DADA2 and filtering
### Raw Fastq reads are available on ENA accession number: PRJEB58635

## removing primers with cutadapt
for i in `cat group`; do cutadapt --discard-untrimmed -o $i.gyrB.R1.fq -p $i.gyrB.R2.fq -g MGNCCNGSNATGTAYATHGG -G ACNCCRTGNARDCCDCCNGA -e 0.1  -O 20 $i*L001_R1_001.fastq.gz $i*_L001_R2_001.fastq.gz; done

 
# Denoising in R with DADA2
```{r}
library(dada2); packageVersion("dada2")

list.files(path)

fnFs <- sort(list.files(path, pattern="gyrB.R1.fq", full.names = TRUE))

fnRs <- sort(list.files(path, pattern="gyrB.R2.fq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:8]) #select 200 (run1)

plotQualityProfile(fnRs[1:8]) #select 160 (run1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names

names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,160),

                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,

                     compress=TRUE, multithread=TRUE) #

head(out)

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE) #without pooling or pseudo-pooling (no need to detect rare ASV)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

table(nchar(getSequences(seqtab)))
```

## Since gyrB is a protein-coding genes only triplets should be conserved (244-247-250-253-256-259-262-265-268)
```{r}
seqtab244 <- seqtab[,nchar(colnames(seqtab)) %in% 244]

seqtab247 <- seqtab[,nchar(colnames(seqtab)) %in% 247]

seqtab250 <- seqtab[,nchar(colnames(seqtab)) %in% 250]

seqtab253 <- seqtab[,nchar(colnames(seqtab)) %in% 253]

seqtab256 <- seqtab[,nchar(colnames(seqtab)) %in% 256]

seqtab259 <- seqtab[,nchar(colnames(seqtab)) %in% 259]

seqtab262 <- seqtab[,nchar(colnames(seqtab)) %in% 262]

seqtab265 <- seqtab[,nchar(colnames(seqtab)) %in% 265]

seqtab268 <- seqtab[,nchar(colnames(seqtab)) %in% 268]
```

## Merge all files
```{r}
seq.final <- cbind(seqtab244, seqtab247, seqtab250, seqtab253, seqtab256, seqtab259, seqtab262, seqtab265, seqtab268)

dim(seq.final)

sum(seq.final)/sum(seqtab)

#Detect/Remove chimera

seqtab.nochim <- removeBimeraDenovo(seq.final, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seq.final)

## Summary

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)
```




# Figures of the paper

# Figure 1 - Meta-analysis radish samples
```{r}
SV_gyrB<-read.table("Subset3-gyrB-MiSeq_table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(SV_gyrB)
dim(SV_gyrB)
meta2 <- read.table("Metadata_gyrB_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(meta2)
dim(meta2)
SV_gyrB_use1<-merge(meta2,SV_gyrB,by="SampleID")
dim(SV_gyrB_use1)

head(SV_gyrB_use1)

#Subset only radish samples
SV_gyrB_R=subset(SV_gyrB_use1, Plant=="Radish")
head(SV_gyrB_R)
dim(SV_gyrB_R) #295 samples
```
```{r}
matrix<-SV_gyrB_R[c(40:ncol(SV_gyrB_R))]
dim(matrix)
head(matrix)
```

```{r}
##Keeping ASVs with at least 100 reads in the meta-analysis dataset
dim(matrix)
matrix_use<-matrix[,colSums(matrix)>=100]
dim(matrix_use)

#transposing
matrix_uset=t(matrix_use)
head(matrix_uset)
dim(matrix_uset)
```


```{r}
## Calculate prevalence and relative abundance for each ASV
#Code Shade lab: https://github.com/ShadeLab/PAPER_Shade_CurrOpinMicro/blob/master/script/Core_prioritizing_script.R
#presence-absence data
SV_gyrB_PA <- 1*((matrix_uset>0)==1)                                              
# occupancy calculation
SV_gyrB_Prevalence <- rowSums(SV_gyrB_PA)/ncol(SV_gyrB_PA) 
# relative abundance  
library(vegan)
library(dplyr)
SV_gyrB_relative_abundance <- apply(decostand(matrix_uset, method="total", MARGIN=2),1, mean)     

# combining occupancy and relative abundance of each SV_gyrB in a table
SV_gyrBprev_rel <- add_rownames(as.data.frame(cbind(SV_gyrB_Prevalence, SV_gyrB_relative_abundance)),'SV') 
head(SV_gyrBprev_rel)
dim(SV_gyrBprev_rel)
```


```{r}
## Merge prevalence/rel abund data with taxonomic info for each SV
taxo_gyrB<-read.table("Subset1-2_All_studies_merged_gyrB-rep-seqs-FINAL-filtered-taxonomy-final.tsv", header=TRUE, check.names = FALSE, sep = "\t")
SV_gyrBprev_rel_taxo<-merge(SV_gyrBprev_rel,taxo_gyrB,by="SV")
dim(SV_gyrBprev_rel_taxo)

head(SV_gyrBprev_rel_taxo)
```



## Plot Figure 1 Abundance-occupancy graph - gyrB - All phyla
```{r}
library(ggplot2)
tax_colors <-  c('Xanthomonas campestris'='red','Plantibacter sp'='lightblue' ,'Stenotrophomonas rhizophila'='#ffa600', "Erwinia persicina"="#bc5090", "Enterobacter cancerogenus"="purple", "Paenibacillus sp"= "#ff7c43", "Pseudomonas fluorescens 1"="#488f31", "Pseudomonas fluorescens 2"="#b7d7a7", "Pseudomonas fluorescens 3"="#004900", "Pseudomonas fluorescens 4"="#aebd38", "Pantoea agglomerans"="#216fd2", "Pseudomonas viridiflava"= "#f9b0d3", "Other taxa"="lightgrey" )

figure1=ggplot(data=SV_gyrBprev_rel_taxo, aes(x=SV_gyrB_relative_abundance, y=SV_gyrB_Prevalence)) +
    geom_point(aes(shape=Core_radish,color=Taxa), size=2) + xlab("Mean Taxa (ASV) Relative Abundance") + ylab("Taxa (ASV) Prevalence")+scale_x_log10(labels = scales::percent_format(accuracy = 0.1))+ scale_y_log10(labels = scales::percent_format(accuracy = 1))+ theme_classic()+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=9, face="bold"))+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold")) +ggtitle("gyrB gene - Bacteria")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+scale_color_manual(values=tax_colors)
figure1
```


```{r}
###Subset only strains to prepare the table in Figure 1
SV_gyrBprev_rel_taxo_strain=subset(SV_gyrBprev_rel_taxo, Taxa!="Other taxa")
SV_gyrBprev_rel_taxo_strain
dim(SV_gyrBprev_rel_taxo_strain) 
```


# Figure 2 Seedling bacterial biomass (CFU) for the different inoculation treatments
```{r}
pop=read.table("phenotype-pop-single-syncom.txt", header = T, sep = "\t")
head(pop)
dim(pop)
```

```{r}
pop$Condition<-ordered(pop$Condition, levels=c( "Control","Enterobacter cancerogenus","Erwinia persicina","Paenibacillus sp", "Pantoea agglomerans","Plantibacter sp","Pseudomonas fluorescens 1", "Pseudomonas fluorescens 2","Pseudomonas fluorescens 3","Pseudomonas fluorescens 4","Pseudomonas viridiflava","Stenotrophomonas rhizophila","Xanthomonas campestris", "6 strains", "8 strains", "12 strains"))
colors <-  c('Xanthomonas campestris'='red','Plantibacter sp'='lightblue' ,'Stenotrophomonas rhizophila'='#ffa600', "Erwinia persicina"="#bc5090", "Enterobacter cancerogenus"="purple", "Paenibacillus sp"= "#ff7c43", "Pseudomonas fluorescens 1"="#488f31", "Pseudomonas fluorescens 2"="#b7d7a7", "Pseudomonas fluorescens 3"="#004900", "Pseudomonas fluorescens 4"="#aebd38", "Pantoea agglomerans"="#216fd2", "Pseudomonas viridiflava"= "#f9b0d3", "6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")
Figure2=ggplot(pop, aes(x=log_CFU, y=Condition, fill=Condition)) + 
  geom_boxplot(outlier.shape = NA) + xlab("Log CFU by seedling")+ylab("Condition")+ theme_classic()+ theme(axis.title = element_text(color="black", size=12, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))	+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))+scale_fill_manual(values=c(colors))+facet_grid(Type~., scales  = "free", space = "free")+ theme(strip.text.y = element_text(size=8, face = "bold")) 
Figure2
```
## Stats Figure 2

```{r}
###Fit the  model. 
ggplot(pop, aes(log_CFU))+geom_histogram(binwidth=0.1)
```

```{r}
### Run first the full model
#Gaussian or Gamma (continuous data)
library(lme4)
library(MASS)
model.nb<-glm(log_CFU ~ Condition,  data=pop, family = Gamma(link="log"))
summary(model.nb)
#check if model fits the data P>0.05 is good
1- pchisq(summary(model.nb)$deviance, summary(model.nb)$df.residual)
plot(model.nb)
```

```{r}
###Estimate P values for main effects and interactions
library(car)
Anova(model.nb)
```


```{r}
###Estimate P values for main effects and interactions
fit <- aov(log_CFU ~ Condition, pop)# analyse de variance
summary(fit)
```

```{r}
###Posthoc tests
library(multcompView)
library(emmeans)
warp.emm <- emmeans(model.nb, ~Condition)
multcomp::cld(warp.emm)
```



# Metabarcoding figures

## Import gyrB dataset
```{r}
meta2 <- read.table("metadata_syncom_Anne_gyrB.txt", header=TRUE, check.names = FALSE, sep = "\t")
SV_16S<-read.table("SynCom_Anne_gyrB_SV_tax_filtered.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(SV_16S)
head(meta2)
dim(meta2)
dim(SV_16S)
SV_16S_use1<-merge(meta2,SV_16S,by="Sample_ID")
dim(SV_16S_use1)

head(SV_16S_use1)
``` 



```{r}
##Subset - remove neg controls and low samples (<1000 reads/sample)
SV_16S_use1=subset(SV_16S_use1, Final=="Yes")
dim(SV_16S_use1)
```

```{r}
### Make matrix of just SVs without metadata (needed to make distance matrix for NMDS and stats). Absolutely no metadata should be present, select just the columns with SV names.
matrix<-SV_16S_use1[c(11:ncol(SV_16S_use1))]

dim(matrix)
head(matrix)
```

```{r}
##Check rarefaction curves
#determine what is the lowest sequence count in a sample (if not already rarefied). 
#First column could be sampleID to be displayed on the graph.
library(vegan)
raremax <- min(rowSums(matrix))
raremax
rar=rarecurve(matrix, step=50, sample = raremax, cex = 0.6, ylab = "ASV richness", label = FALSE)
```
##Rarefy dataset to the lowest number of sequences observed in a sample (raremax) or replace "raremax" by the level you want to use
```{r}
#create rarefied matrix
matrix2=rrarefy(matrix, raremax)
#verify that now we have the same number of sequence per sample (Sample size) across the dataset
rar=rarecurve(matrix2, step=20, sample = raremax, cex = 0.6)
```


```{r}
### Remove singletons and SVs with abundance of 0 in all samples - OPTIONAL if not already done
dim(matrix2)
matrix_use<-matrix2[,colSums(matrix2)>=1]
dim(matrix_use)
```

```{r warning=FALSE, results='hide'}
# 1. Ordination: NMDS with Bray-Curtis distances
##Prepare an ordination (NMDS) based on Bray-Curtis dissimilarity matrix to analyze the community structure
library(vegan)
NMDS <- metaMDS(matrix_use, distance = "bray", trymax = 100)
```
```{r}
# Check if the NMDS is representing correctly the data in 2 dimensions: Stress must be inferior to 0.2 and R2 fit on the stressplot should be high >0.95
NMDS
stressplot(NMDS)
```

```{r}
##Extract scores for sites (samples) and add these results as new columns in SV_use1 dataframe with the metadata
NMDSsites=scores(NMDS, display="sites")
SV_16S_use1=cbind(SV_16S_use1,NMDSsites)
```

# Figure 3 panels A 
```{r warning=FALSE}
library(ggplot2)
tax_colors_16S <-  c('Control'='#ffbb94',"6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")
SV_16S_use1$Condition<-ordered(SV_16S_use1$Condition, levels=c("Control", "6 strains", "8 strains", "12 strains"))

p1=ggplot(data=SV_16S_use1, aes(NMDS1, NMDS2,
color=Condition, shape=Sample_type))+geom_point(size=2.5)+theme_classic()+xlab("NMDS1")+ylab("NMDS2")+scale_shape_manual(values=c(8, 16, 0))+ theme(legend.text = element_text(color="black", size=10, face="bold"))+scale_color_manual(values=tax_colors_16S)+ labs(shape = "Sample Type", color = "Condition")+ theme(legend.title = element_text(color="black", size=12, face="bold"))	+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=8, face="bold"))
p1
```

# Figure 3 panels D
```{r warning=FALSE}
library(ggplot2)
#Plot without controls
SV_16S_use1_nocontrols=subset(SV_16S_use1, Condition!="Control")
dim(SV_16S_use1_nocontrols)
matrix_nocontrols<-SV_16S_use1_nocontrols[c(11:4943)]
dim(matrix_nocontrols)
library(vegan)
matrix3=rrarefy(matrix_nocontrols, 14116)
matrix_use2<-matrix3[,colSums(matrix3)>=1]
dim(matrix_use2)
NMDS2 <- metaMDS(matrix_use2, distance = "bray", trymax = 100)
NMDS2
stressplot(NMDS2)
NMDSsites=scores(NMDS2, display="sites")
SV_16S_use1_nocontrols=cbind(SV_16S_use1_nocontrols,NMDSsites)
```
## Permanova Fig 3D
```{r warning=FALSE}
# Permanova - Distance based multivariate analysis of variance
##Effect of SynCom and sample type
Alltreatments=SV_16S_use1_nocontrols[c(2:10)] 
head(Alltreatments)
#Adonis function in Vegan
adonis2 <- adonis(matrix_use2 ~ Condition*Sample_type, method = "bray", permutations = 999, data=Alltreatments)
adonis2

#looking at difference between conditions
library(pairwiseAdonis)
pairwise.adonis(matrix_use2, Alltreatments$Condition)
```
## plot Figure 3D
```{r warning=FALSE}
p1=ggplot(data=SV_16S_use1_nocontrols, aes(NMDS1, NMDS2,
color=Condition))+geom_point(size=3)+xlab("NMDS1")+ylab("NMDS2")+theme_classic()+facet_wrap(~Sample_type, scales ="free", ncol = 3)+scale_color_manual(values=tax_colors_16S)+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))	+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=8, face="bold"))+ labs(shape = "Sample Type", color = "Condition")+ theme(strip.text.x = element_text(size=10, face = "bold")) 
p1
```

```{r}
# Permanova - Distance based multivariate analysis of variance
##Effect of SynCom and sample type
Alltreatments=SV_16S_use1[c(2:10)] 
head(Alltreatments)
#Adonis function in Vegan
adonis <- adonis(matrix_use ~ Condition*Sample_type, method = "bray", permutations = 999, data=Alltreatments)
adonis

#looking at difference between conditions
library(pairwiseAdonis)
pairwise.adonis(matrix_use, Alltreatments$Condition)
```


# Figure 3B Analysis beta dispersion on seeds and seedlings separately
```{r}
#Seeds
SV_16S_use1_seeds=subset(SV_16S_use1, Sample_type=="Seed")
Alltreatments=SV_16S_use1_seeds[c(2:10)] 
head(Alltreatments)

matrix_seed<-SV_16S_use1_seeds[c(11:ncol(SV_16S_use1_seeds))]
matrix_seed<-matrix_seed[,colSums(matrix_seed)>=1]
## Calculate multivariate dispersions
dis <- vegdist(matrix_seed)
mod <- betadisper(dis, Alltreatments$Condition)
mod
## Perform test
anova(mod)

## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
mod.HSD <- TukeyHSD(mod)
plot(mod.HSD)

## Plot the groups and distances to centroids 
plot(mod)

df <- data.frame(Distance_to_centroid=mod$distances,Group=mod$group)
groups <- mod$group

tax_colors_16S <-  c('Control'='#ffbb94',"6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")
p<- ggplot(data=df,aes(x=Group,y=Distance_to_centroid,fill=groups))
p <-p + geom_boxplot()
p <-p +theme_classic()+scale_fill_manual(values=tax_colors_16S)+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))	+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+ labs( color = "Condition")+ xlab("Condition")+ylab("Beta-dispersion: Distance to centroid")+theme(legend.position = "none") +ggtitle("Seeds")+theme(plot.title = element_text(hjust = 0.5, face="bold"))
p
```
# Figure 3C
```{r}
#Seedlings
SV_16S_use1_seedling=subset(SV_16S_use1, Sample_type=="Seedling")
Alltreatments=SV_16S_use1_seedling[c(2:10)] 
head(Alltreatments)

matrix_seedling<-SV_16S_use1_seedling[c(11:ncol(SV_16S_use1_seedling))]
matrix_seedling<-matrix_seedling[,colSums(matrix_seedling)>=1]
## Calculate multivariate dispersions
dis <- vegdist(matrix_seedling)
mod <- betadisper(dis, Alltreatments$Condition)
mod
## Perform test
anova(mod)

## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
mod.HSD <- TukeyHSD(mod)
plot(mod.HSD)

## Plot the groups and distances to centroids 
plot(mod)

df <- data.frame(Distance_to_centroid=mod$distances,Group=mod$group)
groups <- mod$group

tax_colors_16S <-  c('Control'='#ffbb94',"6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")
p<- ggplot(data=df,aes(x=Group,y=Distance_to_centroid,fill=groups))
p <-p + geom_boxplot()
p <-p +theme_classic()+scale_fill_manual(values=tax_colors_16S)+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))	+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+ labs( color = "Condition")+ xlab("Condition")+ylab("Beta-dispersion: Distance to centroid")+theme(legend.position = "none") +ggtitle("Seedlings")+theme(plot.title = element_text(hjust = 0.5, face="bold"))
p
```

# Figure 4 
```{r}
library(microbiome)
##Convert data as phyloseq object
OTU2 = otu_table(matrix_use, taxa_are_rows = FALSE)
taxo_gyrB2 <- read.table(file="SynCom_Anne_gyrB_SV_Taxo.txt", sep='\t', header=TRUE,check.names=FALSE,row.names=1) 
dim(taxo_gyrB2)
taxo_gyrB2=as.matrix(taxo_gyrB2)
TAXO_gyrB = tax_table(taxo_gyrB2)
meta=SV_16S_use1[c(1:10)]
META=sample_data(meta)
physeq_gyrB = phyloseq(OTU2,META,TAXO_gyrB)
physeq_gyrB
```

#### ASV of SynComs 
```{r}
library(tidyverse)
dataframe_taxo_gyrB2 <- as.data.frame(taxo_gyrB2, raw.names = FALSE)
dataframe_taxo_gyrB2 <- rownames_to_column(dataframe_taxo_gyrB2, var = "ASV")%>%as_tibble()
ASV_SC <- subset(dataframe_taxo_gyrB2, SynCom == "yes")
```

#### Colors used in graphs
```{r}
tax_colors <-  c(
                 "Erwinia cancerogenus"="purple",
                 "Erwinia persicina"="#bc5090", 
                 "Paenibacillus sp"= "#ff7c43", 
                 "Pantoea agglomerans"="#216fd2",
                 'Plantibacter sp'='lightblue' ,
                 "Pseudomonas fluorescens 1"="#488f31", 
                 "Pseudomonas fluorescens 2"="#b7d7a7", 
                 "Pseudomonas fluorescens 3"="#004900", 
                 "Pseudomonas fluorescens 4"="#aebd38", 
                 "Pseudomonas viridiflava"= "#f9b0d3",
                 'Stenotrophomonas rhizophila'='#ffa600',
                 'Xanthomonas campestris'='red',
                 "Other taxa"="lightgrey" 
                 )


treatment_colors <-  c('Control'='#ffbb94',"6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")
```


# Figure 4B Taxonomic profile of the inoculum, seeds and seedlings
```{r}
library(ggplot2)
#relative abundance
count_to_prop <- function(x) {return( x / sum(x))}
physeq_gyrB_abrel <- transform_sample_counts(physeq_gyrB, count_to_prop)
sample_sums(physeq_gyrB_abrel)

#phyloseq to data.frame 
otutable <- otu_table(physeq_gyrB_abrel)
gyrB_abrel <- as.data.frame(otutable)
gyrB_abrel$sample <- sample_data(physeq_gyrB_abrel)$Sample_ID
gyrB_abrel <- gyrB_abrel %>% pivot_longer(-sample) 
colnames(gyrB_abrel) <- c("Sample_ID", "ASV", "relative.abundance")
head(gyrB_abrel)

data_SC=merge(x=gyrB_abrel,y=meta,by.x="Sample_ID",by.y="Sample_ID",all.x = TRUE,all.y = TRUE)
data_SC=merge(x=data_SC,y=dataframe_taxo_gyrB2,by.x="ASV",by.y="ASV",all.x = TRUE,all.y = TRUE, na.rm=TRUE)

data_SC$Strains <- factor(data_SC$Strains, labels =  c("Other taxa", "Erwinia cancerogenus" ,       "Erwinia persicina"   ,        "Paenibacillus sp"        ,    "Pantoea agglomerans"       ,  "Plantibacter sp"         ,    "Pseudomonas fluorescens 1",   "Pseudomonas fluorescens 2" , "Pseudomonas fluorescens 3" ,  "Pseudomonas fluorescens 4"  , "Pseudomonas viridiflava"  ,   "Stenotrophomonas rhizophila", "Xanthomonas campestris"))


#rename sample-type
data_SC$Sample_type <- factor(data_SC$Sample_type, labels = c("Inoculum", "Seeds", "Seedlings"))

#reorder Strains factors alphabetically 
old.lvl <- levels(data_SC$Strains)
data_SC$Strains <- factor(data_SC$Strains, levels=c("Other taxa", sort(old.lvl[old.lvl!="Other taxa"], decreasing=F)))

#subset for graph
data_SC6 <- subset(data_SC, Condition == "6 strains")
data_SC8 <- subset(data_SC, Condition == "8 strains")
data_SC12 <- subset(data_SC, Condition == "12 strains")
data_control <- subset(data_SC, Condition == "Control")

#graph
plot_SC6<- ggplot(data_SC6, aes(x=Sample_ID, y = round(relative.abundance, 4)))+
  geom_col(aes(fill=Strains))+ 
  scale_fill_manual(values = tax_colors)+ 
  #facet_grid(~Sample_type, scales = "free", space = "free")+
  facet_grid(~Sample_type, scales = "free_x")+
  theme_classic()+
  ggtitle("6 Strains") +  
  xlab("") + 
  ylab("Relative Abundance")+ 
  labs(fill="")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.title.y = element_text(size=25, face = "bold"),
        plot.title = element_text(size=25, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
       # legend.text = element_blank(),
        strip.text =  element_text(size = 25),
       strip.background = element_blank(),
       legend.text = element_text(size = 24, face = "italic"))


plot_SC8<- ggplot(data_SC8, aes(x=Sample_ID, y = round(relative.abundance, 4)))+
  geom_col(aes(fill=Strains))+ 
  scale_fill_manual(values = tax_colors)+ 
  #facet_grid(~Sample_type, scales = "free", space = "free")+
  facet_grid(~Sample_type, scales = "free_x")+
  theme_classic()+
  ggtitle("8 Strains") +  
  xlab("") + 
  ylab("Relative Abundance")+ 
  labs(fill="")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.title.y = element_text(size=25, face = "bold"),
        plot.title = element_text(size=25, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
       # legend.text = element_blank(),
        strip.text =  element_text(size = 25),
       strip.background = element_blank(),
       legend.text = element_text(size = 24, face = "italic"))

plot_SC12<- ggplot(data_SC12, aes(x=Sample_ID, y = round(relative.abundance, 4)))+
  geom_col(aes(fill=Strains))+ 
  scale_fill_manual(values = tax_colors)+ 
  #facet_grid(~Sample_type, scales = "free", space = "free")+
  facet_grid(~Sample_type, scales = "free_x")+
  theme_classic()+
  ggtitle("12 Strains") +  
  xlab("") + 
  ylab("Relative Abundance")+ 
  labs(fill="")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.title.y = element_text(size=25, face = "bold"),
        plot.title = element_text(size=25, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
       # legend.text = element_blank(),
        strip.text = element_text(size = 25),
       strip.background = element_blank(),
       legend.text = element_text(size = 24, face = "italic"))

data_control <- add_row(data_control)
data_control[18817,5] <- "Inoculum"
data_control[18817,6] <- "Control"
plot_control<- ggplot(data_control, aes(x=Sample_ID, y = round(relative.abundance, 4)))+
  geom_col(aes(fill=Strains))+ 
  scale_fill_manual(values = tax_colors)+ 
  #facet_grid(~Sample_type, scales = "free", space = "free")+
  facet_grid(~Sample_type, scales = "free_x")+
  theme_classic()+
  ggtitle("Control") +  
  xlab("") + 
  ylab("Relative Abundance")+ 
  labs(fill="")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.title.y = element_text(size=25, face = "bold"),
        plot.title = element_text(size=25, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
       # legend.text = element_blank(),
        strip.text = element_text(size = 25),
       strip.background = element_blank(),
       legend.text = element_text(size = 24, face = "italic"))


library(ggpubr)

plot <- ggarrange(plot_control, plot_SC6, plot_SC8, plot_SC12, ncol=1, nrow=4, common.legend = TRUE, legend = "right")
ggsave("plot.png" , width = 40, height = 40, units = "cm")
```


# Figure 4A - Alpha-diversity (Observed richness)
```{r, echo=FALSE, warning=FALSE, message=FALSE, results="hide"}
library(phyloseq)
physeq_SC <- subset_taxa(physeq_gyrB, SynCom == "yes")

#Table of alpha-diversity estimators
table_rgyrB <- estimate_richness(physeq_SC, split = TRUE, measures=c("Observed", "Shannon"))
#Add evenness
HgyrB <- table_rgyrB$Shannon
S1gyrB <- table_rgyrB$Observed
SgyrB <- log(S1gyrB)
eve_gyrB <- HgyrB/SgyrB
table_rgyrB$Evenness <- eve_gyrB
#Bind sampledata + table of alpha_div (here R1 only)
datargyrB <- cbind(sample_data(physeq_SC), table_rgyrB)
#Reorder manually strains variables
datargyrB$Condition <- factor(datargyrB$Condition, levels=c("Control", "6 strains", "8 strains", "12 strains"))

#Plot

library(ggplot2); packageVersion("ggplot2")

treatment_colors <-  c('Control'='#ffbb94',"6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")

p_gyrB_ric <- ggplot(data=datargyrB, 
                             aes_string(x='Condition',y='Observed', fill='Condition')) + 
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~Sample_type, scales = "free_y") +
  geom_jitter(aes(colour = Condition), size=2, alpha=0.3) +
  scale_color_manual(values = treatment_colors)+
  scale_fill_manual(values = treatment_colors)+ 
  theme_classic()+
  #theme(axis.text.x = element_text(angle = 90)) +
  theme(strip.background = element_blank(), strip.text =  element_text(face = "bold", size = 12))+
  ggtitle("A - Richness")+
  xlab("") + 
  ylab("ASV Richness")+
  scale_y_continuous(limits = c(0,12.5),breaks = seq(0,12,1))+
  theme(axis.text.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold", size = 11))	+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))
p_gyrB_ric


ggsave("richness.png" , width = 30, height = 10, units = "cm")
```

# Figure 5 - Transmission Seed-seedling

```{r}
library(microbiome)
##Convert data as phyloseq object
OTU2 = otu_table(matrix_use, taxa_are_rows = FALSE)
meta=SV_16S_use1[c(1:10)]
META=sample_data(meta)
physeq_OTU = phyloseq(OTU2,META)
physeq_OTU
```

```{r, echo=FALSE, message=FALSE, warning=FALSE}
#Transform into relative abundance
physeq_OTU_rel <- transform_sample_counts(physeq_OTU, function(x) x / sum(x) )
#Extract OTU table
physeq_OTU_rel_table <- data.frame(otu_table(physeq_OTU_rel))
physeq_OTU_rel_table <- cbind(Sample = rownames(physeq_OTU_rel_table), physeq_OTU_rel_table)
row.names(physeq_OTU_rel_table) <- NULL
physeq_OTU_rel_table$Sample <- as.factor(physeq_OTU_rel_table$Sample)
head(physeq_OTU_rel_table)
#Extract metadata
META2 <- data.frame(sample_data(physeq_OTU_rel))
META2 <- cbind(Sample = rownames(META2), META2)
row.names(META2) <- NULL
META2$Sample_ID <- as.factor(META2$Sample_ID)
```


```{r}
### Transform table matrix into long table format
library(reshape2)
table_rel_abund_long=setNames(melt(physeq_OTU_rel_table), c('Sample', 'ASV', 'Relative_Abundance'))
head(table_rel_abund_long)
dim(table_rel_abund_long)
## merge with metadata
table_rel_abund_long_meta<-merge(META2,table_rel_abund_long,by="Sample")
dim(table_rel_abund_long_meta)
head(table_rel_abund_long_meta)
##merge with taxonomic info
taxo<-read.table("SynCom_Anne_gyrB_SV_Taxo.txt", header=TRUE, check.names = FALSE, sep = "\t")
table_rel_abund_long_meta<-merge(taxo,table_rel_abund_long_meta,by="ASV")
dim(table_rel_abund_long_meta)
head(table_rel_abund_long_meta)


###filter rows with Relative_Abundance= 0 
library(dplyr)
table_rel_abund_long_meta2 <- filter(table_rel_abund_long_meta, Relative_Abundance > 0)
head(table_rel_abund_long_meta2)
dim(table_rel_abund_long_meta2)
```


# Figure 5: Transmission
```{r}
##Keep only inoculated ASVs
table_rel_abund_long_meta2_strain=subset(table_rel_abund_long_meta2, SynCom=="yes")
dim(table_rel_abund_long_meta2_strain)

#Average rel abund by condition and compartment
library(Rmisc)
library(ggplot2)
grouped_stat_slope_rich=summarySE(table_rel_abund_long_meta2_strain, measurevar=c("Relative_Abundance"),groupvars=c("Condition","Sample_type", "Strains"), na.rm = TRUE)
```
# Figure 5C
```{r}
#slop plot by strains - seed vs seedling - without controls - with standard errors
seedVSseedlings_rich=grouped_stat_slope_rich[grouped_stat_slope_rich$Sample_type!="Inoculum",]
seedVSseedlings_rich=seedVSseedlings_rich[seedVSseedlings_rich$Condition!="Control",]

#remove rows of strains that have not been inoculated in SynComs
seedVSseedlings_rich2seedling <- subset(seedVSseedlings_rich, Sample_type=="Seedling" & N>16)
seedVSseedlings_rich2seed <- subset(seedVSseedlings_rich, Sample_type=="Seed")
seedVSseedlings_rich_final=rbind(seedVSseedlings_rich2seed,seedVSseedlings_rich2seedling)

#ordered by response type
seedVSseedlings_rich_final$Strains<-ordered(seedVSseedlings_rich_final$Strains, levels=c('Stenotrophomonas rhizophila',"Pseudomonas viridiflava","Paenibacillus sp", "Pseudomonas fluorescens 1","Pantoea agglomerans", "Pseudomonas fluorescens 4","Erwinia persicina","Enterobacter cancerogenus", "Pseudomonas fluorescens 3",'Xanthomonas campestris','Plantibacter sp' , "Pseudomonas fluorescens 2"))
tax_colors_16S <-  c('Control'='#ffbb94',"6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")
seedVSseedlings_rich_final$Condition<-ordered(seedVSseedlings_rich_final$Condition, levels=c("6 strains", "8 strains", "12 strains"))
plot_strains2=ggplot(seedVSseedlings_rich_final, aes(x = Sample_type, y = Relative_Abundance, color = Condition,group=Condition)) + geom_point() + geom_line() +  facet_wrap(~Strains)+theme_classic()+ylab("Relative abundance (%)")+xlab("Sample Type")+scale_color_manual(values=tax_colors_16S)+ scale_y_log10(labels = scales::percent_format(accuracy = 0.1))+ theme(strip.text.x = element_text(size=10,  face = "bold")) + theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))	+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold")) + geom_errorbar(data=seedVSseedlings_rich_final, aes(x=Sample_type, ymin=Relative_Abundance-se, ymax=Relative_Abundance+se, color=Condition), width=.1, alpha=0.7) 
plot_strains2
```

# Figure 5A - slope inoculum - seed
```{r}
#slop plot by strains - seed vs seedling - without controls - with standard errors
inocVSseed_rich=grouped_stat_slope_rich[grouped_stat_slope_rich$Sample_type!="Seedling",]
inocVSseed_rich=inocVSseed_rich[inocVSseed_rich$Condition!="Control",]

#ordered by response type
inocVSseed_rich$Strains<-ordered(inocVSseed_rich$Strains, levels=c('Stenotrophomonas rhizophila',"Pseudomonas viridiflava","Paenibacillus sp", "Pseudomonas fluorescens 1","Pantoea agglomerans", "Pseudomonas fluorescens 4","Erwinia persicina","Enterobacter cancerogenus", "Pseudomonas fluorescens 3",'Xanthomonas campestris','Plantibacter sp' , "Pseudomonas fluorescens 2"))
tax_colors_16S <-  c('Control'='#ffbb94',"6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")
plot_strains2=ggplot(inocVSseed_rich, aes(x = Sample_type, y = Relative_Abundance, color = Condition,group=Condition)) + geom_point() + geom_line() +  facet_wrap(~Strains)+theme_classic()+ylab("Relative abundance (%)")+xlab("Sample Type")+scale_color_manual(values=tax_colors_16S)+ scale_y_log10(labels = scales::percent_format(accuracy = 0.1))+ theme(strip.text.x = element_text(size=10,  face = "bold")) + theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))	+ theme(axis.title = element_text(color="black", size=11, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold")) + geom_errorbar(data=inocVSseed_rich, aes(x=Sample_type, ymin=Relative_Abundance-se, ymax=Relative_Abundance+se, color=Condition), width=.1, alpha=0.7) 
plot_strains2
```
# Figure 5B Ratio abundance Seed / Inoculum
```{r, echo=FALSE, message=FALSE, warning=FALSE}
PS_run2_rich_no0=transform_sample_counts(physeq_OTU, function(x) x+1 )
PS_run2_rich_no0relative=transform_sample_counts(PS_run2_rich_no0, function(x) x/sum(x)*100)
data_fitness=as.data.frame(otu_table(PS_run2_rich_no0relative))
data_fitness2 <- tibble::rownames_to_column(data_fitness, "Sample_ID")
data_fitness3=pivot_longer(data_fitness2,cols=colnames(data_fitness),names_to = "ASV")
colnames(data_fitness3)=c("Sample_ID","ASV","Relative_abundance")
data_fitness3=as.data.frame(data_fitness3)
data_fitness4=left_join(data_fitness3, taxo,by="ASV")
#data_fitness4$Strain_Figure_name = data_fitness4$Strain_Figure_name %>% replace_na('Native strain')
#keep only strains of SynComs
data_fitness4b <- subset(data_fitness4, SynCom == "yes")
sample_data_fitness=as.data.frame(sample_data(physeq_OTU))
sample_data_fitness$Sample_ID=row.names(sample_data_fitness)
data_fitness5=left_join(data_fitness4b,sample_data_fitness,by="Sample_ID")
#####
data_fitness5_inoc=data_fitness5[data_fitness5$Sample_type=="Inoculum",]
names(data_fitness5_inoc)[names(data_fitness5_inoc) == "Relative_abundance"] <- "Relative_abundance_inoc"
data_fitness5_inoc_stat=summarySE(data_fitness5_inoc,measurevar=c("Relative_abundance_inoc"), groupvars=c("ASV","Condition","Strains"), na.rm = TRUE)
data_fitness5_seed=data_fitness5[data_fitness5$Sample_type=="Seed",]
names(data_fitness5_seed)[names(data_fitness5_seed) == "Relative_abundance"] <- "Relative_abundance_seed"
data_fitness5_seed_stat=summarySE(data_fitness5_seed,measurevar=c("Relative_abundance_seed"), groupvars=c("ASV","Strains","Condition"), na.rm = TRUE)
ALL_data_seed_inocs_rich=full_join(data_fitness5_inoc_stat,data_fitness5_seed_stat,by=c("ASV","Condition"))

ALL_data_seed_inocs_rich$fitness=ALL_data_seed_inocs_rich$Relative_abundance_seed/ALL_data_seed_inocs_rich$Relative_abundance_inoc
names(ALL_data_seed_inocs_rich)[names(ALL_data_seed_inocs_rich) == "Strains.x"] <- "Strains"
ALL_data_seed_inocs_rich_nocontrol <- subset(ALL_data_seed_inocs_rich, Condition != "Control")
#remove rows of strains that have not been inoculated in SynComs
ALL_data_seed_inocs_rich_nocontrol <- subset(ALL_data_seed_inocs_rich_nocontrol, sd.y != 0)
#ALL_data_seed_inocs_rich_meta=merge(taxo,ALL_data_seed_inocs_rich,by="Strains",all.x = T)
treatment_colors <-  c('Control'='#ffbb94',"6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")
ALL_data_seed_inocs_rich_nocontrol$Condition<-ordered(ALL_data_seed_inocs_rich_nocontrol$Condition, levels=c("6 strains", "8 strains", "12 strains"))

plot_strains_fitness=ggplot(ALL_data_seed_inocs_rich_nocontrol,aes(x=log10(fitness),y=Strains,color=Condition))+geom_point(size=5)+scale_color_manual(values=treatment_colors)+theme_classic()+theme(text = element_text(size = 20))+guides(fill=guide_legend(ncol=1))+ylab("Strains")+xlab("log10(Relative abundance ratio Seed / Inoculum)")+geom_vline(xintercept = 0)
plot_strains_fitness
#write.csv(ALL_data_seed_inocs_rich_nocontrol, "SynCom_radis_strain_fitness_inoc.csv")
```


# Figure 5D Ratio abundance Seedling / Seed
```{r, echo=FALSE, message=FALSE, warning=FALSE}
PS_run2_rich_no0=transform_sample_counts(physeq_OTU, function(x) x+1 )
PS_run2_rich_no0relative=transform_sample_counts(PS_run2_rich_no0, function(x) x/sum(x)*100)
data_fitness=as.data.frame(otu_table(PS_run2_rich_no0relative))
data_fitness2 <- tibble::rownames_to_column(data_fitness, "Sample_ID")
data_fitness3=pivot_longer(data_fitness2,cols=colnames(data_fitness),names_to = "ASV")
colnames(data_fitness3)=c("Sample_ID","ASV","Relative_abundance")
data_fitness3=as.data.frame(data_fitness3)
data_fitness4=left_join(data_fitness3, taxo,by="ASV")
#data_fitness4$Strain_Figure_name = data_fitness4$Strain_Figure_name %>% replace_na('Native strain')
#keep only strains of SynComs
data_fitness4b <- subset(data_fitness4, SynCom == "yes")
sample_data_fitness=as.data.frame(sample_data(physeq_OTU))
sample_data_fitness$Sample_ID=row.names(sample_data_fitness)
data_fitness5=left_join(data_fitness4b,sample_data_fitness,by="Sample_ID")
#####
data_fitness5_seedling=data_fitness5[data_fitness5$Sample_type=="Seedling",]
names(data_fitness5_seedling)[names(data_fitness5_seedling) == "Relative_abundance"] <- "Relative_abundance_seedling"
data_fitness5_seedling_stat=summarySE(data_fitness5_seedling,measurevar=c("Relative_abundance_seedling"), groupvars=c("ASV","Condition","Strains"), na.rm = TRUE)
data_fitness5_seed=data_fitness5[data_fitness5$Sample_type=="Seed",]
names(data_fitness5_seed)[names(data_fitness5_seed) == "Relative_abundance"] <- "Relative_abundance_seed"
data_fitness5_seed_stat=summarySE(data_fitness5_seed,measurevar=c("Relative_abundance_seed"), groupvars=c("ASV","Strains","Condition"), na.rm = TRUE)
ALL_data_seed_seedlings_rich=full_join(data_fitness5_seedling_stat,data_fitness5_seed_stat,by=c("ASV","Condition"))

ALL_data_seed_seedlings_rich$fitness=ALL_data_seed_seedlings_rich$Relative_abundance_seedling/ALL_data_seed_seedlings_rich$Relative_abundance_seed
names(ALL_data_seed_seedlings_rich)[names(ALL_data_seed_seedlings_rich) == "Strains.x"] <- "Strains"
ALL_data_seed_seedlings_rich_nocontrol <- subset(ALL_data_seed_seedlings_rich, Condition != "Control")

#remove rows of strains that have not been inoculated in SynComs
ALL_data_seed_seedlings_rich_nocontrol <- subset(ALL_data_seed_seedlings_rich_nocontrol, sd.y != 0)

#ALL_data_seed_seedlings_rich_meta=merge(taxo,ALL_data_seed_seedlings_rich,by="Strains",all.x = T)
treatment_colors <-  c('Control'='#ffbb94',"6 strains"="ivory3", "8 strains"="ivory4", "12 strains"="black", "Control"="white")
ALL_data_seed_seedlings_rich_nocontrol$Condition<-ordered(ALL_data_seed_seedlings_rich_nocontrol$Condition, levels=c("6 strains", "8 strains", "12 strains"))

plot_strains_fitness=ggplot(ALL_data_seed_seedlings_rich_nocontrol,aes(x=log10(fitness),y=Strains,color=Condition))+geom_point(size=5)+scale_color_manual(values=treatment_colors)+theme_classic()+theme(text = element_text(size = 20))+guides(fill=guide_legend(ncol=1))+ylab("Strains")+xlab("log10(Relative abundance ratio Seedling / Seed)")+geom_vline(xintercept = 0)
plot_strains_fitness
#write.csv(ALL_data_seed_seedlings_rich_nocontrol, "SynCom_radis_strain_fitness.csv")
```


# Figure 6 Effect on germination and seedling phenotypes
```{r}
pheno=read.table("Phenotype_Syncom_strain.txt", header = T, sep = "\t")
head(pheno)
dim(pheno)

```
## Graph Figure 6 phenotype
```{r}
pheno$Phenotype<-ordered(pheno$Phenotype, levels=c("Non Germinated","Abnormal", "Normal"))

pheno$Condition<-ordered(pheno$Condition, levels=c("6 strains", "8 strains", "12 strains","Enterobacter cancerogenus","Erwinia persicina","Paenibacillus sp", "Pantoea agglomerans","Plantibacter sp","Pseudomonas fluorescens 1", "Pseudomonas fluorescens 2","Pseudomonas fluorescens 3","Pseudomonas fluorescens 4","Pseudomonas viridiflava","Stenotrophomonas rhizophila","Xanthomonas campestris", "Control"))
	
library(ggplot2)
Figure6=ggplot(pheno, aes(x=Condition, y=Frequency, fill=Phenotype)) + 
  geom_bar(position = "stack", stat="identity") + ylab("Frequency of phenotype")+xlab("Condition")+scale_fill_manual(values=c("darkred","#E69F00","#41AB5D"))+ theme(axis.title = element_text(color="black", size=12, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+ scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ theme_classic()	+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))+ coord_flip()+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+facet_grid(Type~., scales  = "free", space = "free")+ theme(strip.text.y = element_text(size=9, face = "bold")) + theme(legend.position = c(-0.26, 0.09)) + theme(legend.background = element_rect(fill="white", size=0.5,linetype="solid", colour ="black"))
Figure6
```

# Figure 7- SynCom effect on phenotype

# Figure 7 effect of phenotype on microbiota structure
```{r warning=FALSE}
##subset seedling to look at seedling type effect
SV_16S_use1_seedling=subset(SV_16S_use1, Sample_type=="Seedling")

#without controls
SV_16S_use1_seedlings_nocontrols=subset(SV_16S_use1_seedling, Condition!="Control")
SV_16S_use1_seedlings_nocontrols$Seedling_type<-ordered(SV_16S_use1_seedlings_nocontrols$Seedling_type, levels=c("PN", "PA"))

matrix3<-SV_16S_use1_seedlings_nocontrols[c(11:4943)]
dim(matrix3)
head(matrix3)

#create rarefied matrix
matrix4=rrarefy(matrix3, raremax)
#verify that now we have the same number of sequence per sample (Sample size) across the dataset
dim(matrix4)
matrix_use2<-matrix4[,colSums(matrix4)>=1]
dim(matrix_use2)

##Prepare an ordination (NMDS) based on Bray-Curtis dissimilarity matrix to analyze the community structure
NMDS2 <- metaMDS(matrix_use2, distance = "bray", trymax = 100)
NMDS2
stressplot(NMDS2)

NMDSsites2=scores(NMDS2, display="sites")
SV_16S_use1_seedlings_nocontrols=cbind(SV_16S_use1_seedlings_nocontrols,NMDSsites2)
```
# Figure 7A
```{r warning=FALSE}
library(ggplot2)
SV_16S_use1_seedlings_nocontrols$Condition<-ordered(SV_16S_use1_seedlings_nocontrols$Condition, levels=c( "6 strains", "8 strains", "12 strains"))
p1=ggplot(data=SV_16S_use1_seedlings_nocontrols, aes(NMDS1, NMDS2,
color=Seedling_type))+geom_point(size=3)+xlab("NMDS1")+ylab("NMDS2")+facet_wrap(~Condition, scales ="free", ncol = 2)+theme_classic()+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))	+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=8, face="bold"))+ labs(color = "Seedling Phenotype")+ theme(strip.text.x = element_text(size=10, face = "bold")) +
  scale_color_manual(labels = c("Normal", "Abnormal"), values = c("#41AB5D", "#E69F00"))+ theme(legend.position = c(0.65, 0.3))
p1
```

```{r}
# Permanova - Distance based multivariate analysis of variance
##Effect of SynCom and sample type
Alltreatments=SV_16S_use1_seedlings_nocontrols[c(2:10)] 
head(Alltreatments)
#Adonis function in Vegan
adonis <- adonis(matrix_use2 ~ Condition*Seedling_type, method = "bray", permutations = 999, data=Alltreatments)
adonis

#looking at difference between conditions
library(pairwiseAdonis)
pairwise.adonis(matrix_use2, Alltreatments$Condition)

```

# Figure 7B rel abund strains in function of seedling type
```{r warning=FALSE}
table_rel_abund_long_meta2_seedling=subset(table_rel_abund_long_meta2, Sample_type=="Seedling")
table_rel_abund_long_meta2_seedling=table_rel_abund_long_meta2_seedling[table_rel_abund_long_meta2_seedling$Condition!="Control",]
library(Rmisc)
library(ggplot2)
Diversity_stat <- summarySE(table_rel_abund_long_meta2_seedling, measurevar="Relative_Abundance", groupvars=c("Seedling_type", "Condition", "Strains"), na.rm = TRUE)
```

## Figure 7B rel abund strains in function of seedling type - just 12 strains SynCom
```{r warning=FALSE}
table_rel_abund_long_meta2_seedling12strains=table_rel_abund_long_meta2_seedling[table_rel_abund_long_meta2_seedling$Condition=="12 strains",]
dim(table_rel_abund_long_meta2_seedling12strains)

##Remove empty rows qith ASVs other than the inoculated strains
table_rel_abund_long_meta2_seedling12strains <- table_rel_abund_long_meta2_seedling12strains[-which(table_rel_abund_long_meta2_seedling12strains$Strains == ""), ]
dim(table_rel_abund_long_meta2_seedling12strains)

table_rel_abund_long_meta2_seedling12strains$Seedling_type<-ordered(table_rel_abund_long_meta2_seedling12strains$Seedling_type, levels=c("PN","PA"))

p3=ggplot(data=table_rel_abund_long_meta2_seedling12strains, aes(x=Strains, y=Relative_Abundance, fill=Seedling_type)) +geom_boxplot(outlier.shape = NA) + theme_classic()+xlab("Strains")+ylab("ASV Relative Abundance")+ theme(axis.title = element_text(color="black", size=12, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+ scale_y_log10(labels = scales::percent_format(accuracy = 0.1))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))	+ theme(legend.text = element_text(color="black", size=8, face="bold"))+ theme(legend.title = element_text(color="black", size=10, face="bold"))+ scale_fill_manual(name = "Seedling Phenotype", labels = c("Normal", "Abnormal"), values = c("#41AB5D", "#E69F00"))
p3
```


```{r}
##stats rel abund in different phenotypes
#choose the strain
strain=subset(table_rel_abund_long_meta2_seedling12strains, Strains=="Xanthomonas campestris")
dim(strain)
ggplot(strain, aes(Relative_Abundance))+geom_histogram(binwidth=0.001)
```

```{r}
#Gaussian or Gamma (continuous data)
library(lme4)
library(MASS)
model.nb<-glm(Relative_Abundance ~ Seedling_type,  data=strain, family = Gamma(link="log"))
summary(model.nb)
#check if model fits the data P>0.05 is good
1- pchisq(summary(model.nb)$deviance, summary(model.nb)$df.residual)
plot(model.nb)
```

```{r}
###Estimate P values for main effects and interactions
library(car)
Anova(model.nb)

fit <- aov(Relative_Abundance ~ Seedling_type, strain)# analyse de variance
summary(fit)
#posthoc
library(multcompView)
library(emmeans)
warp.emm <- emmeans(model.nb, ~Seedling_type)
multcomp::cld(warp.emm)
```

```{r}
#stats on all strains
ggplot(table_rel_abund_long_meta2_seedling12strains, aes(Relative_Abundance))+geom_histogram(binwidth=0.005)
```

```{r}
#Gaussian or Gamma (continuous data)
library(lme4)
library(MASS)
model.nb<-glm(Relative_Abundance ~ Strains*Seedling_type,  data=table_rel_abund_long_meta2_seedling12strains, family = Gamma(link="log"))
summary(model.nb)
#check if model fits the data P>0.05 is good
1- pchisq(summary(model.nb)$deviance, summary(model.nb)$df.residual)
plot(model.nb)
```

```{r}
###Estimate P values for main effects and interactions
library(car)
Anova(model.nb)

fit <- aov(Relative_Abundance ~ Strains*Seedling_type, table_rel_abund_long_meta2_seedling12strains)# analyse de variance
summary(fit)
#posthoc
library(multcompView)
library(emmeans)
warp.emm <- emmeans(model.nb, ~Seedling_type|Strains)
multcomp::cld(warp.emm)
```


