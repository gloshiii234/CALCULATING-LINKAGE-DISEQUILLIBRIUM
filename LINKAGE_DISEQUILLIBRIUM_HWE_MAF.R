
##THE REQUIREMENT IS TO IDENTIFY THE LINKAGE DISEQUILLIBRIUM BETWEEN SNPS IN ANKR GENE


##LOADING REQUIRED LIBRARIES

library(genetics)
library(haplo.stats)
library(dplyr)
library(snpStats)
library(tidyverse)
require(RColorBrewer)
library(tidyr)
library(magrittr)


##STEP 0NE -LOADING THE FAMUSS DATASET

famuss <- read.csv("C:/Users/msaim/Downloads/famuss(1).csv")
##viewing the first 10 rows of the dataset

head(famuss,10)
##This is a partial dataset , i will be downloading the full dataset below

fam <- read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt",header = T, 

                  ##Viewing the first 5 rows and columns of the dataset                 sep = "\t")
fam[1:5,1:5]

dim(fam)

##selecting out the ankr gene from the dataset
 ankr_sel <- fam[, grepl("^ankr",names(fam))]
 head(ankr_sel)

 
##Creating a genotype object from the selected genes
 
 geno_ankr <- lapply(ankr_sel,genotype, sep="")
 geno_ankr

 ##Removing NA's from each element of the list
 
 
geno_ankr_fil <- lapply(geno_ankr, function(x) x[complete.cases(x)])
 
 # Print the modified list
 head(geno_ankr_fil)

 
 
 ##CONVERTING OUR GENOTYPE OBJECT TO A DATAFRAME
 
 geno_df <- as.data.frame(geno_ankr)
 head(geno_df)
 
 ##CALCULATING D'
 matrix_ankr = LD(geno_df)$`D'`
 head(matrix_ankr)
 
 ##Generating a heatmap to show LD of SNPS
 library(pheatmap)
 pheatmap(matrix_ankr, cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE)
 
 
 ##Calculating mean LD
 
 mean_ld_ankr = mean(matrix_ankr,na.rm = TRUE)
 print(mean_ld_ankr)
 
 
 
 
 ##DETERMINING HARDY-WEINBERG EQUILLIBRIUM
 
 #Chisq test
 chisq_pvalues<- vector()
 for (col in geno_df){
   chisq <- HWE.chisq(col)
   chisq_pvalues=c(chisq_pvalues,chisq$p.value)
 }
 sort(chisq_pvalues)
 
 #naming p-values based on the respective SNPs
 names(chisq_pvalues)=colnames(geno_df)
 sort(chisq_pvalues)
 
 sum(chisq_pvalues<0.05)

 names(chisq_pvalues[chisq_pvalues<0.05])

 
 #Fisher’s Exact test
 
 
 
 exact_pvalues<- vector()
 for (col in geno_df){
   exact <- HWE.exact(col)
   exact_pvalues=c(exact_pvalues,exact$p.value)
 }
 names(exact_pvalues)=colnames(geno_df)
 sort(exact_pvalues)
 
 sum(exact_pvalues<0.05)
 names(exact_pvalues[exact_pvalues<0.05])
 
 ##ADJUSTING VALUES
 
 #set.seed(1000)
 #Adj_pval <- p.adjust(chisq_pvalues, method="bonferroni")
 #sort(Adj_pval)
 
 set.seed(42)
 adj_pval <- p.adjust(chisq_pvalues, method="bonferroni")
 sort(adj_pval)
 
 sum(adj_pval<0.05)
names(adj_pval[adj_pval<0.05]) 



##Identifying missing values

missing_geno <- as.data.frame(summary(is.na(geno_df))[3,]) #Indexing the third Row for missing data
names(missing_geno) <- c("Missing_Individuals") #Asigning the columnn name to the data
missing_geno

##We take it that our allele frequencies will remain the same despite the missing data


##CALCULATING MINOR ALLELE FREQUENCY

round(summary(geno_df$ankrd6_t233m)$"allele.freq",2)
round(summary(geno_df$ankrd6_a550t)$"allele.freq",2)
round(summary(geno_df$ankrd6_q122e)$"allele.freq",2)
##MAF gives us an idea about frequency of variants in the population






 
   