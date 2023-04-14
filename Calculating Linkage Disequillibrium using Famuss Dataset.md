This code is an example of calculation of LD betwen 4 SNPs in the Famuss dataset including ANKRD6 (Ankyrin repeat domain 6), ACTN3 (Actinin, alpha 3), LEPR (Leptin receptor) and RETN (Resistin). These genes  are among the 17genes which were using in a study of various genetic loci associated with muscle size and strength at baseline and in response to resistance training.



#FIRST STEP IS READING IN THE FAMUSS DATASET
fm <- read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt",header = T, 
                 sep = "\t")

#THE PACKAGES TO BE DOWNLOADED
Packages to be loaded
library(genetics)
library(dplyr)
library(LDheatmap)
library(snpStats)
library(pheatmap)
library(knitr)
library(purrr)
library(RColorBrewer)
library(haplo.stats)
library(janitor)

THIRD STEP- PULLING OUT THE COLUMNS WITH THE SELECTED GENES
sel_genecols <- fm[, grepl("^ankr|^actn3|^lepr|^resistin", names(fm))]  %>% drop_na()
  
#EXPLAINING THE ABOVE CODE
fm is our dataset, grepl is a function used to match column names that start with ankr,actn3...,the names function is used to retrieve the names of columns in fm so that grepl() can select those chosen. 
The %>% is a pipe operator that tells R to take the first argument and pass it to the second function,ie after grepl selects the columns, then na values should be dropped from those columns.

attach(sel_genecols)-This attaches our new dataframe into the search path.

ncol(sel_genecols)- This finds out the number of columns in the new dataframe-23

colnames(sel_genecols)-This will give you the names of the columns in the new dataframe

STEP FOUR- CALCULATING THE NUMBER OF SNPS IN THE SELECTED GENES

-The first option is identifying the unique SNPS in all the selected genes and this can be achieved as below;
##combining the selected columns into a single vector
single_vector <- unlist(sel_genecols),

##then count the number of unique SNPS
no_SNPS <- length(unique(single_vector))
no_SNPS = 12

##The result can be printed as follows
cat("The total number of SNPS is:",no_SNPS)

-The second option is counting all the SNPS in the selected columns(with duplicates) and it can be achieved as follows:

no_snps <- sum(sapply(sel_genecols, function(x) length(x)))
no_snps = 8901
cat("The total number of SNPS in all columns is:",no_snps)

##The code above applies the length function to all columns in the selected dataframe and counts all none na values and then the sum function adds up all the values obtained from all the columns. 



STEP FIVE- CONVERTING THE SELECTED GENOTYPE DATAFRAME INTO A GENOTYPE MATRIX

##Converting the dataframe of selected genes into genotype matrix

#First, apply the genotype function to all columns
geno = lapply(sel_genecols, genotype, sep"")
head(geno)
#The above code computes the genotype of each selected column and concartenates the result together without separation
because of the empty sep string.

#Second,convert geno to dataframe
geno_converted <- as.dataframe(geno)
head(geno_converted)

#Creating a matrix
matrix <- LD(geno_converted)$`D'`
#The above code is calculating the D' measure of LD which is the pairwise calculation and returns the scores in a matrix called matrix, other LD measures include r and r^2

STEP SIX-GENERATING A HEAT MAP TO VISUALIZE LD BETWEEN SNPS IN THE SELECTED GENES
##First, change the data to genotype format for easy analysis of LD
genotypes <- makeGenotypes(sel_gencols,sep"",method=as.genotype)
#Here, the function makeGenotypes from the GenABEL package is used to convert genotypes data into a format suitable for LD analysis and stores it in variable genotypes. The sep="" specifies that there is are no delimeters between the SNPs in the genotype database.

##Second is creating a color palette
color = palette(brewer.pal(n =8, name = "Reds"))
#The above code creates a color palette using the "Reds color scheme from the RColorBrewer package.



STEP SEVEN- SETTING PARAMETERS FOR THE LD HEATMAP

LDheatmap(matrix,LDmeasure="D'",color=rev(color),SNP.name=names(genotypes),text=TRUE,
title= "LD heatmap for SNPS in the four selected genes")

#The above code is generating a heatmap where matrix contains the LD scores of measure D' for all the genes and reversing the color so that it is easy to understand, with SNP.name as column names in the dataset genotypes and returning the values of LD measure on the heatmap with the given title.


STEP EIGHT- CONSTRUCTION OF HEATMAP SHOWING LD ON ANY ONE SELECTED GENE

#SELECTED GENE IN THIS EXAMPLE IS LEPR

lep_SNPs <- sel_genecols[, grepl("^lep", names(sel_genecols))]
lep_SNPs_geno <- lapply(lep_SNPs, genotype, sep="")
lep_SNPs_geno_df = as.data.frame(lep_SNPs_geno)
lep_matrix <- LD(lep_SNPs_geno_df)$`D'`
LDheatmap(lep_matrix, LDmeasure = "D'", color = rev(color), SNP.name = names(lep_SNPs_geno_df), text = TRUE,
          title = "LD heatmap for the SNPs in the Leptin receptor gene")
