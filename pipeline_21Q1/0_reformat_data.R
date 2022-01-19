##############
#### 21Q1 ####
##############

##############
# reformat dependency and mutation data
# 1) change KO gene names
# 2) subset the mutation data to CERES cell lines

##############
# load library
library(readr)
library(dplyr)
library(stringr)
##############
# set wd and define directories
setwd("C:/Users/pty01/OneDrive - Princeton University/000_ragr/000_github/SuperDendrix/pipeline_21Q1")

rawdir = "../data/21Q1/raw/"
featuredir = "../data/21Q1/features/"
profiledir = "../data/21Q1/profiles/"
rawprocesseddir = "../data/21Q1/raw_preprocessed/"

#############
# helpfer?
'%!in%' <- function(x,y)!('%in%'(x,y))


#############
# load dependency data
Achilles_gene_effect <- read_csv(paste(rawdir, "Achilles_gene_effect.csv", sep=""))
colnames(Achilles_gene_effect)[1] <- "DepMap_ID"
# change space in gene name to underscore
for (i in 1:length(colnames(Achilles_gene_effect))){
  colnames(Achilles_gene_effect)[i] <- str_replace(colnames(Achilles_gene_effect)[i], " ", "_")
}

CERES_samples <- Achilles_gene_effect$DepMap_ID

#############
# load mutation data
CCLE_mutations <- read_csv(paste(rawdir, "CCLE_mutations.csv", sep=""))
CERES_mutations <- CCLE_mutations[CCLE_mutations$DepMap_ID %in% CERES_samples,]

# change colname for OncoKB
ind <- match("Protein_Change", colnames(CERES_mutations))
colnames(CERES_mutations)[ind] <- "HGVSp_Short"


#############
# load sample data
sample_info <- read_csv(paste(rawdir, "sample_info.csv", sep=""))
CERES_sample_info <- sample_info[sample_info$DepMap_ID %in% CERES_samples,]



############
# write output files
output.file <- file(paste(rawprocesseddir, "CERES_scores.tsv", sep=""), "wb")
write.table(Achilles_gene_effect,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)

output.file <- file(paste(rawprocesseddir, "CERES_mutations.tsv", sep=""), "wb")
write.table(CERES_mutations,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)

output.file <- file(paste(rawprocesseddir, "CERES_sample_info.tsv", sep=""), "wb")
write.table(CERES_sample_info,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)