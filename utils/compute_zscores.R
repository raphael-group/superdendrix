# Compute CERES zscores before fitting with TMM or GMM

## Parse input arguments
args = commandArgs(trailingOnly=TRUE)
dependency_scores_filename = args[1]
zscores_filename = args[2]

## load library
library(readr)

## load data files
dependency_scores <- read_csv(dependency_scores_filename)
colnames(dependency_scores)[1] <- "Sample"
scores <- dependency_scores
samples <- dependency_scores$Sample
scores$Sample <- NULL

## compute the z-scores according to Meyers et al., 2017
means <- apply(scores,2,function(x){mean(x, na.rm = TRUE)})
meandiff <- sweep(scores,2,c(means),"-")
stds <- apply(scores,2,function(x){sd(x, na.rm = TRUE)})
means2 <- apply(meandiff,1,function(x){mean(x, na.rm = TRUE)})
stds2 <- apply(meandiff,1,function(x){sd(x, na.rm = TRUE)})
zscores <- sweep(meandiff,1,means2,"-")
zscores <- sweep(zscores,1,stds2,"/")

zscores_df <- cbind(Sample = samples, zscores, stringsAsFactors = FALSE)

## write output files
output.file <- file(zscores_filename,"wb")
write.table(zscores_df,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t",
            file = output.file,
            quote = FALSE)
close(output.file)
