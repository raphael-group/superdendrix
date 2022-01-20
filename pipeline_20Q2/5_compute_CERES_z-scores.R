##############
#### 21Q1 ####
##############

##############
### Compute CERES zscores before fitting with TMM or GMM


##############
# load library
library(readr)
##############
# set wd and define directories
setwd("C:/Users/pty01/OneDrive - Princeton University/000_ragr/000_github/SuperDendrix/pipeline_21Q1")

rawdir = "../data/21Q1/raw/"
featuredir = "../data/21Q1/features/"
profiledir = "../data/21Q1/profiles/"
rawprocesseddir = "../data/21Q1/raw_preprocessed/"

### load data files
Achilles_gene_effect <- read_delim(paste(rawprocesseddir,"CERES_scores.tsv",sep=""),
                         "\t", escape_double = FALSE, trim_ws = FALSE)
scores <- Achilles_gene_effect
samples <- scores$DepMap_ID
scores$DepMap_ID <- NULL

### compute CERES z-scores according to Meyers et al., 2017
means <- apply(scores,2,function(x){mean(x, na.rm = TRUE)})
meandiff <- sweep(scores,2,c(means),"-")
stds <- apply(scores,2,function(x){sd(x, na.rm = TRUE)})
means2 <- apply(meandiff,1,function(x){mean(x, na.rm = TRUE)})
stds2 <- apply(meandiff,1,function(x){sd(x, na.rm = TRUE)})
zscores <- sweep(meandiff,1,means2,"-")
zscores <- sweep(zscores,1,stds2,"/")

CERES_zscores <- cbind(DepMap_ID = samples, zscores, stringsAsFactors = FALSE)

### record sixsigma  summary ####

sixsigmaL <- data.frame(apply(zscores,2,function(x){x < -6}))
sixsigmaR <- data.frame(apply(zscores,2,function(x){x > 6}))
sixsigmaLcount <- apply(sixsigmaL,2,function(x){sum(x, na.rm = TRUE)})
sixsigmaRcount <- apply(sixsigmaR,2,function(x){sum(x, na.rm = TRUE)})
minimums <- apply(zscores,2,function(x){min(x, na.rm = TRUE)})
maximums <- apply(zscores,2,function(x){max(x, na.rm = TRUE)})

sixsigma_summary <- data.frame(gene=colnames(zscores),sixsigmaLcount=sixsigmaLcount,sixsigmaRcount=sixsigmaRcount, min = minimums, max = maximums,stringsAsFactors = FALSE)
sixsigma_summary$sixsigmaTotalcount <- with(sixsigma_summary, sixsigmaLcount + sixsigmaRcount)
df_sixsigma_summary <- sixsigma_summary[,c(1,2,3,6,4,5)]


### CERES scores of six-sigma genes (store for fitting with gaussian mixture)
sixsigma_genes <- df_sixsigma_summary[df_sixsigma_summary$sixsigmaTotalcount > 0,]$gene

### write output files
output.file <- file(paste(profiledir,"CERES_zscores.tsv",sep=""),"wb")
write.table(CERES_zscores,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t",
            file = output.file,
            quote = FALSE)
close(output.file)

output.file <- file(paste(profiledir,"CERES_zscores_sixsigma_genes.tsv",sep=""),"wb")
write.table(sixsigma_genes,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t",
            file = output.file,
            quote = FALSE)
close(output.file)

output.file <- file(paste(profiledir,"CERES_zscores_sixsigma_summary.tsv",sep=""),"wb")
write.table(df_sixsigma_summary,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)