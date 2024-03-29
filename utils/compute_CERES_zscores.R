### fit with one and two-component mixture of t-distributions

rawdir = "data/raw/"
featuredir = "data/features/"
profiledir = "data/profiles/"


rawdir = "data/20Q2/raw/"
featuredir = "data/20Q2/features/"
profiledir = "data/20Q2/profiles/"


### load library
library(readr)

### load data files
gene_effect_cor <- read_csv(paste(rawdir,"gene_effect.csv",sep=""))

colnames(gene_effect_cor)[1] <- "Sample"

scores <- gene_effect_cor
samples <- scores$Sample
scores$Sample <- NULL

### remove columns with NA entries
scores <- scores[ , colSums(is.na(scores)) == 0]


### compute CERES z-scores according to Meyers et al., 2017
means <- apply(scores,2,function(x){mean(x)})
meandiff <- sweep(scores,2,c(means),"-")
stds <- apply(scores,2,function(x){sd(x)})
means2 <- apply(meandiff,1,function(x){mean(x)})
stds2 <- apply(meandiff,1,function(x){sd(x)})
zscores <- sweep(meandiff,1,means2,"-")
zscores <- sweep(zscores,1,stds2,"/")

CERES_zscores <- cbind(Sample = samples, zscores, stringsAsFactors = FALSE)

### sixsigma  summary ####

sixsigmaL <- data.frame(apply(zscores,2,function(x){x < -6}))
sixsigmaR <- data.frame(apply(zscores,2,function(x){x > 6}))
sixsigmaLcount <- apply(sixsigmaL,2,function(x){sum(x)})
sixsigmaRcount <- apply(sixsigmaR,2,function(x){sum(x)})
minimums <- apply(zscores,2,function(x){min(x)})
maximums <- apply(zscores,2,function(x){max(x)})

sixsigma_summary <- data.frame(gene=colnames(zscores),sixsigmaLcount=sixsigmaLcount,sixsigmaRcount=sixsigmaRcount, min = minimums, max = maximums,stringsAsFactors = FALSE)
sixsigma_summary$sixsigmaTotalcount <- with(sixsigma_summary, sixsigmaLcount + sixsigmaRcount)
df_sixsigma_summary <- sixsigma_summary[,c(1,2,3,6,4,5)]

### CERES scores of six-sigma genes (store for fitting with gaussian mixture)
sixsigma_genes <- df_sixsigma_summary[df_sixsigma_summary$sixsigmaTotalcount > 0,]$gene

#gene_effect_cor_sixsigma <- gene_effect_cor[,sixsigma_genes]
gene_effect_cor_sixsigma <- cbind(Sample = samples, gene_effect_cor[,sixsigma_genes], stringsAsFactors = FALSE)



### write output files

output.file <- file(paste(profiledir,"CERES_z-scores.tsv",sep=""),"wb")
write.table(CERES_zscores,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t",
            file = output.file,
            quote = FALSE)
close(output.file)

output.file <- file(paste(profiledir,"CERES_scores_sixsigma.tsv",sep=""),"wb")
write.table(gene_effect_cor_sixsigma,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t",
            file = output.file,
            quote = FALSE)
close(output.file)

output.file <- file(paste(profiledir,"CERES_z-scores_sixsigma_summary.tsv",sep=""),"wb")
write.table(df_sixsigma_summary,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)
