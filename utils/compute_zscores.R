### fit with one and two-component mixture of t-distributions

rawdir = "../data/raw/"
featuredir = "../data/features/"
profiledir = "../data/profiles/"

### load library
library(readr)

### load data files
gene_effect_cor <- read_csv(paste(rawdir,"gene_effect_corrected.csv",sep=""))

colnames(gene_effect_cor)[1] <- "Broad_ID"

scores <- gene_effect_cor
samples <- scores$Broad_ID
scores$Broad_ID <- NULL



### compute CERES z-scores according to Meyers et al., 2017
means <- apply(scores,2,function(x){mean(x)})
meandiff <- sweep(scores,2,c(means),"-")
stds <- apply(scores,2,function(x){sd(x)})
means2 <- apply(meandiff,1,function(x){mean(x)})
stds2 <- apply(meandiff,1,function(x){sd(x)})
zscores <- sweep(meandiff,1,means2,"-")
zscores <- sweep(zscores,1,stds2,"/")

CERES_zscores <- cbind(Sample = samples, zscores, stringsAsFactors = FALSE)

output.file <- file(paste(profiledir,"CERES_differential_dependencies.tsv",sep=""),"wb")
write.table(CERES_zscores,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t",
            file = output.file,
            quote = FALSE)
close(output.file)


#### six_sigma differential dependency summary ####

sixsigmaL <- data.frame(apply(zscores,2,function(x){x < -6}))
sixsigmaR <- data.frame(apply(zscores,2,function(x){x > 6}))
sixsigmaLcount <- apply(sixsigmaL,2,function(x){sum(x)})
sixsigmaRcount <- apply(sixsigmaR,2,function(x){sum(x)})
minimums <- apply(zscores,2,function(x){min(x)})
maximums <- apply(zscores,2,function(x){max(x)})

sixsigma_summary <- data.frame(gene=colnames(zscores),sixsigmaLcount=sixsigmaLcount,sixsigmaRcount=sixsigmaRcount, min = minimums, max = maximums)
sixsigma_summary$sixsigmaTotalcount <- with(sixsigma_summary, sixsigmaLcount + sixsigmaRcount)
df_sixsigma_summary <- sixsigma_summary[,c(1,2,3,6,4,5)]


output.file <- file(paste(profiledir,"CERES_zscores_sixsigma_summary.tsv",sep=""),"wb")
write.table(df_sixsigma_summary,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)