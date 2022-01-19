# fit with one and two-component mixture of t-distributions

rawdir <- "../data/20Q2/raw/"
featuredir <- "../data/20Q2/features/"
profiledir <- "../data/20Q2/profiles/"

# integrated
#rawdir = "../data/integrated/raw/"
#featuredir = "../data/integrated/features/"
#profiledir = "../data/integrated/profiles/"

#install.packages("EMMIXskew_1.0.3 (2).tar.gz", repos = NULL, type ='source')

rand.seed <- 2019
### load library
library(readr)
library(EMMIXskew)
library(dplyr)
#library(rgr)
### load data files
#CERES_scores_sixsigma <- read_delim(paste(profiledir,"CERES_scores_sixsigma.tsv",sep=""), 
                      #"\t", escape_double = FALSE, trim_ws = TRUE)



CERES_z <- read_delim(paste(profiledir,"CERES_zscores.tsv",sep=""), 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
CERES_z$Sample <- NULL

CERES_zscores_sixsigma_summary <- read_delim(paste(profiledir,"CERES_zscores_sixsigma_summary.tsv",sep=""), 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)

### select six-sigma genes
#CERES_zscores_sixsigma_summary$sixsigmaTotalcount

sixsigma_only <- FALSE

if (sixsigma_only){
  genes <- CERES_zscores_sixsigma_summary[CERES_zscores_sixsigma_summary$sixsigmaTotalcount > 0,]$gene
}else{
  genes <- colnames(CERES_z)
}


tmm_summary <- data.frame(gene=character(0), diffBIC = numeric(0), mean1 = numeric(0), dof1 = numeric(0), w1 = numeric(0), n1 = numeric(0),
                          mean2 = numeric(0),dof2 = numeric(0), w2 = numeric(0),
                          n2 = numeric(0), k1ll = numeric(0),k2ll = numeric(0),  k1bic = numeric(0), k2bic = numeric(0), direction = character(0), stringsAsFactors = FALSE)

### default initialization parameters
initobj1 <- list()
initobj1$pro <- c(0.5, 0.5)
initobj1$mu <- c(-1, 1)
initobj1$sigma <- c(1,1)
initobj1$dof <- c(100,100)
initobj1$delta <- c(0,0)

#na.omit(datacollected) 

### fit with mixture models
for (i in 1:length(genes)){
#for (i in 1:2){
  print(i)
  set.seed(rand.seed)
  gene = genes[i]
  #rm_mat <- remove.na(CERES_z[[gene]]) # rgr library
  #mat <- rm_mat$x
  mat <- na.omit(CERES_z[[gene]])

  k2 <- EmSkew(mat, g=2, init=initobj1, distr="mvt")
  k1 <- EmSkew(mat, g=1, distr="mvt")
  diffBIC <- k1$bic - k2$bic
  if (k2$mu[2] < k2$mu[1]){
    mu1 <- k2$mu[2]
    mu2 <- k2$mu[1]
    
    dof1 <- k2$dof[2]
    dof2 <- k2$dof[1]
    
    w1 <- k2$pro[2]
    w2 <- k2$pro[1]
    
    df <- data.frame(k2$tau)
    df$LR <- with(df,log(X2/X1))
    df$N1 <- with(df,(LR > 0)+0)
    df$N2 <- with(df,(LR < 0)+0)
    n1 <- sum(df$N1)
    n2 <- sum(df$N2)
  }
  if (k2$mu[1] < k2$mu[2]){
    mu2 <- k2$mu[2]
    mu1 <- k2$mu[1]
    
    dof2 <- k2$dof[2]
    dof1 <- k2$dof[1]
    
    w2 <- k2$pro[2]
    w1 <- k2$pro[1]
    
    df <- data.frame(k2$tau)
    df$LR <- with(df,log(X1/X2))
    df$N1 <- with(df,(LR > 0)+0)
    df$N2 <- with(df,(LR < 0)+0)
    n1 <- sum(df$N1)
    n2 <- sum(df$N2)
  }
  
  if (n1 < n2){
    direction <- "increased dependency"
  }else {
    direction <- "decreased dependency"
  }
  
  currrow <- list(gene,diffBIC, mu1, dof1, w1, n1, mu2, dof2, w2, n2, k1$loglik, k2$loglik, k1$bic, k2$bic, direction)
  tmm_summary[nrow(tmm_summary) + 1,] = currrow
}

tmm_summary <- tmm_summary %>% arrange(desc(diffBIC))


output.file <- file(paste(profiledir,"tmm_summary.tsv",sep=""), "wb")
write.table(tmm_summary,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)

#differential_dependencies <- tmm_summary[tmm_summary$diffBIC>0,]$gene
#CERES_scores_2C <- CERES_scores_sixsigma[,differential_dependencies]

#output.file <- file(paste(profiledir,"CERES_scores_2C.tsv",sep=""), "wb")
#write.table(CERES_scores_2C,
            #row.names = FALSE,
            #col.names = TRUE,
            #file = output.file,
            #sep = "\t",
            #quote = FALSE)
#close(output.file)


