# fit with one and two-component mixture of t-distributions

## Parse input arguments
args = commandArgs(trailingOnly=TRUE)
scores_filename = args[1]
fit_results_filename = args[2]

rand.seed <- 2019

## load library
library(readr)
library(EMMIXskew)
library(dplyr)

scores <- read_delim(scores_filename, 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
scores$Sample <- NULL
genes <- colnames(scores)
num_genes <- length(genes)

## initialize dataframe
tmm_summary <- data.frame(gene=character(num_genes), diffBIC = numeric(num_genes), 
                          mean1 = numeric(num_genes), dof1 = numeric(num_genes), w1 = numeric(num_genes), n1 = numeric(num_genes),
                          mean2 = numeric(num_genes),dof2 = numeric(num_genes), w2 = numeric(num_genes), n2 = numeric(num_genes), 
                          k1ll = numeric(num_genes),k2ll = numeric(num_genes),  
                          k1bic = numeric(num_genes), k2bic = numeric(num_genes), direction = character(num_genes), 
                          stringsAsFactors = FALSE)

## default initialization parameters
initobj1 <- list()
initobj1$pro <- c(0.5, 0.5)
initobj1$mu <- c(-1, 1)
initobj1$sigma <- c(1,1)
initobj1$dof <- c(100,100)
initobj1$delta <- c(0,0)

## fit with mixture models
for (i in 1:num_genes){
  set.seed(rand.seed)
  gene = genes[i]
  mat <- na.omit(scores[[gene]])

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
  tmm_summary[i,] = currrow
}

tmm_summary <- tmm_summary %>% arrange(desc(diffBIC))

# Save results
output.file <- file(fit_results_filename, "wb")
write.table(tmm_summary,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)