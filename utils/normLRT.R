library(MASS)
library(sn)
library(readr)
library(dplyr)

#rawscore <- read_delim("Q1_gene_effect.tsv", 
                       #"\t", escape_double = FALSE, trim_ws = TRUE)
#X2 <- matrix(1, nrow = 558, ncol = 1)
#normal <- fitdistr(rawscore$PIGR, 'normal')
#skewt <- st.mple(x=X2, y = rawscore$PIGR)

#normal$loglik
#skewt$logL

##sample_info <- read_csv(paste(rawdir,"sample_info.csv",sep=""))
CERES_FC_processed <- read_csv("../data/20Q2/raw/Achilles_gene_effect.csv", na="")
CERES_FC_processed$DepMap_ID <- NULL

#CERES_FC_processed <- read_delim("../data/integrated/raw/CERES_FC_processed.tsv", 
#                                 "\t", escape_double = FALSE, trim_ws = TRUE)

#CERES_FC_processed$Sample <- NULL
genes <- colnames(CERES_FC_processed)
#head(CERES_FC_processed)



# compute normLRT score: 2*[ln(likelihood for Skewed-t) - ln(likelihood for Gaussian)]

results <- data.frame(gene = genes, LRT = 0)
#X2 <- matrix(1, nrow = dim(CERES_FC_processed)[1], ncol = 1)

for (i in 1:length(genes)){
  gene <- genes[i]
  curr <- as.vector(na.omit(CERES_FC_processed[,i][[1]]))
  X2 <- matrix(1,nrow = length(curr),ncol=1)
#  normal <- fitdistr(CERES_FC_processed[,i][[1]], 'normal')
#  skewt <- st.mple(x=X2, y = CERES_FC_processed[,i][[1]])

  normal <- fitdistr(curr, 'normal')
  skewt <- st.mple(x=X2, y = curr)
  
  
  LRT <- 2*(skewt$logL - normal$loglik)
  
  results$LRT[i] <- LRT
  print(i)
  print(length(curr))
}

sorted_results <- results %>% arrange(desc(LRT))

output.file <- file("../data/20Q2/normLRT.tsv", "wb")
write.table(sorted_results,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)
