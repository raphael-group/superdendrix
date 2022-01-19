# load library
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

version = "20Q2/"
rawdir = paste("../data/",version,"raw/",sep="")
featuredir = paste("../data/",version,"features/",sep="")
profiledir = paste("../data/",version,"profiles/",sep="")

# load maf
maf <- read_delim(paste(featuredir,"mutations_oncokb.csv",sep=""), 
                     ",", escape_double = FALSE, trim_ws = TRUE)



# load oncoKB
oncokb <- read_delim(paste(featuredir,"oncoKB_list.tsv",sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
oncogenes <- oncokb[oncokb$type == "oncogene",]$Gene
tsg <- oncokb[oncokb$type == "tsg",]$Gene

# load CCLE maf
ccle <- read_delim(paste(rawdir,"depmap_19Q1_mutation_calls_v2.csv",sep=""), 
                   ",", escape_double = FALSE, trim_ws = TRUE)
ccle$X1 <- NULL

# load CERES cell lines
sample_info <- read_csv(paste(rawdir,"cell_line_info.csv",sep=""))

# restrict to cell lines with CERES scores
ccle <- ccle %>% filter(DepMap_ID %in% sample_info$DepMap_ID)


#### exclude non-damaging alterations

tokeep <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Del", "Splice_Site", "De_novo_Start_OutOfFrame",
            "Frame_Shift_Ins", "Start_Codon_Del", "Stop_Codon_Del", "Start_Codon_Ins", "Stop_Codon_Ins")

ccle <- ccle[ccle$Variant_Classification %in% tokeep,]

#### CCLE inactivating alterations

inactivating <- c("Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Del",
           "Frame_Shift_Ins", "Start_Codon_Del", "Stop_Codon_Del", "Start_Codon_Ins", "Stop_Codon_Ins", "De_novo_Start_OutOfFrame", "Splice_Site")


#### annotate each mutation in the maf according to OncoKB
output <- character (nrow(ccle))
mutname <- character (nrow(ccle))
condition0 <- ccle$Hugo_Symbol %in% oncokb$Gene
condition1 <- !is.na(ccle$Protein_Change)
condition2 <- ccle$Hugo_Symbol %in% tsg
condition3 <- ccle$Hugo_Symbol %in% oncogenes


for (i in (1:nrow(ccle))[condition0]){
  gene = ccle$Hugo_Symbol[i]
  if (condition1[i]){
    aachn = substr(ccle$Protein_Change[i],3,nchar(ccle$Protein_Change[i]))
  }
  mut_class = ccle$Variant_Classification[i]
  if (condition2[i]){   # gene is tsg
    if (mut_class %in% inactivating){ # matching variant classification
      output[i] <- "inactivating"
      mutname[i] <- paste(gene,"_I_MUT",sep="")
    }else if (mut_class == "Missense_Mutation" & grepl(aachn,oncokb[oncokb$Gene == gene,2])){  # or matching aa change
      output[i] <- "inactivating"
      mutname[i] <- paste(gene,"_I_MUT",sep="")
    }else{
      output[i] <- "other"
      mutname[i] <- paste(gene,"_O_MUT",sep="")
    }
  }else if (condition3[i]){   # gene is oncogene
    if (mut_class == "Missense_Mutation" & grepl(aachn, oncokb[oncokb$Gene == gene,2])){  # matching variant classification and aa change
      output[i] <- "activating"
      mutname[i] <- paste(gene,"_A_MUT",sep="")
    }else{
      output[i] <- "other"
      mutname[i] <- paste(gene,"_O_MUT",sep="")
    }
  }else{
    output[i] <- "other"
    mutname[i] <- paste(gene,"_O_MUT",sep="")
  }
}

for (i in (1:nrow(ccle))[!condition0]){
  gene = ccle$Hugo_Symbol[i]
  output[i] <- "noONCOKB"
  mutname[i] <- paste(gene,"_NA_MUT",sep="")
}

ccle$"Variant annotation" <- output
ccle$"Alteration" <- mutname

ccle_oncokb <- ccle[!grepl("noONCOKB", ccle$"Variant annotation"),]


#output.file <- file("../data/features/depmap_19Q1_mutation_calls_oncoKB.tsv", "wb")
#write.table(ccle,
            #row.names = FALSE,
            #col.names = TRUE,
            #file = output.file,
            #sep = "\t",
            #quote = FALSE)
#close(output.file)


##############################################
#### change to sample - alteration format ####
##############################################


parse_maf <- function(input_df, outdir, cancertypes,sample_data){
  agg <- aggregate(Tumor_Sample_Barcode ~ Alteration, data = input_df,paste, collapse = ",")
  
  df <- agg %>% mutate(Tumor_Sample_Barcode=str_split(Tumor_Sample_Barcode, ","))
  df <- df %>% unnest(Tumor_Sample_Barcode) %>% unique()  # remove duplicate entries
  
  agg <- aggregate(Tumor_Sample_Barcode ~ Alteration, data = df,paste, collapse = ",")
  agg$numCelllines <- with(agg,str_count(agg$Tumor_Sample_Barcode,",")+1)
  

  # filter for minimum number of cell lines
  agg2 <- agg[!((agg$numCelllines < 10) & grepl("_O",agg$Alteration)), ]
  agg2 <- agg2[!(agg2$numCelllines < 10 & grepl("_NA",agg2$Alteration)), ]

  ### optional output number of cell lines for each alteration
  df <- agg2 %>% mutate(Tumor_Sample_Barcode=str_split(Tumor_Sample_Barcode, ","))
  df <- df %>% unnest(Tumor_Sample_Barcode) %>% unique()
  newdf<- df
  newdf$numCelllines <- NULL
  newdf <- aggregate(Alteration ~ Tumor_Sample_Barcode, data = newdf, paste, collapse = "\t")
  
  newdf$numGenes <- with(newdf,str_count(newdf$Alteration,"\t")+1)
  
  ### optional output number of alterations for each cell lines
  newdf$numGenes <- NULL
  colnames(newdf) <- c("#Sample","Events")
  
  
  ## optionally add cancer-type features
  if (cancertypes){
    cancertypes <- sample_data[,c("DepMap_ID","primary_tissue")]
    colnames(cancertypes) <- c("#Sample","primary_tissue")
    newdf <- merge(cancertypes,newdf,by="#Sample")
    newdf$mut <- with(newdf, paste(newdf$primary_tissue,newdf$Events,sep = "\t"))
    newdf <- newdf[,c("#Sample","mut")]
    colnames(newdf) <- c("#Sample","Events")
  }
  return(newdf)
}

features_full <- parse_maf(ccle,featuredir, FALSE,sample_info)
features_full_ct <- parse_maf(ccle, featuredir, TRUE,sample_info)
features_oncokb <- parse_maf(ccle_oncokb, featuredir, FALSE,sample_info)
features_oncokb_ct <- parse_maf(ccle_oncokb, featuredir, TRUE,sample_info)

feature_list <- features_oncokb_ct %>% mutate(Events=str_split(Events, "\t"))
feature_list <- feature_list %>% unnest(Events)
features <- data.frame(Alteration=unique(feature_list$Events)
colnames(features) <- "#alteration"

# write output files

output.file <- file(paste(featuredir,"Q1_CERES_full_MUT.tsv",sep=""), "wb")
write.table(features_full,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(featuredir,"Q1_CERES_full_CT.tsv",sep=""), "wb")
write.table(features_full_ct,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(featuredir,"Q1_CERES_oncoKB_MUT.tsv",sep=""), "wb")
write.table(features_oncokb,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(featuredir,"Q1_CERES_oncoKB_CT.tsv",sep=""), "wb")
write.table(features_oncokb_ct,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)


output.file <- file(paste(featuredir,"Q1_CERES_oncoKB_feature_list.csv",sep=""), "wb")
write.table(features,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = ",",
            quote=FALSE)
close(output.file)



