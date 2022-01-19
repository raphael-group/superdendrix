# load library
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# directory names
rawdir = "../data/integrated/raw/"
featuredir = "../data/integrated/features/"
profiledir = "../data/integrated/profiles/"

# load maf
ccle <- read_delim("../data/integrated/features/CCLE_mutations_integrated_oncoKB.tsv", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)
#colnames(ccle)
#nrow(ccle) # 1293716

# load CERES
CERES <- read_delim("../data/integrated/raw/CERES_FC_processed.tsv",
                  "\t", escape_double = FALSE, trim_ws = TRUE)


# load sample info
sample_info <- read_csv(paste(rawdir,"sample_info.csv",sep=""))
sample_info <- sample_info %>% filter(BROAD_ID %in% CERES$Sample)


# restrict to cell lines with CERES scores
ccle_CERES <- ccle %>% filter(SAMPLE_ID %in% CERES$DepMap_ID)

length(unique(ccle_CERES$SAMPLE_ID))

length(unique(ccle2$SAMPLE_ID))
unique(ccle2$Variant_Classification)
unique(ccle2$GENE_IN_ONCOKB)
unique(ccle2$VARIANT_IN_ONCOKB)
unique(ccle2$MUTATION_EFFECT)
unique(ccle2$ONCOGENIC)

# preprocessing to 
# exclude non-damaging alterations
tokeep <- c("Missense_Mutation", 
            "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "De_novo_Start_OutOfFrame",
            "Frame_Shift_Ins", "Frame_Shift_Del",
            "Start_Codon_SNP", "Start_Codon_Del",  "Start_Codon_Ins",
            "Stop_Codon_Del", "Stop_Codon_Ins",
            "In_Frame_Del", "In_Frame_Ins")

inactivating <- c("Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "De_novo_Start_OutOfFrame",
            "Frame_Shift_Ins", "Frame_Shift_Del",
            "Start_Codon_SNP", "Start_Codon_Del",  "Start_Codon_Ins",
            "Stop_Codon_Del", "Stop_Codon_Ins",
            "In_Frame_Del", "In_Frame_Ins")

ccle_full <- ccle_CERES %>% filter(Variant_Classification %in% tokeep)


##### construct oncokb features
ccle_full <- ccle_CERES %>% filter(Variant_Classification %in% tokeep)

oncogenic <- c("Likely Oncogenic", "Oncogenic", "Predicted Oncogenic")
lof <- c("Likely Loss-of-function", "Loss-of-function")
gof <- c("Likely Gain-of-function", "Gain-of-function")


output <- character (nrow(ccle_full))
mutname <- character (nrow(ccle_full))

condition0 <- ccle_full$GENE_IN_ONCOKB == "TRUE"
condition1 <- ccle_full$ONCOGENIC %in% oncogenic
#condition2 <- ccle_full$Variant_Classification %in% inactivating
condition2 <- ccle_full$MUTATION_EFFECT %in% lof
condition3 <- ccle_full$MUTATION_EFFECT %in% gof

# Annotated in OncoKB
for (i in (1:nrow(ccle_full))[condition0]){
  gene = ccle_full$Hugo_Symbol[i]
  if (condition1[i]){
    # oncogenic
    if (condition2[i]){
      # LOF
      output[i] <- "inactivating"
      mutname[i] <- paste(gene,"_I_MUT",sep="")
    }
    else if (condition3[i]){
      # GOF
      output[i] <- "activating"
      mutname[i] <- paste(gene,"_A_MUT",sep="")
    }
    else {
      output[i] <- "other"
      mutname[i] <- paste(gene,"_O_MUT",sep="")
    }
  }
  else{
    # not oncogenic
    output[i] <- "other"
    mutname[i] <- paste(gene,"_O_MUT",sep="")
    # print("here444")
  }
}

# Not in OncoKB
for (i in (1:nrow(ccle_full))[!condition0]){
  gene = ccle_full$Hugo_Symbol[i]
  output[i] <- "noONCOKB"
  mutname[i] <- paste(gene,"_NA_MUT",sep="")
}

ccle_full$"Variant annotation" <- output
ccle_full$"Alteration" <- mutname


###### filter by frequency (exclude O and NA with < 10 cell lines)
O_NA <- ccle_full %>% filter(`Variant annotation` %in% c("other", "noONCOKB"))
length(unique(O_NA$Alteration))
O_NA2 <- O_NA %>% group_by(Alteration) %>% summarise(numCelllines = n())
to_exclude <- O_NA2 %>% filter(numCelllines < 10)

# nrow ccle_full 433354 
ccle_full_2 <- ccle_full %>% filter(!(Alteration %in% to_exclude$Alteration))
# nrow ccle_full_2 403713

ccle_oncokb <- ccle_full_2[!grepl("noONCOKB", ccle_full_2$"Variant annotation"),]

length(unique(ccle_oncokb$Hugo_Symbol))
length(unique(ccle_oncokb$Alteration))

genes <- data.frame(Gene = unique(ccle_oncokb$Hugo_Symbol))
mutations <- data.frame(Mutation = unique(ccle_oncokb$Alteration))


# save files
output.file <- file("../data/20Q4/features/CCLE_mutations_20Q4_oncoKB_features.tsv", "wb")
write.table(ccle_oncokb,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)

output.file <- file("../data/20Q4/features/CCLE_mutations_20Q4_full_features.tsv", "wb")
write.table(ccle_full_2,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)

output.file <- file("../data/20Q4/features/OncoKB_gene_list.tsv", "wb")
write.table(genes,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)

output.file <- file("../data/20Q4/features/OncoKB_mutation_list.tsv", "wb")
write.table(mutations,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)


### optional output number of cell lines for each alteration
####
library(readr)
CCLE_mutations_20Q4_oncoKB_features <- read_delim("C:/Users/pty01/OneDrive - Princeton University/000_ragr/000_github/SuperDendrix/data/20Q4/features/CCLE_mutations_20Q4_oncoKB_features.tsv", 
                                                  "\t", escape_double = FALSE, trim_ws = TRUE)

input_df <- CCLE_mutations_20Q4_oncoKB_features
agg <- aggregate(SAMPLE_ID ~ Alteration, data = input_df,paste, collapse = ",")

df <- agg %>% mutate(SAMPLE_ID=str_split(SAMPLE_ID, ","))
df <- df %>% unnest(SAMPLE_ID) %>% unique()  # remove duplicate entries

agg <- aggregate(SAMPLE_ID ~ Alteration, data = df,paste, collapse = ",")
agg$numCelllines <- with(agg,str_count(agg$SAMPLE_ID,",")+1)

output.file <- file(paste(featuredir,"number_of_celllines_per_oncokb_feature.tsv",sep=""), "wb")
write.table(agg,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

##############################################
#### change to sample - alteration format ####
##############################################

parse_maf <- function(input_df, outdir, cancertypes,sample_data,metadata=FALSE){
  agg <- aggregate(SAMPLE_ID ~ Alteration, data = input_df,paste, collapse = ",")
  
  df <- agg %>% mutate(SAMPLE_ID=str_split(SAMPLE_ID, ","))
  df <- df %>% unnest(SAMPLE_ID) %>% unique()  # remove duplicate entries
  
  ### optional output number of cell lines for each alteration
  agg <- aggregate(SAMPLE_ID ~ Alteration, data = df,paste, collapse = ",")
  agg$numCelllines <- with(agg,str_count(agg$SAMPLE_ID,",")+1)
  #####
  
  df <- agg %>% mutate(SAMPLE_ID=str_split(SAMPLE_ID, ","))
  df <- df %>% unnest(SAMPLE_ID) %>% unique()
  newdf<- df
  newdf$numCelllines <- NULL
  newdf <- aggregate(Alteration ~ SAMPLE_ID, data = newdf, paste, collapse = "\t")
  
  newdf$numGenes <- with(newdf,str_count(newdf$Alteration,"\t")+1)
  
  ### optional output number of alterations for each cell lines
  if (metadata){
    if (cancertypes){
      ct <- "_cancertypes"
      fn = file(paste(outdir,"features_per_sample",ct,".tsv",sep=""),"wb")
    }
    nummutpersample <- data.frame(sample=newdf$SAMPLE_ID,numAlterations=newdf$numGenes)
    
    fn = file(paste(outdir,"features_per_sample_oncokb",".tsv",sep=""),"wb")
    write.table(nummutpersample,file=fn,row.names=FALSE,col.names=TRUE,sep="\t")
    close(fn)
  }
  newdf$numGenes <- NULL
  colnames(newdf) <- c("#Sample","Events")
  
  
  ## optionally add cancer-type features
  if (cancertypes){
    cancertypes <- sample_data[,c("DepMap_ID","primary_disease")]
    colnames(cancertypes) <- c("#Sample","primary_disease")
    newdf <- merge(cancertypes,newdf,by="#Sample")
    newdf$mut <- with(newdf, paste(newdf$primary_disease,newdf$Events,sep = "\t"))
    newdf <- newdf[,c("#Sample","mut")]
    colnames(newdf) <- c("#Sample","Events")
  }
  return(newdf)
}

features_full <- parse_maf(ccle_full_2, featuredir, FALSE, sample_info, TRUE)
features_oncokb <- parse_maf(ccle_oncokb, featuredir, FALSE, sample_info, TRUE)

features_full_ct <- parse_maf(ccle_full_2, featuredir, TRUE, sample_info, FALSE)
features_oncokb_ct <- parse_maf(ccle_oncokb, featuredir, TRUE, sample_info, FALSE)

# write output files

output.file <- file(paste(featuredir,"CERES_full_MUT.tsv",sep=""), "wb")
write.table(features_full,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(featuredir,"CERES_full_CT.tsv",sep=""), "wb")
write.table(features_full_ct,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(featuredir,"CERES_oncoKB_MUT.tsv",sep=""), "wb")
write.table(features_oncokb,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(featuredir,"CERES_oncoKB_CT.tsv",sep=""), "wb")
write.table(features_oncokb_ct,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)
