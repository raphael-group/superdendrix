##############
#### 21Q1 ####
##############

##############
# feature selection for OncoKB-annotated mutation
# 1) select non-synonymous mutations (exclude silent)
# 2) construct OncoKB features
# 3) filter based on mutation frequency (keep mutations with >= 10)
# 4) convert to sample-mutation format

##############
# load library
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
##############
# set wd and define directories
setwd("C:/Users/pty01/OneDrive - Princeton University/000_ragr/000_github/SuperDendrix/pipeline_21Q1")

rawdir = "../data/21Q1/raw/"
featuredir = "../data/21Q1/features/"
profiledir = "../data/21Q1/profiles/"
rawprocesseddir = "../data/21Q1/raw_preprocessed/"

#############
# load maf
ccle_CERES <- read_delim(paste(featuredir,"CERES_mutations_OncoKB.tsv",sep=""),
                   "\t", escape_double = FALSE, trim_ws = FALSE, col_types = cols(ONCOGENIC = col_character()))

#str(ccle)
colnames(ccle_CERES)
nrow(ccle_CERES) # 589,442 

# load sample info
sample_info <- read_delim(paste(rawprocesseddir,"CERES_sample_info_cancer_type.tsv",sep=""),
           "\t", escape_double = FALSE, trim_ws = TRUE)

unique(ccle_CERES$Variant_Classification)

# 1) exclude non-damaging alterations
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

oncogenic <- c("Likely Oncogenic", "Oncogenic", "Predicted Oncogenic")
lof <- c("Likely Loss-of-function", "Loss-of-function")
gof <- c("Likely Gain-of-function", "Gain-of-function")



### construct OncoKB features
output <- character (nrow(ccle_full))
mutname <- character (nrow(ccle_full))

condition0 <- ccle_full$GENE_IN_ONCOKB == "TRUE"
condition1 <- ccle_full$ONCOGENIC %in% oncogenic
#condition2 <- ccle_full$Variant_Classification %in% inactivating
condition2 <- ccle_full$MUTATION_EFFECT %in% lof
condition3 <- ccle_full$MUTATION_EFFECT %in% gof

# for genes Annotated in OncoKB
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

# for those Not in OncoKB
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
output.file <- file(paste(featuredir,"CERES_mutations_OncoKB_features.tsv",sep=""), "wb")
write.table(ccle_oncokb,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)


output.file <- file(paste(featuredir,"CERES_mutations_full_features.tsv",sep=""), "wb")
write.table(ccle_full_2,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)

output.file <- file(paste(featuredir,"OncoKB_gene_list.tsv",sep=""), "wb")
write.table(genes,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)

output.file <- file(paste(featuredir,"OncoKB_mutation_list.tsv",sep=""), "wb")
write.table(mutations,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote = FALSE)
close(output.file)


### record number of cell lines for each alteration
library(readr)
input_df <- ccle_oncokb
agg <- aggregate(DepMap_ID ~ Alteration, data = input_df,paste, collapse = ",")
df <- agg %>% mutate(DepMap_ID=str_split(DepMap_ID, ","))
df <- df %>% unnest(DepMap_ID) %>% unique()  # remove duplicate entries
agg <- aggregate(DepMap_ID ~ Alteration, data = df,paste, collapse = ",")
agg$numCelllines <- with(agg,str_count(agg$DepMap_ID,",")+1)

output.file <- file(paste(featuredir,"number_of_cell_lines_per_OncoKB_feature.tsv",sep=""), "wb")
write.table(agg,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)
###############

input_df <- ccle_full_2
agg <- aggregate(DepMap_ID ~ Alteration, data = input_df,paste, collapse = ",")
df <- agg %>% mutate(DepMap_ID=str_split(DepMap_ID, ","))
df <- df %>% unnest(DepMap_ID) %>% unique()  # remove duplicate entries
agg <- aggregate(DepMap_ID ~ Alteration, data = df,paste, collapse = ",")
agg$numCelllines <- with(agg,str_count(agg$DepMap_ID,",")+1)

output.file <- file(paste(featuredir,"number_of_cell_lines_per_full_feature.tsv",sep=""), "wb")
write.table(agg,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)


### record number of mutations for each cell lines
####
input_df <- ccle_oncokb
agg <- aggregate(Alteration ~ DepMap_ID, data = input_df,paste, collapse = ",")
df <- agg %>% mutate(Alteration=str_split(Alteration, ","))
df <- df %>% unnest(Alteration) %>% unique()  # remove duplicate entries
agg <- aggregate(Alteration ~ DepMap_ID, data = df,paste, collapse = ",")
agg$numCelllines <- with(agg,str_count(agg$Alteration,",")+1)

output.file <- file(paste(featuredir,"number_of_OncoKB_features_per_cell_line.tsv",sep=""), "wb")
write.table(agg,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

####
input_df <- ccle_full_2
agg <- aggregate(Alteration ~ DepMap_ID, data = input_df,paste, collapse = ",")
df <- agg %>% mutate(Alteration=str_split(Alteration, ","))
df <- df %>% unnest(Alteration) %>% unique()  # remove duplicate entries
agg <- aggregate(Alteration ~ DepMap_ID, data = df,paste, collapse = ",")
agg$numCelllines <- with(agg,str_count(agg$Alteration,",")+1)

output.file <- file(paste(featuredir,"number_of_full_features_per_cell_line.tsv",sep=""), "wb")
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

parse_maf <- function(input_df, outdir, cancertypes, sample_data, metadata=FALSE, full_or_OncoKB="OncoKB"){
  agg <- aggregate(DepMap_ID ~ Alteration, data = input_df,paste, collapse = ",")
  
  df <- agg %>% mutate(DepMap_ID=str_split(DepMap_ID, ","))
  df <- df %>% unnest(DepMap_ID) %>% unique()  # remove duplicate entries
  
  ### optional output number of cell lines for each alteration
  agg <- aggregate(DepMap_ID ~ Alteration, data = df,paste, collapse = ",")
  agg$numCelllines <- with(agg,str_count(agg$DepMap_ID,",")+1)
  #####
  
  df <- agg %>% mutate(DepMap_ID=str_split(DepMap_ID, ","))
  df <- df %>% unnest(DepMap_ID) %>% unique()
  newdf<- df
  newdf$numCelllines <- NULL
  newdf <- aggregate(Alteration ~ DepMap_ID, data = newdf, paste, collapse = "\t")
  
  newdf$numGenes <- with(newdf,str_count(newdf$Alteration,"\t")+1)
  
  ### optional output number of alterations for each cell lines
  if (metadata){
    if (cancertypes){
      fn = file(paste(outdir,"number_of_",full_or_OncoKB,"_features_per_cell_line", "_cancertypes",".tsv",sep=""),"wb")
    }
    nummutpersample <- data.frame(sample=newdf$DepMap_ID,numAlterations=newdf$numGenes)
    
    fn = file(paste(outdir,"number_of_",full_or_OncoKB,"_features_per_cell_line.tsv",sep=""),"wb")
    write.table(nummutpersample,file=fn,row.names=FALSE,col.names=TRUE,sep="\t")
    close(fn)
  }
  newdf$numGenes <- NULL
  colnames(newdf) <- c("#Sample","Events")
  
  
  ## optionally add cancer-type features
  if (cancertypes){
    cancertypes <- sample_data[,c("DepMap_ID","Cancer_type")]
    colnames(cancertypes) <- c("#Sample","primary_disease")
    newdf <- merge(cancertypes,newdf,by="#Sample")
    newdf$mut <- with(newdf, paste(newdf$primary_disease,newdf$Events,sep = "\t"))
    newdf <- newdf[,c("#Sample","mut")]
    colnames(newdf) <- c("#Sample","Events")
  }
  return(newdf)
}

features_full <- parse_maf(ccle_full_2, featuredir, FALSE, sample_info, FALSE)
features_oncokb <- parse_maf(ccle_oncokb, featuredir, FALSE, sample_info, FALSE)

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

output.file <- file(paste(featuredir,"CERES_OncoKB_MUT.tsv",sep=""), "wb")
write.table(features_oncokb,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(featuredir,"CERES_OncoKB_CT.tsv",sep=""), "wb")
write.table(features_oncokb_ct,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)
