# load library
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Parse input arguments
args = commandArgs(trailingOnly=TRUE)

mutations_filename = args[1]
sample_info_filename = args[2]
output_dir = args[3]

# load mutations
mutations <- read_delim(mutations_filename,
                   "\t", escape_double = FALSE, trim_ws = FALSE, col_types = cols(ONCOGENIC = col_character(), 
                                                                                  Chromosome = col_character()))

# load sample info
sample_info <- read_delim(sample_info_filename,
           "\t", escape_double = FALSE, trim_ws = TRUE)

# restrict to cell lines with CERES scores
mutations_avana <- mutations %>% filter(Sample_ID %in% sample_info$DepMap_ID)
# 20Q2 :767 cell lines only 

# exclude non-damaging alterations
mutaion_classes <- c("Missense_Mutation", 
            "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "De_novo_Start_OutOfFrame",
            "Frame_Shift_Ins", "Frame_Shift_Del",
            "Start_Codon_SNP", "Start_Codon_Del",  "Start_Codon_Ins",
            "Stop_Codon_Del", "Stop_Codon_Ins",
            "In_Frame_Del", "In_Frame_Ins")

mutations_all <- mutations_avana %>% filter(Variant_Classification %in% mutaion_classes)

##### construct oncokb features
oncogenic <- c("Likely Oncogenic", "Oncogenic", "Predicted Oncogenic")
lof <- c("Likely Loss-of-function", "Loss-of-function")
gof <- c("Likely Gain-of-function", "Gain-of-function")

output <- character(nrow(mutations_all))
mutname <- character(nrow(mutations_all))

condition0 <- mutations_all$GENE_IN_ONCOKB == "TRUE"
condition1 <- mutations_all$ONCOGENIC %in% oncogenic
condition2 <- mutations_all$MUTATION_EFFECT %in% lof
condition3 <- mutations_all$MUTATION_EFFECT %in% gof

# Annotated in OncoKB
for (i in (1:nrow(mutations_all))[condition0]){
  gene = mutations_all$Hugo_Symbol[i]
  if (condition1[i]){ # oncogenic
    if (condition2[i]){ # LOF
      output[i] <- "inactivating"
      mutname[i] <- paste(gene,"_I_MUT",sep="")
    }
    else if (condition3[i]){ # GOF
      output[i] <- "activating"
      mutname[i] <- paste(gene,"_A_MUT",sep="")
    }
    else { # Other
      output[i] <- "other"
      mutname[i] <- paste(gene,"_O_MUT",sep="")
    }
  }
  else{ # not oncogenic
    output[i] <- "other"
    mutname[i] <- paste(gene,"_O_MUT",sep="")
  }
}

# Not in OncoKB
for (i in (1:nrow(mutations_all))[!condition0]){
  gene = mutations_all$Hugo_Symbol[i]
  output[i] <- "noONCOKB"
  mutname[i] <- paste(gene,"_NA_MUT",sep="")
}

mutations_all$"Variant annotation" <- output
mutations_all$"Alteration" <- mutname

###### filter by frequency (exclude O and NA with < 10 cell lines)
O_NA <- mutations_all %>% filter(`Variant annotation` %in% c("other", "noONCOKB"))
O_NA2 <- O_NA %>% group_by(Alteration) %>% summarise(numCelllines = n())
to_exclude <- O_NA2 %>% filter(numCelllines < 10)
mutations_all_freq <- mutations_all %>% filter(!(Alteration %in% to_exclude$Alteration))

mutations_oncokb <- mutations_all_freq[!grepl("noONCOKB", mutations_all_freq$"Variant annotation"),]


##############################################
#### change to sample - alteration format ####
##############################################

parse_mutations <- function(input_df, cancertypes,sample_data){
  agg <- aggregate(Sample_ID ~ Alteration, data = input_df,paste, collapse = ",")
  
  df <- agg %>% mutate(Sample_ID=str_split(Sample_ID, ","))
  df <- df %>% unnest(Sample_ID) %>% unique()  # remove duplicate entries
  
  ### optional output number of cell lines for each alteration
  agg <- aggregate(Sample_ID ~ Alteration, data = df,paste, collapse = ",")
  agg$numCelllines <- with(agg,str_count(agg$Sample_ID,",")+1)
  #####
  
  df <- agg %>% mutate(Sample_ID=str_split(Sample_ID, ","))
  df <- df %>% unnest(Sample_ID) %>% unique()
  newdf<- df
  newdf$numCelllines <- NULL
  newdf <- aggregate(Alteration ~ Sample_ID, data = newdf, paste, collapse = "\t")
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

# Construct mutation features and mutation+cancer_type features
# including 1) all mutations and 2) oncokb features only.

features_all_mut <- parse_mutations(mutations_all_freq, FALSE, sample_info)
features_all_ct <- parse_mutations(mutations_all_freq, TRUE, sample_info)

features_oncokb_mut <- parse_mutations(mutations_oncokb, FALSE, sample_info)
features_oncokb_ct <- parse_mutations(mutations_oncokb, TRUE, sample_info)

# write output files
output.file <- file(paste(output_dir,"features_all_MUT.tsv",sep=""), "wb")
write.table(features_all_mut,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(output_dir,"features_all_CT.tsv",sep=""), "wb")
write.table(features_all_ct,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(output_dir,"features_oncokb_MUT.tsv",sep=""), "wb")
write.table(features_oncokb_mut,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)

output.file <- file(paste(output_dir,"features_oncokb_CT.tsv",sep=""), "wb")
write.table(features_oncokb_ct,
            row.names = FALSE,
            col.names = TRUE,
            file = output.file,
            sep = "\t",
            quote=FALSE)
close(output.file)
