# SCRIPT TO CREATE THE 3-WAY INTERACTION MODEL
# MUT * CNA * STAGE



#####----- CHOOSING A COMPUTER -----#####
## CNIO:
# setwd("/home/llolle/Downloads/B_CANCER_SAMPLES/MSK_MET_DATASET/")
## HP LAPTOP: 
setwd("C:/Users/laiut/Documents/MASTER/PRACTICAS_TFM/B_CANCER_SAMPLES/MSK_MET_DATASET/")




#####----- LOADING LIBRARIES -----#####
library("stringr")
library("reshape2")
library("dplyr")
library("tibble")




#####----- MODEL IMPLEMENTATION -----#####

## THRESHOLD SETTNG ##
FREQ <- "MUTFREQ1-CNVFREQ5"
freqmut_threshold <- 0.01
freqcnv_threshold <- 0.05
# FREQ <- "MUTFREQ2-CNVFREQ5"
# freqmut_threshold <- 0.02
# freqcnv_threshold <- 0.05

## DATASET SPLITTING ##
# SPLITMOD <- "TISSUE"
# SPLITMOD <- "SUBTYPE"
SPLITMOD <- "TISSUE-STAGE"

## CHECKING DIRECTION ##
CHECKDIR <- "CHECKING-DIRECTION"

# FDR CORRECTION ##
FDR <- "FDR10"
# FDR <- "FDR20"
# FDR <- "FDR30"
# FDR <- "FDR40"




#####---- LOADING INPUT FILES ----#####
setwd("./MERGED_MATs_SP/")
matnames <- list.files(pattern="*Mut_CNA_matrix_input.txt")
cancernames <- unlist(lapply(matnames, function(x) strsplit(x,"_")[[1]][1]))
mutcnv_mats <- lapply(matnames,function(x) read.csv(x, sep = "\t", header = TRUE))
mutcnv_mats <- setNames(mutcnv_mats,cancernames)
setwd("..")

cancerGeneList <- read.csv("cancerGeneList.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Storing cancerGeneDF with genes that have a "Yes" in the OncoKB.Annotated column. Keeping OG and TSG information
cancgenedf <- cancerGeneList[cancerGeneList$OncoKB.Annotated == "Yes",][,c(1,9,10)]
cancgenedf <- cancgenedf[-which(cancgenedf$Is.Oncogene == "No" & cancgenedf$Is.Tumor.Suppressor.Gene == "No"),]
colnames(cancgenedf) <- c("Gene","Is.OG","Is.TSG")

setwd(sprintf("./PSEUDOCOUNTS_new/%s/%s/ANALYSIS_2WAY/SIGNIFICANT-GENES_2WAY_%s",SPLITMOD,FREQ,CHECKDIR))
# Retrieving sigificant genes from 2-way interaction model
gain_file_names <- list.files(pattern=sprintf("*_Gain_%s_%s.tsv",FDR,tolower(FREQ)))
gain_genes_2way <- lapply(gain_file_names,function(x) read.csv(x, sep="\t", header=TRUE, stringsAsFactors=FALSE)$Gene)
# Merging primary+metastasis significant 2-way genes
if(SPLITMOD=="TISSUE-STAGE"){
  gain_genes_2way <- lapply(seq(1,length(gain_genes_2way),by=2),function(ind){
    return(unique(c(gain_genes_2way[[ind]],gain_genes_2way[[ind+1]])))
  })
}
loss_file_names <- list.files(pattern=sprintf("*_Loss_%s_%s.tsv",FDR,tolower(FREQ)))
loss_genes_2way <- lapply(loss_file_names,function(x) read.csv(x, sep ="\t", header=TRUE, stringsAsFactors=FALSE)$Gene)
if(SPLITMOD=="TISSUE-STAGE"){
  loss_genes_2way <- lapply(seq(1,length(loss_genes_2way),by=2),function(ind){
    return(unique(c(loss_genes_2way[[ind]],loss_genes_2way[[ind+1]])))
  })
}

setwd("../..")




####---- FUNCTION DECLARATION ----#####

# NOT IN FUNCTION
"%!in%" <- Negate("%in%")

# FUNCTION TO FILTER LIST OF DATASETS BY 2WAY SIGNIFICANT GENES
filtering_datasets <- function(list_datasets, list_genes){
  # Returns a list with the filtered datasets
  new_list_datasets <- lapply(1:length(list_datasets),function(dataset_index){
    # Stores dataset
    dataset <- list_datasets[[dataset_index]]
    # Stores gene list
    genes_2way <- list_genes[[dataset_index]]
    # Converts stage column into binary
    stage <-ifelse(dataset$SAMPLE_TYPE=="Primary",0,1)
    # Retrieves column indexes for the genes of the dataframe that appear in the gene list
    gene_column_indexes <- sort(unlist(lapply(genes_2way, function(x) grep(paste0(x,"_"),colnames(dataset)))))
    # Filters dataset by column index
    dataset <- dataset[gene_column_indexes]
    # Adds stage column in binary
    if (ncol(dataset)!=0){dataset$STAGE <- stage}
    return(dataset)
  })
  new_list_datasets <- setNames(new_list_datasets, names(list_datasets))
}

# FUNCTION TO CREATE MODEL INPUT FROM BINARY MATRIX COLUMNS (MUT/LOSS/GAIN)
model_input <- function(df,freqmut_threshold,freqcnv_threshold){
  if(ncol(df)==0){
    inputs_table <- data.frame(matrix(nrow=0,ncol=19))
    colnames(inputs_table) <- c("MutFreq","LossFreq","GainFreq",
                                "NoMutLossPrim","NoMutLossMeta","MutLossPrim","MutLossMeta",
                                "NoMutWTPrim","NoMutWTMeta","MutWTPrim","MutWTMeta",
                                "NoMutGainPrim","NoMutGainMeta","MutGainPrim","MutGainMeta",
                                "Gene","Mutation_filter","Loss_filter","Gain_filter")
  }else{
  # vector with indexes that correspond to each mutation gene column without last index (STAGE column)
  muts <- seq(1,ncol(df)-1,by=3)
  # vector with genes that appear in the dataframe
  genes <- unlist(lapply(names(df)[muts], function(x) strsplit(x,"_")[[1]][1]))
  inputs_table <- as.data.frame(t(sapply(muts,function(mutcol){
    # Takes stage vector
    st <- df$STAGE
    # Takes mutation vector
    mut <- df[,mutcol]
    # Creates vector with loss vector and gain vector
    cnvs <- as.factor(paste0(df[,mutcol+1],df[,mutcol+2]))
    # Eliminates row where both gain/loss are 1/1, and eliminates that row from the mutation and the stage vector also
    if (nlevels(cnvs)==4){
      elimin_index <- which(cnvs=="11")
      st <- st[-elimin_index]
      mut <- mut[-elimin_index]
      cnvs <- cnvs[-elimin_index]}
    # Converts 10 into -1 (loss), 00 into o (WT), 01 into 1 (gain)
    cnvs <- as.character(cnvs)
    cnvs <- replace(cnvs, cnvs=="10",-1)
    cnvs <- replace(cnvs,cnvs=="00",0)
    cnvs <- replace(cnvs,cnvs=="01",1)
    # Vector with the possible occurrences
    possib_occur <- c("00-1","10-1","01-1","11-1","000","100","010","110","001","101","011","111")
    # Count of each occurrence
    N <- unlist(lapply(1:12,function(x) sum(paste0(st,mut,cnvs)==possib_occur[x])))
    # Adding pseudocounts
    N <- N + 1
    # Adding mutation filter column (boolean)
    mutfreq <- sum(mut)/length(mut)
    # Adding CNA loss filter column (boolean)
    lossfreq <- length(which(cnvs==-1))/length(cnvs)
    # Adding CNA gain filter column (boolean)
    gainfreq <- length(which(cnvs==1))/length(cnvs)
    return(c(mutfreq,lossfreq,gainfreq,N))})))
  colnames(inputs_table) <- c("MutFreq","LossFreq","GainFreq",
                              "NoMutLossPrim","NoMutLossMeta","MutLossPrim","MutLossMeta",
                              "NoMutWTPrim","NoMutWTMeta","MutWTPrim","MutWTMeta",
                              "NoMutGainPrim","NoMutGainMeta","MutGainPrim","MutGainMeta")
  # Adding gene names
  inputs_table$Gene <- genes
  # Adding mutation filter column (boolean)
  inputs_table$Mutation_filter <- inputs_table$MutFreq >= freqmut_threshold
  # Adding CNA loss filter column (boolean)
  inputs_table$Loss_filter <- inputs_table$LossFreq >= freqcnv_threshold
  # Adding CNA gain filter column (boolean)
  inputs_table$Gain_filter <- inputs_table$GainFreq >= freqcnv_threshold}
  return(inputs_table)}

# FUNCTION TO GENERATE REGRESSION MODEL
# Returns vector with "Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
# "NoMutWTPrim","NoMutWTMeta","MutWTPrim","MutWTMeta","NoMutCNVPrim","NoMutCNVMeta","MutCNVPrim","MutCNVMeta" columns
# The vector contains NAs for genes not trespassing the mutation frequency or cnv frequency threshold
reg_model <- function(dataf){
  # Stablishing model formula (3WAY)
  form <- "N ~ st*mut*cnv"
  # Creating model with input occurrences
  mod <- glm(form, family = poisson(link = "log"), data = dataf)
  result <- c(unname(summary(mod)$coef[8,]), 
                summary(mod)$null.deviance, 
                summary(mod)$deviance,
                dataf$N)
  return(result)}

# FUNCTION TO CREATE RESULT .TSV FILE
retrieving_results <- function(df,cnv_model) {
  if (nrow(df)==0){
    resulttable <- data.frame(matrix(nrow=0,ncol=15))
    return(resulttable)
  }else{
    if (cnv_model=="Loss"){
      # Only keeping genes that pass the mutation frequency threshold
      mut_df <- filter(df,Mutation_filter==TRUE)
      # For the loss model, keeping genes that pass the CNA loss frequency threshold
      loss_df <- filter(mut_df,Loss_filter==TRUE)
      # Calculate model results for each gene
      loss_regmod_df <- as.data.frame(t(sapply(loss_df$Gene, function(gene){
        # Retrieve gene row with occurrences
        datarow <- loss_df[loss_df$Gene==gene,]
        # Store occurrences in vector
        N <- c(datarow$NoMutWTPrim,datarow$NoMutWTMeta,datarow$MutWTPrim,datarow$MutWTMeta,datarow$NoMutLossPrim,datarow$NoMutLossMeta,datarow$MutLossPrim,datarow$MutLossMeta)
        # Put vector in df as input for the model
        dataf <- data.frame(st=c(0,1,0,1,0,1,0,1),mut=c(0,0,1,1,0,0,1,1),cnv=c(0,0,0,0,1,1,1,1), N=N)
        # Retrieving model results
        result <- reg_model(dataf)
        return(result)
      })))
      colnames(loss_regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                                    "NoMutWTPrim","NoMutWTMeta","MutWTPrim","MutWTMeta","NoMutLossPrim","NoMutLossMeta","MutLossPrim","MutLossMeta")
      # Adding gene names
      loss_regmod_df$Gene <- rownames(loss_regmod_df)
      # Adding information on function (Is.OG, Is.TSG)
      loss_regmod_df <- merge(loss_regmod_df,cancgenedf,by="Gene")
      return(loss_regmod_df)
    }else if(cnv_model=="Gain"){
      # Only keeping genes that pass the mutation frequency threshold
      mut_df <- filter(df,Mutation_filter==TRUE)
      # For the gain model, keeping genes that pass the CNA loss frequency threshold
      gain_df <- filter(mut_df,Gain_filter==TRUE)
      # Calculate model results for each gene
      gain_regmod_df <- as.data.frame(t(sapply(gain_df$Gene, function(gene){
        # Retrieve gene row with occurrences
        datarow <- gain_df[gain_df$Gene==gene,]
        # Store occurrences in vector
        N <- c(datarow$NoMutWTPrim,datarow$NoMutWTMeta,datarow$MutWTPrim,datarow$MutWTMeta,datarow$NoMutGainPrim,datarow$NoMutGainMeta,datarow$MutGainPrim,datarow$MutGainMeta)
        # Put vector in df as input for the model
        dataf <- data.frame(st=c(0,1,0,1,0,1,0,1),mut=c(0,0,1,1,0,0,1,1),cnv=c(0,0,0,0,1,1,1,1), N=N)
        # Retrieving model results
        result <- reg_model(dataf)
        return(result)
      })))
      if (ncol(gain_regmod_df)==0){gain_regmod_df<-data.frame(matrix(nrow=0,ncol=14))}
      colnames(gain_regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                                    "NoMutWTPrim","NoMutWTMeta","MutWTPrim","MutWTMeta","NoMutGainPrim","NoMutGainMeta","MutGainPrim","MutGainMeta")
      # Adding gene names
      gain_regmod_df$Gene <- rownames(gain_regmod_df)
      # Adding information on function (Is.OG, Is.TSG)
      gain_regmod_df <- merge(gain_regmod_df,cancgenedf,by="Gene")
      return(gain_regmod_df)
    }
  }
}




####---- MODEL IMPLEMENTATION ----####

# MODIFYING MUT-CNV MATRIXES ACCORDING TO THE MODEL (TISSUE/SUBTYPE/CANCER STAGE)
ready_datasets <- lapply(1:length(mutcnv_mats),function(matrix_index) {
  mutcnv_matrix <- mutcnv_mats[[matrix_index]]
  # Subsetting matrixes
  if (SPLITMOD == "TISSUE"){
    sub_list <- split(mutcnv_matrix,mutcnv_matrix$ORGAN_SYSTEM)
  }else if (SPLITMOD == "SUBTYPE"){
    sub_list <- split(mutcnv_matrix,mutcnv_matrix$SUBTYPE)
  }else if (SPLITMOD == "TISSUE-STAGE"){
    sub_list <- split(mutcnv_matrix,mutcnv_matrix$ORGAN_SYSTEM)
  }
  sub_list <- setNames(sub_list,lapply(names(sub_list),function(x) str_replace_all(x," ","-")))
  return(sub_list)
})
ready_datasets <- setNames(ready_datasets,cancernames)
ready_datasets <- unlist(ready_datasets,recursive=FALSE)

# Filtering datasets by 2way significant genes for each CNV 
# Significant genes for the Gain model and the Loss model differently
gain_datasets <- filtering_datasets(ready_datasets,gain_genes_2way)
loss_datasets <- filtering_datasets(ready_datasets,loss_genes_2way)

gain_inputs <- lapply(gain_datasets,model_input,freqmut_threshold,freqcnv_threshold)
loss_inputs <- lapply(loss_datasets,model_input,freqmut_threshold,freqcnv_threshold)

gain_results <- lapply(1:length(gain_inputs),function(df_index){
  print(df_index)
  df <- gain_inputs[[df_index]]
  results <- retrieving_results(df,"Gain")
  if(nrow(results)!=0){
    results$Tissue <- strsplit(names(gain_datasets)[[df_index]],"\\.")[[1]][1]
    results$Subtype <- strsplit(names(gain_datasets)[[df_index]],"\\.")[[1]][2]
    results$FDR_2way_signif <- FDR}
  return(results)})
gain_results <- setNames(gain_results,names(gain_datasets))
loss_results <- lapply(1:length(loss_inputs),function(df_index){
  df <- loss_inputs[[df_index]]
  results <- retrieving_results(df,"Loss")
  if(nrow(results)!=0){
    results$Tissue <- strsplit(names(loss_datasets)[[df_index]],"\\.")[[1]][1]
    results$Subtype <- strsplit(names(loss_datasets)[[df_index]],"\\.")[[1]][2]
    results$FDR_2way_signif <- FDR}
  return(results)})
loss_results <- setNames(loss_results,names(loss_datasets))


####---- STORING RESULTS ----####


folders <- c(sprintf("RESULTS_3WAY_%s_%s",CHECKDIR,FDR))
for (folder in folders){
  if (file.exists(folder)) {
    cat("The folder already exists\n")
  } else {
    dir.create(folder)
  }
  setwd(folder)}

sapply(1:length(gain_results),function(results_index){
  results <- gain_results[[results_index]]
  write.table(results,
              sprintf("%s_Gain_3way_%s_results.tsv",names(gain_results)[results_index],FDR),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)
  })

sapply(1:length(loss_results),function(results_index){
  results <- loss_results[[results_index]]
  write.table(results,
              sprintf("%s_Loss_3way_%s_results.tsv",names(loss_results)[results_index],FDR),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)})

setwd("..")

folders <- c("INPUTS_3WAY","Gain_inputs")
for (folder in folders){
  if (file.exists(folder)) {
    cat("The folder already exists\n")
  } else {
    dir.create(folder)
  }
  setwd(folder)}

sapply(1:length(gain_inputs), function(x){
  inputs <- gain_inputs[[x]]
  write.table(inputs,
              sprintf("GAIN_%s_%s_%s_%s_results.tsv",names(gain_inputs)[x],FDR,tolower(SPLITMOD),tolower(FREQ)),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)})

setwd("..")

folders <- c("Loss_inputs")
for (folder in folders){
  if (file.exists(folder)) {
    cat("The folder already exists\n")
  } else {
    dir.create(folder)
  }
  setwd(folder)}

sapply(1:length(loss_inputs), function(x){
  inputs <- loss_inputs[[x]]
  write.table(inputs,
              sprintf("LOSS_%s_%s_%s_%s_results.tsv",names(loss_inputs)[x],FDR,tolower(SPLITMOD),tolower(FREQ)),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)})
