# SCRIPT TO CREATE THE 2-WAY INTERACTION MODEL




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




#####----- LOADING INPUT FILES -----#####
setwd("./MERGED_MATs_SP/")
matnames <- list.files(pattern="*Mut_CNA_matrix_input.txt")
cancernames <- unlist(lapply(matnames, function(x) strsplit(x,"_")[[1]][1]))
mutcnv_mats <- lapply(matnames,function(x) read.csv(x, sep = "\t", header = TRUE))
mutcnv_mats <- setNames(mutcnv_mats,cancernames)
setwd("..")

cancerGeneList <- read.csv("cancerGeneList.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Storing cancerGeneDF with genes that have a "Yes" in the OncoKB.Annotated column. Keeping OG and TSG information
cancgenedf <- filter(cancerGeneList, OncoKB.Annotated == "Yes")[,c(1,9,10)]
cancgenedf <- cancgenedf[-which(cancgenedf$Is.Oncogene == "No" & cancgenedf$Is.Tumor.Suppressor.Gene == "No"),]
colnames(cancgenedf) <- c("Gene","Is.OG","Is.TSG")




#####----- FUNCTION DECLARATION -----#####

# NOT IN FUNCTION
"%!in%" <- Negate("%in%")

# FUNCTION TO CREATE MODEL INPUT FROM BINARY MATRIX COLUMNS (MUT/LOSS/GAIN)
model_input <- function(df,freqmut_threshold,freqcnv_threshold){
  # vector with indexes that correspond to each mutation gene column without last index (STAGE column)
  muts <- seq(1,ncol(df),by=3)
  # vector with genes that appear in the dataframe
  genes <- unlist(lapply(names(df)[muts], function(x) strsplit(x,"_")[[1]][1]))
  inputs_table <- as.data.frame(t(sapply(muts,function(mutcol){
    # Takes mutation vector
    mut <- df[,mutcol]
    # Creates vector with loss vector and gain vector
    cnvs <- as.factor(paste0(df[,mutcol+1],df[,mutcol+2]))
    # Eliminates row where both gain/loss are 1/1, and eliminates that row from the mutation vector also
    if (nlevels(cnvs)==4){
      elimin_index <- which(cnvs=="11")
      mut <- mut[-elimin_index]
      cnvs <- cnvs[-elimin_index]}
    # Converts 10 into -1 (loss), 00 into o (WT), 01 into 1 (gain)
    cnvs <- as.character(cnvs)
    cnvs <- replace(cnvs, cnvs=="10",-1)
    cnvs <- replace(cnvs,cnvs=="00",0)
    cnvs <- replace(cnvs,cnvs=="01",1)
    # Vector with the possible occurrences
    possib_occur <- c("0-1","1-1","00","10","01","11")
    # Count of each occurrence
    N <- unlist(lapply(1:6,function(x) sum(paste0(mut,cnvs)==possib_occur[x])))
    # Adding pseudocounts
    N <- N + 1
    # Calculating mutation frequency from binary vector
    mutfreq <- sum(mut)/length(mut)
    # Calculating CNA loss frequency from binary vector
    lossfreq <- length(which(cnvs==-1))/length(cnvs)
    # Calculating CNA gain frequency from binary vector
    gainfreq <- length(which(cnvs==1))/length(cnvs)
    return(c(mutfreq,lossfreq,gainfreq,N))})))
  colnames(inputs_table) <- c("MutFreq","LossFreq","GainFreq","NoMutLoss","MutLoss","NoMutWT","MutWT","NoMutGain","MutGain")
  # Adding gene names
  inputs_table$Gene <- genes
  # Adding mutation filter column (boolean)
  inputs_table$Mutation_filter <- inputs_table$MutFreq >= freqmut_threshold
  # Adding CNA loss filter column (boolean)
  inputs_table$Loss_filter <- inputs_table$LossFreq >= freqcnv_threshold
  # Adding CNA gain filter column (boolean)
  inputs_table$Gain_filter <- inputs_table$GainFreq >= freqcnv_threshold
  return(inputs_table)}

# FUNCTION TO GENERATE REGRESSION MODEL
# Returns vector with "Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance","NoMutWT","MutWT","NoMutCNV","MutCNV" columns
# The vector contains NAs for genes not trespassing the mutation frequency or cnv frequency threshold
reg_model <- function(dataf){
  # Stablishing model formula (2WAY)
  form <- "N ~ mut*cnv"
  # Creating model with input occurrences
  mod <- glm(form, family = poisson(link = "log"), data = dataf)
  result <- c(unname(summary(mod)$coef[4,]), 
              summary(mod)$null.deviance, 
              summary(mod)$deviance,
              dataf$N)
  return(result)}

# FUNCTION TO CREATE RESULT .TSV FILE
retrieve_results <- function(df_index) {
  # Retrieving dataset frame
  df <- model_inputs[[df_index]]
  # Only keeping genes that pass the mutation frequency threshold
  mut_df <- filter(df,Mutation_filter==TRUE)
  # For the loss model, keeping genes that pass the CNA loss frequency threshold
  loss_df <- filter(mut_df,Loss_filter==TRUE)
  # If the dataframe is empty, return an empty df
  if(nrow(loss_df)==0){
    loss_regmod_df <- data.frame(matrix(nrow=0,ncol=13))
    colnames(loss_regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                                  "NoMutWT","MutWT","NoMutCNV","MutCNV","Gene","Is.OG","Is.TSG")
  }else{
  # If the dataframe is not empty, calculate model results for each gene
  loss_regmod_df <- as.data.frame(t(sapply(loss_df$Gene, function(gene){
    # Retrieve gene row with occurrences
    datarow <- loss_df[loss_df$Gene==gene,]
    # Store occurrences in vector
    N <- c(datarow$NoMutWT,datarow$MutWT,datarow$NoMutLoss,datarow$MutLoss)
    # Put vector in df as input for the model
    dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,-1,-1), N=N)
    # Retrieving model results
    result <- reg_model(dataf)
    return(result)
  })))
  colnames(loss_regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                             "NoMutWT","MutWT","NoMutCNV","MutCNV")
  # Adding gene names
  loss_regmod_df$Gene <- rownames(loss_regmod_df)
  # Adding information on function (Is.OG, Is.TSG)
  loss_regmod_df <- merge(loss_regmod_df,cancgenedf,by="Gene")}
  # For the gain model, keeping genes that pass the CNA gain frequency threshold
  gain_df <- filter(mut_df,Gain_filter==TRUE)
  # If the dataframe is empty, return an empty df
  if(nrow(gain_df)==0){
    gain_regmod_df <- data.frame(matrix(nrow=0,ncol=13))
    colnames(gain_regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                                  "NoMutWT","MutWT","NoMutCNV","MutCNV","Gene","Is.OG","Is.TSG")
  }else{
  # If the dataframe is not empty, calculate model results for each gene
  gain_regmod_df <- as.data.frame(t(sapply(gain_df$Gene, function(gene){
    # Retrieve gene row with occurrences
    datarow <- gain_df[gain_df$Gene==gene,]
    # Store occurrences in vector
    N <- c(datarow$NoMutWT,datarow$MutWT,datarow$NoMutGain,datarow$MutGain)
    # Put vector in df as input for the model
    dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,1,1), N=N)
    # Retrieving model results
    result <- reg_model(dataf)
    return(result)
  })))
  colnames(gain_regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                                "NoMutWT","MutWT","NoMutCNV","MutCNV")
  # Adding gene names
  gain_regmod_df$Gene <- rownames(gain_regmod_df)
  # Adding information on function (Is.OG, Is.TSG)
  gain_regmod_df <- merge(gain_regmod_df,cancgenedf,by="Gene")}
  return(list(Loss=loss_regmod_df,Gain=gain_regmod_df))
  }




#####----- DATA PREPARATION -----#####

# MODIFYING MUT-MAT MATRIX BY FLTERING BY GENES THAT ARE ANNOTATED IN THE ONCOKB DB
oncokb_filtered_datasets  <- lapply(mutcnv_mats,function(mutcnv_matrix){
  # Storing mut-mat gene names (3x each (gene-mut, gene-gain, gene-loss))
  mat_genes <- unlist(lapply(names(mutcnv_matrix),function(x) strsplit(x,"_")[[1]][1]))
  # Array with indexes of the genes that DO NOT appear in OncoKB database
  notfound_index <- which(mat_genes %!in% cancgenedf$Gene)
  # Eliminating genes (COLUMNS) that do not appear in the OncoKB database
  filtered_matrix <- mutcnv_matrix[-notfound_index]
  # Appending tissue, subtype, stage,  to dataframe
  filtered_matrix$TISSUE <- mutcnv_matrix$ORGAN_SYSTEM
  filtered_matrix$SUBTYPE <- mutcnv_matrix$SUBTYPE
  filtered_matrix$STAGE <- mutcnv_matrix$SAMPLE_TYPE
  return(filtered_matrix)})




#####----- MODEL IMPLEMENTATION -----#####

## DATASET SPLITTING ##
SPLITMOD <- "TISSUE"
# SPLITMOD <- "SUBTYPE"
# SPLITMOD <- "TISSUE-STAGE"
# SPLITMOD <- "SUBTYPE-STAGE"


## THRESHOLD SETTNG ##
FREQ <- "MUTFREQ1-CNVFREQ5"
freqmut_threshold <- 0.01
freqcnv_threshold <- 0.05

# FREQ <- "MUTFREQ2-CNVFREQ5"
# freqmut_threshold <- 0.02
# freqcnv_threshold <- 0.05


ready_datasets <- unlist(lapply(oncokb_filtered_datasets,function(filtered_matrix) {
  finalcol <- ncol(filtered_matrix)-3
  # Subsetting matrixes
  if (SPLITMOD == "TISSUE"){
    sub_list <- split(filtered_matrix,filtered_matrix$TISSUE)
  }else if (SPLITMOD == "SUBTYPE"){
    sub_list <- split(filtered_matrix,filtered_matrix$SUBTYPE)
  }else if (SPLITMOD == "TISSUE-STAGE"){
    sub_list <- split(filtered_matrix,filtered_matrix$STAGE)
  }else if (SPLITMOD == "SUBTYPE-STAGE"){
    sub_list <- split(filtered_matrix,list(filtered_matrix$SUBTYPE,filtered_matrix$STAGE))}
  # Cuttling last rows of the original dataframe, which contain
  # ORGAN_SYSTEM, SUBTYPE, SAMPLE_TYPE
  sub_list <- lapply(sub_list, function(d){d <- d[,1:finalcol]})
  sub_list <- setNames(sub_list,lapply(names(sub_list),function(x) str_replace_all(x," ","-")))
}),recursive=FALSE)

model_inputs <- lapply(ready_datasets,model_input,freqmut_threshold,freqcnv_threshold)

model_results <- lapply(1:length(model_inputs),retrieve_results)
model_results <- setNames(model_results,names(model_inputs))




#####----- MODEL STORING -----#####

folders <- c("PSEUDOCOUNTS_new",SPLITMOD,FREQ,"RESULTS_2WAY")
for (folder in folders){
  if (file.exists(folder)) {
    cat("The folder already exists\n")
  } else {
    dir.create(folder)
  }
  setwd(folder)}

sapply(1:length(model_results), function(x){
  for (cnv in 1:2){
    if (cnv == 1){cna = "Loss"}else if (cnv == 2){cna = "Gain"}
    results <- model_results[[x]][[cnv]]
    nameslist <- strsplit(names(model_results)[x],"\\.")[[1]]
    if (nrow(results)==0){
      newcols <- data.frame(matrix(ncol=4,nrow=0))
      colnames(newcols) <- c("Tssue","Subtype","Stage","CNA_type")
      results <- cbind(results,newcols)
    }else{
      results$Tissue <- nameslist[1]
      if (nchar(SPLITMOD) == 6){
        results$Subtype <- NA
        results$Stage <- NA
      }else if(nchar(SPLITMOD) == 7){
        results$Subtype <- nameslist[2]
        results$Stage <- NA
      }else if(nchar(SPLITMOD) == 12){
        results$Stage <- nameslist[2]
        results$Subtype <- NA
      }else if(nchar(SPLITMOD) == 13){
        results$Subtype <- nameslist[2]
        results$Stage <- nameslist[3]
      }
      results$CNA_type <- cna}
    write.table(results,
                sprintf("%s_%s_%s_%s_results.tsv",names(model_results)[x],cna,tolower(SPLITMOD),tolower(FREQ)),
                sep="\t",
                quote=FALSE,
                row.names = FALSE)}})


Poss_genes <- sapply(ready_datasets, function(x) ncol(x)/3)
Loss_genes <- sapply(model_results, function(x) nrow(x$Loss))
Gain_genes <- sapply(model_results, function(x) nrow(x$Gain))

ngenes_table <- data.frame(Model = names(Poss_genes),
                           N_possible_tested_genes = Poss_genes,
                           Loss_tested_genes = Loss_genes,
                           Gain_tested_genes = Gain_genes)

write.table(ngenes_table,
            sprintf("number_genes_%s_%s.tsv",tolower(SPLITMOD),tolower(FREQ)),
            sep="\t",
            quote=FALSE,
            row.names = FALSE)

setwd("..")

folders <- c("INPUTS_2WAY")
for (folder in folders){
  if (file.exists(folder)) {
    cat("The folder already exists\n")
  } else {
    dir.create(folder)
  }
  setwd(folder)}

sapply(1:length(model_inputs), function(x){
  inputs <- model_inputs[[x]]
  write.table(inputs,
              sprintf("%s_%s_%s_results.tsv",names(model_inputs)[x],tolower(SPLITMOD),tolower(FREQ)),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)})
write.table(bind_rows(model_inputs),
            sprintf("TOTAL-INPUTS_%s_%s_results.tsv",tolower(SPLITMOD),tolower(FREQ)),
            sep="\t",
            quote=FALSE,
            row.names = FALSE)

folders <- c("Gain_inputs")
for (folder in folders){
  if (file.exists(folder)) {
    cat("The folder already exists\n")
  } else {
    dir.create(folder)
  }
  setwd(folder)}

sapply(1:length(model_inputs), function(x){
  inputs <- model_inputs[[x]]
  inputs <- filter(inputs, Mutation_filter==TRUE & Gain_filter == TRUE)
  write.table(inputs,
              sprintf("GAIN_%s_%s_%s_results.tsv",names(model_inputs)[x],tolower(SPLITMOD),tolower(FREQ)),
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

sapply(1:length(model_inputs), function(x){
  inputs <- model_inputs[[x]]
  inputs <- filter(inputs, Mutation_filter==TRUE & Loss_filter == TRUE)
  write.table(inputs,
              sprintf("LOSS_%s_%s_%s_results.tsv",names(model_inputs)[x],tolower(SPLITMOD),tolower(FREQ)),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)})
