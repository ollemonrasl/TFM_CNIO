# SCRIPT TO CREATE THE 3-WAY INTERACTION MODEL 
# Gene A MUT * Gene A CNA * Gene B MUT



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
SPLITMOD2 <- c("TISSUE-STAGE","tissue")
SPLITMOD1 <- c("TISSUE","tissue")
# SPLITMOD2 <- c("SUBTYPE-STAGE","subtype")
# SPLITMOD1 <- c("SUBTYPE","subtype")

## CHECKING DIRECTION ##
CHECKDIR <- "CHECKING-DIRECTION"

# FDR CORRECTION ##
FDR <- "FDR10"
# FDR <- "FDR20"
# FDR <- "FDR30"
# FDR <- "FDR40"




#####---- LOADING INPUT FILES ----#####
setwd(sprintf("./TWOGENES/%s",SPLITMOD1[2]))
matnames <- list.files(pattern="*_two_genes.txt")
cancernames <- str_replace_all(unlist(lapply(matnames,function(x) strsplit(x,"_3")[[1]][1]))," ","-")
twogene_mats <- lapply(matnames,function(x) read.csv(x, sep = "\t", header = TRUE))
twogene_mats <- setNames(twogene_mats,cancernames)
# Separating twogene mats into gain and loss according to the model
if(SPLITMOD2[2]=="tissue"){
  twogene_mats_gain <- twogene_mats[1:20]
  twogene_mats_loss <- twogene_mats[21:40]
}else if(SPLITMOD2[2]=="subtype"){
  twogene_mats_gain <- twogene_mats[1:100]
  twogene_mats_loss <- twogene_mats[101:200]
}
setwd("../..")

cancerGeneList <- read.csv("cancerGeneList.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Storing cancerGeneDF with genes that have a "Yes" in the OncoKB.Annotated column. Keeping OG and TSG information
cancgenedf <- filter(cancerGeneList, OncoKB.Annotated == "Yes")[,c(1,9,10)]
cancgenedf <- cancgenedf[-which(cancgenedf$Is.Oncogene == "No" & cancgenedf$Is.Tumor.Suppressor.Gene == "No"),]
colnames(cancgenedf) <- c("Gene","Is.OG","Is.TSG")


setwd(sprintf("./PSEUDOCOUNTS_new/%s/%s/ANALYSIS_2WAY/SIGNIFICANT-GENES_2WAY_%s",SPLITMOD1[1],FREQ,CHECKDIR))
# Retrieving sigificant genes from 2-way interaction model (tissue AND cancer stage subsetting)
# No cancer stage subsetting significant genes
gain_file_names_nostage <- list.files(pattern=sprintf("*_Gain_%s_%s.tsv",FDR,tolower(FREQ)))
signif_gain_names_nostage <- str_replace_all(unlist(lapply(gain_file_names_nostage,function(x){
  cnv <- strsplit(x,"_")[[1]][3]
  if (SPLITMOD1[1]=="TISSUE-STAGE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  }else if(SPLITMOD1[1]=="SUBTYPE-STAGE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  }else if(SPLITMOD1[1]=="TISSUE"){
    splitmod <- strsplit(x,"_")[[1]][2]
  }else if(SPLITMOD1[1]=="SUBTYPE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  }
  return(paste0(cnv,"_",splitmod))})
)," ","-")
gain_genes_2way_nostage <- lapply(gain_file_names_nostage,function(x) read.csv(x, sep="\t", header=TRUE, stringsAsFactors=FALSE)$Gene)
gain_genes_2way_nostage <- setNames(gain_genes_2way_nostage,signif_gain_names_nostage)
loss_file_names_nostage <- list.files(pattern=sprintf("*_Loss_%s_%s.tsv",FDR,tolower(FREQ)))
signif_loss_names_nostage <- str_replace_all(unlist(lapply(loss_file_names_nostage,function(x){
  cnv <- strsplit(x,"_")[[1]][3]
  if (SPLITMOD1[1]=="TISSUE-STAGE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  }else if(SPLITMOD1[1]=="SUBTYPE-STAGE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  } else if(SPLITMOD1[1]=="TISSUE"){
    splitmod <- strsplit(x,"_")[[1]][2]
  }else if(SPLITMOD1[1]=="SUBTYPE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  }
  return(paste0(cnv,"_",splitmod))})
)," ","-")
loss_genes_2way_nostage <- lapply(loss_file_names_nostage,function(x) read.csv(x, sep ="\t", header=TRUE, stringsAsFactors=FALSE)$Gene)
loss_genes_2way_nostage <- setNames(loss_genes_2way_nostage,signif_loss_names_nostage)

setwd("../../../../")

setwd(sprintf("./%s/%s/ANALYSIS_2WAY/SIGNIFICANT-GENES_2WAY_%s",SPLITMOD2[1],FREQ,CHECKDIR))
# Retrieving sigificant genes from 2-way interaction model (tissue AND cancer stage subsetting)
# Cancer stage subsetting significant genes
gain_file_names_stage <- list.files(pattern=sprintf("*_Gain_%s_%s.tsv",FDR,tolower(FREQ)))
signif_gain_names_stage <- str_replace_all(unlist(lapply(gain_file_names_stage,function(x){
  cnv <- strsplit(x,"_")[[1]][3]
  if (SPLITMOD2[1]=="TISSUE-STAGE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  }else if(SPLITMOD2[1]=="SUBTYPE-STAGE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  }else if(SPLITMOD2[1]=="TISSUE"){
    splitmod <- strsplit(x,"_")[[1]][2]
  }
  return(paste0(cnv,"_",splitmod))})
)," ","-")
gain_genes_2way_stage <- lapply(gain_file_names_stage,function(x) read.csv(x, sep="\t", header=TRUE, stringsAsFactors=FALSE)$Gene)
gain_genes_2way_stage <- setNames(gain_genes_2way_stage,signif_gain_names_stage)
loss_file_names_stage <- list.files(pattern=sprintf("*_Loss_%s_%s.tsv",FDR,tolower(FREQ)))
signif_loss_names_stage <- str_replace_all(unlist(lapply(loss_file_names_stage,function(x){
  cnv <- strsplit(x,"_")[[1]][3]
  if (SPLITMOD2[1]=="TISSUE-STAGE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  }else if(SPLITMOD2[1]=="SUBTYPE-STAGE"){
    splitmod <- str_replace_all(strsplit(x,"_")[[1]][2],"\\.","_")
  } else if(SPLITMOD2[1]=="TISSUE"){
    splitmod <- strsplit(x,"_")[[1]][2]
  }
  return(paste0(cnv,"_",splitmod))})
)," ","-")
loss_genes_2way_stage <- lapply(loss_file_names_stage,function(x) read.csv(x, sep ="\t", header=TRUE, stringsAsFactors=FALSE)$Gene)
loss_genes_2way_stage <- setNames(loss_genes_2way_stage,signif_loss_names_stage)

setwd("../../../../../")




####---- FUNCTION DECLARATION ----#####

# NOT IN FUNCTION
"%!in%" <- Negate("%in%")

# FUNCTION TO MERGE SIGNIFICANT GENES
merging_genes <- function(genes_2way_nostage, genes_2way_stage){
  # Returns list with significant genes for each model (tissue-cna)
  genes_2way <- lapply(1:length(genes_2way_nostage),function(index){
    # Stores genes from non-cancer stage separated 2way model
    stage_genes <- genes_2way_nostage[[index]]
    # Stogres genes from metastatic model
    met_genes <- genes_2way_stage[[index*2-1]]
    # Stores genes from primary model
    prim_genes <- genes_2way_stage[[index*2]]
    # Returns merged genes
    totalgenes <- unique(c(stage_genes,met_genes,prim_genes))
    return(totalgenes)
  })
  genes_2way <- setNames(genes_2way,names(genes_2way_nostage))
}

# FUNCTION TO FILTER 2GENE TABLES BY 2WAY SIGNIFICANT GENES
filtering_datasets <- function(twogene_mats,genes_2way,cnvmod){
  # Replicates gene sets (in order to have one for primary twogene mat and one for metastatic twogene mat)
  genes_2way <- rep(genes_2way,each=2)
  # Returns a list with the filtered datasets
  filtered_dataset <- lapply(1:length(twogene_mats),function(index){
    # Stores twogene matrix
    twogene_mat <- twogene_mats[[index]]
    # Stores gene list
    significant_2way_genes <- genes_2way[[index]]
    # Stores indexes of genes A in the twogene matrix that appear in the significant genes list
    indexes_siggenes <- sort(unlist(lapply(significant_2way_genes,function(gene){
      grep(paste0(gene,"_"),twogene_mat$Gene)})))
    # Filteres twogene matrix by indexes of significant genes
    filtered_mat <- twogene_mat[indexes_siggenes,]
    # Filters twogene matrix by mutation threshold for genes B
    filtered_mat <- filter(filtered_mat,filtered_mat$Freq_MutB >= freqmut_threshold*100)
    # Adds gene names column
    filtered_mat$Gene <- as.character(filtered_mat$Gene)
    # Vector with Genes B that appear in the filtered twogene matrix
    genesB <- unlist(lapply(filtered_mat$Gene, function(x) strsplit(x,"_")[[1]][2]))
    # Filtering for genes that appear in the OncoKB database
    filtered_mat <- filtered_mat[genesB %in% cancgenedf$Gene,]
    # Returns two-elements list (1 - significant 2way genes (A), 2- filtered twogene matrix)
    return(list(SIG_GENES = significant_2way_genes, FILT_MAT = filtered_mat))
  })
  if (cnvmod=="Gain"){
    # Filtering datasets by CNA gain threshold
    really_filtered_dataset <- lapply(filtered_dataset,function(list){
      filtered_dataframe <- list[[2]]
      filtered_dataframe <- filter(filtered_dataframe,filtered_dataframe$Freq_GainB >= freqcnv_threshold*100)
      return(list(SIG_GENES = list[[1]], FILT_MAT = filtered_dataframe))
    })
  }else if (cnvmod=="Loss"){
    # Filtering dataset by CNA loss threshold
    really_filtered_dataset <- lapply(filtered_dataset,function(list){
      filtered_dataframe <- list[[2]]
      filtered_dataframe <- filter(filtered_dataframe,filtered_dataframe$Freq_LossB >= freqcnv_threshold*100)
      return(list(SIG_GENES = list[[1]], FILT_MAT = filtered_dataframe))
    })
  }
  really_filtered_dataset <- setNames(really_filtered_dataset,names(twogene_mats))
  return(really_filtered_dataset)
}

# FUNCTION TO GENERATE REGRESSION MODEL
# Returns vector with "Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
# "MutACNVMutB","NoMutACNVMutB","MutAWTMutB","NoMutAWTMutB","MutACNVNoMutB","NoMutACNVNoMutB","MutAWTNoMutB","NoMutAWTNoMutB" columns
reg_model <- function(dataf){
  # Stablishing model formula (3WAY)
  form <- "N ~ mutA*cnvA*mutB"
  # Creating model with input occurrences
  mod <- glm(form, family = poisson(link = "log"), data = dataf)
  result <- c(unname(summary(mod)$coef[8,]), 
                summary(mod)$null.deviance, 
                summary(mod)$deviance,
                dataf$N)
  return(result)}

# FUNCTION TO CREATE RESULT .TSV FILE
retrieving_results <- function(filtered_mats_twogenes) {
  # Storing list of significant genes (A)
  signif_genes <- filtered_mats_twogenes[[1]]
  # Storing matrix 
  filtered_mat <- filtered_mats_twogenes[[2]]
  if (length(signif_genes)==0 | nrow(filtered_mat)==0){
    regmod_df <- data.frame(matrix(nrow=0,ncol=20))
    return(regmod_df)
  }else{
    # Calculate model results for each gene
    regmod_df <- as.data.frame(t(sapply(unique(filtered_mat$Gene),function(gene){
      # Retrieve gene row with occurrences
      datarows <- filtered_mat[filtered_mat$Gene==gene,]
      # Store occurrences in vector
      N <- melt(t(datarows[,c(11:14)]))$value
      # Add pseudocounts
      N <- N + 1
      # Put vector in df as input for the model
      dataf <-  data.frame(mutA=c(1,0,1,0,1,0,1,0),
                           cnvA=c(1,1,0,0,1,1,0,0),
                           mutB=c(1,1,1,1,0,0,0,0), 
                           N=N)
      # Retrieving model results
      results <- reg_model(dataf)
      return(results)
    })))
    colnames(regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                               "MutACNVMutB","NoMutACNVMutB","MutAWTMutB","NoMutAWTMutB","MutACNVNoMutB","NoMutACNVNoMutB","MutAWTNoMutB","NoMutAWTNoMutB")
    # Adding gene pair names
    regmod_df$Genes <- rownames(regmod_df)
    # Adding gene A names
    regmod_df$Gene <- unlist(lapply(regmod_df$Genes,function(x) strsplit(x,"_")[[1]][1]))
    # Adding information on function (Is.OG, Is.TSG) for gene A
    regmod_df <- merge(regmod_df,cancgenedf,by="Gene")
    # Renaming columns
    colnames(regmod_df)[c(1,17,18)] <- c("GeneA","Is.GeneA.OG","Is.GeneA.TSG")
    # Relocating columns
    regmod_df <- regmod_df[c(2:16,1,17,18)]
    # Adding gene B names
    regmod_df$Gene <- unlist(lapply(regmod_df$Genes,function(x) strsplit(x,"_")[[1]][2]))
    # Adding information on function (Is.OG, Is.TSG) for gene B
    regmod_df <- merge(regmod_df,cancgenedf,by="Gene")
    # Renaming columns
    colnames(regmod_df)[c(1,20,21)] <- c("GeneB","Is.GeneB.OG","Is.GeneB.TSG")
    # Relocating columns
    regmod_df <- regmod_df[c(2:19,1,20,21)]
  return(regmod_df)}
  }

####---- MODEL IMPLEMENTATION ----####

gain_genes_2way <- merging_genes(gain_genes_2way_nostage,gain_genes_2way_stage)
loss_genes_2way <- merging_genes(loss_genes_2way_nostage,loss_genes_2way_stage)

twogenes_gain_mats_filtered <- filtering_datasets(twogene_mats_gain,gain_genes_2way,"Gain")
twogenes_loss_mats_filtered <- filtering_datasets(twogene_mats_loss,loss_genes_2way,"Loss")

twogenes_gain_results <- lapply(twogenes_gain_mats_filtered,retrieving_results)
twogenes_loss_results <- lapply(twogenes_loss_mats_filtered,retrieving_results)


####---- STORING RESULTS ----####

setwd("./TWOGENES_new/")

folders <- c(sprintf("3WAYTWOGENES_from%s_%s_%s",FDR,SPLITMOD1[1],tolower(FREQ)),"INPUTS")
for (folder in folders){
  if (file.exists(folder)) {
    cat("The folder already exists\n")
  } else {
    dir.create(folder)
  }
  setwd(folder)}

sapply(1:length(twogenes_gain_mats_filtered),function(mat_index){
  write.table(twogenes_gain_mats_filtered[[mat_index]][[2]],
              sprintf("%s_TWOGENES_%s_%s_input.tsv",names(twogenes_gain_mats_filtered)[mat_index],FDR,tolower(SPLITMOD1[[1]])),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)
})
sapply(1:length(twogenes_loss_mats_filtered),function(mat_index){
  write.table(twogenes_loss_mats_filtered[[mat_index]][[2]],
              sprintf("%s_TWOGENES_%s_%s_input.tsv",names(twogenes_loss_mats_filtered)[mat_index],FDR,tolower(SPLITMOD1[[1]])),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)
})
setwd("..")

folders <- c("RESULTS")
for (folder in folders){
  if (file.exists(folder)) {
    cat("The folder already exists\n")
  } else {
    dir.create(folder)
  }
  setwd(folder)}

sapply(1:length(twogenes_gain_results),function(results_index){
  write.table(twogenes_gain_results[[results_index]],
              sprintf("%s_TWOGENES_%s_%s_results.tsv",names(twogenes_gain_results)[results_index],FDR,tolower(SPLITMOD1[[1]])),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)
  })
sapply(1:length(twogenes_loss_results),function(results_index){
  write.table(twogenes_loss_results[[results_index]],
              sprintf("%s_TWOGENES_%s_%s_results.tsv",names(twogenes_loss_results)[results_index],FDR,tolower(SPLITMOD1[[1]])),
              sep="\t",
              quote=FALSE,
              row.names = FALSE)
})

