##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# Metadata generating
# Perform correction
# Evaluate correction
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################
BiocManager::install("ProteoMM")
BiocManager::install("waveICA")
BiocManager::install("sva")
BiocManager::install("pmp")
BiocManager::install("statTarget")

# setup environment
library(data.table)
library(dplyr)
library(stringr)
setwd("")

##############################################################################################################################################################
# Metadata generating
##############################################################################################################################################################

# load peak table in csv
# load in format: 1st column: 1. qc1(s78_1) b1 QC1(KMS 53).X --- ro. nameN(nameN_NREPEAT) BatchIndex QCN(Label N)
# "150. qc30 b6 QC30.lcd" 
# "167. s64_2 b7 MS 45.lcd"
#wd_1_pos <- c("C:/Users/Plyush/Yura/NAVUKA/rats_dir/ipo/qc_pos") # folder with files for IPO optimization process
wd_2_pos <- c("") # folder with all study samples

#wd_1_neg <- c("C:/Users/Plyush/Yura/NAVUKA/rats_dir/ipo/qc_neg") # folder with files for IPO optimization process
wd_2_neg <- c("") # folder with all study samples

#files_QC_pos <- list.files(wd_1_pos, recursive = TRUE, full.names = TRUE, pattern = ".CDF") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.
#files_QC_neg <- list.files(wd_1_neg, recursive = TRUE, full.names = TRUE, pattern = ".CDF") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.


files_all_pos <- list.files(wd_2_pos, recursive = TRUE, full.names = TRUE, pattern = ".mzXML") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.
files_all_neg <- list.files(wd_2_neg, recursive = TRUE, full.names = TRUE, pattern = ".mzXML") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.


dsr <-as.data.frame(fread(input = "hops_recalc MVI rf2_pos.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]

#rname <- rownames(dsr)
#rname <- str_remove(rname,"_hops_Seg1Ev2.mzXML")

#rname <- as.data.frame(cbind(rname,c(1:length(rname))))

pd <- data.frame(sample_name = sub(basename(files_all_neg), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE), stringsAsFactors = FALSE) # download filenames

rname <- pd
all_id <- sapply(1:nrow(rname), function(y) unlist(str_split(rname[y,], "_"))) # split info from rownames
ro_id <- unlist(lapply(1:length(all_id), function(y) return(all_id[[y]][1])))
ro_id <- str_remove(ro_id,"ro")

batch_id <- unlist(lapply(1:length(all_id), function(y) return(all_id[[y]][2])))
batch_id <- str_remove(batch_id,"bn")

type_id <- unlist(lapply(1:length(all_id), function(y) return(all_id[[y]][3])))

info <- cbind(ro_id,batch_id)
info <- cbind(info,type_id)
dsm <- cbind(info,dsr)

##### We are actually using this for veryfication
###############################################################################
########################################## QC-GB (xgboost/catboost)
###############################################################################

# generate batch data
batch <- dsm$batch_id
#batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
class <- as.character(dsm$type_id)

# generate run order data
order <- as.numeric(dsm$ro_id)

# Sample data
s_d <- data.frame(cbind(sample = rownames(dsm), batch, class, order))
rownames(s_d) <- rownames(dsm)

# Feature data
f_d <- data.frame(dsm[,-c(1:3)])

# Parameters for catboost
fit_params <- list(
  iterations = 100,
  loss_function = 'RMSE',
  border_count = 32,
  depth = 2,
  learning_rate = 0.03,
  l2_leaf_reg = 3.5,
  train_dir = 'train_dir')

################################################## Single mode
QC.SM.GB <- function(int_data, order, class, qc_label, model = "xgboost", max.depth = 2, nrounds = 100, params) {
  
  library(pbapply)
  
  qc_id <-  grep(qc_label, class)
  int_data_ro <- cbind(order, int_data)
  int_data_ro_qc <- int_data_ro[qc_id,]
  ro_qc <- int_data_ro[qc_id,1]
  ro <- order
  print("Start Correction")
  
  if (model == "catboost") {
    library(catboost)
    catb_train <- pblapply(2:ncol(f_d_ro_qc), function(t) catboost.load_pool(data = as.matrix(as.numeric(ro_qc)), label = as.matrix(as.numeric(f_d_ro_qc[,t]))))
    catb_test <- pblapply(2:ncol(f_d_ro), function(t) catboost.load_pool(data = as.matrix(as.numeric(ro)), label = as.matrix(as.numeric(f_d_ro[,t]))))
    
    fit_catb <- pblapply(1:length(catb_train), function(t) catboost.train(catb_train[[t]], params = params))
    predict_gb <- pblapply(1:length(fit_catb), function(t) catboost.predict(fit_catb[[t]], catb_test[[t]], verbose = F)) 
  } else {
    library(xgboost)
    xgb_train <- pblapply(2:ncol(int_data_ro_qc), function(t) xgb.DMatrix(data = as.matrix(as.numeric(ro_qc)), label = as.matrix(as.numeric(int_data_ro_qc[,t]))))
    xgb_test <- pblapply(2:ncol(int_data_ro), function(t) xgb.DMatrix(data = as.matrix(as.numeric(ro)), label = as.matrix(as.numeric(int_data_ro[,t]))))
    
    fit_xgb <- pblapply(1:length(xgb_train), function(t) xgboost(data = xgb_train[[t]], max.depth = max.depth, nrounds = nrounds, verbose = 0))
    predict_gb <- pblapply(1:length(fit_xgb), function(t) predict(fit_xgb[[t]], xgb_test[[t]]))
  } 
  
  val_sm <- pblapply(1:ncol(int_data), function(t) int_data[,t]/predict_gb[[t]])
  print("Done Correction")
  res_sm <- as.data.frame(t(do.call(rbind, val_sm)))
  rownames(res_sm) <- rownames(int_data)
  colnames(res_sm) <- colnames(int_data)
  res_sm <- res_sm*1000
  return(res_sm)
}

######################### perform xgboost in single mode
# Parameters of QC.SM.GB function:  
# int_data -> numeric feature table, order -> numeric vector of the run order, class -> sample group variable, 
# qc_label -> label for QC samples in group, 
# model -> type of the model ("xgboost" for xgboost, "catboost" - catboost), 
# max_depth -> maximum depth of a tree in xgboost model, nrounds	-> max number of boosting iterations in xgboost, params -> parameters for the catboost model.
qc_sm <- QC.SM.GB(int_data = f_d, order = order, class = class, qc_label = "qc", model = "xgboost", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_sm, "xcms after pos hops_recalc MVI QC-XGB.csv", row.names = T)
