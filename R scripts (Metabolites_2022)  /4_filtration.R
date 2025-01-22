##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
library(data.table)
library(dplyr)
library(stringr)

setwd("/4_filtration/")

##############################################################################################################################################################
# Metadata generating & repeated measurements calculation
##############################################################################################################################################################

dsr <-as.data.frame(fread(input = "xcms after pos hops_recalc MVI QC-XGB.csv", head = T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]


# metadata generating
rname_all <- rownames(dsr) # obtain all info from rownames
qc_id <- grep(pattern = "qc", x = rname_all) # find by pattern in info from rownames
dsr_qc <- dsr[qc_id,] # select only observations by pattern from peak table
rname <- rownames(dsr_qc) # obtain all info from rownames
rname <- str_remove(rname, ".mzXML") # remove some pattern from vendor-specific format
all_id <- as.data.frame(t(sapply(1:length(rname), function(y) unlist(str_split(rname[y], "_"))))) # split info from rownames
all_id$V2 <- str_remove_all(all_id$V2,"bn")
b_id <- all_id$V2 # obtain batch ID (3 row)

ds_qc <- cbind(b_id, dsr_qc)
ds <- ds_qc
colnames(ds)[1] <- "batch"


tn <- ds$batch # batch variables
tnu <- unique(tn)
vb <- list()
for (i in (1:length(tnu))){
  vb[[i]] <- which(tn == tnu[i])}

# RSD calculation by batches
RSD_by_batch_results <- list()
RSD_by_batch <- function (x) { 
  for (i in (1:length(vb))){
    RSD_by_batch_results[[i]] <- apply(x[vb[[i]],], 2, function(y) sd(y, na.rm = T)/mean(y, na.rm = T)*100)}
  RSD_by_batch_results[[(length(vb)+1)]] <- apply(x, 2, function(y) sd(y, na.rm = T)/mean(y, na.rm = T)*100) # all batches
  for (i in (1:length(RSD_by_batch_results))){
    RSD_by_batch_results[[i]] <- data.frame(colnames(x),RSD_by_batch_results[[i]])}
  for (i in (1:length(RSD_by_batch_results))){
    colnames(RSD_by_batch_results[[i]]) <- c("name", "rsd")}
  n <- c(paste(c(1:length(vb)), "batch"), "all batches")
  names(RSD_by_batch_results) <- n
  RSD_by_batch_results
}

n <- c(paste(c(1:length(vb)), "batch"), "all batches") # names
ds1 <- ds[,-1]# sample in row only metabolites variables
rsd <- RSD_by_batch(ds1)
rsd_all_b <- rsd[[length(rsd)]]

# cutoff by RSD in QC
cutoff <- 30
cutoff_all <- lapply(1:length(rsd), function(y) nrow(dplyr::filter(rsd[[y]], rsd < cutoff)))
cutoff_f <- filter(rsd_all_b, rsd <= cutoff)
nrow(cutoff_f) # number of retained features
round(nrow(cutoff_f)/ncol(dsr)*100, 0) # percent of retained features
nan_l <- length(which(is.nan(rsd_all_b$rsd))) # number of NaN (when 0 values) 
dsr_QC_RSD <- dsr[,cutoff_f$name]

# save
fwrite(dsr_QC_RSD, "xcms after IPO MVI XGBOOST filter by QC_pos.csv", row.names = T)
