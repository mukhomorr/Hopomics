##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# IPO for XCMS params optimization
# XCMS with the best params
# Description of peaks table
  
##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
setwd("")

# folder with all study samples
wd_2_pos <- c("C:/Users/Plyush/Yura/NAVUKA/Hop and beer metabolomics/recacl_hops_dir/POS all samples no bl") 
wd_2_neg <- c("C:/Users/Plyush/Yura/NAVUKA/Hop and beer metabolomics/recacl_hops_dir/NEG all samples no bl") 

files_QC_pos <- list.files(wd_1_pos, recursive = TRUE, full.names = TRUE, pattern = ".mzXML") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.
files_QC_neg <- list.files(wd_1_neg, recursive = TRUE, full.names = TRUE, pattern = ".mzXML") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.


files_all_pos <- list.files(wd_2_pos, recursive = TRUE, full.names = TRUE, pattern = ".mzXML") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.
files_all_neg <- list.files(wd_2_neg, recursive = TRUE, full.names = TRUE, pattern = ".mzXML") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.

##############################################################################################################################################################
# IPO for XCMS params optimization
##############################################################################################################################################################

library(IPO)
library(xcms)
library(doParallel)

nCore=detectCores()-1

#PeakPickingParameters
peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
peakpickingParameters$noise=c(100,1000)
peakpickingParameters$value_of_prefilter=c(3,800)
peakpickingParameters$min_peakwidth<- c(3,15)
peakpickingParameters$max_peakwidth<- c(30,40)
peakpickingParameters$ppm<- c(15,35)
param=SnowParam(workers = nCore)
nSlaves=1 # or nCore or more
###### POS #######
resultPeakpicking_pos <-
  optimizeXcmsSet(files = files_QC_pos,
                  params = peakpickingParameters,
                  BPPARAM = param,
                  nSlaves = nSlaves,
                  subdir = NULL,
                  plot = TRUE)

optimizedXcmsSetObject_pos <- resultPeakpicking_pos$best_settings$xset
save(resultPeakpicking_pos, file = "IPO_optimiz_xcms_mzXML_5QC_pos.RData")
####### NEG #######

resultPeakpicking_neg <-
  optimizeXcmsSet(files = files_QC_neg,
                  params = peakpickingParameters,
                  BPPARAM = param,
                  nSlaves = nSlaves,
                  subdir = NULL,
                  plot = TRUE)

optimizedXcmsSetObject_neg <- resultPeakpicking_neg$best_settings$xset
save(resultPeakpicking_neg, file = "IPO_optimiz_xcms_mzXML_5QC_neg.RData")

# Retention Time Alignment Optimization
retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- c(0.33,1)
retcorGroupParameters$gapExtend <- c(2.0,3.0)
retcorGroupParameters$minfrac=c(0.05,0.5)
retcorGroupParameters$response=c(9.0,18.0)
retcorGroupParameters$gapInit=c(0.10, 0.50)
retcorGroupParameters$mzwid=c(0.020, 0.050)
BiocParallel::register(BiocParallel::SerialParam())
nSlaves=1 # or nCore or more

###### POS #######
resultRetcorGroup_pos <-
  optimizeRetGroup(xset = optimizedXcmsSetObject_pos,
                   params = retcorGroupParameters,
                   nSlaves = nSlaves,
                   subdir = NULL,
                   plot = TRUE)

resultRetcorGroup_pos$best_settings
save(resultRetcorGroup_pos, file = "IPO_optimiz_ret_align_mzXML_5QC_rats_pos.RData")

###### NEG #######

resultRetcorGroup_neg <-
  optimizeRetGroup(xset = optimizedXcmsSetObject_neg,
                   params = retcorGroupParameters,
                   nSlaves = nSlaves,
                   subdir = NULL,
                   plot = TRUE)

resultRetcorGroup_neg$best_settings
save(resultRetcorGroup_neg, file = "IPO_optimiz_ret_align_mzXML_5QC_rats_neg.RData")

# Get final results
###### POS #######

writeRScript(resultPeakpicking_pos$best_settings$parameters, resultRetcorGroup_pos$best_settings)
param_pos <- c(resultPeakpicking_pos$best_settings$parameters, resultRetcorGroup_pos$best_settings)
save(param_pos,file = 'all params IPO mzXML 5QC  pos.RData')

funs_params_pos <- capture.output(writeRScript(resultPeakpicking_pos$best_settings$parameters, resultRetcorGroup_pos$best_settings), type = "message")
save(funs_params_pos,file = 'funs params mzXML 5QC neg .RData')

###### NEG #######

writeRScript(resultPeakpicking_neg$best_settings$parameters, resultRetcorGroup_neg$best_settings)
param_neg <- c(resultPeakpicking_neg$best_settings$parameters, resultRetcorGroup_neg$best_settings)
save(param_neg,file = 'all params IPO 5QC  neg.RData')

funs_params_neg <- capture.output(writeRScript(resultPeakpicking_neg$best_settings$parameters, resultRetcorGroup_neg$best_settings), type = "message")
save(funs_params_neg,file = 'funs params IPO 5QC rats.RData')

##############################################################################################################################################################
# XCMS with the best params
##############################################################################################################################################################

library(xcms)
library(data.table)
library(dplyr)
library(stringr)

# load best parameters from IPO for both polarities
# load("IPO_optimiz_xcms_CDF_7QC.RData")
# load("IPO_optimiz_ret_align_CDF_7QC.RData")

###### 

# Create a phenodata data.frame
pd <- data.frame(sample_name = sub(basename(files_all_pos), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE), stringsAsFactors = FALSE) # download filenames

rname <- pd
all_id <- sapply(1:nrow(rname), function(y) unlist(str_split(rname[y,], "_"))) # split info from rownames
ro_id <- unlist(lapply(1:length(all_id), function(y) return(all_id[[y]][1])))
ro_id <- str_remove(ro_id,"ro")
type_id <- unlist(lapply(1:length(all_id), function(y) return(all_id[[y]][3])))
info <- cbind(ro_id,type_id)
n_gr_t <- cbind(pd,info)
rownames(n_gr_t) <- n_gr_t$sample_name
n_gr_t$group_id <- n_gr_t$type_id[which(type_id != "qc")] <- "sample"

vec_gr <- as.numeric(as.factor(data_l[,2]))
sample_gr <- unique(as.numeric(as.factor(data_l[,2])))
n_gr <- sapply(1:length(sample_gr), function(y) length(vec_gr[vec_gr == y]))
min_frac_man <- min(round(n_gr/length(vec_gr), 1)) # calculate min_frac manually

# download files
raw_data <- readMSData(files = files_all_pos, pdata = NULL, mode = "onDisk") # or use pd only as: pdata = new("NAnnotatedDataFrame", pd) or pdata = new("NAnnotatedDataFrame", n_gr_t)
raw_data$polarity <-1L
raw_data <- filterRt(raw_data, c(0,3600)) # time range in sec for unified rt range

# parallel processing
cores = detectCores()-1
register(bpstart(SnowParam(cores)))
BiocParallel::register(BiocParallel::SerialParam())
##### rename your results
resultPeakpicking <- resultPeakpicking_pos
resultRetcorGroup <- resultRetcorGroup_pos
# feature detection
cwp <- CentWaveParam(ppm = resultPeakpicking$best_settings$parameters$ppm, 
                     peakwidth = c(resultPeakpicking$best_settings$parameters$min_peakwidth, resultPeakpicking$best_settings$parameters$max_peakwidth),
                     snthresh = resultPeakpicking$best_settings$parameters$snthresh,
                     prefilter = c(resultPeakpicking$best_settings$parameters$prefilter, resultPeakpicking$best_settings$parameters$value_of_prefilter),
                     mzCenterFun = resultPeakpicking$best_settings$parameters$mzCenterFun,
                     integrate = resultPeakpicking$best_settings$parameters$integrate,
                     mzdiff = resultPeakpicking$best_settings$parameters$mzdiff,
                     fitgauss = resultPeakpicking$best_settings$parameters$fitgauss,
                     noise = resultPeakpicking$best_settings$parameters$noise)

feat_det <- findChromPeaks(raw_data, param = cwp)

# Merge peaks by rt, mz, proportion, see ?MergeNeighboringPeaksParam
data_merge <- refineChromPeaks(feat_det, MergeNeighboringPeaksParam(minProp = 0.05, expandRt = 4, expandMz = 0, ppm = 50)) # adjust to your data

save(data_merge, file = "xcms obj feat_det_hopsrecalc_pos_reduced_epxmz0.RData")

# retention time correction
BiocParallel::register(BiocParallel::SerialParam())
ret_cor <- adjustRtime(data_merge, param = ObiwarpParam(
  binSize = resultRetcorGroup$best_settings$profStep,
  center = resultRetcorGroup$best_settings$center,
  response = resultRetcorGroup$best_settings$response,
  distFun = resultRetcorGroup$best_settings$distFunc,
  gapInit = resultRetcorGroup$best_settings$gapInit,
  gapExtend = resultRetcorGroup$best_settings$gapExtend,
  factorDiag = resultRetcorGroup$best_settings$factorDiag,
  factorGap = resultRetcorGroup$best_settings$factorGap,
  localAlignment = ifelse(resultRetcorGroup$best_settings$localAlignment==0, F,T)))

save(ret_cor, file = "xcms obj ret_cor_hopsrecalc_pos_reduced_expmz0.RData")

# peak grouping
pgp <- PeakDensityParam(sampleGroups = n_gr_t$group_id, # rep(1, length(fileNames(feat_det))) or as.numeric(as.factor(n_gr_t))
                        bw = 0.88, #resultRetcorGroup$best_settings$bw
                        minFraction =  0.5, # or resultRetcorGroup$best_settings$minfrac
                        minSamples = resultRetcorGroup$best_settings$minsamp, 
                        binSize = resultRetcorGroup$best_settings$mzwid,
                        maxFeatures = resultRetcorGroup$best_settings$max) 

pk_gr <- groupChromPeaks(ret_cor, param = pgp)

save(pk_gr, file = "hops_recalc pos xcms obj pk_gr bw088 mf05.RData")
# peak filling
pk_fil <- fillChromPeaks(pk_gr)
save(pk_fil, file = "hops_recalc pos  xcms obj pk_fil bw088 mf05.RData")
# final feature table
ft_tbl <- featureValues(pk_fil, value = "into")

# final peak info table
ft_inf <- featureDefinitions(pk_fil)

# join peak table and save
ft_tbl_f <- data.frame(t(ft_tbl))
colnames(ft_tbl_f) <- paste(ft_inf$mzmed, ft_inf$rtmed, sep = " / ")
fwrite(ft_tbl_f, "xcms hops_recalc_cor_pos_bw10_mf_05_reduced_expmz0.csv", row.names = T)

# save all xcms objects
save(feat_det, file = "xcms obj feat_det_hopsrecalc_pos.RData")
save(ret_cor, file = "xcms obj ret_cor_hopsrecalc_pos.RData")
save(pk_gr, file = "xcms obj pk_gr.RData")
save(pk_fil, file = "xcms obj pk_fil.RData")


##############################################################################################################################################################
# Description of peaks table
##############################################################################################################################################################

# dataset
ds_pt <- ft_tbl_f

# Missing value %
tn <- nrow(ds_pt)*ncol(ds_pt)
mv_c <- sum(is.na(ds_pt)) # or length(which(is.a(...)))
pr_mv <- round(mv_c/tn*100)
pr_mv

# Number of peaks
ncol(ds_pt) # or nrow

