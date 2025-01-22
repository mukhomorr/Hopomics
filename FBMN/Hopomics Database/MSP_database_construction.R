##### .msp database construction (for Thermo .raw data) #####
setwd()

metadata <- read.csv("...") #load table with metadata for your compounds 


#Convert your raw data to .mgf to extract target scans
library(MSnbase)


mgf <-readMgfData("your_mgf.mgf")
#Substract scan number from scan name
library(stringr)
all_scans <- str_sub(mgf@featureData@data$TITLE,start = -7)

all_scans <- gsub("[^0-9.-]", "", all_scans)


inds <- which(as.numeric(all_scans) %in% scans)
subset_mgf <- mgf[inds]
subset_mgf@featureData@data$TITLE
selected_scans <- str_sub(subset_mgf@featureData@data$TITLE,start = -7)
selected_scans <- as.numeric(gsub("[^0-9.-]", "", selected_scans))

scans <- metadata$`scan number`
metadata <- metadata[match(selected_scans,scans),]
#save subseted .mgf
writeMgfData(subset_mgf)
# Load it and add metadata
library(Spectra)
library(MsBackendMsp)
library(MsBackendMgf)
spectra_mgf <-Spectra("experiment.mgf", source = MsBackendMgf())
spectra_mgf@backend@spectraData@listData$rtime <- spectra_mgf@backend@spectraData@listData$rtime/60
spectra_mgf@backend@spectraData@listData$NAME <- metadata$Compound
spectra_mgf@backend@spectraData@listData$SMILES <- metadata$smiles
spectra_mgf@backend@spectraData@listData$FORMULA <- metadata$Formula

tmpf <- tempfile()

export(spectra_mgf, backend = MsBackendMsp(), file = paste0(getwd(),"/your_database.msp"),
       mapping = spectraVariableMapping(MsBackendMsp()))
head(readLines(paste0(getwd(),"/your_database.msp")))
