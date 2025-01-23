# Okey Ive lost everything, but at least i have models saved
setwd("")

remotes::install_github("PaoloBnn/Retip")
library(Retip)
library(data.table)
library(rcdk)

## load your data
for_prediction <- read.csv("Compounds_article - For_prediction.csv",sep = ";")[,-1] # data for prediction
colnames(for_prediction) <- c("NAME","InchiKey","SMILES")  # rename columns in correct way
fp_desc <- getCD(for_prediction) # calculate descriptors for it
xgb <- readRDS("xgb_model_CHON_no_cutoff.rds") # load model if you don't want to compute it yourself
CDesk <- read.csv("desc.csv") # load precalculated descriptors for SMRT dataset
# Remove logP associated descriptors
library(stringi)
library(stringr)
CDesk <- CDesk[,-str_detect(colnames(CDesk),"LogP")]
# Keep only CHON containing molecules
chon_vec <- grepl(c("Br|Cl|S|F|I|P|B|Se|Si"),CDesk[,4])
Cdesk_CHON <- CDesk[which(chon_vec == F),]
Cdesk <- Cdesk_CHON

#> Clean dataset from NA and low variance value
db_rt <- CDesk
db_rt$RT <- db_rt$RT/60
preProc <- cesc(db_rt) #Build a model to use for center and scale a dataframe 
db_rt_cs <- predict(preProc,db_rt) # use the above created model for center and scale dataframe
# Split for train/test purposes
inTraining <- caret::createDataPartition(db_rt_cs$nB, p = .8, list = FALSE)
training <- db_rt_cs[ inTraining,]
testing  <- db_rt_cs[-inTraining,]
# Compute model
xgb <- fit.xgboost(training)
# Calculate times
pred_xgb <- RT.spell(testing,fp_desc, model = xgb, cesc = preProc)
