# setup environment
library(data.table)
setwd("C:/Users/Plyush/Yura/NAVUKA/Hop and beer metabolomics/recacl_hops_dir/9_Final_Stats/unsupervised classification")

# PEAK TABLE
# dataset with intensities, annotation and label column
ds1 <- as.data.frame(fread(input = "neg genetic_n.csv", header=T))
rownames(ds1) <- ds1[,1]
#ds1 <- ds1[,-c(1,3:5)]
ds1 <- ds1[,-1]
colnames(ds1)[-1] <- paste0("n",colnames(ds1)[-1])
ds1[,-1] <- sapply(ds1[,-1], as.numeric)
ds1$Label <- as.factor(ds1$Label)
ds2 <- as.data.frame(fread(input = "pos genetic_n.csv", header=T))

rownames(ds2) <- ds2[,1]
#ds2 <- ds2[,-c(1,3:5)]
ds2 <- ds2[,-1]
colnames(ds2)[-1] <- paste0("p",colnames(ds2)[-1])
ds2[,-1] <- sapply(ds2[,-1], as.numeric)
ds2$Label <- as.factor(ds2$Label)
ds <- cbind(ds1,ds2[,-1])


####new labels load

test_ds <- fread("Êîïèÿ pos_gen_w_kmeans_lab_c_and_S_Lloyd.csv")
ds$Label <- test_ds$`PCA on HCA sum (1 -NA 2 - EU)`
#### filter signals by time
library(stringr)
signals <- colnames(ds)
split_signals <- as.data.frame(t(as.data.frame(str_split(signals[-1],"\\.\\.\\."))))
rownames(split_signals) <- str_remove_all(rownames(split_signals),"c\\.\\.")
split_signals$V2 <- as.numeric(split_signals$V2)
retained_molecules <- rownames(split_signals[which(split_signals$V2>180),])
retained_id <- which(split_signals$V2>180)
labs <- as.character(ds$Label)                             
#ds <- ds[-1,]
ds <- ds[,retained_id+1]
ds <- cbind(as.factor(labs),ds)
colnames(ds)[1] <- "Label"




####centerin n scaling

ds_norm_sc_c <- data.frame(scale(ds[,-1], center = T, scale = apply(ds[,-1], 2, sd, na.rm = T))) # for centering: center = T
ds_norm_sc_c <- cbind(ds[,1],ds_norm_sc_c)
dsz <- ds ####not scaled data
ds <- ds_norm_sc_c ####scaled data
colnames(ds)[1] <- "Label"
###### get only qc intensities for your data to filter by it
df1 <- fread("hops_recalc MVI rf2_neg.csv")
rownames(df1) <- df1$V1
#df1 <- df1[,-1]
colnames(df1) <- paste0("n",colnames(df1))

df2 <- fread("hops_recalc MVI rf2_pos.csv")
rownames(df2) <- df2$V1
#df2 <- df2[,-1]
colnames(df2) <- paste0("p",colnames(df2))

df <- cbind(df1[,-1],df2[,-1])
rownames(df) <- df1$nV1

dq <- df[which(grepl("qc",rownames(df))),]

mean_df <- as.data.frame(lapply(dq,mean))

#### some information on your signals
min(mean_df)
max(mean_df)
mean(as.numeric(mean_df[1,]))
length(which(as.numeric(mean_df[1,])>100000,T))

df_id <- as.numeric(df[,which(as.numeric(mean_df[1,])>100000)])
high <- colnames(mean_df[,df_id])

####filter by int
gg <- intersect(high, str_remove(colnames(ds),"X"))
gg <- str_replace(gg,"n","nX")
gg <- str_replace(gg,"p","pX")
gg <- c("Label",gg)
dsf <- ds[,gg]

ds <- dsf


  ############
#####stats######
  ############

############################################### 
############################################### MODERATED T-TEST
###############################################
####### T-test part (moderated) and fold changw calculation

library(limma)
library(dplyr)

mdl_mtrx <- model.matrix(~Label, ds) # adjust to your data
lmf <- lmFit(t(ds[,-c(1,2)]), method = "ls", design = mdl_mtrx, maxit = 1000) # "robust" or "ls"
efit <- eBayes(lmf)
####your top table
tableTop <- topTable(efit, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method

cutoff <- 0.05 # set cutoff
lim_pval <- rownames(filter(tableTop, as.numeric(adj.P.Val) <= cutoff)) # features

####### Fold change part (actually we are doing limma based moderated t-test, but with logg dataframe)

# transform data into log2 base
ds_log <- as.data.frame(log2(ds[,-c(1,2)]))
ds_log <- cbind(ds[,1], ds_log)

colnames(ds_log)[1] <- "Label"
mdl_mtrx_log <- model.matrix(~Label, ds_log) # adjust to your data
lmf_log <- lmFit(t(ds_log[,-1]), method = "ls", design = mdl_mtrx_log, maxit = 1000) # "robust" or "ls"
efit_log <- eBayes(lmf_log)
####your top table
tableTop_log <- topTable(efit_log, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method

FC_cutoff <- 1
ii_fc <- rownames(filter(tableTop_log, abs(as.numeric(logFC)) >= FC_cutoff)) # features
####check FC and t test intersect
combine_simple <- Reduce(intersect, list(ii_fc, lim_pval))
length(combine_simple)
###############################################
############################################### VIP FROM PLS MODELS
############################################### 

library(ropls)

pls <- opls(ds[,-1], ds[,1], orthoI = 0, predI = NA, crossvalI = 10, permI = 100) # orthoI -> 0/1 PLS/OPLS, predI -> No of components, crossvalI -> cross-validation, permI -> No of permutations
# plot(pls)
vip <- as.data.frame(getVipVn(pls))
vip <- cbind(name = rownames(vip), VIP = vip) 
colnames(vip)[2] <- "VIP"
vip <- vip[order(vip$VIP, decreasing = T),]
th_vip <- 1 # set value for filtration
vip_th <- subset(vip, vip$VIP > th_vip)
ii_vip_pls<- rownames(vip_th) # features

combine <- Reduce(intersect, list(ii_vip_pls,ii_fc, lim_pval))

########## Get table with metrics #########

t_test_table <- tableTop
t_test_table <- t_test_table[order(row.names(t_test_table)),]

FC_table <- tableTop_log
FC_table <- FC_table[order(row.names(FC_table)),]

vip_table <- vip
vip_table <- vip_table[order(row.names(vip_table)),]

combined_table <- cbind(vip_table[-length(vip_table$name),],t_test_table$adj.P.Val,FC_table$logFC)


############################################################
combined_table_select <- combined_table[combine,]
combined_table_select <- combined_table_select[-which(is.na(combined_table_select$`FC_table$logFC`)),]
combined_table_select <- combined_table_select[order(combined_table_select$VIP,decreasing = T),]
top50 <- combined_table_select[1:50,]

qua <- grepl("n",top50$name)
negative <- length(which(qua,qua == T))

negative
50-negative


#########Nested cross validation#########

############## my 1st Nested CV Feature Selection

#function for stat tests

library(limma)
library(dplyr)
library(ropls)

perform_tests <- function(inner_indecies,dataframe){
  ds <- dataframe[inner_indecies,]
  mdl_mtrx <- model.matrix(~Label, ds) # adjust to your data
  lmf <- lmFit(t(ds[,-1]), method = "ls", design = mdl_mtrx, maxit = 1000) # "robust" or "ls"
  efit <- eBayes(lmf)
  ####your top table
  tableTop <- topTable(efit, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method
  
  cutoff <- 0.05 # set cutoff
  lim_pval <- rownames(filter(tableTop, as.numeric(adj.P.Val) <= cutoff)) # features
  ####### Fold change part (actually we are doing limma based moderated t-test, but with logg dataframe)
  
  # transform data into log2 base
  ds_log <- as.data.frame(log2(ds[,-c(1,2)]))
  ds_log <- cbind(ds[,1], ds_log)
  
  colnames(ds_log)[1] <- "Label"
  mdl_mtrx_log <- model.matrix(~Label, ds_log) # adjust to your data
  lmf_log <- lmFit(t(ds_log[,-1]), method = "ls", design = mdl_mtrx_log, maxit = 1000) # "robust" or "ls"
  efit_log <- eBayes(lmf_log)
  ####your top table
  tableTop_log <- topTable(efit_log, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method
  
  FC_cutoff <- 1
  ii_fc <- rownames(filter(tableTop_log, abs(as.numeric(logFC)) >= FC_cutoff)) # features
  #combine_simple <- Reduce(intersect, list(ii_fc, lim_pval))
  ####VIP calc part
  pls <- opls(ds[,-1], ds[,1], orthoI = 0, predI = NA, crossvalI = 10, permI = 100,printL = F,plotL = F,
             fig.pdfC = "none",
              info.txtC = "none") # orthoI -> 0/1 PLS/OPLS, predI -> No of components, crossvalI -> cross-validation, permI -> No of permutations
   plot(pls)
  vip <- as.data.frame(getVipVn(pls))
  vip <- cbind(name = rownames(vip), VIP = vip) 
  colnames(vip)[2] <- "VIP"
  vip <- vip[order(vip$VIP, decreasing = T),]
  th_vip <- 1 # set value for filtration
  vip_th <- subset(vip, vip$VIP > th_vip)
  ii_vip_pls<- rownames(vip_th) # features
  
  combine <- Reduce(intersect, list(ii_vip_pls,ii_fc, lim_pval))
return(combine)
}

# create outer folds (this actually creates test groups)
outer <- createFolds(ds$Label, 10)
l_r <- list()

for (i in 1:length(outer)) { #for each outer fold do 
  inner <- createResample(ds[-outer[[i]],]$Label, 10) # do bootstrap for each outer train (outer train = inner)
    is <- lapply(inner, perform_tests,ds[-outer[[i]],]) ### do statistic test part for each inner
    ifs <- Reduce(intersect, is) # table(unlist(is)) or Reduce(intersect, is)
    l_r[[i]] <- ifs
}

ifs_sum <- Reduce(intersect, l_r)

####function to work with rsample bootstrap
library(rsample)
library(caret)
get_data <- function(x){
  return(x[["in_id"]])
}

set.seed(20140102)
outer <- createFolds(ds$Label, 10)
inner_func <- function(outr){
  inner_ <- bootstraps(ds[-outr,],5,strata = Label) # do bootstrap for each outer train (outer train = inner)
  inner_ <- inner_[[1]]
  inner <- lapply(inner_,get_data)
  is <- lapply(inner, perform_tests,ds[-outr,]) ### do statistic test part for each inner
  ifs <- Reduce(intersect, is) # table(unlist(is)) or Reduce(intersect, is)
  return(ifs)
  
}

all <- lapply(outer,inner_func)
test_sum <- Reduce(intersect, all)
test_sum

qua <- grepl("n",test_sum)
negative <- length(which(qua,qua == T))

negative
length(test_sum)-negative

test_table <- str_split(test_sum,"\\.\\.\\.")
test_table <- lapply(test_table,str_remove,"X")


fwrite(test_table,"res_160322_pls_easier_cutoff_100000.csv")

#### cross validation

outer

#### table for nest selected data with original (whole dataset tests) metrics

combined_table_nest <- combined_table[test_sum,]
combined_table_nest <- combined_table_nest[order(combined_table_nest$VIP,decreasing = T),]

fwrite(combined_table_nest,"res_met_160322_pls_easier_cutoff_100000.csv")

`#### PCA and HCA with nest selected vars

library(factoextra)
library(FactoMineR)
library(dendextend)
library(rafalib)
library(RSEIS)
library(ggsci)
library(pheatmap)
library(Rtsne)
library(NbClust)
library(clustertend)
library(mclust)
library(clValid)
library(fpc)
library(pvclust)
library(parallel)
library(doParallel)

dx <- ds[,c("Label",test_sum)]


base1 <- dx # dataset
mtrx1 <- dx[,-1] # numeric data
grp1 <- as.character(base1[,1])

res.pca <- PCA(mtrx1, ncp = 10, graph = F)

pca <- fviz_pca_ind(res.pca,
                    
                    title = "",
                    geom.ind = "point", # show points only 
                    col.ind = grp1, # color by groups
                    #palette = pallete_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                    addEllipses = T, # Concentration ellipses
                    legend.title = "Groups",
                    label = "ind")
pca


# number of groups
k <- length(unique(grp1)) # groups in HC

# color
Cols = function(vec, ord){
  cols = pal_lancet(palette = c("lanonc"), alpha = 1)(length(unique(vec)))
  return(cols[as.fumeric(vec)[ord]])}

# grey
#Cols = function(vec, ord){
# cols = JGRAY(length(unique(vec)))
# return(cols[as.fumeric(vec)[ord]])}

mtrx <- mtrx1
#mtrx <- data.frame(scale(mtrx1, center = T, scale = T))
rownames(mtrx) = row.names(ds)
res.dist1 <- dist(mtrx1, method = "manhattan")
res.hc1 <- hclust(d = res.dist1, method = "ward.D2")
dend1 <- as.dendrogram(res.hc1, hang =70)
dend1 <- color_branches(dend1, k=k, groupLabels = F, col =unique(Cols(grp1,res.hc1$order))) 
labels_colors(dend1) <-Cols(grp1,res.hc1$order)
dend1 <- assign_values_to_leaves_nodePar(dend1, 0.6, "lab.cex")
dend1 <- dendextend::set(dend1, "branches_lwd", 2.5)
plot(dend1)
#dend1 <- rect.dendrogram(dend1, k=k, border = 1, lty = 1, lwd = 1, col=rgb(0.1, 0.2, 0.4, 0.1))
legend("topright", legend = unique(grp1), fill = pal_lancet(palette = c("lanonc"), alpha = 1)(length(unique(grp1))))


###############################################
############################################### PLOT RESULTS
###############################################

dx <- ds[,c("Label",test_sum)]
library(stringr)
library(reshape2)
library(ggplot2)
colnames(dx) <- str_remove(colnames(dx),"pX")
colnames(dx) <- str_remove(colnames(dx),"nX")
labs <- colnames(dx)
labs <- str_split(labs,"\\.\\.\\.")
labbs <- do.call(rbind,labs)
labbs <- as.data.frame(labbs[-1,])
labbs <- as.data.frame(sapply(labbs,as.numeric))
mz <- round(labbs$V1,4)
rt <- round(labbs$V2,2)
cols <- paste0(mz,"@",rt)
colnames(dx)[-1] <- cols

df.m <- melt(dx, id.var = "Label") # reshape data frame

p <- ggplot(data = df.m, aes(x=variable, y=value)) + xlab("") + ylab("") +
  geom_boxplot(aes(fill=Label)) + theme(legend.position="bottom") + theme_classic() + scale_x_discrete(labels=c("")) 

pp <- p + facet_wrap( ~ variable, scales="free") + theme_classic   () + theme(legend.position="bottom") 

pp
pp$data
