setwd("")

library(data.table)

ds <- as.data.frame(fread(input = "", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
class1 <- which(ds$Label == 1)
class2 <- which(ds$Label != 1)
ds$Label[class1] <- "Class_1"
ds$Label[class2] <- "Class_2"
ds$Label <- as.factor(ds$Label)


##### Simple feature selection #####

library(limma)
# Moderated t-test
mdl_mtrx <- model.matrix(~Label, ds) # adjust to your data
lmf <- lmFit(t(ds[,-1]), method = "ls", design = mdl_mtrx, maxit = 1000) # "robust" or "ls"
efit <- eBayes(lmf)
tableTop <- topTable(efit, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method
lim_pval <- rownames(dplyr::filter(tableTop, as.numeric(adj.P.Val) <= 0.05))


# Fold Change
ds_log <- as.data.frame(log2(ds[,-1]))
ds_log <- cbind(ds[,1], ds_log)
colnames(ds_log)[1] <- "Label"
mdl_mtrx_log <- model.matrix(~Label, ds_log) # adjust to your data
lmf_log <- lmFit(t(ds_log[,-1]), method = "ls", design = mdl_mtrx_log, maxit = 1000) # "robust" or "ls"
efit_log <- eBayes(lmf_log)
tableTop_log <- topTable(efit_log, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method
fc <- rownames(dplyr::filter(tableTop_log, abs(as.numeric(logFC)) >= 1))

library(ropls)
# VIP value from PLS
pls <- opls(ds[,-1], ds[,1], orthoI = 0, predI = NA, crossvalI = 10, permI = 100,printL = F,plotL = F, fig.pdfC = "none", info.txtC = "none") # orthoI -> 0/1 PLS/OPLS, predI -> No of components, crossvalI -> cross-validation, permI -> No of permutations
vip <- as.data.frame(getVipVn(pls))
vip <- cbind(name = rownames(vip), VIP = vip)
colnames(vip)[2] <- "VIP"
vip <- vip[order(vip$VIP, decreasing = T),]
vip_th <- subset(vip, vip$VIP > 1)
vip_pls<- rownames(vip_th)

combine <- Reduce(intersect, list(vip_pls, fc, lim_pval))

########## Get table with metrics #########

t_test_table <- tableTop
t_test_table <- t_test_table[order(row.names(t_test_table)),]

FC_table <- tableTop_log
FC_table <- FC_table[order(row.names(FC_table)),]

vip_table <- vip
vip_table <- vip_table[order(row.names(vip_table)),]

combined_table <- cbind(vip_table,t_test_table$adj.P.Val,FC_table$logFC)

write.csv(combined_table,"Metrics for all features.csv")

combined_table_FS <- combined_table[combine,]

write.csv(combined_table_FS,"Metrics for simple selected features.csv")

##Let's look at the HCA and PCA

#### PCA and HCA with nest selected vars

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

dx <- ds[,c("Label",combine)]


base1 <- dx # dataset
mtrx1 <- dx[,-1] # numeric data
grp1 <- as.character(base1[,1])

res.pca <- PCA(mtrx1, ncp = 10, graph = F)

pca <- fviz_pca_ind(res.pca,
                    
                    title = "Simple Feature Selection (by FC, adj. p-value, VIP from PLS cut offs)",
                    geom.ind = "text", # show points only 
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


mtrx <- mtrx1

rownames(mtrx) = row.names(ds)
res.dist1 <- dist(mtrx1, method = "manhattan")
res.hc1 <- hclust(d = res.dist1, method = "ward.D2")
dend1 <- as.dendrogram(res.hc1, hang =70)
dend1 <- color_branches(dend1, k=k, groupLabels = F, col =unique(Cols(grp1,res.hc1$order))) 
labels_colors(dend1) <-Cols(grp1,res.hc1$order)
dend1 <- assign_values_to_leaves_nodePar(dend1, 0.6, "lab.cex")
dend1 <- dendextend::set(dend1, "branches_lwd", 2.5)
plot(dend1)

legend("topright", legend = unique(grp1), fill = pal_lancet(palette = c("lanonc"), alpha = 1)(length(unique(grp1))))


##### Selection with 10-fold resampling and validation on the keeped outer folds (no bootstrap) ##### 

# wrapper function with feature selection methods (in this example: intersection of 2 algorithms (p-value in moderated t-test, Fold-Change, VIP in PLS))
perform_tests <- function(inner_indecies,dataframe, cutoff_tt = 0.05, cutoff_fc = 1, cutoff_vip = 1){
  ds <- dataframe[inner_indecies,]
  
  # Moderated t-test
  mdl_mtrx <- model.matrix(~Label, ds) # adjust to your data
  lmf <- lmFit(t(ds[,-1]), method = "ls", design = mdl_mtrx, maxit = 1000) # "robust" or "ls"
  efit <- eBayes(lmf)
  tableTop <- topTable(efit, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method
  lim_pval <- rownames(dplyr::filter(tableTop, as.numeric(adj.P.Val) <= cutoff_tt))
  
  # Fold Change
  ds_log <- as.data.frame(log2(ds[,-1]))
  ds_log <- cbind(ds[,1], ds_log)
  colnames(ds_log)[1] <- "Label"
  mdl_mtrx_log <- model.matrix(~Label, ds_log) # adjust to your data
  lmf_log <- lmFit(t(ds_log[,-1]), method = "ls", design = mdl_mtrx_log, maxit = 1000) # "robust" or "ls"
  efit_log <- eBayes(lmf_log)
  tableTop_log <- topTable(efit_log, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method
  fc <- rownames(dplyr::filter(tableTop_log, abs(as.numeric(logFC)) >= cutoff_fc))
  
  # VIP value from PLS
  pls <- opls(ds[,-1], ds[,1], orthoI = 0, predI = NA, crossvalI = 10, permI = 100,printL = F,plotL = F, fig.pdfC = "none", info.txtC = "none") # orthoI -> 0/1 PLS/OPLS, predI -> No of components, crossvalI -> cross-validation, permI -> No of permutations
  vip <- as.data.frame(getVipVn(pls))
  vip <- cbind(name = rownames(vip), VIP = vip)
  colnames(vip)[2] <- "VIP"
  vip <- vip[order(vip$VIP, decreasing = T),]
  vip_th <- subset(vip, vip$VIP > cutoff_vip)
  vip_pls<- rownames(vip_th)
  
  combine <- Reduce(intersect, list(vip_pls, fc, lim_pval))
  return(combine)
}
library(rsample)
# Set parameters
cutoff_tt <- 0.05 # set p-value for moderated t-test
cutoff_fc <- 1 # set cutoff for Fold Change
cutoff_vip <- 1 # set cutoff for VIP in PLS
outer <- createFolds(ds$Label, 10) # Type of resampling in outer loop, adjust to your data
tc_boot <- trainControl(method = "cv", number = 10, savePredictions = T, classProbs = T) # Type of resampling in outer loop for machine learning, adjust to your data
l_r <- list()
results <- data.frame(fold = numeric(), par = numeric(), acc_test = numeric()) 

# Perform
for (i in 1:length(outer)) { 
  inner <- c(1:54)[-outer[[i]]] # Type of resampling in inner loop, adjust to your data
  is <- perform_tests(inner,ds)
  l_r[[i]] <- is
  
  
}

nfs <- Reduce(intersect, l_r) # Selected features (intersect from all inner and outer loops)

combined_table_resample_selection <- combined_table[colnames(ds_loop),]

write.csv(combined_table_resample_selection,"Metrics for resample selected features.csv")

# Lets look at the PCA and HCA

dx <- ds[,c("Label",nfs)]


base1 <- dx # dataset
mtrx1 <- dx[,-1] # numeric data
grp1 <- as.character(base1[,1])

res.pca <- PCA(mtrx1, ncp = 10, graph = F)

pca <- fviz_pca_ind(res.pca,
                    title = "Feature Selection with 10-fold based resampling (by FC, adj. p-value, VIP from PLS cut offs)",
                    geom.ind = "text", # show points only 
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

### Scheme we used in our work



set.seed(1234)
# Perform
for (i in 1:length(outer)) { 
  inner <- createResample(ds[-outer[[i]],]$Label, 5) # Type of resampling in inner loop, adjust to your data
  is <- lapply(1:length(inner), function(y) perform_tests(inner[[y]], ds[-outer[[i]],], cutoff_tt = cutoff_tt, cutoff_fc = cutoff_fc, cutoff_vip = cutoff_vip)) # adjust cutoffs to your data
  ifs <- Reduce(intersect, is) # table(unlist(is)) or Reduce(intersect, is)
  l_r[[i]] <- ifs
  
  # PLS (as example) see http://topepo.github.io/caret/train-models-by-tag.html for other regression algorithm
  #m_inner <- train(x = ds[-outer[[i]], l_r[[i]]], y = ds[-outer[[i]],]$Label, method = "glm", trControl = tc_boot) 
  #acc <- confusionMatrix(predict(m_inner, ds[outer[[i]], l_r[[i]]]), ds[outer[[i]],]$Label) 
  #results <- dplyr::add_row(results, fold = i, par = as.numeric(m_inner$bestTune), acc_test = acc$overall["Accuracy"]) # adjust ti your data
}

nfs_b_recalc_1234 <- Reduce(intersect, l_r) # Selected features (intersect from all inner and outer loops)

article_nfs <- read.csv("",header = F) ### features from article
article_nfs <- article_nfs$V1
article_nfs[1] <- str_remove_all(article_nfs[1],"﻿")

dx <- ds[,c("Label",article_nfs)]


base1 <- dx # dataset
mtrx1 <- dx[,-1] # numeric data
grp1 <- as.character(base1[,1])

res.pca <- PCA(mtrx1, ncp = 10, graph = F)

pca <- fviz_pca_ind(res.pca,
                    
                    title = "Feature selection scheme used in the article",
                    geom.ind = "text", # show points only 
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



####ROPLS TEST
ds_simp <- ds[,combine]
ds_loop <- ds[,nfs]
ds_boot <- ds[,nfs_b]
ds_boot_5 <- ds[,nfs_b_5]
ds_article <- ds[,article_nfs]

pls_simple_fs <- opls(ds_simp, ds_simp[,1], orthoI = 0, predI = NA, crossvalI = 10, permI = 100,printL = F,plotL = F, fig.pdfC = "none", info.txtC = "none") # orthoI -> 0/1 PLS/OPLS, predI -> No of components, crossvalI -> cross-validation, permI -> No of permutations
 pls_loop_fs <- opls(ds_loop[,-1], ds_loop[,1], orthoI = 0, predI = NA, crossvalI = 10, permI = 100,printL = F,plotL = F, fig.pdfC = "none", info.txtC = "none") # orthoI -> 0/1 PLS/OPLS, predI -> No of components, crossvalI -> cross-validation, permI -> No of permutations

pls_article <- opls(ds_article[,-1], ds_article[,1], orthoI = 0, predI = NA, crossvalI = 10, permI = 100,printL = F,plotL = F, fig.pdfC = "none", info.txtC = "none") # orthoI -> 0/1 PLS/OPLS, predI -> No of components, crossvalI -> cross-validation, permI -> No of permutations
