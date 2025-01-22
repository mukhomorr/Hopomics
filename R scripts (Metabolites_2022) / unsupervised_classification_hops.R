# setup environment
library(data.table)
setwd("D:/NAVUKA/Hops article/recacl_hops_dir/9_Final_Stats/unsupervised classification")

# PEAK TABLE
# dataset with intensities, annotation and label column
ds1 <- as.data.frame(fread(input = "neg genetic_n.csv", header=T))
rownames(ds1) <- ds1[,1]
#ds1 <- ds1[,-c(1,3:5)]
ds1 <- ds1[,-1]
ds1[,-1] <- sapply(ds1[,-1], as.numeric)
ds1$Label <- as.factor(ds1$Label)
ds2 <- as.data.frame(fread(input = "pos genetic_n.csv", header=T))

rownames(ds2) <- ds2[,1]
#ds2 <- ds2[,-c(1,3:5)]
ds2 <- ds2[,-1]
ds2[,-1] <- sapply(ds2[,-1], as.numeric)
ds2$Label <- as.factor(ds2$Label)
ds <- cbind(ds1,ds2[,-1])

# Settings for projection
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

# dataset
base1 <- ds # dataset
mtrx1 <- ds[,-1] # numeric data
grp1 <- as.character(base1[,1]) # label of dataset or try base1[,1]


###############################################
############################################### k-MEANS CLUSTERING
############################################### 
ds_norm_sc_c <- data.frame(scale(ds[,-1], center = T, scale = apply(ds[,-1], 2, sd, na.rm = T))) # for centering: center = T
ds_norm_sc_c <- cbind(ds[,1],ds_norm_sc_c)
dsz <- ds
ds <- ds_norm_sc_c
colnames(ds)[1] <- "Label"

base1 <- ds # dataset
mtrx1 <- ds[,-1] # numeric data
grp1 <- as.character(base1[,1])

k <- length(unique(grp1)) # groups in KM
km.res1 <- kmeans(mtrx1, centers = 2, nstart = 2,iter.max = 100,algorithm = "Lloyd")
fviz_cluster(list(data = mtrx1, cluster = km.res1$cluster), repel = T,
             ellipse.type = "manhattan", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic())

####check cluster
kmeans_lab <- km.res1$cluster
ds_new <- cbind(kmeans_lab,ds)
ds_new[,1:2]

write.csv(ds_new,"pos_gen_w_kmeans_lab_c_and_S_Lloyd.csv")

#############################
##### HCA on PCA matrix #####
#############################

base1 <- ds # dataset
mtrx1 <- ds[,-1] # numeric data
grp1 <- as.character(base1[,1])

mtrx_old <- mtrx1


res.pca <- PCA(mtrx_old, ncp = 2, graph = F)
pca_mat <- as.data.frame(res.pca$ind$coord)
pca_mat <- as.data.frame(cbind(as.character(ds$Label),pca_mat))

base1 <- pca_mat # dataset
mtrx1 <- pca_mat[,-1] # numeric data
grp1 <- as.character(base1[,1]) # label of dataset or try base1[,1]


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

mtrx1_1 <- mtrx1
#mtrx1_1 <- data.frame(scale(mtrx1, center = T, scale = T))
rownames(mtrx1_1) = make.names(grp1, unique=TRUE)
res.dist1 <- dist(mtrx1, method = "manhattan") #{euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
res.hc1 <- hclust(d = res.dist1, method = "ward.D2") #{ward (ward.D), (ward.D2)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}
hca <- fviz_dend(res.hc1, k = k, # Cut in k groups
                 type = "rectangle",
                 cex = 0.7, # label size 0.3/0.7
                 
                 k_colors = unique(Cols(grp1,res.hc1$order)), # "lancet" color "jco" gray JGRAY(k_hc)
                 color_labels_by_k = F, # color labels by groups
                 label_cols = Cols(grp1,res.hc1$order),#Cols(pca_mat[,1])[res.hc1$order], #as.fumeric(pca_mat[,1])[res.hc1$order]
                 rect = T, # Add rectangle around groups
                 rect_fill = T,
                 rect_border = unique(Cols(grp1,res.hc1$order)), #"lancet"# color "jco" gray JGRAY(k_hc)
                 horiz = F,
                 lwd = 0.3, # lines size 0.3/0.7
                 show_labels = T,
                 main = "HCA on PCA matrix results",
                 ylab = "")
hca


pca <- fviz_pca_ind(res.pca,
                    
                    title = "",
                    #geom.ind = "point", # show points only 
                    col.ind = grp1, # color by groups
                    palette = c("#008000", "#FFBF00"), # color "jco" gray JGRAY(length(unique(grp1)))
                    addEllipses = F, # Concentration ellipses
                    legend.title = "Groups",
                    label = "ind")
pca

####load summary table (with group information)

#extract different label vectors
test_ds <- fread("Êîïèÿ pos_gen_w_kmeans_lab_c_and_S_Lloyd.csv")
original_lab <- test_ds$`ds[, 1]`
hca_on_pca_lab <- as.character(test_ds$`PCA on HCA sum (1 -North America 2 - EU)`)
kmeans_lloyd_lab <- as.character(test_ds$`kmeans_lab lloyd`)
Geographic <- as.character(test_ds$`geographic origin`)

#fviz way with ellipses and custom color, but no ability to ad label dependent shape
pca_f <- fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (#text for sample names)
             col.ind =hca_on_pca_lab, repel = T,addEllipses=T,title = "PCA matrix based HCA classification",
             palette = c("#008000", "#FFBF00")) #use htnl color codes for customization

pca_f

#fviz + ggplot way
pca_t <- fviz_pca_ind(res.pca, label = "none") #simple pca plot with no aesthetics
                      

ggplot(cbind(pca_f$data,hca_on_pca_lab,geo_lab),
       aes(x=x,y=y, shape=geo_lab, fill=hca_on_pca_lab)) +#shape and color for different levels
  geom_point(size = 2.5) + theme_minimal() +
  scale_shape_manual(values = c(21, 22, 24)) +# chose shapes
  scale_fill_manual(values = c("#008000", "#FFBF00"))+ #chose points color
  scale_color_manual(values = c("#008000", "#FFBF00"))+ #chose elipse color
  xlim(-100,100)+ylim(-100,100)+ #set axes limits to display
  # set elipses to draw
  stat_ellipse(aes(group=hca_on_pca_lab,color=hca_on_pca_lab),type = "norm", size=0.1,alpha=0.1,geom = "polygon")+
  geom_vline(xintercept=0, linetype="dotted")+ #add axes
  geom_hline(yintercept=0, linetype="dotted")+
  labs(shape = "Geographic origin:",
       fill = "HCA based classification:",
       y = "PC2 (10.5 %)",
       x = "PC1 (25.7 %)")+### correct legend
  guides(size = "legend", colour = "none")
  

