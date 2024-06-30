#knn,spark,BIC
rm(list=ls())
library(DR.SC)
library(SingleCellExperiment)
library(mvtnorm)
library(GiRaF)
library(dplyr)
library(SC.MEB)
library(orclus)
library(Seurat)
library(igraph)

source("source-easy.R")

simulDRcluster <- DR.SC:::simulDRcluster
selectClustNumber <- DR.SC:::selectClustNumber
find_neighbors <- DR.SC:::find_neighbors

name_ID <- as.character(c(151507, 151508, 151509, 151510, 151669, 151670,151671, 151672, 151673, 151674, 151675, 151676))
trueK_set <- c(rep(7,4), rep(5,4), rep(7,4)) #contains K for each sample
n_ID <- length(name_ID)
num_cut <- 2000
p=num_cut
q=15
final=matrix(NA,12,14)
KMat=matrix(NA,12,7)

address="DLPFC/"
data=list()
d2=list()
label=list()
position=list()
n=rep(NA,12)
for (i in 1:12){
  dlpfc <- readRDS(paste(address,paste0(name_ID[i], ".rds"), sep = ""))
  ## use  top 2000 genes from SPARK
  load(paste(address,paste0("brain_", name_ID[i],"_spark.Rdata"),sep=""))
  adjPval <- PvalDF[,2]
  names(adjPval) <- row.names(PvalDF)
  sort_adjPval <- sort(adjPval)
  if(sum(sort_adjPval<0.05)<= num_cut){
    sp_sig_genes <- row.names(PvalDF)[PvalDF[,2] < 0.05]
  }else{
    sp_sig_genes <- names(sort_adjPval)[1:num_cut]
  }
  logCount <- assay(dlpfc, "logcounts") #p*n
  sp_logCount <- logCount[sp_sig_genes, ]
  
  X <- as.matrix(t(sp_logCount)) # obtain data, n*p
  pos=cbind(dlpfc$row, dlpfc$col) # spatial coordinates
  y <- dlpfc$layer_guess_reordered # true label
  na_index=which(is.na(y))
  y=y[-na_index]
  
  d2[[i]] <- as.matrix(t(logCount))[-na_index,]#original p
  data[[i]]=X[-na_index,] #reduced p
  position[[i]]=pos[-na_index,]
  label[[i]]=y
  n[i]=nrow(X[-na_index,])
}

method="BIC"
neighbor="knn"
zeta=0.01
seed=2023
usage="spark"
est_label=list()


for (i in 1:12){
  print(i)
  Xst=as.matrix(data[[i]]) #used the reduced p
  source=do.call(rbind, d2[-i]) #use the original p first
  Xsc=source[,colnames(Xst)]
  pos=as.matrix(position[[i]])
  z_st=as.numeric(factor(label[[i]], levels = unique(label[[i]])))
  z_temp=do.call(c, label[-i])
  z_sc=as.numeric(factor(z_temp, levels = unique(z_temp)))
  K=trueK_set[i]

  ##For TransST
  all=scale(rbind(Xst,Xsc))
  n_st=nrow(Xst)
  n_all=nrow(all)
  Xst_new=all[1:n_st,]
  Xsc_new=all[(n_st+1):n_all,]

  raw=cbind(Xst_new,pos) #
  rownames(raw) <- paste0("spot", seq_len(nrow(Xst_new)))
  hZ=approxPCA(Xst, q=q)$PCs

  ##Creat Seurat object
  counts=t(Xst)
  pp=ncol(Xst)
  rownames(counts) <- paste0("gene", seq_len(pp))
  colnames(counts) <- paste0("spot", seq_len(nrow(Xst)))
  ## Make array coordinates - filled rectangle
  cdata <- list()
  cdata$row <- pos[,1]
  cdata$col <- pos[,2]
  cdata <- as.data.frame(do.call(cbind, cdata))
  cdata$imagerow <- cdata$row
  cdata$imagecol <- cdata$col
  rownames(cdata) <- colnames(counts)
  seu <-  CreateSeuratObject(counts= counts, meta.data=cdata) #
  seu <- NormalizeData(seu,verbose=FALSE, normalization.method = "LogNormalize")
  seu <- FindVariableFeatures(seu, nfeatures = p,verbose=FALSE)

  ### DR.SC -----------------------------------------------------------------
  set.seed(seed+i)
  res_drsc <- DR.SC(seu, q=q,K=K, platform = "Visium",  verbose=FALSE)
  res_drsc_K <- DR.SC(seu, q=q,K=2:10, platform = "Visium",  verbose=FALSE)
  KMat[i, 2] =length(unique(res_drsc_K$spatial.drsc.cluster))
  #
  ### Seurat -----------------------------------------------------------------
  brain <- SCTransform(seu, verbose = FALSE)
  brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
  brain <- FindNeighbors(brain, reduction = "pca", dims = 1:q)
  brain <- FindClusters(brain, verbose = FALSE)
  KMat[i, 3] =length(unique(brain$seurat_clusters))

  # k-means -----------------------------------------------------------------
  y_kmeans = kmeans(hZ, centers=K)$cluster

  ### GMM
  res_gmm0=Mclust(hZ, G=K, modelNames='EEE',verbose=FALSE)
  res_gmm0_K=Mclust(hZ, G=NULL, modelNames='EEE',verbose=FALSE)
  KMat[i, 4] =res_gmm0_K$G

  # ### SC-MEB
  sce <- SingleCellExperiment(assays=list(logcounts=counts), colData=cdata,reducedDims=SimpleList(PCA=hZ))
  Adj_sp <- find_neighbors2(sce, platform = "Visium")
  set.seed(seed+i)
  scmeb_o = SC.MEB(hZ, Adj_sp, beta_grid = seq(0.5, 5, by=0.5), K_set= K, parallel=TRUE, num_core = 3, PX = TRUE)
  out1 = SC.MEB::selectK(scmeb_o, K_set = K,criterion = "BIC")
  set.seed(seed+i)
  scmeb_K = SC.MEB(hZ, Adj_sp, beta_grid = seq(0.5, 5, by=0.5), K_set= 2:10, parallel=TRUE, num_core = 3, PX = TRUE)
  out2 = SC.MEB::selectK(scmeb_K, K_set = 2:10,criterion = "BIC")
  KMat[i, 5] =out2$best_K_BIC

  # Louvain
  g.jaccard = scran::buildSNNGraph(sce, use.dimred="PCA", type="jaccard")
  y_louvain <- igraph::cluster_louvain(g.jaccard)$membership
  KMat[i, 6] =length(unique(y_louvain))

  ## Leiden
  ldc <- cluster_leiden(g.jaccard, resolution_parameter = quantile(strength(g.jaccard))[2] / (gorder(g.jaccard) - 1))
  KMat[i, 7] =length(unique(ldc$membership))

  ####select K
  res=trans_k_bic(q=q,K_set=2:10,Ux=cbind(hZ,pos),zeta=zeta,method="HBIC",verbose=FALSE)
  K_ours=res$K

  KMat[i, 1] =K_ours

  ######
  res6=spGMM(K=K,U=cbind(hZ,pos),nS=n_st,maxIter=100,zeta=zeta,neighbor=neighbor,verbose=FALSE)
  set.seed(seed+i)
  res2=TransST3(Xsc=Xsc_new,z_sc=z_sc,Xst=raw,q=q,K=K,nS=n_st,maxIter=100,zeta=zeta,combine=FALSE,neighbor=neighbor,verbose=FALSE)

  if(K_ours==K){
    res2_K=res2
    res6_K=res6
  }else{
    res6_K=spGMM(K=K_ours,U=cbind(hZ,pos),nS=n_st,maxIter=100,zeta=zeta,neighbor=neighbor,verbose=FALSE)
    set.seed(seed+i)
    res2_K=TransST3(Xsc=Xsc_new,z_sc=z_sc,Xst=raw,q=q,K=K_ours,nS=n_st,maxIter=100,zeta=zeta,combine=FALSE,neighbor=neighbor,verbose=FALSE)
  }

  #########
  KMat[i,]=abs(KMat[i,]-K)

  final[i,]=c(adjustedRandIndex(z_st, y_kmeans),adjustedRandIndex(res_gmm0$classification, z_st),
              adjustedRandIndex(out1$best_K_label, z_st),adjustedRandIndex(brain$seurat_clusters, z_st),
              adjustedRandIndex(res2$z, z_st),adjustedRandIndex(res6$z, z_st),adjustedRandIndex(res_drsc$spatial.drsc.cluster, z_st),
              #for determine K
              adjustedRandIndex(res_gmm0_K$classification, z_st),adjustedRandIndex(z_st, y_louvain),
              adjustedRandIndex(z_st, ldc$membership),adjustedRandIndex(res2_K$z, z_st),adjustedRandIndex(res6_K$z, z_st),
              adjustedRandIndex(res_drsc_K$spatial.drsc.cluster, z_st),adjustedRandIndex(out2$best_K_label, z_st))

  if (i>1){
    print(colMeans(final[1:i,],na.rm=TRUE))
    print(colMeans(KMat[1:i,],na.rm=TRUE))
  }
  est_label[[i]]=cbind(y_kmeans,res_gmm0$classification,out1$best_K_label,brain$seurat_clusters,
                       res2$z,res6$z,res_drsc$spatial.drsc.cluster,res_gmm0_K$classification,y_louvain,
                       ldc$membership,res2_K$z,res6_K$z,res_drsc_K$spatial.drsc.cluster,out2$best_K_label,pos,z_st)
}


