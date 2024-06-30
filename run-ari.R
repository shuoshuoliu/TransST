rm(list=ls())
library(DR.SC)
#library(SingleCellExperiment)
library(mvtnorm)
library(GiRaF)
library(dplyr)
library(SC.MEB)
library(orclus)
library(Seurat)
library(igraph)

source("help/source-easy.R")

simulDRcluster <- DR.SC:::simulDRcluster
selectClustNumber <- DR.SC:::selectClustNumber
find_neighbors <- DR.SC:::find_neighbors

n=50
p=1000
q=10
beta=1
K=4

iter=50
final=matrix(NA,iter,14)
KMat=matrix(NA,iter,7)
error=2 #DRSC needs this to be small
umax=10 #large for Trans and DRSC,bad for others

set.seed(2023)
W_q1=diag(1,q,q)
W_q2=matrix(0,q,p-q)
W=cbind(W_q1,W_q2)+matrix(runif(q*p,-0.1,0.1),q,p)
diag(W[1:q,1:q])=1


for (i in 1:iter){
  print(i)
  set.seed(2023+i)
  seu0 <- gendata_ST(height=n, width=n,W=W,p=p, umax=umax,platform="ST",K=K,beta=beta,error=error,tau=5)#for ST data
  y <- seu0$true_clusters
  counts.st=t(as.matrix(seu0@assays$RNA@counts))
  X=counts.st#log(counts.st+1)#scale(counts.st) #dim:n*p
  
  ind1=which(y==1);n1=length(ind1)
  ind2=which(y==2);n2=length(ind2)
  ind3=which(y==3);n3=length(ind3)
  ind4=which(y==4);n4=length(ind4)
  ind_st=c(ind1[1:round(n1/4)],ind2[1:round(n2/4)],ind3[1:round(n3/4)],ind4[1:round(n4/4)])
  Xsc=X[-ind_st,]
  z_sc=y[-ind_st]
  Xst=X[ind_st,]
  z_st=y[ind_st]
  pos=cbind(seu0@meta.data$row,seu0@meta.data$col)[ind_st,]

  ###Creat Seurat object
  counts=t(Xst)
  rownames(counts) <- paste0("gene", seq_len(p))
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
  
  ### Extract data
  Xst=log(Xst+1)
  raw=cbind(Xst,pos) #
  rownames(raw) <- paste0("spot", seq_len(nrow(Xst)))
  hZ=approxPCA(Xst, q=q)$PCs
  
  ### DR.SC -----------------------------------------------------------------
  set.seed(2023+i)
  res_drsc <- DR.SC(seu, q=q,K=K, platform = "Visium",  verbose=FALSE)
  res_drsc_K <- DR.SC(seu, q=q,K=2:10, platform = "Visium",  verbose=FALSE)
  KMat[i, 2] =length(unique(res_drsc_K$spatial.drsc.cluster))
  # ij <- find_neighbors(pos, platform="Visium")
  # library(Matrix)
  # Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1)
  # res_drsc <- simulDRcluster(Xst, Adj_sp=Adj_sp,  q=10, K=K, verbose=F)
  
  ### Seurat -----------------------------------------------------------------
  brain <- SCTransform(seu, verbose = FALSE)
  brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
  brain <- FindNeighbors(brain, reduction = "pca", dims = 1:10)
  brain <- FindClusters(brain, verbose = FALSE)
  KMat[i, 3] =length(unique(brain$seurat_clusters))

  # k-means -----------------------------------------------------------------
  y_kmeans = kmeans(hZ, centers=K)$cluster
  
  ### GMM
  res_gmm0=Mclust(hZ, G=K, modelNames='EEE',verbose=FALSE)
  res_gmm0_K=Mclust(hZ, G=NULL, modelNames='EEE',verbose=FALSE)
  KMat[i, 4] =res_gmm0_K$G

  # ### SC-MEB
  sce <- SingleCellExperiment(colData=cdata,reducedDims=SimpleList(PCA=hZ))
  Adj_sp <- find_neighbors2(sce, platform = "Visium")
  scmeb = SC.MEB(hZ, Adj_sp, K_set= K, parallel=TRUE, num_core = 3, PX = TRUE)
  out1 = SC.MEB::selectK(scmeb, K_set = 4,criterion = "BIC")
  
  scmeb_K = SC.MEB(hZ, Adj_sp, K_set= 2:10, parallel=TRUE, num_core = 3, PX = TRUE)
  out2 = SC.MEB::selectK(scmeb_K, K_set = 2:10,criterion = "BIC")
  KMat[i, 5] =out2$best_K_BIC

  # Louvain
  g.jaccard = scran::buildSNNGraph(sce, use.dimred="PCA", type="jaccard")
  y_louvain <- igraph::cluster_louvain(g.jaccard)$membership
  KMat[i, 6] =length(unique(y_louvain))

  ## Leiden
  resolution=quantile(strength(g.jaccard))[2] / (gorder(g.jaccard) - 1)
  ldc <- cluster_leiden(g.jaccard, resolution_parameter = resolution)
  KMat[i, 7] =length(unique(ldc$membership))
  
  ##Proposed
  res=trans_k_bic(q=q,K_set=2:10,Ux=cbind(hZ,pos),zeta=0.00001,method="HBIC",verbose=FALSE)
  K_ours=res$K
  
  res6=spGMM(K=K,U=cbind(hZ,pos),nS=nrow(hZ),maxIter=100,zeta=0.00001,neighbor=neighbor,verbose=FALSE)
  set.seed(2023+i)
  res2=TransST3(Xsc=Xsc,z_sc=z_sc,Xst=raw,q=q,K=K,nS=nrow(raw),maxIter=100,zeta=0.00001,combine=FALSE,neighbor=neighbor,verbose=FALSE)
  
  if(K_ours==K){
    res2_K=res2
    res6_K=res6
  }else{
    res6_K=spGMM(K=K_ours,U=cbind(hZ,pos),nS=nrow(hZ),maxIter=100,zeta=0.00001,neighbor=neighbor,verbose=FALSE)
    set.seed(2023+i)
    res2_K=TransST3(Xsc=Xsc,z_sc=z_sc,Xst=raw,q=q,K=K_ours,nS=nrow(raw),maxIter=100,zeta=0.00001,combine=FALSE,neighbor=neighbor,verbose=FALSE)
  }
  KMat[i, 1] =K_ours
  
  
  KMat[i,]=abs(KMat[i,]-K)

  #5,11 TransST; 6,12 spGMM
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
  
}

colMeans(final)
colMeans(KMat)



