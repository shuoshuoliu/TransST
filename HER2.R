rm(list=ls())
library(Seurat)
library(ggplot2)
library(dbscan)
library(MASS)
library(preprocessCore)
library(readr)

source("source-easy.R")

name_ID=c("A1","B1","C1","D1","E1","F1","G2","H1")
s_len=length(name_ID)

q=15
p=2000

data=list()
label=list()
position=list()
n=rep(NA,s_len)
for (i in 1:s_len){
  Xst<-read.delim(paste0(name_ID[i], ".tsv"))
  rownames(Xst)=Xst$X
  Xst$X=NULL
  meta=read.delim(paste0(name_ID[i], "_labeled_coordinates.tsv"))
  common_id=paste(round(as.numeric(meta$x)), "x",round(as.numeric(meta$y)), sep = "")
  if ("NAxNA"%in%common_id){
    na_id=which(common_id=="NAxNA")
    meta=meta[-na_id,]
    common_id=common_id[common_id!="NAxNA"]
    rownames(meta)=common_id
    meta=meta[common_id,]
  }else{
    rownames(meta)=common_id
    meta=meta[common_id,]
  }
  meta$Row.names=NULL
  meta$x=NULL
  meta$y=NULL
  if (nrow(meta)>nrow(Xst)){
    meta=meta[rownames(Xst),]
  }else if(nrow(meta)<nrow(Xst)){
    Xst=Xst[rownames(meta),]
  }else{Xst=Xst[rownames(meta),]}
  
  if(i%in%c(1,3)){
    un_index=which(meta$label%in%c("undetermined","immune infiltrate"))
  }else{
    un_index=which(meta$label=="undetermined")
  }
  
  if(length(un_index)==0){
    Xst=Xst
    meta0=meta
  }else{
    Xst=Xst[-un_index,]
    meta0=meta[-un_index,]
  }
  y=meta0$label
  
  # Divide each cell by the total number of molecules measured in the cell
  # Multiply that number by a scaling factor (i.e. 10000
  # Add 1, and take a natural log
  
  data[[i]]=t(apply(t(Xst), 2, process_column))
  #data[[i]]=log(Xst+1)
  position[[i]]=meta0[,1:2] # spatial coordinates
  label[[i]]=y
  n[i]=nrow(Xst)
  print(length(unique(y)))
}

d2=data

common_columns <- Reduce(intersect, lapply(d2, colnames))
subM <- lapply(d2, function(mat) mat[, common_columns, drop = FALSE])

method="BIC"
neighbor="knn"
zeta=0.001
seed=2023
est_label=list()

final=matrix(NA,s_len,14)
KMat=matrix(NA,s_len,7)
trueK_set=rep(NA,s_len)
for (i in 1:s_len){
  trueK_set[i]=length(unique(label[[i]]))
}

for (i in 1:s_len){
  print(i)
  tryCatch({
    Xst=as.matrix(subM[[i]])
    #top_columns <- colnames(Xst)[order(apply(Xst, 2, var), decreasing = TRUE)[1:p]]            #By variance
    top_columns <- colnames(Xst)[order(apply(Xst, 2, process_column2), decreasing = TRUE)[1:p]] #By CV
    Xst <- Xst[, top_columns]
    source=do.call(rbind, subM[-i])
    Xsc=source[,top_columns]
    
    pos=as.matrix(position[[i]])
    z_st=as.numeric(factor(label[[i]], levels = unique(label[[i]])))
    z_temp=do.call(c, label[-i])
    z_sc=as.numeric(factor(z_temp, levels = unique(z_temp)))
    K=trueK_set[i]
    n2=n[i]
    
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
    K_ours=trans_k_bic(q=q,K_set=2:10,Ux=cbind(hZ,pos),zeta=zeta,method="HBIC",verbose=FALSE)$K
    
    KMat[i, 1] =K_ours
    
    ######
    res6=spGMM(K=K,U=cbind(hZ,pos),nS=n_st,maxIter=100,zeta=zeta)
    set.seed(seed+i)
    res2=TransST3(Xsc=as.matrix(Xsc_new),z_sc=z_sc,Xst=as.matrix(raw),q=q,K=K,nS=n_st,zeta=zeta,weighted=FALSE)
    
    if(K_ours==K){
      res2_K=res2
      res6_K=res6
    }else{
      res6_K=spGMM(K=K_ours,U=cbind(hZ,pos),nS=n_st,zeta=zeta)
      set.seed(seed+i)
      res2_K=TransST3(Xsc=as.matrix(Xsc_new),z_sc=z_sc,Xst=as.matrix(raw),q=q,K=K_ours,nS=n_st,zeta=zeta,weighted=FALSE)
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
  }, error=function(e){})
}




