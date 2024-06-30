library(dbscan)
#library(irlba)
library(mclust)
library(pracma)
library(Rcpp)
library(parallel)
library(bigstatsr)

find_neighbors <- DR.SC:::find_neighbors

sourceCpp('help/utilSimulDRcluster.cpp')
sourceCpp('help/wpca.cpp')
#source("help/help.R")

#svd issue appreas in approxPCA
approxPCA <- function(X, q){ ## speed the computation for initial values.
  n <- nrow(X)
  svdX  <- svd(X, nu=q,nv = q)
  PCs <- svdX$u %*% diag(svdX$d[1:q])
  loadings <- svdX$v
  # dX <- PCs %*% t(loadings) - X
  # Lam_vec <- colSums(dX^2)/n
  return(list(PCs = PCs, loadings = loadings))
}

wpca <- function(X, q, weighted=TRUE){
  if(!is.matrix(X)) stop("wpca: X must be a matrix!")
  if(q< 1) stop("wpca: q must be a positive integer!")
  X <- scale(X, scale=F) # centralize
  out <- wpcaCpp(X, q, weighted)
  return(out)
}

mycluster <- function(Z, G,init="GMM", int.model='EEE'){ #Mclust package: for each sample
  #set.seed(111)
  mclus2 <- Mclust(Z, G=G, modelNames=int.model,verbose=FALSE)
  return(mclus2)
}

parfun_int <- function(k, Z, int.model='EEE'){ 
  #for each sample
  #Z: low-dimensional representation
  #k: # of clusters
  mclus2 <- mycluster(Z, k, int.model)
  yveck <- mclus2$classification
  Mu0k <- mclus2$parameters$mean #q*K matrix
  Sigma0k <- mclus2$parameters$variance$sigma #array[,,k]
  return(list(yveck=yveck, Mu0k = Mu0k, Sigma0k=Sigma0k,K=mclus2$G))
}

fast_init=function(K,U,q){
  # needs to determine q first
  mu <- matrix(NA, K,q)
  Sigma <- array(NA, dim = c(q, q, K))
  #princ <- wpca(U, q, weighted = TRUE) #approxPCA(U,q) #use PCA for dim redu
  princ <- approxPCA(U,q)
  W0 <- princ$loadings
  hZ <- princ$PCs
  res=parfun_int(K,hZ)
  
  z=res$yveck #labels
  P=matrix(1/K,K,K)
  Lambda=runif(ncol(U),0.5,2) #is a vector
  for (k in 1:K){
    mu[k,]=res$Mu0k[,k]
    Sigma[,,k]=res$Sigma0k[,,k] 
  }
  return(list(z=z,W=W0,P=P,Lambda=Lambda,mu=mu,Sigma=Sigma,K=res$K,hZ=hZ))
}


TransST3=function(Xsc,z_sc,Xst,q,K,nS,maxIter=100,zeta=0.00001,neighbor="knn",combine=FALSE,weighted=FALSE,lambda_len=20,verbose=FALSE){
  ################
  #pLDR with Xsc and Xst by EM+update W+spGMM
  ################
  message("Running Algorithm TransST3")
  
  Xst1=Xst[,1:(ncol(Xst)-2)]
  Xst2=Xst[,(ncol(Xst)-1):ncol(Xst)]
  
  if(combine==FALSE){
    res=pLDR(q,K,Xsc,z_sc,maxIter=100,standardize=FALSE,weighted=weighted,zeta=zeta,verbose=FALSE)
    W0=res$W
    z_st=NULL
  }else if(combine==TRUE){
    if (is.null(Xsc)|is.null(z_sc)){
      princ <- approxPCA(Xst1,q)
      W0 <- princ$loadings
      z_st=NULL
    }else{
      res_pLDR=pLDR3(q,K,Xsc,z_sc,Xst1,maxIter=maxIter,standardize=FALSE,zeta=zeta,verbose=FALSE)
      W0=res_pLDR$W
      z_st=res_pLDR$z_st
    }
  }
  
  #Step 2: estimate ST v by PCA
  lambda_list=c(0,0.01*(1:lambda_len))
  lambda_error=rep(NA,length(lambda_list))
  for (lambda in lambda_list){
    num_rows <- nrow(Xst1)
    shuffled_indices <- sample(1:num_rows)
    half_rows <- num_rows %/% 2
    Xst1_1 <- Xst1[shuffled_indices[1:half_rows], ]
    Xst1_2 <- Xst1[shuffled_indices[(half_rows + 1):num_rows], ]
    res_W1=AdaptiveW(Xst1_1,W0,lambda,q,maxIter=50)
    res_W2=AdaptiveW(Xst1_2,W0,lambda,q,maxIter=50)
    error=(norm(res_W1$v%*%t(res_W2$W)-Xst1_1,type="F"))^2+(norm(res_W2$v%*%t(res_W1$W)-Xst1_2,type="F"))^2
    lambda_error[which(lambda_list==lambda)]=error
  }
  lambda=lambda_list[which.min(lambda_error)]
  ad_res=AdaptiveW(Xst1,W0,lambda,q,maxIter=50)
  v=ad_res$v
  #v=v=Xst1%*%W0%*%solve(t(W0)%*%W0)
  v1=cbind(v,Xst2)
  
  #Step 3: clustering ST data based on v
  res=spGMM(K,v1,nS,maxIter=maxIter,zeta=zeta,neighbor=neighbor,verbose=FALSE)
  
  #z: estimated label; v: low-d; z_st: initially estimated label
  return(list(z=res$z,v=v,z_st=z_st))
}


pLDR <- function(q,K,U,z,maxIter=100,standardize=FALSE,weighted=FALSE,zeta=zeta,verbose=FALSE) { 
  #use source only for W
  
  K=length(unique(z))
  
  mu <- matrix(NA, K,q)
  Sigma <- array(NA, dim = c(q, q, K))
  if(weighted==FALSE){
    princ <- approxPCA(U,q)
  }else if(weighted==TRUE){
    princ <- wpca(U, q, weighted = TRUE)
  }
  
  W <- princ$loadings
  hZ <- princ$PCs
  Lambda=runif(ncol(U),0.5,2) #is a vector
  
  for (k in 1:K){
    mu[k,]=colMeans(hZ[z==k,])
    Sigma[,,k]=cov(hZ[z==k,])
  }
  
  iter <- 1
  loglik <- rep(NA,maxIter)
  loglik[1] <- -1000000
  p=ncol(U)
  n=nrow(U)
  
  Rik=matrix(0,n,K)
  indices <- cbind(1:n, z)  # Creating row and column indices for the 1 values
  Rik[indices] <- 1
  N <- colSums(Rik)
  Iden=0.00001*diag(1,q,q)
  
  for (iter in 2:maxIter) { #iteration
    res=runICM_sp2(U, W,Lambda,mu, Sigma, Iden)
    Ex_temp=res$Ez
    Cki_ara=res$Cki_ara
    
    #update mu
    for (k in 1:K) {mu[k,] <- t(Rik[, k])%*%Ex_temp[,,k]/N[k]} #mu:K*q matrix
    #update Sigma
    Sigma = update_Sigma0(R=Rik, Ez=Ex_temp, Ci=Cki_ara, Mu=mu, N=N)
    #update W
    W = update_W0(X=U, R=Rik, Ez=Ex_temp, Ci=Cki_ara, N=N)
    #update Lambda
    Lambda = as.vector(update_Lam(R=Rik, X=U, W=W, Ez=Ex_temp, Ci=Cki_ara))
    loglik[iter]=res$loglik
    if (abs((loglik[iter] - loglik[iter - 1]) / loglik[iter - 1]) < zeta) {break} #used 0.001 in the simulation
  }
  a=loglik[-1]
  a=a[!is.na(a)]
  
  return(list(loglike=a,q=q,K=K,W=W,Lambda=Lambda))
}


pLDR3 <- function(q,K,Xsc,z_sc,Xst,maxIter=100,standardize=FALSE,zeta=zeta,verbose=FALSE) {
  #combine Xsc and Xst, z_st by model in nature
  K=max(K,length(unique(z_sc)))

  U=rbind(Xsc,Xst)
  # int_values=fast_init(K,U,q)
  # W=int_values$W
  # Lambda=int_values$Lambda
  # hZ=int_values$hZ
  # mu=int_values$mu #mu, Sigma
  # Sigma=int_values$Sigma
  
  mu <- matrix(NA, K,q)
  Sigma <- array(NA, dim = c(q, q, K))
  princ <- approxPCA(U,q)
  W <- princ$loadings
  hZ <- princ$PCs
  Lambda=runif(ncol(U),0.5,2) #is a vector
  
  for (k in 1:K){
    mu[k,]=colMeans(hZ[z_sc==k,])
    Sigma[,,k]=cov(hZ[z_sc==k,]) 
  }
  
  iter <- 1
  loglik <- rep(NA,maxIter)
  loglik[1] <- -1000000

  p=ncol(U)
  n=nrow(U)
  n1=nrow(Xsc)
  n2=nrow(Xst)
  
  Ux=matrix(NA,n,K)
  Cki_ara=array(NA, dim = c(q, q, K))
  Ex_temp=array(NA, dim = c(n, q, K))
  #R=matrix(NA,n1,length(unique(z_sc)))
  
  for (iter in 2:maxIter) { #iteration
    lambda_inverse=pinv(diag(Lambda)) 
    WtLW <- t(W) %*% lambda_inverse%*%W #q*q
    XLW <- t(W) %*%lambda_inverse%*%t(U) #W'*Lambda*U: q*n
    for (k in 1:K){
      C_k=WtLW+pinv(Sigma[,,k]) #q*q
      Ci=pinv(C_k)#inverse of C_k: q*q
      Cki_ara[,,k]=Ci
      res=multi_det_SkCpp2(X=U, Lam_vec0=Lambda, W0=W, Ck=C_k, Muk=mu[k,], Sigmak=Sigma[,,k])
      Ux[, k] <- -0.5 * res$logdSk + 0.5 * res$mSk # negative loglike
      Ex_temp[,,k] <- t(Ci %*% (XLW + matrix(rep(solve(Sigma[,,k])%*%mu[k,], n), ncol = n, byrow = FALSE))) #<xi>: n*q
    }
    maxA1 <- apply(-Ux, 1, max)
    Ux <- (-Ux - matrix(rep(maxA1, each = K), nrow = nrow(Ux), ncol = K,byrow=TRUE))
    loglik_more_vec <- rowSums(exp(Ux))
    LL <- sum(log(loglik_more_vec) + maxA1)
    Rik <- exp(Ux) / matrix(rep(loglik_more_vec, each = K), nrow = nrow(Ux), ncol = K,byrow=TRUE)
    z_st <- apply(Ux, 1, which.max)[(n1+1):n]  # label
    indices <- cbind(1:n1, z_sc)  # Creating row and column indices for the 1 values
    Rik[indices] <- 1
    N <- colSums(Rik)
    #N2=colSums(Rik[1:n1,])
    
    #update mu
    for (k in 1:K) {mu[k,] <- t(Rik[, k])%*%Ex_temp[,,k]/N[k]} #mu:K*q matrix
    #update Sigma
    Sigma = update_Sigma0(R=Rik[,], Ez=Ex_temp[,,], Ci=Cki_ara, Mu=mu, N=N)
    #update W
    W = update_W0(X=U, R=Rik, Ez=Ex_temp, Ci=Cki_ara, N=N)
    #W = update_W0(X=U, R=R, Ez=Ex_temp, Ci=Cki_ara, N=N)
    #update Lambda
    Lambda = as.vector(update_Lam(R=Rik[,], X=U[,], W=W, Ez=Ex_temp[,,], Ci=Cki_ara))
    loglik[iter]=LL
    if (abs((loglik[iter] - loglik[iter - 1]) / loglik[iter - 1]) < zeta) {break} #used 0.001 in the simulation
  }
  a=loglik[-1]
  a=a[!is.na(a)]
  
  return(list(loglike=a,q=q,K=K,W=W,Lambda=Lambda,z_st=z_st))
}


AdaptiveW=function(Xst1,W0,lambda,q,maxIter=50){
  iter=1
  loss <- rep(NA,maxIter)
  loss[1] <- 1e9
  W=W0
  #v=Xst1%*%W0%*%solve(t(W0)%*%W0)
  for (iter in 2:maxIter) {
    v=Xst1%*%W%*%solve(t(W)%*%W)
    W=(lambda*W0+t(Xst1)%*%v)%*%solve(diag(lambda,q,q)+t(v)%*%v)
    loss[iter]=(norm(Xst1-v%*%t(W),"F"))^2#+lambda*(norm(W-W0,"F"))^2
    if (abs((loss[iter] - loss[iter - 1]) / loss[iter - 1]) < 1e-3) {break} #used 0.001 in the simulation
  }
  return(list(W=W,v=v))
}


spGMM <- function(K,Ux,nS,maxIter=100,zeta=zeta,neighbor="knn",verbose=FALSE) {
  # This function fits model v_i|z=k; z~MRF
  # nS: a vector containing sizes of diff samples
  
  U=Ux[,1:(ncol(Ux)-2)] #if spatial, then take the first q
  q=ncol(U)
  n=nrow(U)
  pos=Ux[,(ncol(Ux)-1):ncol(Ux)]
  
  res=parfun_int(K,U)
  z=res$yveck #labels
  mu <- matrix(NA, K,q)
  Sigma <- array(NA, dim = c(q, q, K))
  for (k in 1:K){
    mu[k,]=res$Mu0k[,k]
    Sigma[,,k]=res$Sigma0k[,,k] 
  }
  
  iter <- 1
  loglik <- rep(NA,maxIter)
  loglik[1] <- -10000
  Cki_ara=array(0, dim = c(q, q, K))
  Ex_temp=array(U, dim = c(dim(U), K))
  U2=matrix(NA,n,K)
  
  S_len=length(nS)
  if (neighbor=="knn"){
    if (S_len==1){
      id_matrix=kNN(pos, k = 5)$id
      temp_mt=matrix(0,nrow(pos),nrow(pos))
      temp=match_matrices(temp_mt,id_matrix)
      Adj <- as(temp, "sparseMatrix")
    }else{
      Adj=matrix(0,n,n)
      S_ind=cumsum(nS)
      for (s in 1:S_len){
        temp_mt=matrix(0,nS[s],nS[s])
        if(s==1){first=1}else{first=S_ind[s-1]+1}
        id_matrix <- kNN(pos[first:S_ind[s],], k = 5)$id
        temp=match_matrices(temp_mt,id_matrix)
        Adj[first:S_ind[s],first:S_ind[s]]=temp #Matrix::sparseMatrix(temp_mt)#Adj_temp#getAdj(pos[first:S_ind[s],], platform="ST")
      }
      Adj=as(Adj, "sparseMatrix")
    }
  }else if (neighbor=="l1"){
    if (S_len==1){
      #Adj <- getAdj(pos, platform="ST")
      ij <- find_neighbors(pos, platform="Visium")
      Adj <- Matrix::sparseMatrix(ij[,1], ij[,2], x = 1)
    }else{
      Adj=matrix(0,n,n)
      S_ind=cumsum(nS)
      for (s in 1:S_len){
        temp_mt=matrix(0,nS[s],nS[s])
        if(s==1){first=1}else{first=S_ind[s-1]+1}
        ij <- find_neighbors(pos[first:S_ind[s],], platform="Visium")
        temp_mt[ij]=1
        Adj[first:S_ind[s],first:S_ind[s]][ij]=temp_mt #Matrix::sparseMatrix(temp_mt)#Adj_temp#getAdj(pos[first:S_ind[s],], platform="ST")
      }
      Adj=as(Adj, "sparseMatrix")
    }
  }
  
  for (iter in 2:maxIter) { #iteration for whole
    for (k in 1:K){
      for (i in 1:n){
        U2[i,k] = 0.5*log(det(Sigma[,,k]))+ 0.5 * t(U[i,] - mu[k,])%*%pinv(Sigma[,,k])%*%(as.matrix(U[i,] - mu[k,]))
      }
    }
    res=runICM_sp(U,U2, z, mu, Sigma, Adj,alpha=rep(0,K), beta_grid=seq(0,1,0.1), beta=2, maxIter_ICM=50)
    Rik=res$R
    LL=res$loglik
    N <- colSums(Rik)
    z=res$y

    #update mu
    for (k in 1:K) {mu[k,] <- t(Rik[, k])%*%U/N[k]} #mu:K*q matrix
    #update Sigma
    Sigma = update_Sigma0(R=Rik, Ez=Ex_temp, Ci=Cki_ara, Mu=mu, N=N)
    loglik[iter]=res$loglik
    
    if (is.na(loglik[iter])){break}
    if (abs((loglik[iter] - loglik[iter - 1]) / loglik[iter - 1]) < zeta) {break} #used 0.001 in the simulation
  }
  
  #LL is the likelihood
  return(list(z=z,LL=LL,q=q,K=K))
}

trans_k_bic <- function(q,K_set=2:15,Ux,zeta,method="HBIC",verbose=FALSE) {
  n=nrow(Ux)
  num_cores  <- parallel::detectCores()
  cl  <- parallel::makeCluster(num_cores-2)
  registerDoParallel(cl)
  
  k_min=min(K_set)-1
  
  #Run forloop in Parallel
  x<- foreach(k=K_set, .combine=c) %dopar% {
    source("help/source-easy.R")
    res=spGMM(K=k,Ux=Ux,nS=n,maxIter=100,zeta=zeta,neighbor="knn",verbose=FALSE)
    LL=res$LL
    
    if(n<1000){
      const=0.5
    }else if(n>=1000){
      const=2
    }#old is 2
    
    pp=k*(q+q*(q+1)/2)+1 #number of parameters
    if (method=="BIC"){
      bic <-  -2*LL + const*pp* log(n)# BIC
    }else if (method=="MBIC"){
      bic <-  -2*LL + const*pp* log(n)* log(log(q+n))# mBIC
    }else if (method=="HBIC"){
      bic <-  -2*LL + const*pp* log(log(n)+pp)# HBIC
    }
    c(bic)
  }
  stopCluster(cl) #Stop cluster
  
  return(list(K=which.min(x)+k_min))
}


#### Generate Spatial data with ST platform
gendata_ST <- function(height=30, width=30, platform="ST", W=NULL,mu=NULL, p =100, q=10, K=7, umax=2,beta,G=4,homo=FALSE, tau=5, error=2,const=0,const_w=0){
  if(q <2) stop("error:gendata_sp::q must be greater than 2!")
  
  n <- height * width # # of cell in each indviduals 
  
  if(platform=="ST"){
    beta= beta
  }else if(platform=='SC'){
    beta = 0
  }
  ## generate deterministic parameters, fixed after generation
  if (homo==TRUE){
    Lambda <- rep(2,p)
  }else{
    Lambda <- runif(p,0.1,error)
  }
  
  if (is.null(W)){
    W1 <- matrix(rnorm(p*q), p, q)
    W <- t(qr.Q(qr(W1)))
  }
  W=W+const_w*matrix(rnorm(p*q),q,p)
  
  if (is.null(mu)){
    mu <- matrix(0, q,  K)
    if(q > K){
      q1 <- floor(K/2)+1
      for(j in 1:q1){
        mu[j,j] <- tau
      }
      mu[(q1+1):q, K] <- tau
    }else if(q <= K){
      for(k in 1:K)
        mu[,k] <- rep(tau/8 *k, q) #
    }
  }
  
  diagmat = array(0, dim = c(q, q, K))
  #Sigma <- outer(1:q, 1:q, function(x,y){0.5^(abs(x-y))})
  for(k in 1:K){
    tmp2  <- rep(1, q)
    if(k <= K/2){
      tmp2=runif(q,0.1,umax)#[q] <- 10
    }
    #diagmat[,,k]=matrix(0.1,q,q)
    #diagmat[,,k]=Sigma
    diag(diagmat[,,k]) <- tmp2
    eps <- rnorm(q,sd=0.1)
    diagmat[,,k] <-diagmat[,,k] + eps %*% t(eps)
  }
  
  Mu <- t(mu)
  Mu=Mu+const*matrix(runif(K*q),K,q)
  Sigma <- diagmat
  
  y <- sampler.mrf(iter = n, sampler = "Gibbs", h = height, w = width, ncolors = K, nei = G, param = beta,initialise = FALSE) #label
  y <- c(y) + 1
  
  Z <- matrix(0, n, q)
  for(k in 1:K){
    nk <- sum(y==k)
    Z[y==k, ] <- MASS::mvrnorm(nk, Mu[k,], Sigma[,,k])
  }
  #print(Z[1:10,])
  Ez <- colMeans(Z)
  #Mu <- Mu - matrix(Ez, K, q, byrow=T) # center Z
  X <- Z %*% W + MASS::mvrnorm(n, mu=rep(0,p), Sigma=diag(Lambda))
  #print(Z[1:5,1:5])
  # make position
  pos <- cbind(rep(1:height, width), rep(1:height, each=width))
  
  counts <- t(X) - min(X)
  p <- ncol(X); n <- nrow(X)
  rownames(counts) <- paste0("gene", seq_len(p))
  colnames(counts) <- paste0("spot", seq_len(n))
  counts <- as.matrix(exp(counts)-1)
  ## Make array coordinates - filled rectangle
  
  if(platform=="ST"){
    cdata <- list()
    cdata$row <- pos[,1]
    cdata$col <- pos[,2]
    cdata <- as.data.frame(do.call(cbind, cdata))
    cdata$imagerow <- cdata$row
    cdata$imagecol <- cdata$col 
    cdata$Z=Z
    row.names(cdata) <- colnames(counts)
    seu <-  CreateSeuratObject(counts= counts, meta.data=cdata) #
  }else if(platform=='SC'){
    seu <-  CreateSeuratObject(counts= counts)
  }else{
    stop("gendata_RNAExp: Unsupported platform \"", platform, "\".")
  }
  
  seu$true_clusters <- y
  return(seu)
}


match_matrices <- function(A, B) {
  result_matrix <- A  # Create a copy of matrix A to avoid modifying the original matrix
  for (i in 1:nrow(B)) {
    result_matrix[i, B[i,]] <- 1
  }
  return(result_matrix)
}

