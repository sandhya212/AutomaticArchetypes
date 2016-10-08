###
### R code for Prabhakaran, Sandhya, Sudhir Raman, Julia E. Vogt, and Volker Roth. "Automatic model selection in archetype analysis." In Joint DAGM (German Association for Pattern Recognition) and OAGM Symposium, pp. 458-467. Springer Berlin Heidelberg, 2012.
###
### Code authors - Prof. Dr. Volker Roth and Sandhya Prabhakaran,
###                Department of Mathematics and Computer Science
###                University of Basel, Switzerland
###
### February 2012.
###


rm(list=ls())

library(MCMCpack)

do.simulate <- TRUE
if(do.simulate){
  n.archet <- 3  ## number of true archetypes
  n <- 1000 #1000001   ## number of points
  d <- 100 ## dimension
  X <- matrix( rnorm(n.archet*d),nrow = n.archet, ncol = d)
  X <- scale(X, scale = FALSE)
  Z <- matrix(0,nrow=d,ncol=n)
  Z[,(1:n.archet)] <- t(X) 
 
  for(j in (n.archet+1):n){
    rd <- t(rdirichlet(1,rep(2,n.archet)))
    
    Z[,j] <- t(X) %*% rd + 0.5*rnorm(d) ## convex combinations + noise
  }
  
  X <- t(Z) ## the first n.archet rows in X are the true archetypes, the other rows are convex  combinations
}


X <- scale(X, scale = FALSE) ## center X

eig <- eigen(t(X) %*% X) ## eigenvalue decomposition needed for PCA



max.n.archet <- 20  ## max number of archetypes considered
proj.dim <- max.n.archet ## no need to consider larger PCA-subspace... 

evecs <- (eig$vectors)[,1:proj.dim]


X <- X%*%evecs ## project onto first proj.dim eigenvectors


Z.orig <- X[1:n.archet,] ## the true archetypes (stored for plotting)


n <- nrow(X)
d <- ncol(X)

Z <- X[( n.archet+1):(n.archet+1+max.n.archet),] ## initialize archetypes to some points (which do NOT include the true archetypes)

X.pro <- X[,1:2] ## just for plotting...





do.convh <- TRUE ## sample points near convex hull by combining points on convex hull in 2D projections
if(do.convh){
  require("geometry")
  max.conv.dim <- 5
  convh.idx <- NULL
  if(max.conv.dim > ncol(X)){max.conv.dim <- ncol(X)}
  for(i in 1:(max.conv.dim-1)){
    for(j in (i+1):max.conv.dim){
      conv.tri <- convhulln(X[,i:j],"QJ") 
      convh.idx <- unique(c(convh.idx,unlist(apply(conv.tri,2,unique))))
      #points(X.pro[convh.idx,1:2], col = "red" , pch = 15 )
    }
  }

  
  
  X <- X[convh.idx,]
  X.pro <- X.pro[convh.idx,]
  
  n <- nrow(X)
  d <- ncol(X)
  
}



plot(X.pro)
points(Z.orig[,1:2],pch = 2, col = "magenta",cex=2)


## for C++
n <- nrow(X)

A <- matrix(0,nrow= n,ncol = nrow(Z))
B <- t(A)

NNlasso <- function(y,X,kappa, epsilon, XtX =  t(X) %*% X){

  d <- ncol(X)
  n <- length(y)
  r <- y
  beta <- rep(0,d)
  y.hat <- rep(0,n)
 

  nonzero <- 0
  l.1 <- 0
  j <- -1
  while(l.1 < kappa-epsilon & nonzero < min(d,n)+1 & sum(r*r) > 1e-10){
   
    if(j<0){
      correl <- t(X) %*% r
    }else{
                                       
      correl <- correl-epsilon*XtX[,j]
    }
   
   
    if(max(correl) <1e-10){
      break
    }
    j <- which.max(correl)
    beta[j] <-  beta[j] + epsilon
    r <- r -epsilon*X[,j]
    y.hat <- y.hat + epsilon*X[,j]
    l.1 <- l.1 +  epsilon
    nonzero <- length(which(beta>1e-10))
    
  }

  list(beta = beta, l.1 = l.1, nzn = nonzero,r=r, y.hat= y.hat)

}

NNlasso.large <- function(y,X,kappa, epsilon){

  d <- ncol(X)
  n <- length(y)
  r <- y
  beta <- rep(0,d)
  y.hat <- rep(0,n)
 

  nonzero <- 0
  l.1 <- 0
  j <- -1
  while(l.1 < kappa-epsilon & nonzero < min(d,n)+1 & sum(r*r) > 1e-10){
   
   
    correl <- t(X) %*% r
  
   
    if(max(correl) <1e-10){
      break
    }
    j <- which.max(correl)
    beta[j] <-  beta[j] + epsilon
    r <- r -epsilon*X[,j]
    y.hat <- y.hat + epsilon*X[,j]
    l.1 <- l.1 +  epsilon
    nonzero <- length(which(beta>1e-10))
    
  }

  list(beta = beta, l.1 = l.1, nzn = nonzero,r=r, y.hat= y.hat)

}


dyn.load("MultiTaskLassoForR.so")

lambda <- 1e-1
init.idx <- c(0,1)


X.work <- X

for(iter in 1:10){


  n <- nrow(X.work)
  
   
  
  for(i in 1:n){
    
    NNL <- NNlasso(X.work[i,],t(Z),1, 5e-2)
    beta <- NNL$beta
                                        #cat(NNL$l.1, sum((NNL$r)^2),NNL$nzn,'\n')
    A[i,] <- beta
    
  }
  
  
  
  
  if(iter < 2){
    Z.tilde <- solve(t(A) %*% A + lambda *diag(1,ncol(A))) %*% t(A) %*% X.work ## start with one "standard" estimate (without the lasso)
  }else{
    
    Amat <- A
    for(i in 2:ncol(X)){
      Amat <- rbind(Amat,A)
    }
    Amat <- cbind(Amat,rep(1e-10,nrow(Amat)))
    
    
    
    N <- nrow(Amat) 
    n.tasks <- ncol(X) ## number of tasks
    dim <- ncol(Amat) ## dimension
    task.labels <- rep(0:(n.tasks-1), each = nrow(X.work))
    y <- as.vector(X.work) ## the targets
    beta.path <- matrix(0, nrow = dim, ncol= n.tasks)
    LOGREG <- 0
    infty <- 0
    kappa <- 25
    errtol <- 3e-2
    stepsize <- 2e-4
    iterN <- 7
    modulo <- 500
    update.groups <- 2
    max.iter <- 10000
    max.groups <- nrow(Z)+1
    
    
    MTL <- .C("MultiTaskLasso",   as.integer(init.idx), as.integer(N),  as.integer(n.tasks),  as.integer(dim),  as.double(t(Amat)),  as.integer(task.labels), as.double(y),  beta=as.double(beta.path), as.integer( LOGREG), as.integer( infty),  as.double( kappa), as.double( errtol), as.double( stepsize), as.integer( iterN), as.integer( modulo), as.integer( update.groups), as.integer(max.iter), as.integer(  max.groups ))
    
    beta.path <- matrix(MTL[["beta"]], nrow = dim, ncol= n.tasks, byrow=TRUE)
    Z.tilde <- beta.path[1:(dim-1),]
    L.1.order <- order(apply(abs(Z.tilde), 1, sum), decreasing = TRUE)
    init.idx <- L.1.order[1:2]-1 ## select best 2 initial groups for next run of the group-lasso 
  }
  
  plot(X.pro)
  
  points(Z.tilde[,1:2], col = "cyan" , pch = 11 ,cex=2)
  
  
  
  for(i in 1:nrow(Z.tilde)){
    
    NNL <- NNlasso.large(Z.tilde[i,], t(X.work), 1, 2e-2)
    beta <- NNL$beta
    
    B[i,] <- beta
    
    Z[i,]  <-  NNL$y.hat
  }
  

  
  points(Z[,1:2], col = "red" , pch = 15 )
  
  B.colsum <- apply(B,2,sum) 
  w.colsum <- which(B.colsum > 5e-5) ## the "active" points
  points(X.work[w.colsum,1:2], col = "green"  )
  
  
  points(Z.orig[,1:2],pch = 2, col = "magenta",cex=2  )
  
  X.work <- X.work[w.colsum,] 
  A <- A[w.colsum,]
  B <- t(A)
  
  
  
}

title(main = "Automatic Archetype analysis")
legend("topright", legend=c("X", "X.work","Z.tilde","Z","Z.orig"), col=c("black","green","cyan","red","magenta"), pch=c(1,1,11,15,2), cex=0.8)
