
#step Input
delta <- 1e-5;
mu_temp <- 0.99*sqrt(sum(X^2)); #0.99*norm(X);

mu_bar <- delta*mu_temp;
eta <- 0.9;
tol <- 0.000001*norm(X,"F");
stop_crit <- FALSE;

n <- nrow(X)
d <- ncol(X)

#step 1
L_k <- matrix(0,nrow=n,ncol=d); # low-rank output matrix
C_k <- matrix(0,nrow=n,ncol=d);	# column-sparse output matrix
L_new <- matrix(0,nrow=n,ncol=d);
C_new <- matrix(0,nrow=n,ncol=d);
L_k_1 <- matrix(0,nrow=n,ncol=d); #L_{k-1}
C_k_1 <- matrix(0,nrow=n,ncol=d); #C_{k-1}
t_k <- 1;
t_k_1 <- 1;
k <- 0;
lambda <- 1


#step 2
while (stop_crit == FALSE){
	
	Y_L <- L_k + ((t_k_1 -1)/t_k) * (L_k - L_k_1); #step 3
	Y_C <- C_k + ((t_k_1 -1)/t_k) * (C_k - C_k_1); #step 3	
    	X_diff <- Y_L + Y_C - X;
	G_L <- Y_L - 0.5 * (X_diff); #step 4
	G_C <- Y_C - 0.5 * (X_diff); #step 4	

	#step 5
   	s <- svd(G_L);
	S <- diag(s$d);
	for ( i in 1:nrow(S)){  # soft thresholding
		if (S[i,i] > mu_temp/2){
			S[i,i] = S[i,i] - mu_temp/2;
		}else if(S[i,i] < (mu_temp/2)){ # the -ve sign
			S[i,i] = S[i,i] + mu_temp/2;	
		}else{ 
			S[i,i] = 0;
		}#end if
		#if (S[i,i] <= mu_temp/2){
		#	S[i,i] = 0;
		#}else{ 
		#	S[i,i] = S[i,i] - sign(S[i,i]) * mu_temp/2;	
		#}
	
	}#end for

	L_new <- s$u %*% S %*% t(s$v);

	#step 6
	eps <- 0.5 * (mu_temp * lambda);
	for ( i in 1:ncol(G_C)){
		temp_col <- G_C[,i];
		norm_temp_col <- sqrt(sum(temp_col^2)); # Euclidean norm of column
		if (norm_temp_col > eps){
   			temp_col = temp_col - temp_col * eps/norm_temp_col;
	  	}else{
			temp_col = matrix(0,nrow=nrow(G_C),1);
		}
		
		C_new[,i] <- temp_col;	
	}
	
	t_new <- 0.5 * (1 + sqrt(4*t_k^2 + 1)); #step 7
	mu_new <- max(eta * mu_temp, mu_bar); #step 7	

	#step 8 stopping criterion
	S_L=2 * (Y_L-L_new)+(L_new+C_new-Y_L-Y_C);
  	S_C=2 * (Y_C-C_new)+(L_new+C_new-Y_L-Y_C);
	if (norm(S_L, "F")^2 + norm(S_C, "F")^2 <= tol^2){
		stop_crit = TRUE;
	}else{
    		L_k_1 <- L_k;
		L_k <- L_new;
		C_k_1 <- C_k;
		C_k <- C_new;
		t_k_1 <- t_k;
		t_k <- t_new;
		mu_temp <- mu_new;
		k=k+1
		end  	
	}

}# end while

X <- L_new;
eig <- eigen(t(X) %*% X) ## eigenvalue decomposition needed for PCA
proj.dim <- max.n.archet ## no need to consider larger PCA-subspace... 
evecs <- (eig$vectors)[,1:proj.dim]

X <- X%*%evecs;

