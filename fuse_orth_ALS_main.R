#  #######################################################################
#       File-Name:      fused-orth_ALS_main.R
#       Date:           Wed Apr 27 17:36:23 2022
#       Author:         JCW
#       Purpose:        
#       Input Files:    NONE
#       Output Files:   NONE
#       Data Output:    NONE
#       Dependencies:   NONE
#       Status:         In Progress
#  #######################################################################
library(rTensor)
library(genlasso)
library(gtools)

fuse_orth_als_graph <- function(tnsr, r, max_iter, lambda,knn_size=2, tol = 1e-4){
  # get dimension along each mode
  d1 = tnsr@modes[1]
  d2 = tnsr@modes[2]
  d3 = tnsr@modes[3]
  
  # get initialization for each factor matrix
  A = matrix_init(d1,r)
  B = matrix_init(d2,r)
  C = matrix_init(d3,r)
  
  # matricization of tensor along different mode
  tnsr_unfold1 = unfold(tnsr,row_idx=1,col_idx=c(2,3))@data
  tnsr_unfold2 = unfold(tnsr,row_idx=2,col_idx=c(1,3))@data
  tnsr_unfold3 = unfold(tnsr,row_idx=3,col_idx=c(1,2))@data
  
  # fuse distance matrix 
  # D = fused_lasso_distance(d3)
  C_tmp = cp(tnsr, num_components = 2)$U[[3]]
  D_tmp = genlasso_distance(nrow(C_tmp))
  w_tmp = kernel_weights(C_tmp)
  wsparse = knn_weights(w_tmp,knn_size,d3)
  Dsparse = wsparse*D_tmp
  D_graph = Dsparse[apply(Dsparse,1,function(x)!all(x==0)),]
  
  
  # ALS update (without orthogonalization)
  iter = 0
  err = NULL
  while(iter<max_iter){
    A.old = A
    B.old = B
    C.old = C
    A_bar = qr_orth(A.old)
    B_bar = qr_orth(B.old)
    C_bar = qr_orth(C.old)
    A_est = tnsr_unfold1 %*% khatri_rao(C_bar, B_bar)
    B_est = tnsr_unfold2 %*% khatri_rao(C_bar, A_bar) 
    C_est = tnsr_unfold3 %*% khatri_rao(B_bar, A_bar) 
    C_est_fuse = general_fuse(C_est, lambda, D_graph)
    A = matrix_col_normalize(A_est)$matrix
    B = matrix_col_normalize(B_est)$matrix
    C = matrix_col_normalize(C_est_fuse)$matrix
    iter = iter + 1
    tnsr_hat = array(0, c(d1,d2,d3))
    for (i in 1:r){
      wi = Tabc(tnsr@data, A[,i], B[,i], C[,i])
      tnsr_hat = tnsr_hat + wi*outer(outer(A[,i],B[,i]),C[,i])
    }
    if (sum(norm(A.old-A,type ='2')+norm(B.old-B,type ='2')+norm(C.old-C,type ='2'))<tol){break}
    new = sqrt(sum((tnsr_hat -  tnsr@data)^2))/ sqrt(sum((tnsr@data)^2))
    err <- append(err, new)
    message(paste(iter,"-th  iteration -- relative error is", new," -----------------"))
  }
  tnsr_hat = array(0, c(d1,d2,d3))
  w = c(rep(NA,r))
  for (i in 1:r){
    w[i] = Tabc(tnsr@data, A[,i], B[,i], C[,i])
    tnsr_hat = tnsr_hat + w[i]*outer(outer(A[,i],B[,i]),C[,i])
  }
  return(list(w = w,
              A = A, 
              B = B,
              C = C, 
              D_graph = D_graph,
              err = err))
}

# dependent functions
qr_orth <- function(a){
  
  row_size <- nrow(a)
  col_size <- ncol(a)
  for (i in 1:(col_size-1)){
    a[,i] = a[,i]/norm(a[,i],'2')
    for (j in (i+1):col_size){
      a[,j] = a[,j] - as.numeric(a[,i]%*%a[,j])*a[,i]
    }
  }
  return(a)
}
genlasso_distance <-function(n){
  dist  <- matrix(rep(0, (n-1)*n/2*n), ncol = n)
  diff = gtools::combinations(n,2)
  for (i in 1:nrow(diff)){
    id1 = diff[i,][1]
    id2 = diff[i,][2]
    dist[i,id1] = 1
    dist[i,id2] = -1
  }
  return(dist = dist)
}
kernel_weights <- function(matrix, phi = 0.5){
  diff = gtools::combinations(nrow(matrix),2)
  w = c(rep(NA,nrow(diff)))
  for (i in 1:nrow(diff)){
    w[i] = exp(-phi*norm(matrix[diff[i,1],] - matrix[diff[i,2],],type = '2'))
  }
  return(w = w)
}

# find k nearest neighbors for each row and calculate the sparse distance.
tri2vec <- function(i,j,n) {
  return(n*(i-1) - i*(i-1)/2 + j -i)
}
knn_weights <- function(w,k,n) {
  i <- 1
  neighbors <- tri2vec(i,(i+1):n,n)
  keep <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  for (i in 2:(n-1)) {
    group_A <- tri2vec(i,(i+1):n,n)
    group_B <- tri2vec(1:(i-1),i,n)
    neighbors <- c(group_A,group_B)
    knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
    keep <- union(knn,keep)
  }
  i <- n
  neighbors <- tri2vec(1:(i-1),i,n)
  knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  keep <- union(knn,keep)
  w[-keep] <- 0
  return(w)
}

# Function to get the random initialization with unit form column
matrix_init <- function(d,r){
  init = matrix(rep(NA,d*r),ncol = r)
  for (i in 1:r){
    init[,i] = rnorm(d)
    init[,i] = init[,i]/norm(init[,i],'2')
  }
  return(init)
}
# Function to normalize the column of matrix to be norm 1
matrix_col_normalize <- function(a){
  w = c(rep(NA,ncol(a)))
  for (i in 1:ncol(a)){
    norm_col = norm(a[,i],'2')
    w[i] = norm_col
    a[,i] = a[,i]/norm_col
  }
  return(list(matrix = a,
              w = w))
}
Tabc = function(T, a, b, c){
  #compute tensor vector product T(a,b,c) for tensor T and vectors, a,b,c with different dimensions.
  
  d1 = dim(T)[1]	
  d2 = dim(T)[2]
  d3 = dim(T)[3]
  tmp = 0
  for(i in 1:d1){
    for(j in 1:d2){	
      for(k in 1:d3){		
        tmp = tmp + a[i]*b[j]*c[k]*T[i,j,k]
      }
    }
  }
  tmp
}

general_fuse <- function(mat, lambda, D){
  mat_fuse = matrix(rep(NA,nrow(mat)*ncol(mat)),nrow = nrow(mat))
  for (i in 1:ncol(mat)){
    output = genlasso(mat[,i], diag(nrow(mat)), D)
    mat_fuse[,i] = coef(output, lambda = lambda)$beta
  }
  return(mat_fuse = mat_fuse)
}
