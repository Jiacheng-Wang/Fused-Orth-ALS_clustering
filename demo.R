#  #######################################################################
#       File-Name:      demo.R
#       Date:           Wed Apr 27 18:09:59 2022
#       Author:         JCW
#       Purpose:        
#       Input Files:    NONE
#       Output Files:   NONE
#       Data Output:    NONE
#       Dependencies:   NONE
#       Status:         In Progress
#  #######################################################################
set.seed(123)
# generate a random tensor and specify the 4 clusters for the third mode
generate_tensor_cp_randind <- function(d, mu, n){
  beta_11 <- c(mu, -mu, 0.5*mu, -0.5*mu, rep(0, d-4))/norm(c(mu, -mu, 0.5*mu, -0.5*mu, rep(0, d-4)), type = "2")
  beta_21 <- beta_11
  beta_12 <- c(0,0,0,0, mu, -mu, 0.5*mu, -0.5*mu, rep(0, d-8))/norm(c(0,0,0,0, mu, -mu, 0.5*mu, -0.5*mu, rep(0, d-8)), type = "2")
  beta_22 <- beta_12
  index = split(1:n, 1:4)
  beta_31 = c(rep(NA,n))
  beta_32 = c(rep(NA,n))
  beta_31[index[[1]]] = mu; beta_32[index[[1]]] = mu;
  beta_31[index[[2]]] = -mu; beta_32[index[[2]]] = mu;
  beta_31[index[[3]]] = mu; beta_32[index[[3]]] = -mu;
  beta_31[index[[4]]] = -mu; beta_32[index[[4]]] = -mu;
  
  w_1 = norm(c(mu, -mu, 0.5*mu, -0.5*mu, rep(0, d-4)), type = "2")^2 * norm(beta_31, type = "2")
  w_2 = norm(c(0,0,0,0, mu, -mu, 0.5*mu, -0.5*mu, rep(0, d-8)), type = "2")^2 * norm(beta_32, type = "2")
  
  
  beta_31 = beta_31/norm(beta_31,type = "2")
  beta_32 = beta_32/norm(beta_32,type = "2")
  
  # generate tensor example using the singular vectors via tensor CP decomposition structure
  tensor = as.tensor(w_1 * (beta_11 %o% beta_21 %o% beta_31) + w_2 * (beta_12 %o% beta_22 %o% beta_32))
  return(list(beta = list(beta_11, beta_21, beta_31, beta_12, beta_22, beta_32),
              w = c(w_1, w_2),
              index = index,
              tensor = tensor))
}

# generate a random tensor
tnsr <- generate_tensor_cp_randind(10,1.2,20)
# perform the fused orthogonalized ALS algorithm
res <- fuse_orth_als_graph(tnsr$tensor, 2, max_iter=100, lambda=0.1,knn_size=2, tol = 1e-4)

# define the clustering arruracy evaluation
clus.dist = function(clus1, clus2) {
  ## Input two vectors of membership, output the disagreement between them. if one of them is truth, output is the clustering error.
  if (length(clus1)!=length(clus2)) return("cluster sizes don't match!")
  n=length(clus1)
  s=0
  for ( i in 2:n ) {
    for ( j in 1:(i-1)) {
      s=s+as.numeric(clus1[i]==clus1[j] & clus2[i]!=clus2[j])
      s=s+as.numeric(clus1[i]!=clus1[j] & clus2[i]==clus2[j])
    }
  }
  return(2*s/(n*(n-1)))
}

# transfer the true clustering assignment into a vector
clus.true <- function(index, group_size){
  label = c(rep(NA, length(unlist(index))))
  for (i in 1:group_size){
    ind = index[[i]]
    for (j in 1:length(ind)){
      label[ind[j]] = i
    }
  }
  return(label)
}

# get the cluster assignment estimates from the algorithm
clus.est <- function(matrix, group_size,iter.max = 100, nstart = 100){
  label = kmeans(matrix, centers = group_size, nstart = nstart,iter.max = 100)$cluster
  return(label)
}

# compare the estimate and true clustering assignment
label_true = clus.true(tnsr$index, group_size = 4)
label_est = clus.est(res$C, group_size = 4)
clus.dist(label_true, label_est) # clustering error is 0
