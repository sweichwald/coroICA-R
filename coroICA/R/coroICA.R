##' Estimates the unmixing and confounded sources of the coroICA model
##' X=A(S+H).
##'
##' For further details see the references.
##' @title coroICA
##' @param X data matrix. Each column corresponds to one predictor
##'   variable.
##' @param group_index vector coding to which group each sample
##'   belongs, with length(\code{group_index})=nrow(\code{X}). If no
##'   group index is provided a rigid grid with \code{groupsize}
##'   samples per group is used (which defaults to all samples if
##'   \code{groupsize} was not set).
##' @param partition_index vector coding to which partition each
##'   sample belongs, with
##'   length(\code{partition_index})=nrow(\code{X}). If no partition
##'   index is provided a rigid grid with \code{partitionsize} samples
##'   per partition is used.
##' @param n_components number of components to extract. If NA is
##'   passed, the same number of components as the input has
##'   dimensions is used.
##' @param n_components_uwedge number of components to extract during
##'   uwedge approximate joint diagonalization of the matrices. If NA
##'   is passed, the same number of components as the input has
##'   dimensions is used.
##' @param rank_components boolean, optional. When TRUE, the
##'   components will be ordered in decreasing stability.
##' @param pairing either 'complement', 'neighbouring' or
##'   'allpairs'. If 'allpairs' the difference matrices are computed
##'   for all pairs of partition covariance matrices, if 'complement'
##'   a one-vs-complement scheme is used and if 'neighbouring'
##'   differences with the right neighbour parition are used.
##' @param max_matrices float or 'no_partitions', optional
##'   (default=1).  The fraction of (lagged) covariance matrices to
##'   use during training or, if 'no_partitions', at most as many
##'   covariance matrices are used as there are partitions.
##' @param groupsize int, optional. Approximate number of samples in
##'   each group when using a rigid grid as groups. If NA is passed,
##'   all samples will be in one group unless group_index is passed
##'   during fitting in which case the provided group index is used
##'   (the latter is the advised and preferred way).
##' @param partitionsize int or vector of ints, optional. Approximate
##'   number of samples in each partition when using a rigid grid as
##'   partition. If NA is passed, a (hopefully sane) default is used,
##'   again, unless partition_index is passed during fitting in which
##'   case the provided partition index is used. If a vector is
##'   passed, each element is used to construct a grid and all
##'   resulting partitions are used.
##' @param timelags vector of ints, optional. Specifies which timelags
##'   should be included. 0 correpsonds to covariance matrix.
##' @param instantcov boolean, default TRUE. Specifies whether to
##'   include covariance matrix when timelags are used.
##' @param max_iter int, optional. Maximum number of iterations for
##'   the uwedge approximate joint diagonalisation during fitting.
##' @param tol float, optional. Tolerance for terminating the uwedge
##'   approximate joint diagonalisation during fitting.
##' @param minimize_loss boolean, optional. Parameter is passed to
##'   uwedge and specifies whether to compute loss function in each
##'   iteration step of uwedge.
##' @param condition_threshold float, optional. Parameter is passed to
##'   uwedge and specifies whether and at which threshold to terminate
##'   uwedge iteration depending on the condition number of the
##'   unmixing matrix.
##' @param silent boolean whether to supress status outputs.
##' 
##' @return object of class 'CoroICA' consisting of the following
##'   elements
##'
##' \item{V}{the unmixing matrix.}
##' 
##' \item{coverged}{boolean indicating whether the approximate joint
##' diagonalisation converged due to tol.}
##'
##' \item{n_iter}{number of iterations of the approximate joint
##' diagonalisation.}
##'
##' \item{meanoffdiag}{mean absolute value of the off-diagonal values
##' of the to be jointly diagonalised matrices, i.e., a proxy of the
##' approximate joint diagonalisation objective function.}
##' 
##' @export
##'
##' @import stats utils MASS
##'
##' @author Niklas Pfister and Sebastian Weichwald
##'
##' @references
##' Pfister, N., S. Weichwald, P. Bühlmann and B. Schölkopf (2018).
##' Robustifying Independent Component Analysis by Adjusting for Group-Wise Stationary Noise
##' ArXiv e-prints (arXiv:1806.01094).
##'
##' Project website (https://sweichwald.de/coroICA/)
##'
##' @seealso The function \code{\link{uwedge}} allows to perform to
##'   perform an approximate joint matrix diagonalization.
##'
##' @examples
##' ## Example
##' set.seed(1)
##' 
##' # Generate data from a block-wise variance model
##' d <- 2
##' m <- 10
##' n <- 5000
##' group_index <- rep(c(1,2), each=n)
##' partition_index <- rep(rep(1:m, each=n/m), 2)
##' S <- matrix(NA, 2*n, d)
##' H <- matrix(NA, 2*n, d)
##' for(i in unique(group_index)){
##'   varH <- abs(rnorm(d))/4
##'   H[group_index==i, ] <- matrix(rnorm(d*n)*rep(varH, each=n), n, d)
##'   for(j in unique(partition_index[group_index==i])){
##'     varS <- abs(rnorm(d))
##'     index <- partition_index==j & group_index==i
##'     S[index,] <- matrix(rnorm(d*n/m)*rep(varS, each=n/m),
##'                                                      n/m, d)
##'   }
##' }
##' A <- matrix(rnorm(d^2), d, d)
##' A <- A%*%t(A)
##' X <- t(A%*%t(S+H))
##' 
##' # Apply coroICA
##' res <- coroICA(X, group_index, partition_index, pairing="neighbouring", rank_components=TRUE)
##' 
##' # Compare results
##' par(mfrow=c(2,2))
##' plot((S+H)[,1], type="l", main="true source 1", ylab="S+H")
##' plot(res$Shat[,1], type="l", main="estimated source 1", ylab="Shat")
##' plot((S+H)[,2], type="l", main="true source 2", ylab="S+H")
##' plot(res$Shat[,2], type="l", main="estimated source 2", ylab="Shat")
##' cor(res$Shat, S+H)
                

coroICA <- function(X,
                    group_index=NA,
                    partition_index=NA,
                    n_components=NA,
                    n_components_uwedge=NA,
                    rank_components=FALSE,
                    pairing='complement',
                    max_matrices=1,
                    groupsize=1,
                    partitionsize=NA,
                    timelags=NA,
                    instantcov=TRUE,
                    max_iter=1000,
                    tol=1e-12,
                    minimize_loss=FALSE,
                    condition_threshold=NA,
                    silent=TRUE){

  d <- dim(X)[2]
  n <- dim(X)[1]

  # check timelags consistency
  if(is.na(timelags)){
    if(!instantcov){
      stop("No timelags and instantcov = FALSE. Change settings.")
    }
    else{
      timelags <- 0
    }
  }
  if(0 %in% timelags){
    if(!instantcov){
      warning("intantcov and timelags are inconsistent. Including timelag 0 in the model")
    }
  }
  else{
    if(instantcov){
      timelags <- c(0, timelags)
    }
  }
  no_timelags <- length(timelags)

  # generate group index as needed
  if(!is.numeric(group_index) & is.na(groupsize)){
    group_index <- rep(0, n)
  }
  else if(!is.numeric(group_index)){
    group_index <- rigidgroup(n, groupsize)
  }
  no_groups <- length(unique(group_index))

  # generate partiton indices as needed
  if(!is.numeric(partition_index) & is.na(partitionsize)){
    smallest_group <- min(unique(group_index, return_counts=TRUE)$counts)
    partition_indices <- list(rigidpartition(group_index,
                                           max(c(d, floor(smallest_group/2)))))
  }
  else if(!is.numeric(partition_index)){
    partition_indices <- lapply(partitionsize,
                                function(x) rigidpartition(group_index, x))
  }
  else{
    partition_indices <- list(partition_index)
  }
  
  
  if(!silent){
    print("coroICA: computing covmats")
  }
  
  # estimate covariances (depending on pairing)
  if(pairing == "complement"){
    if(max_matrices == "no_partitions"){
      max_matrices <- 1
    }
    no_pairs <- 0
    for(partition_index in partition_indices){
      for(env in unique(group_index)){
        if(length(unique(partition_index[group_index == env]))>1){
          no_pairs <- no_pairs+length(unique(partition_index[group_index == env]))
        }
      }
    }
    covmats <- vector("list", no_pairs*no_timelags)
    idx <- 1
    for(partition_index in partition_indices){
      for(env in unique(group_index)){
        unique_partitions <- unique(partition_index[group_index == env])
        unique_partitions <- sample(unique_partitions,
                                    ceiling(max_matrices*length(unique_partitions)))
        if(length(unique_partitions)==1){
          warning(paste("Removing group", toString(env),
                        "since the partition is trivial, i.e., contains only one set"))
        }
        else{
          for(subenv in unique_partitions){
            ind1 <- ((partition_index == subenv) &
                       (group_index == env))
            ind2 <- ((partition_index != subenv) &
                       (group_index == env))
            for(timelag in timelags){
              covmats[[idx]] <- autocov(X[ind1,], timelag) - autocov(X[ind2,], timelag)
              idx <- idx + 1
            }
          }
        }
      }
    }
    covmats <- covmats[1:(idx-1)]
  }
  else if(pairing == "allpairs"){
    no_pairs <-  0
    for(part_ind in 1:length(partition_indices)){
      partition_index <- partition_indices[[part_ind]]
      subvec_list[[part_ind]] <- rep(0, no_groups)
      for(i in 1:no_groups){
        env <- unique(group_index)[i]
        subvec_list[[part_ind]][i] <- length(unique(partition_index[group_index == env]))
        no_pairs <- no_pairs + subvec_list[[part_ind]][i]*(subvec_list[[part_ind]][i]-1)/2
      }
    }
    covmats <- vector("list", no_pairs*no_timelags)
    idx <- 1
    for(part_ind in 1:length(partition_indices)){
      partition_index <- partition_indices[[part_ind]]
      for(count in 1:no_groups){
        env <- unique(group_index)[count]
        unique_subs <- unique(partition_index[group_index == env])
        if(subvec_list[[part_ind]][count] == 1){
          warning(paste("Removing group", toString(env),
                        "since the partition is trivial, i.e., contains only one set"))
        }
        else{
          for(i in 1:(subvec_list[[part_ind]][count]-1)){
            for(j in (i+1):subvec_list[[part_ind]][count]){
              ind1 <- ((partition_index == unique_subs[i]) &
                         (group_index == env))
              ind2 <- ((partition_index == unique_subs[j]) &
                         (group_index == env))
              for(timelag in timelags){
                covmats[[idx]] <- autocov(X[ind1,], timelag) - autocov(X[ind2,], timelag)
                idx <- idx + 1
              }
            }
          }
        }
      }
    }
  }
  else if(pairing == "neighbouring"){
    if(max_matrices == "no_partitions"){
      max_matrices <- 1
    }
    no_pairs <- 0
    for(partition_index in partition_indices){
      for(env in unique(group_index)){
        if(length(unique(partition_index[group_index == env]))>1){
          no_pairs <- no_pairs+length(unique(partition_index[group_index == env])) - 1
        }
      }
    }
    covmats <- vector("list", no_pairs*no_timelags)
    idx <- 1
    for(partition_index in partition_indices){
      for(env in unique(group_index)){
        unique_partitions <- unique(partition_index[group_index == env])
        unique_partitions <- sort(sample(unique_partitions,
                                         ceiling(max_matrices*length(unique_partitions))))
        if(length(unique_partitions)==1){
          warning(paste("Removing group", toString(env),
                        "since the partition is trivial, i.e., contains only one set"))
        }
        else{
          for(subenv_ind in 1:(length(unique_partitions)-1)){
            subenv1 <- unique_partitions[subenv_ind]
            subenv2 <- unique_partitions[subenv_ind+1]
            ind1 <- ((partition_index == subenv1) &
                       (group_index == env))
            ind2 <- ((partition_index == subenv2) &
                       (group_index == env))
            for(timelag in timelags){
              covmats[[idx]] <- autocov(X[ind1,], timelag) - autocov(X[ind2,], timelag)
              idx <- idx + 1
            }
          }
        }
      }
    }
  }
  else{
    stop('no appropriate pairing specified')
  }
  
  # check if there are sufficiently many covariance matrices
  if(length(covmats)<=0){
    stop("Not sufficiently many covariance matrices.")
  }

  # compute total observational covariance for normalization
  Rx0 <- cov(X)

  if(!silent){
    print('coroICA: computed cov matrices')
  }

  # joint diagonalisation
  adj_res <- uwedge(covmats,
                    init=NA,
                    Rx0=Rx0,
                    return_diag=FALSE,
                    tol=tol,
                    max_iter=max_iter,
                    n_components=n_components_uwedge,
                    minimize_loss=minimize_loss,
                    condition_threshold=condition_threshold,
                    silent=silent)
  V <- adj_res$V
  
  if(!silent){
    print('coroICA: finished uwedge ajd')
  }

  # normalise V
  normaliser <- diag(V %*% Rx0 %*% t(V))
  V <- V / matrix((sign(normaliser)*sqrt(abs(normaliser))), nrow(V), d)

  # rank components
  if(rank_components | !is.na(n_components)){
    A <- ginv(V)
    colcorrs <- rep(0, dim(V)[1])
    # running average
    for(k in 1:length(covmats)){
      colcorrs <- colcorrs + diag(abs(cor(A, covmats[[k]] %*% t(V)))) / length(covmats)
    }
    sorting <- order(colcorrs, decreasing=FALSE)
    V <- V[sorting,]
  }
  if(!is.na(n_components)){
    V <- V[1:n_components,]
  }

  # source estimation
  Shat <- t(V %*% t(X))


  return(list(V=V,
              Shat=Shat,
              converged=adj_res$converged,
              n_iter=adj_res$iteration,
              meanoffdiag=adj_res$meanoffdiag))
}



rigidpartition <- function(group, nosamples){
  partition <- rep(0, length(group))
  for(e in unique(group)){
    partition[group == e] <- rigidgroup(sum(group == e),
                                               nosamples) + max(partition)
  }
  return(partition)
}


rigidgroup <- function(len, nosamples){
  groups <- floor(len / nosamples)
  changepoints <- round(seq(1, len, length.out=groups+1))
  index <- rep(0, len)
  for(i in 1:(length(changepoints)-1)){
    index[changepoints[i]:changepoints[i+1]] <- i
  }
  return(index)
}

center_rowmeans <- function(x){
  xcenter <- rowMeans(x)
  return(x - rep(xcenter, ncol(x)))
}

autocov <- function(X, lag=0){
  if(lag==0){
    return(cov(X))
  }
  else{
    n <- ncol(X)
    A <- center_rowmeans(X[, lag:n])
    B <- center_rowmeans(X[, 1:(n-lag)])
    new <- ncol(A) - 1
    return((A %*% t(B))/new)
  }
}

