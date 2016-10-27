#' @title Permute a matrix
#' @description This function generates a null distribution of matrices by maintaining row sums
#' @param mat a binary matrix
#' @param iter the number of null matrices to create. Defaults to 100
#' @export
#' @family null models

permutes_r <- function(mat, iter = 100){
  
  pattern1 <- matrix(c(0,1), nrow = 1, ncol = 2)
  pattern2 <- matrix(c(1,0), nrow = 1, ncol = 2)
  count <- 0
  
  mat.list <- list()
  
  while(count < iter){
    srow <- sample(1:nrow(mat), 1)
    scol <- sample(1:ncol(mat), 2)
    
    test <- mat[srow, scol]
    
    if(sum(test == pattern1) == 2){
      count <- count + 1
      mat[srow, scol] <- pattern2
      mat.list[[count]] <- mat
      
      next
    } else if(sum(test == pattern2) == 2){
      count <- count + 1
      mat[srow, scol] <- pattern1
      mat.list[[count]] <- mat
      
      next
    } else {next}
  }
  
  matrices <- lapply(mat.list, as.matrix)
  return(permuted.matrices = matrices)
}


#' @title permutes_c
#' @description This function generates a null distribution of matrices by maintaining column sums
#' @param mat a binary matrix
#' @param iter the number of null matrices to create. Defaults to 100
#' @export
#' @family null models

permutes_c <- function(mat, iter = 100){
  
  pattern1 <- matrix(c(0,1), nrow = 2, ncol = 1)
  pattern2 <- matrix(c(1,0), nrow = 2, ncol = 1)
  count <- 0
  
  mat.list <- list()
  
  while(count < iter){
    srow <- sample(1:nrow(mat), 2)
    scol <- sample(1:ncol(mat), 1)
    
    test <- mat[srow, scol]
    
    if(sum(test == pattern1) == 2){
      count <- count + 1
      mat[srow, scol] <- pattern2
      mat.list[[count]] <- mat
      
      next
    } else if(sum(test == pattern2) == 2){
      count <- count + 1
      mat[srow, scol] <- pattern1
      mat.list[[count]] <- mat
      
      next
    } else {next}
  }
  
  matrices <- lapply(mat.list, as.matrix)
  return(permuted.matrices = matrices)
}


#' @title permutes_rc
#' @description This function generates a null distribution of matrices by maintaining row and column sums
#' @param mat a binary matrix
#' @param iter the number of null matrices to create. Defaults to 100
#' @export
#' @family null models

permutes_rc <- function(mat, iter = 100){
 
  pattern1 <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)
  pattern2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)
  count <- 0
  
  mat.list <- list()
  
  while(count < iter){
    srow <- sample(1:nrow(mat), 2)
    scol <- sample(1:ncol(mat), 2)
    
    test <- mat[srow, scol]
    
    if(sum(test == pattern1) == 4){
      count <- count + 1
      mat[srow, scol] <- pattern2
      mat.list[[count]] <- mat
      
      next
    } else if(sum(test == pattern2) == 4){
      count <- count + 1
      mat[srow, scol] <- pattern1
      mat.list[[count]] <- mat
      
      next
    } else {next}
  }
  
  matrices <- lapply(mat.list, as.matrix)
  return(permuted.matrices = matrices)
}


#' @title permutes_rcd
#' @description This function generates a null distribution of matrices by maintaining row and column sums as well as the diagonal
#' @param mat a binary matrix
#' @param iter the number of null matrices to create. Defaults to 100
#' @export
#' @family null models

permutes_rcd <- function(mat, iter = 100){
  
  pattern1 <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)
  pattern2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)
  count <- 0
  
  mat.list <- list()
  
  while(count < iter){
    
    s <- sample(1:nrow(mat), 4)
      
    test <- mat[s[1:2], s[3:4]]
    
    if(sum(test == pattern1) == 4){
      count <- count + 1
      mat[s[1:2], s[3:4]] <- pattern2
      mat.list[[count]] <- mat
      
      next
    } else if(sum(test == pattern2) == 4){
      count <- count + 1
      mat[s[1:2], s[3:4]] <- pattern1
      mat.list[[count]] <- mat
      
      next
    } else {next}
  }
  
  matrices <- lapply(mat.list, as.matrix)
  return(permuted.matrices = matrices)
}


#' @title permutes_prob
#' @description This function generates a null distribution of matrices by maintaining row and column probabilities
#' @param mat a binary matrix
#' @param iter the number of null matrices to create. Defaults to 100
#' @export
#' @family null models

permutes_prob <- function(mat, iter = 100){
  cS <- colSums(mat)/nrow(mat)
  rS <- rowSums(mat)/ncol(mat)
  
  pmat <- matrix(nrow =nrow(mat), ncol = ncol(mat))
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      pmat[i,j] <- sum((rS[i]+cS[j])/2)
    }
  }
  
  mat.list <- list()
  for(q in 1:iter){
    mat2 <- matrix(nrow =nrow(mat), ncol = ncol(mat))
    for(i in 1:nrow(mat)){
      for(j in 1:ncol(mat)){
        mat2[i,j] <- rbinom(1,1,prob = pmat[i,j])
      }
    }
    mat.list[[q]] <- mat2
  }
  
  return(permuted.matrices = mat.list)
}


#' @title make_null
#' @description This function generates a null distribution of matrices using one of 5 methods
#' @param mat a binary matrix
#' @param iter the number of null matrices to create. Defaults to 100
#' @param method which null model you want to use. Defaults to method = "rc
#' @export

make_null <- function(mat, iter = 100, method = "rc"){
  if(method == "r"){p <- permutes_r(mat, iter)}
  
  if(method == "c"){p <- permutes_c(mat, iter)}
  
  if(method == "rc"){p <- permutes_rc(mat, iter)}
  
  if(method == "rcd"){p <- permutes_rcd(mat, iter)}
  
  if(method == "rcp"){p <- permutes_prob(mat, iter)}
  
  return(p)
}