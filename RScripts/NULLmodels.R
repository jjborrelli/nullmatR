library(igraph)
library(bipartite)




###
### CURVE-BALL ALGORITHM
###
# Strona, G. et al. 2014. A fast and unbiased procedure to randomize 
# ecological binary matrices with fixed row and column totals.
# -Nat. Comm. 5: 4114. [doi: 10.1038/ncomms5114]
# (http://www.nature.com/ncomms/2014/140611/ncomms5114/full/ncomms5114.html)

curve_ball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}

curving <- function(adjmat, n){
  newmat <- lapply(1:n, function(x) matrix(0, nrow(adjmat), ncol(adjmat)))
  newmat[[1]] <- adjmat
  
  for(i in 2:n){
    newmat[[i]] <- curve_ball(newmat[[i-1]])
  }
  return(newmat)
}




testcurve <- function(N, C, n){
  rg <- erdos.renyi.game(N, C, "gnp")
  rm <- get.adjacency(rg, sparse = F)
  
  cb <- curving(rm, n)
  csc <- sapply(cb, C.score)
  p1 <- sum(csc <= csc[[1]])/n
  p2 <- sum(csc >= csc[[1]])/n
  
  return(c(p1, p2))
}


tc1 <- sapply(1:100, function(x) testcurve(20, .1, 100000))
tc2 <- sapply(1:100, function(x) testcurve(20, .1, 100000))
tc3 <- sapply(1:100, function(x) testcurve(20, .1, 100000))
tc4 <- sapply(1:100, function(x) testcurve(20, .1, 100000))
tc5 <- sapply(1:100, function(x) testcurve(20, .1, 100000))