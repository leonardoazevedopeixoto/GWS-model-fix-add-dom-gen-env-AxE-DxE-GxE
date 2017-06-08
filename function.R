

freq.loci <- function(X){ (sum(X==2) + sum(X==1)*.5)/length(X) }

scale.marker <- function(X,freq){
  if(X==0){
    X = -2*(1-freq)^2
  } else if(X==1){
    X = 2*freq*(1-freq)
  } else if(X==2){
    X = -2*freq^2
  }
  return(X)
}

get.all.loci <- function(X){
  freq = freq.loci(X)
  apply(as.matrix(X,ncol=1),1,scale.marker,freq)
}

scale.dom <- function(X){
  apply(X,2,get.all.loci)
}

design.matrix <- function(X, Nids){
  fix  = as.factor(X)
  lev  = unique(X)
  nfix = length(lev)
  fix1 = Matrix(0, nrow = Nids, ncol = nfix, sparse = TRUE)
  for(i in 1:nfix){
    idx = X == lev[i]
    fix1[idx,i] = 1
  }
  rownames(fix1) = 1:Nids
  return(fix1)
}

recode.genotype <- function(M){
  ##
  geno1 = M #matrix(0,Nrow,Ncol)
  ##
  idx = M == 1
  geno1[idx] = 3
  ##
  idx = M == 0
  geno1[idx] = 1
  ##
  colnames(geno1) = colnames(M)
  rownames(geno1) = rownames(M)
  return(geno1)
}


