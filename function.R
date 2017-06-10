

model <- function(model){
  if (model=="FIXED"){model=model}
}

effect.pred <- function(X, a){
  effect<-matrix(X%*%rowMeans(a))
  return(effect)
}

ind.pred <- function (X, G, A, D, GE, AE, DE){
  genetic_GxE_value<-data.frame(pheno[,1:3],matrix((gen1%*%rowMeans(genetic_effects)+add%*%rowMeans(additive_effects)+
                                                      dom%*%rowMeans(dominance_effects)+GxE%*%rowMeans(GeneticxEnvironment_effects)+
                                                      AxE%*%rowMeans(AdditivexEnvironment_effects)+DxE%*%rowMeans(DominancexEnvironment_effects))))
  colnames(genetic_GxE_value)<-c("Env", "Gen", "Blo", "GV")
  genetic_GxE_value1<-genetic_GxE_value[order(genetic_GxE_value$Env),]
}

family.pred <- function (X, Z, W, k, r, add, dom){
  if (add = TRUE & dom = FALSE){genotypic_value<-data.frame(X[,1:3],effect.pred(X=A, a=k))
  colnames(genotypic_value)<-c("Env", "Gen", "Blo", "GV")}
  else if (add = TRUE & dom = TRUE){genotypic_value<-data.frame(X[,1:3],effect.pred(X=A, a=k), effect.pred(X=W, a=r))
  colnames(genotypic_value)<-c("Env", "Gen", "Blo", "GV")}
  else if (add = FALSE & dom = TRUE){genotypic_value<-data.frame(X[,1:3],effect.pred(X=W, a=r))
  colnames(genotypic_value)<-c("Env", "Gen", "Blo", "GV")}
  genotypic_value1<-genotypic_value[order(genotypic_value$Env),]
  lev  = unique(X[,2])
  NFam = length(lev)
  lev  = unique(X[,1])
  Nemv = length(lev)
  lev  = unique(X[,3])
  Nbloc = length(lev)
  
  family_genotypic_value<-matrix(nrow=ncol(NFam), ncol=1)
  for (x in 1:ncol(Nfam))
  {
    fgv<-matrix(0,nrow=nrow(X), ncol=1)
    for (y in 1:nrow(X))
    {
      if (genotypic_value1[y,2]==x){fgv[y,1]=genotypic_value1[y,4]}
    }
    family_genotypic_value[x,1]<-sum(fgv)/ncol(Nenv)*ncol(Nbloc)
  }
  rownames(family_genotypic_value)<-seq(1, max(X[,2]), 1)
  return(family_genotypic_value)
}

design.interaction <- function(M,X){
  inter<-Matrix(0, nrow=Nids, ncol=(ncol(M)*ncol(X)+ncol(M)), sparse= TRUE)
  for (i in 1:ncol(X)){
    for (j in 1:Nids){
      if(X[j,i]==1){inter[j,(i*ncol(M)+1):((i+1)*ncol(M))]=M[j,1:ncol(M)]}
    }
  }
  internew=inter[,(ncol(M)+1):ncol(inter)]
  return(internew)
}

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
  lev  = unique(X)
  nlev = length(lev)
  incidence = Matrix(0, nrow = Nids, ncol = nlev, sparse = TRUE)
  for(i in 1:nlev){
    idx = X == lev[i]
    incidence[idx,i] = 1
  }
  rownames(incidence) = 1:Nids
  return(incidence)
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


