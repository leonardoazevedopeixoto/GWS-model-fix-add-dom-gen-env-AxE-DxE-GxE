
pred <- function(args, model){
  add = args$add
  dom = args$dom
  gen = args$gen
  env = args$env
  if (model =1) {if (add & env) {prediction= matrix(envtest%*%env_effect[,k] + addtest%*%add_effect[,k])}
    if (dom & env) {prediction= matrix(envtest%*%env_effect[,k] + domtest%*%dom_effect[,k])}
    if (gen & env) {prediction= matrix(envtest%*%env_effect[,k] + gentest%*%gen_effect[,k])}
    if (add & env) {prediction= matrix(envtest%*%env_effect[,k] + addtest%*%add_effect[,k] + AxEtest%*%AxE_effect[,k])}
    if (add & env) {prediction= matrix(envtest%*%env_effect[,k] + domtest%*%dom_effect[,k] + DxEtest%*%DxE_effect[,k])}
    if (add & env) {prediction= matrix(envtest%*%env_effect[,k] + gentest%*%gen_effect[,k] + GxEtest%*%GxE_effect[,k])}
    if (add & gen & env) {prediction= matrix(envtest%*%env_effect[,k] + addtest%*%add_effect[,k] + AxEtest%*%AxE_effect[,k] +
                                               gentest%*%gen_effect[,k] + GxEtest%*%GxE_effect[,k])}
    if (add & don & env) {prediction= matrix(envtest%*%env_effect[,k] + addtest%*%add_effect[,k] + AxEtest%*%AxE_effect[,k] +
                                               domtest%*%dom_effect[,k] + DxEtest%*%DxE_effect[,k])}
    if (gen & don & env) {prediction= matrix(envtest%*%env_effect[,k] + gentest%*%gen_effect[,k] + GxEtest%*%GxE_effect[,k] +
                                               domtest%*%dom_effect[,k] + DxEtest%*%DxE_effect[,k])}
    if (add & gen & don & env) {prediction= matrix(envtest%*%env_effect[,k] + addtest%*%add_effect[,k] + AxEtest%*%AxE_effect[,k] +
                                                     gentest%*%gen_effect[,k] + GxEtest%*%GxE_effect[,k] +
                                                     domtest%*%dom_effect[,k] + DxEtest%*%DxE_effect[,k])}
  }
  else if (model =2) {if (add) {prediction= matrix(addtest%*%add_effect[,k])}
    if (dom) {prediction= matrix(domtest%*%dom_effect[,k])}
    if (gen) {prediction= matrix(gentest%*%gen_effect[,k])}
    if (add & env) {prediction= matrix(addtest%*%add_effect[,k] + AxEtest%*%AxE_effect[,k])}
    if (add & env) {prediction= matrix(domtest%*%dom_effect[,k] + DxEtest%*%DxE_effect[,k])}
    if (add & env) {prediction= matrix(gentest%*%gen_effect[,k] + GxEtest%*%GxE_effect[,k])}
    if (add & gen & env) {prediction= matrix(addtest%*%add_effect[,k] + AxEtest%*%AxE_effect[,k] +
                                               gentest%*%gen_effect[,k] + GxEtest%*%GxE_effect[,k])}
    if (add & don & env) {prediction= matrix(addtest%*%add_effect[,k] + AxEtest%*%AxE_effect[,k] +
                                               domtest%*%dom_effect[,k] + DxEtest%*%DxE_effect[,k])}
    if (gen & don & env) {prediction= matrix(gentest%*%gen_effect[,k] + GxEtest%*%GxE_effect[,k] +
                                               domtest%*%dom_effect[,k] + DxEtest%*%DxE_effect[,k])}
    if (add & gen & don & env) {prediction= matrix(addtest%*%add_effect[,k] + AxEtest%*%AxE_effect[,k] +
                                                     gentest%*%gen_effect[,k] + GxEtest%*%GxE_effect[,k] +
                                                     domtest%*%dom_effect[,k] + DxEtest%*%DxE_effect[,k])}
  }
  else if (model =3) {if (add) {prediction= matrix(addtest%*%add_effect[,k])}
    if (dom) {prediction= matrix(domtest%*%dom_effect[,k])}
    if (gen) {prediction= matrix(gentest%*%gen_effect[,k])}
    if (add & gen) {prediction= matrix(addtest%*%add_effect[,k] + gentest%*%gen_effect[,k])}
    if (add & don) {prediction= matrix(addtest%*%add_effect[,k] + domtest%*%dom_effect[,k])}
    if (gen & don) {prediction= matrix(gentest%*%gen_effect[,k] + domtest%*%dom_effect[,k])}
    if (add & gen & don) {prediction= matrix(addtest%*%add_effect[,k] + gentest%*%gen_effect[,k] + 
                                                     domtest%*%dom_effect[,k])}
  }
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

family.pred <- function (args){
  add = args$add
  dom = args$dom
  
  genotypic_value           <- as.data.frame( X[,1:3] )
  colnames(genotypic_value) <- c( "Env", "Gen", "Blo" )
  
  if ( add ){
    genotypic_value$add <- effect.pred(X = args$Z, a = args$k)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Additive_Value")
  }
  
  if ( dom ){
    genotypic_value$dom <- effect.pred( X = args$W, a = args$r)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Dominance_Value")
  }
  
  if( add & dom ){
    genotypic_value$GV <- genotypic_value$Additive_Value + genotypic_value$Dominance_Value
  }
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
      if (genotypic_value1[y,2]==x & add){fgv[y,1]=genotypic_value1[y,4]}
      if (genotypic_value1[y,2]==x & dom){fgv[y,1]=genotypic_value1[y,4]}
      if (genotypic_value1[y,2]==x & add & dom){fgv[y,1]=genotypic_value1[y,6]}
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


cal.inter <- function(Z, W){
  if (W==1) {X[,(1+ncol(Z)*which(W==1)):(ncol(Z)+ncol(Z)*which(W==1))]=Z}
}

design.interaction <- function(Z,W){
  X<-Matrix(0,nrow=nrow(Z), ncol=(ncol(Z)*ncol(W)+ncol(Z)), sparse= TRUE)
  mapply(cal.inter, split(X, row(X)), split(Z, row(Z)), split(W, row(W)))
  inter<=X[,(ncol(Z)+1):ncol(X)]
  return(inter)
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


