
ETA.pred <-function (args, model){
  if(Fix & Nfix>1) {ETA=list(list(X=args$fixtrain, model = model[1,]))}
  if(Env & Nenv>1) {ETA=list(list(X=args$envtrain, model = model[2,]))}
  if(Gen) {ETA=list(list(X=args$gentrain, model = model[3,]))} 
  if(Add) {ETA=list(list(X=args$addtrain, model = model[4,]))}
  if(Dom) {ETA=list(list(X=args$domtrain, model = model[5,]))}
  if(Fix & Nfix>1 & Env & Nenv>1) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                            list(X=args$envtrain, model = model[2,]))}
  if(Fix & Nfix>1 & Gen) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                   list(X=args$gentrain, model = model[3,]))}
  if(Fix & Nfix>1 & Add) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                   list(X=args$addtrain, model = model[4,]))}
  if(Fix & Nfix>1 & Dom) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                   list(X=args$domtrain, model = model[5,]))}
  if(Env & Nenv>1 & Gen) {ETA=list(list(X=args$envtrain, model = model[2,]),
                                   list(X=args$gentrain, model = model[3,]),
                                   list(X=args$GxEtrain, model = model[6,]))}
  if(Env & Nenv>1 & Add) {ETA=list(list(X=args$envtrain, model = model[2,]),
                                   list(X=args$addtrain, model = model[4,]),
                                   list(X=args$AxEtrain, model = model[7,]))}
  if(Env & Nenv>1 & Dom) {ETA=list(list(X=args$envtrain, model = model[2,]),
                                   list(X=args$domtrain, model = model[5,]),
                                   list(X=args$DxEtrain, model = model[8,]))}
  if(Gen & Add) {ETA=list(list(X=args$gentrain, model = model[3,]),
                          list(X=args$addtrain, model = model[4,]))}
  if(Gen & Dom) {ETA=list(list(X=args$gentrain, model = model[3,]),
                          list(X=args$domtrain, model = model[5,]))}
  if(Dom & Add) {ETA=list(list(X=args$domtrain, model = model[5,]),
                          list(X=args$addtrain, model = model[4,]))} 
  if(Fix & Nfix>1 & Env & Nenv>1 & Gen) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                                  list(X=args$envtrain, model = model[2,]),
                                                  list(X=args$gentrain, model = model[3,]),
                                                  list(X=args$GxEtrain, model = model[6,]))}
  if(Fix & Nfix>1 & Env & Nenv>1 & Add) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                                  list(X=args$envtrain, model = model[2,]),
                                                  list(X=args$addtrain, model = model[4,]),
                                                  list(X=args$AxEtrain, model = model[7,]))}
  if(Fix & Nfix>1 & Env & Nenv>1 & Dom) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                                  list(X=args$envtrain, model = model[2,]),
                                                  list(X=args$domtrain, model = model[5,]),
                                                  list(X=args$DxEtrain, model = model[8,]))}
  if(Fix & Nfix>1 & Gen & Add) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                         list(X=args$gentrain, model = model[3,]),
                                         list(X=args$addtrain, model = model[4,]))}
  if(Fix & Nfix>1 & Gen & Dom) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                         list(X=args$gentrain, model = model[3,]),
                                         list(X=args$domtrain, model = model[5,]))}
  if(Fix & Nfix>1 & Add & Dom) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                         list(X=args$addtrain, model = model[4,]),
                                         list(X=args$domtrain, model = model[5,]))}
  if(Env & Nenv>1 & Gen & Add) {ETA=list(list(X=args$envtrain, model = model[2,]),
                                         list(X=args$gentrain, model = model[3,]),
                                         list(X=args$GxEtrain, model = model[6,]),
                                         list(X=args$addtrain, model = model[4,]),
                                         list(X=args$AxEtrain, model = model[7,]))}
  if(Env & Nenv>1 & Gen & Dom) {ETA=list(list(X=args$envtrain, model = model[2,]),
                                         list(X=args$gentrain, model = model[3,]),
                                         list(X=args$GxEtrain, model = model[6,]),
                                         list(X=args$domtrain, model = model[5,]),
                                         list(X=args$DxEtrain, model = model[8,]))}
  if(Env & Nenv>1 & Add & Dom) {ETA=list(list(X=args$envtrain, model = model[2,]),
                                         list(X=args$addtrain, model = model[4,]),
                                         list(X=args$AxEtrain, model = model[7,]),
                                         list(X=args$domtrain, model = model[5,]),
                                         list(X=args$DxEtrain, model = model[8,]))}
  if(Gen & Add & Dom) {ETA=list(list(X=args$gentrain, model = model[3,]),
                                list(X=args$addtrain, model = model[4,]),
                                list(X=args$domtrain, model = model[5,]))}
  if(Fix & Nfix>1 & Env & Nenv>1 & Gen & Add) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                                        list(X=args$envtrain, model = model[2,]),
                                                        list(X=args$gentrain, model = model[3,]),
                                                        list(X=args$GxEtrain, model = model[6,]),
                                                        list(X=args$addtrain, model = model[4,]),
                                                        list(X=args$AxEtrain, model = model[7,]))}
  if(Fix & Nfix>1 & Env & Nenv>1 & Gen & Dom) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                                        list(X=args$envtrain, model = model[2,]),
                                                        list(X=args$gentrain, model = model[3,]),
                                                        list(X=args$GxEtrain, model = model[6,]),
                                                        list(X=args$domtrain, model = model[5,]),
                                                        list(X=args$DxEtrain, model = model[8,]))}
  if(Fix & Nfix>1 & Env & Nenv>1 & Add & Dom) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                                        list(X=args$envtrain, model = model[2,]),
                                                        list(X=args$addtrain, model = model[4,]),
                                                        list(X=args$AxEtrain, model = model[7,]),
                                                        list(X=args$domtrain, model = model[5,]),
                                                        list(X=args$DxEtrain, model = model[8,]))}
  if(Fix & Nfix>1 & Gen & Add & Dom) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                               list(X=args$gentrain, model = model[3,]),
                                               list(X=args$addtrain, model = model[4,]),
                                               list(X=args$domtrain, model = model[5,]))}
  if(Env & Nenv>1 & Gen & Add & Dom) {ETA=list(list(X=args$envtrain, model = model[2,]),
                                               list(X=args$gentrain, model = model[3,]),
                                               list(X=args$GxEtrain, model = model[6,]),
                                               list(X=args$addtrain, model = model[4,]),
                                               list(X=args$AxEtrain, model = model[7,]),
                                               list(X=args$domtrain, model = model[5,]),
                                               list(X=args$DxEtrain, model = model[8,]))}
  if(Fix & Nfix>1 & Env & Nenv>1 & Gen & Add & Dom) {ETA=list(list(X=args$fixtrain, model = model[1,]),
                                                              list(X=args$envtrain, model = model[2,]),
                                                              list(X=args$gentrain, model = model[3,]),
                                                              list(X=args$GxEtrain, model = model[6,]),
                                                              list(X=args$addtrain, model = model[4,]),
                                                              list(X=args$AxEtrain, model = model[7,]),
                                                              list(X=args$domtrain, model = model[5,]),
                                                              list(X=args$DxEtrain, model = model[8,]))}
  return(ETA)
}

matrix.creation <-function(X, args, train){
  X=as.matrix(X)
  rownames(X)<-seq(1,nrow(X),1)
  if (args$fam==FALSE){
    newM=Matrix(X[train,])
  }
  else if (args$fam==TRUE){
    new=data.frame(agrs$pheno[,1:3], X)
    new1=sort(new[,3], decreasing=F)
    newM=Matrix(new1[,4:ncol(new1)])
  }
  return(newM)
}

pred <- function(args, mod, k){
  add = args$add
  dom = args$dom
  gen = args$gen
  env = args$env
  if (mod==1) {if (add & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$addtest%*%args$add_effect[,k])}
    if (dom & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$domtest%*%args$dom_effect[,k])}
    if (gen & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$gentest%*%args$gen_effect[,k])}
    if (add & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$addtest%*%args$add_effect[,k] + args$AxEtest%*%args$AxE_effect[,k])}
    if (add & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$domtest%*%args$dom_effect[,k] + args$DxEtest%*%args$DxE_effect[,k])}
    if (add & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$gentest%*%args$gen_effect[,k] + args$GxEtest%*%args$GxE_effect[,k])}
    if (add & gen & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$addtest%*%args$add_effect[,k] + args$AxEtest%*%args$AxE_effect[,k] +
                                                        args$gentest%*%args$gen_effect[,k] + args$GxEtest%*%args$GxE_effect[,k])}
    if (add & dom & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$addtest%*%args$add_effect[,k] + args$AxEtest%*%args$AxE_effect[,k] +
                                                        args$domtest%*%args$dom_effect[,k] + args$DxEtest%*%args$DxE_effect[,k])}
    if (gen & dom & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$gentest%*%args$gen_effect[,k] + args$GxEtest%*%args$GxE_effect[,k] +
                                                        args$domtest%*%args$dom_effect[,k] + args$DxEtest%*%args$DxE_effect[,k])}
    if (add & gen & dom & env & Nenv>1) {prediction= matrix(args$envtest%*%args$env_effect[,k] + args$addtest%*%args$add_effect[,k] + args$AxEtest%*%args$AxE_effect[,k] +
                                                              args$gentest%*%args$gen_effect[,k] + args$GxEtest%*%args$GxE_effect[,k] +
                                                              args$domtest%*%args$dom_effect[,k] + args$DxEtest%*%args$DxE_effect[,k])}
  }
  if (mod==2) {if (add & env & Nenv>1) {prediction= matrix(args$addtest%*%args$add_effect[,k] + args$AxEtest%*%args$AxE_effect[,k])}
    if (add & env & Nenv>1) {prediction= matrix(args$domtest%*%args$dom_effect[,k] + args$DxEtest%*%args$DxE_effect[,k])}
    if (add & env & Nenv>1) {prediction= matrix(args$gentest%*%args$gen_effect[,k] + args$GxEtest%*%args$GxE_effect[,k])}
    if (add & gen & env & Nenv>1) {prediction= matrix(args$addtest%*%args$add_effect[,k] + args$AxEtest%*%args$AxE_effect[,k] +
                                                        args$gentest%*%args$gen_effect[,k] + args$GxEtest%*%args$GxE_effect[,k])}
    if (add & dom & env & Nenv>1) {prediction= matrix(args$addtest%*%args$add_effect[,k] + args$AxEtest%*%args$AxE_effect[,k] +
                                                        args$domtest%*%args$dom_effect[,k] + args$DxEtest%*%args$DxE_effect[,k])}
    if (gen & dom & env & Nenv>1) {prediction= matrix(args$gentest%*%args$gen_effect[,k] + args$GxEtest%*%args$GxE_effect[,k] +
                                                        args$domtest%*%args$dom_effect[,k] + args$DxEtest%*%args$DxE_effect[,k])}
    if (add & gen & dom & env & Nenv>1) {prediction= matrix(args$addtest%*%args$add_effect[,k] + args$AxEtest%*%args$AxE_effect[,k] +
                                                              args$gentest%*%args$gen_effect[,k] + args$GxEtest%*%args$GxE_effect[,k] +
                                                              args$domtest%*%args$dom_effect[,k] + args$DxEtest%*%args$DxE_effect[,k])}
  }
  if (mod==3) {if (add) {prediction= matrix(args$addtest%*%args$add_effect[,k])}
    if (dom) {prediction= matrix(args$domtest%*%args$dom_effect[,k])}
    if (gen) {prediction= matrix(args$gentest%*%args$gen_effect[,k])}
    if (add & gen) {prediction= matrix(args$addtest%*%args$add_effect[,k] + args$gentest%*%args$gen_effect[,k])}
    if (add & dom) {prediction= matrix(args$addtest%*%args$add_effect[,k] + args$domtest%*%args$dom_effect[,k])}
    if (gen & dom) {prediction= matrix(args$gentest%*%args$gen_effect[,k] + args$domtest%*%args$dom_effect[,k])}
    if (add & gen & dom) {prediction= matrix(args$addtest%*%args$add_effect[,k] + args$gentest%*%args$gen_effect[,k] + 
                                               args$domtest%*%args$dom_effect[,k])}
  } else {prediction=NULL}
  return(prediction)
}

effect.pred <- function(Z, a){
  effect<-matrix(Z%*%rowMeans(a, na.rm = TRUE))
  return(effect)
}

ind.pred <- function (args){
  add = args$add
  dom = args$dom
  gen = args$gen
  env = args$env
  genotypic_value           <- as.data.frame( X[,1:3] )
  colnames(genotypic_value) <- c( "Env", "Gen", "Blo" )
  if ( gen ){
    genotypic_value$gen <- effect.pred(Z = args$Gen, a = args$genetic_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Genetic_Value")
  }
  if ( add ){
    genotypic_value$add <- effect.pred(Z = args$Add, a = args$additive_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Additive_Value")
  }
  if ( dom ){
    genotypic_value$dom <- effect.pred(Z = args$Dom, a = args$dominance_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Dominance_Value")
  }
  if( gen & add ){
    genotypic_value$GV <- genotypic_value$Genetic_Value + genotypic_value$Additive_Value
  }
  if( gen & dom ){
    genotypic_value$GV <- genotypic_value$Genetic_Value + genotypic_value$Dominance_Value
  }
  if( add & dom ){
    genotypic_value$GV <- genotypic_value$Additive_Value + genotypic_value$Dominance_Value
  }
  if( gen & add & dom ){
    genotypic_value$GV <- genotypic_value$Genetic_Value + genotypic_value$Additive_Value 
    + genotypic_value$Dominance_Value
  }
  if ( gen & env & Nenv>1 ){
    genotypic_value$GV <- effect.pred(Z = args$Gen, a = args$genetic_effects) + effect.pred(Z = args$GxE, a = args$GeneticxEnvironment_effects)
  }
  if ( add & env & Nenv>1 ){
    genotypic_value$GV <- effect.pred(Z = args$Add, a = args$additive_effects) + effect.pred(Z = args$AxE, a = args$AdditivexEnvironment_effects)
  }
  if ( dom & env & Nenv>1 ){
    genotypic_value$GV <- effect.pred(Z = args$Dom, a = args$dominance_effects) + effect.pred(Z = args$DxE, a = args$DominancexEnvironment_effects)
  }
  if( gen & add & env & Nenv>1 ){
    genotypic_value$GV <- genotypic_value$Genetic_Value + genotypic_value$Additive_Value +
      + effect.pred(Z = args$GxE, a = args$GeneticxEnvironment_effects)
    + effect.pred(Z = args$AxE, a = args$AdditivexEnvironment_effects)
  }
  if( gen & dom & env & Nenv>1 ){
    genotypic_value$GV <- genotypic_value$Genetic_Value + genotypic_value$Dominance_Value
    + effect.pred(Z = args$GxE, a = args$GeneticxEnvironment_effects)
    + effect.pred(Z = args$DxE, a = args$DominancexEnvironment_effects)
  }
  if( add & dom & env & Nenv>1 ){
    genotypic_value$GV <- genotypic_value$Additive_Value + genotypic_value$Dominance_Value
    + effect.pred(Z = args$AxE, a = args$AdditivexEnvironment_effects)
    + effect.pred(Z = args$DxE, a = args$DominancexEnvironment_effects)
  }
  if( gen & add & dom ){
    genotypic_value$GV <- genotypic_value$Genetic_Value + genotypic_value$Additive_Value 
    + genotypic_value$Dominance_Value + effect.pred(Z = args$GxE, a = args$GeneticxEnvironment_effects)
    + effect.pred(Z = args$AxE, a = args$AdditivexEnvironment_effects)
    + effect.pred(Z = args$DxE, a = args$DominancexEnvironment_effects)
  }
  genotypic_value1<-genotypic_value[order(genotypic_value$Env),]
  return(genotypic_value1)
}

family.pred <- function (args){
  add = args$add
  dom = args$dom
  
  genotypic_value           <- as.data.frame( X[,1:3] )
  colnames(genotypic_value) <- c( "Env", "Gen", "Blo" )
  
  if ( add ){
    genotypic_value$add <- effect.pred(Z = args$Add, a = args$additive_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Additive_Value")
  }
  
  if ( dom ){
    genotypic_value$dom <- effect.pred(Z = args$Dom, a = args$dominance_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Dominance_Value")
  }
  
  if( add & dom ){
    genotypic_value$GV <- genotypic_value$Additive_Value + genotypic_value$Dominance_Value
  }
  genotypic_value1<-genotypic_value[order(genotypic_value$Env),]
  lev  = unique(X[,2])
  NFam = length(lev)
  lev  = unique(X[,1])
  Nenv = length(lev)
  lev  = unique(X[,3])
  Nbloc = length(lev)
  
  family_genotypic_value<-matrix(nrow=NFam, ncol=1)
  for (x in 1:NFam)
  {
    fgv<-matrix(0,nrow=nrow(X), ncol=1)
    for (y in 1:nrow(X))
    {
      if (genotypic_value1[y,2]==x & add){fgv[y,1]=genotypic_value1[y,4]}
      if (genotypic_value1[y,2]==x & dom){fgv[y,1]=genotypic_value1[y,4]}
      if (genotypic_value1[y,2]==x & add & dom){fgv[y,1]=genotypic_value1[y,6]}
    }
    family_genotypic_value[x,1]<-sum(fgv)/Nenv*Nbloc
  }
  rownames(family_genotypic_value)<-seq(1, max(X[,2]), 1)
  family_genotypic_value1<-data.frame(rownames(family_genotypic_value), family_genotypic_value)
  family_genotypic_value2<-family_genotypic_value1[order(family_genotypic_value1[,2]),]
  colnames(family_genotypic_value2)<-c("ID", "Genotypic value")
  return(family_genotypic_value)
}

## Function to get enviroment level
get.index <- function(X){ which(X==1) }

## Function to compute desired matrix (usually GxE, or AxE, or DxE incidence matrix)
get.matrix <- function( m1, m2, nsnp, nRow, nLev){
  ## Output matrix
  m <- Matrix(data = 0, nrow = nRow, ncol = nsnp*nLev , sparse = TRUE)
  
  ## Arguments
  Index <- apply(m2 , 1 , get.index)
  top   <- Index*nsnp
  bot   <- ifelse(top==nsnp, 1, top-nsnp+1)
  
  ## Creates the matrix in parallel
  system.time(
    foreach(ii = 1:nRow) %dopar% { m[ii,bot[ii]:top[ii]] <- m1[ii,] }
  )
  return(m)
}

##Function to estimate the frequency for each loci
freq.loci <- function(X){ (sum(X==2) + sum(X==1)*.5)/length(X) }

##Function to create the dominance incidence matrix
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

##Function to create the dominance incidence matrix
get.all.loci <- function(X){
  freq = freq.loci(X)
  apply(as.matrix(X,ncol=1),1,scale.marker,freq)
}

##Function to create the dominance incidence matrix
scale.dom <- function(X){
  apply(X,2,get.all.loci)
}

##Function to create the incidence matrix
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

## Function to recode the genotype file aiming to calculate maf, call rate and EWH using the function maf
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
