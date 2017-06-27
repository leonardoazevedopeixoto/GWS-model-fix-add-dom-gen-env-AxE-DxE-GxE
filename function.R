
matrix.creation <-function(X, args){
  if (args$Fam=FALSE){
    newM=X[as.character(args$train),]
    return(newM)
  }
  if (args$Fam=TRUE){
    newM=data.frame(agrs$pheno[,1:3], X)
    newM1=sort(newM[,3], decreasing=F)
    newM2=newM1[,4:ncol(newM1)]
    return(newM2)
  }
}

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
  else if (model =2) {if (add & env) {prediction= matrix(addtest%*%add_effect[,k] + AxEtest%*%AxE_effect[,k])}
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

ind.pred <- function (args){
  add = args$add
  dom = args$dom
  gen = args$gen
  env = args$env
  genotypic_value           <- as.data.frame( X[,1:3] )
  colnames(genotypic_value) <- c( "Env", "Gen", "Blo" )
  if ( gen ){
    genotypic_value$gen <- effect.pred(X = args$gen, a = args$genetic_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Genetic_Value")
  }
  if ( add ){
    genotypic_value$add <- effect.pred(X = args$add, a = args$additive_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Additive_Value")
  }
  if ( dom ){
    genotypic_value$dom <- effect.pred( X = args$dom, a = args$dominance_effects)
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
  if ( gen & env ){
    genotypic_value$gen <- effect.pred(X = args$gen, a = args$genetic_effects) + effect.pred(X = args$GxE, a = args$GeneticxEnvironment_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Genetic_Value")
  }
  if ( add & env ){
    genotypic_value$add <- effect.pred(X = args$add, a = args$additive_effects) + effect.pred(X = args$AxE, a = args$AdditivexEnvironment_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Additive_Value")
  }
  if ( dom & env ){
    genotypic_value$dom <- effect.pred( X = args$dom, a = args$dominance_effects) + effect.pred(X = args$DxE, a = args$DominancexEnvironment_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Dominance_Value")
  }
  if( gen & add & env ){
    genotypic_value$GV <- genotypic_value$Genetic_Value + genotypic_value$Additive_Value +
      + effect.pred(X = args$GxE, a = args$GeneticxEnvironment_effects)
      + effect.pred(X = args$AxE, a = args$AdditivexEnvironment_effects)
  }
  if( gen & dom & env ){
    genotypic_value$GV <- genotypic_value$Genetic_Value + genotypic_value$Dominance_Value
    + effect.pred(X = args$GxE, a = args$GeneticxEnvironment_effects)
    + effect.pred(X = args$DxE, a = args$DominancexEnvironment_effects)
  }
  if( add & dom & env ){
    genotypic_value$GV <- genotypic_value$Additive_Value + genotypic_value$Dominance_Value
    + effect.pred(X = args$AxE, a = args$AdditivexEnvironment_effects)
    + effect.pred(X = args$DxE, a = args$DominancexEnvironment_effects)
  }
  if( gen & add & dom ){
    genotypic_value$GV <- genotypic_value$Genetic_Value + genotypic_value$Additive_Value 
    + genotypic_value$Dominance_Value + effect.pred(X = args$GxE, a = args$GeneticxEnvironment_effects)
    + effect.pred(X = args$AxE, a = args$AdditivexEnvironment_effects)
    + effect.pred(X = args$DxE, a = args$DominancexEnvironment_effects)
  }
  genetic_GxE_value1<-genetic_GxE_value[order(genetic_GxE_value$Env),]
  for(l in 1:ncol(args$env))
  {
    a<-1+((l-1)*ncol(args$fix)*ncol(args$gen))
    b<-l*ncol(args$fix)*ncol(args$gen)
    gv<-data.frame(genetic_GxE_value1[a:b,])
    pheno2<-data.frame(args$pheno[a:b,])
    cat("-----Individual selection - environment " ,l, ": ", "\n")
    rank=gv[order(gv[,4], decreasing = TRUE),]
    print(gv[order(gv[,4], decreasing = TRUE),])
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(gv[,4]))
    write.table(gv[order(gv[,4], decreasing = TRUE),], paste("genetic_value", i, l), quote = FALSE, row.names = FALSE)
    write.table(emp.hpd(gv[,4]), paste("CI_genetic_value", i), quote = FALSE, row.names = FALSE)
    cat("\n")
    plot(gv[,4]~pheno2[,i+3],xlab='Observed',ylab='Predicted',col=2, main=paste('Environment ', l),
         xlim=c(min(pheno1[i,]), max(pheno1[i,])),ylim=c(min(genetic_interaction),max(genetic_interaction))) 
    abline(a=0,b=1,col=4,lwd=2)
    return(rank)
  }
}

family.pred <- function (args){
  add = args$add
  dom = args$dom
  
  genotypic_value           <- as.data.frame( X[,1:3] )
  colnames(genotypic_value) <- c( "Env", "Gen", "Blo" )
  
  if ( add ){
    genotypic_value$add <- effect.pred(X = args$add, a = args$additive_effects)
    colnames(genotypic_value)[ncol(genotypic_value)] <- c("Additive_Value")
  }
  
  if ( dom ){
    genotypic_value$dom <- effect.pred( X = args$dom, a = args$dominance_effects)
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


