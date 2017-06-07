## Script to run bayesian analysis for genomic selection consider GxE interaction
## By Leonardo de Azevedo Peixoto
## Date: April 3rd 2017

####################################################################
## The phenotype file should be display in the way described below: 
## The first colunn should be environment
## The second colunn should be genotypes
## The third colunn should be replication
## Beyond the fourth colunn should be the variables

####################################################################
## Genotype file should contain just the SNPs

####################################################################
setwd("D:\\Leonardo\\POSDOC")
pdf.options(family="Helvetica", height=5, width= 10)
pdf("Full model.pdf", pointsize=6)
nfolds=5
niteraction=2
nIter=nIter
burmIm=burnIn
thin=thin

library(BGLR)
library(boa)
library(TeachingDemos)
library(HapEstXXR)
library(foreach)
library(doParallel)
library(Matrix)

## Reading genotype file
geno<-read.table("RILgeno.txt", h=F)

##calculating alelic frequence 
##Creating a new matrix codified to run HapEstXXR
geno1 = matrix(0,nrow=nrow(geno),ncol=ncol(geno))
#
idx = geno == 2
geno1[idx] = 2
#
idx = geno == 1
geno1[idx] = 3
#
idx = geno == 0
geno1[idx] = 1
colnames(geno1)<-colnames(geno)
newgen=data.frame(maf(geno1, marker.label=colnames(geno1)))
freq<-cbind(((newgen[,1]+(newgen[,2]/2))/newgen[,4]), ((newgen[,3]+(newgen[,2]/2))/newgen[,4]))
colnames(freq)<-c("p", "q")
rm(geno1, newgen)

## Reading phenotype file
pheno<-as.matrix(read.table("RILvalfen.txt", h=F))
pheno1<-as.matrix(pheno[,4:ncol(pheno)])
rownames(pheno1)<-seq(1,nrow(pheno1),1)
nvariable=ncol(pheno)-3


## Creating an incidence matrix for fixed effect 
fix<-data.frame(as.factor(pheno[,3]))
nfix<-table(fix)
fix1<-Matrix(0,nrow=nrow(pheno), ncol=length(nfix), sparse= TRUE)
for(i in 1:nrow(fix))
{
  for (j in 1:length(nfix))
  {
    if (fix[i,]==j){fix1[i,j]=1}
  }
}
rownames(fix1)<-seq(1,nrow(pheno1),1)
rm(fix)

## Creating an incidence matrix for environment effect
env<-data.frame(as.factor(pheno[,1]))
nenv<-table(env)

env1<-Matrix(0,nrow=nrow(pheno), ncol=length(nenv), sparse= TRUE)
for(i in 1:nrow(env))
{
  for (j in 1:length(nenv))
  {
    if (env[i,]==j){env1[i,j]=1}
  }
}
rownames(env1)<-seq(1,nrow(pheno1),1)
rm(env)

## Creating an incidence matrix for genetic effect
gen<-data.frame(as.factor(pheno[,2]))
ngen<-table(gen)
gen1<-Matrix(0,nrow=nrow(pheno), ncol=length(ngen), sparse= TRUE)
for(i in 1:nrow(gen))
{
  for (j in 1:length(ngen))
  {
    if (gen[i,]==j){gen1[i,j]=1}
  }
}
rownames(gen1)<-seq(1,nrow(pheno1),1)
rm(gen)

## Creating an incidence matrix for additive genetic effect
add <- scale(geno,center=T,scale=F)
rownames(add)<-seq(1,nrow(pheno1),1)

## Creating an incidence matrix for dominance genetic effect
dom<-matrix(0,nrow=nrow(geno),ncol=ncol(geno))
for (j in 1:nrow(geno))
{
  for (i in 1:ncol(geno))
  {
    if(geno[j,i] == 2){dom[j,i] = -2*((freq[i,2])^2)}
    if(geno[j,i] == 1){dom[j,i] = 2*(freq[i,2]*freq[i,1])}
    if(geno[j,i] == 0){dom[j,i] = -2*((freq[i,1])^2)}
  }
}
rownames(dom)<-seq(1,nrow(pheno1),1)
rm(geno)

## Creating the incidence matrix of GeneticxEnvironment effect
GxE<-Matrix(0,nrow=nrow(pheno), ncol=length(ngen)*length(nenv), sparse= TRUE)

for(i in 1:nrow(pheno))
{
  GxE[i,]<-Matrix(kronecker(env1[i,], gen1[i,]), sparse= TRUE)
}
rownames(GxE)<-seq(1,nrow(pheno1),1)

## Creating the incidence matrix of AdditivexEnvironment effect
AxE<-Matrix(0,nrow=nrow(pheno), ncol=ncol(add)*length(nenv))
for(i in 1:nrow(pheno))
{
  AxE[i,]<-Matrix(kronecker(env1[i,], add[i,]), sparse= TRUE)
}
rownames(AxE)<-seq(1,nrow(pheno1),1)

## Creating the incidence matrix of DominancexEnvironment effect
DxE<-Matrix(0,nrow=nrow(pheno), ncol=ncol(dom)*length(nenv))

for(i in 1:nrow(pheno))
{
  DxE[i,]<-Matrix(kronecker(env1[i,], dom[i,]))
}
rownames(DxE)<-seq(1,nrow(pheno1),1)

## Defining how many subsets (folds) will be used to run the analysis
subset<-cut(seq(1,nrow(fix1)),breaks=nfolds,labels=FALSE)
subset<-cut(seq(1,nrow(env1)),breaks=nfolds,labels=FALSE)
subset<-cut(seq(1,nrow(gen1)),breaks=nfolds,labels=FALSE)
subset<-cut(seq(1,nrow(add)),breaks=nfolds,labels=FALSE)
subset<-cut(seq(1,nrow(dom)),breaks=nfolds,labels=FALSE)
subset<-cut(seq(1,nrow(GxE)),breaks=nfolds,labels=FALSE)
subset<-cut(seq(1,nrow(AxE)),breaks=nfolds,labels=FALSE)
subset<-cut(seq(1,nrow(DxE)),breaks=nfolds,labels=FALSE)
subset<-cut(seq(1,nrow(pheno1)),breaks=nfolds,labels=FALSE)

## Run bayesian analysis
for (i in 1:nvariable)
{
  genetic_value_accuracy<-matrix(nrow=niteraction, ncol=nfolds)
  prediction_accuracy<-matrix(nrow=niteraction, ncol=nfolds)
  genetic_interaction_accuracy<-matrix(nrow=niteraction, ncol=nfolds)
  residual_variance<-matrix(nrow=niteraction, ncol=nfolds)
  environment_variance<-matrix(nrow=niteraction, ncol=nfolds)
  genetic_variance<-matrix(nrow=niteraction, ncol=nfolds)
  additive_variance<-matrix(nrow=niteraction, ncol=nfolds)
  dominance_variance<-matrix(nrow=niteraction, ncol=nfolds)
  GeneticxEnvironment_variance<-matrix(nrow=niteraction, ncol=nfolds)
  AdditivexEnvironment_variance<-matrix(nrow=niteraction, ncol=nfolds)
  DominancexEnvironment_variance<-matrix(nrow=niteraction, ncol=nfolds)
  fixed_effects<-matrix(nrow=ncol(fix1), ncol=niteraction)
  environment_effects<-matrix(nrow=ncol(env1), ncol=niteraction)
  environment_SD_effects<-matrix(nrow=ncol(env1), ncol=niteraction)
  genetic_effects<-matrix(nrow=ncol(gen1), ncol=niteraction)
  genetic_SD_effects<-matrix(nrow=ncol(gen1), ncol=niteraction)
  additive_effects<-matrix(nrow=ncol(add), ncol=niteraction)
  additive_SD_effects<-matrix(nrow=ncol(add), ncol=niteraction)
  dominance_effects<-matrix(nrow=ncol(dom), ncol=niteraction)
  dominance_SD_effects<-matrix(nrow=ncol(dom), ncol=niteraction)
  GeneticxEnvironment_effects<-matrix(nrow=ncol(GxE), ncol=niteraction)
  GeneticxEnvironment_SD_effects<-matrix(nrow=ncol(GxE), ncol=niteraction)
  AdditivexEnvironment_effects<-matrix(nrow=ncol(AxE), ncol=niteraction)
  AdditivexEnvironment_SD_effects<-matrix(nrow=ncol(AxE), ncol=niteraction)
  DominancexEnvironment_effects<-matrix(nrow=ncol(DxE), ncol=niteraction)
  DominancexEnvironment_SD_effects<-matrix(nrow=ncol(DxE), ncol=niteraction)
  DIC<-matrix(nrow=niteraction, ncol=nfolds)
  for (j in 1:niteraction)
  {
    fix_effect<-matrix(nrow=ncol(fix1), ncol=nfolds)
    env_effect<-matrix(nrow=ncol(env1), ncol=nfolds)
    env_SD_effect<-matrix(nrow=ncol(env1), ncol=nfolds)
    gen_effect<-matrix(nrow=ncol(gen1), ncol=nfolds)
    gen_SD_effect<-matrix(nrow=ncol(gen1), ncol=nfolds)
    add_effect<-matrix(nrow=ncol(add), ncol=nfolds)
    add_SD_effect<-matrix(nrow=ncol(add), ncol=nfolds)
    dom_effect<-matrix(nrow=ncol(dom), ncol=nfolds)
    dom_SD_effect<-matrix(nrow=ncol(dom), ncol=nfolds)
    GxE_effect<-matrix(nrow=ncol(GxE), ncol=nfolds)
    GxE_SD_effect<-matrix(nrow=ncol(GxE), ncol=nfolds)
    AxE_effect<-matrix(nrow=ncol(AxE), ncol=nfolds)
    AxE_SD_effect<-matrix(nrow=ncol(AxE), ncol=nfolds)
    DxE_effect<-matrix(nrow=ncol(DxE), ncol=nfolds)
    DxE_SD_effect<-matrix(nrow=ncol(DxE), ncol=nfolds)
    train<-as.matrix(sample(1:nrow(pheno1), nrow(pheno1), replace=FALSE))
    pheno2<-data.frame(pheno1[as.character(train),i])
    fix2<-fix1[as.character(train),]
    env2<-env1[as.character(train),] 
    gen2<-gen1[as.character(train),]
    add2<-add[as.character(train),]
    dom2<-dom[as.character(train),]
    GxE2<-GxE[as.character(train),]
    AxE2<-AxE[as.character(train),]
    DxE2<-DxE[as.character(train),]
    for (k in 1:nfolds)
    {
      timei<-proc.time() #initial time 
      
      testIndexes<-which(subset==i,arr.ind=TRUE) 
      
      phenotrain<-pheno2[-testIndexes, ] ## fixed effect matrix to train 
      phenotest<-pheno2[testIndexes, ] ## fixed effect matrix to validate 
      
      fixtrain<-fix2[-testIndexes, ] ## fixed effect matrix to train 
      fixtest<-fix2[testIndexes, ] ## fixed effect matrix to validate 
      
      envtrain<-env2[-testIndexes,] ## enviroment effect matrix to train 
      envtest<-env2[testIndexes,]  ## enviroment effect matrix to validate 
      
      gentrain<-gen2[-testIndexes,] ## genetic effect matrix to train 
      gentest<-gen2[testIndexes,] ## genetic effect matrix to validate 
      
      addtrain<-add2[-testIndexes,] ## additive effect matrix to train 
      addtest<-add2[testIndexes,] ## additive effect matrix to validate 
      
      domtrain<-dom2[-testIndexes,] ## dominance effect matrix to train 
      domtest<-dom2[testIndexes,] ## dominance effect matrix to validate 
      
      GxEtrain<-GxE2[-testIndexes,] ## GxE effect matrix to train 
      GxEtest<-GxE2[testIndexes,]  ## GxE effect matrix to validate 
      
      AxEtrain<-AxE2[-testIndexes,] ## GxE effect matrix to train 
      AxEtest<-AxE2[testIndexes,]  ## GxE effect matrix to validate 
      
      DxEtrain<-DxE2[-testIndexes,] ## GxE effect matrix to train 
      DxEtest<-DxE2[testIndexes,]  ## GxE effect matrix to validate 
      
      ETA=list(list(X=fixtrain, model = "FIXED"),
               list(X=envtrain, model = "BRR"),
               list(X=gentrain, model = "BRR"),
               list(X=addtrain, model = "BRR"),
               list(X=domtrain, model = "BRR"),
               list(X=GxEtrain, model = "BRR"),
               list(X=AxEtrain, model = "BRR"),
               list(X=DxEtrain, model = "BRR"))
      results<-BGLR(y=phenotrain, response_type= "gaussian",
                    ETA=ETA,
                    nIter=5000, burnIn=1000, thin=1,
                    saveAt = "results")
      
      ## Genetic Estimate Breeding Value
      residual_variance[j,k]= mean(results$varE) #residual variance
      environment_variance[j,k]= mean(results$ETA[[2]]$varB)
      genetic_variance[j,k]= mean(results$ETA[[3]]$varB)
      additive_variance[j,k]= mean(results$ETA[[4]]$varB)
      dominance_variance[j,k]= mean(results$ETA[[5]]$varB)
      GeneticxEnvironment_variance[j,k]= mean(results$ETA[[6]]$varB)
      AdditivexEnvironment_variance[j,k]= mean(results$ETA[[7]]$varB)
      DominancexEnvironment_variance[j,k]= mean(results$ETA[[8]]$varB)
      fix_effect[,k]<-as.matrix(results$ETA[[1]]$b)
      env_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      gen_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      add_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      dom_effect[,k]<-as.matrix(results$ETA[[5]]$b)
      GxE_effect[,k]<-as.matrix(results$ETA[[6]]$b)
      AxE_effect[,k]<-as.matrix(results$ETA[[7]]$b)
      DxE_effect[,k]<-as.matrix(results$ETA[[8]]$b)
      env_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      gen_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      add_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[5]]$SD.b)
      GxE_SD_effect[,k]<-as.matrix(results$ETA[[6]]$SD.b)
      AxE_SD_effect[,k]<-as.matrix(results$ETA[[7]]$SD.b)
      DxE_SD_effect[,k]<-as.matrix(results$ETA[[8]]$SD.b)
      prediction<-matrix(envtest%*%env_effect[,k] + gentest%*%gen_effect[,k] + 
                               addtest%*%add_effect[,k] + domtest%*%dom_effect[,k] + 
                               GxEtest%*%GxE_effect[,k] + AxEtest%*%AxE_effect[,k] +
                               DxEtest%*%DxE_effect[,k] + fixtest%*%fix_effect[,k])
      genetic_interaction<-matrix(gentest%*%gen_effect[,k] + addtest%*%add_effect[,k] +
                                        domtest%*%dom_effect[,k] + GxEtest%*%GxE_effect[,k] +
                                        AxEtest%*%AxE_effect[,k] + DxEtest%*%DxE_effect[,k])
      genetic_value<-matrix(addtest%*%add_effect[,k])
      statistics<-data.frame(results$fit)
      DIC[j,k]<-statistics[,4]
      
      varE<-scan('resultsvarE.dat')
      plot(varE,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Resildual variance");
      abline(h=results$varE,col=4,lwd=2);
      abline(v=results$burnIn/results$thin,col=4)
      
      varENV<-scan('resultsETA_2_varB.dat')
      plot(varENV,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Environment variance");
      abline(h=results$ETA[[2]]$varB,col=4,lwd=2);
      abline(v=results$burnIn/results$thin,col=4)
      
      varGEN<-scan('resultsETA_3_varB.dat')
      plot(varGEN,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Genetic variance");
      abline(h=results$ETA[[3]]$varB,col=4,lwd=2);
      abline(v=results$burnIn/results$thin,col=4)
      
      varADD<-scan('resultsETA_4_varB.dat')
      plot(varADD,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Additive variance");
      abline(h=results$ETA[[4]]$varB,col=4,lwd=2);
      abline(v=results$burnIn/results$thin,col=4)
      
      varDOM<-scan('resultsETA_5_varB.dat')
      plot(varDOM,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Dominance variance");
      abline(h=results$ETA[[5]]$varB,col=4,lwd=2);
      abline(v=results$burnIn/results$thin,col=4)
      
      varGxE<-scan('resultsETA_6_varB.dat')
      plot(varGxE,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "GxE variance");
      abline(h=results$ETA[[6]]$varB,col=4,lwd=2);
      abline(v=results$burnIn/results$thin,col=4)
      
      varAxE<-scan('resultsETA_7_varB.dat')
      plot(varAxE,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "AxE variance");
      abline(h=results$ETA[[7]]$varB,col=4,lwd=2);
      abline(v=results$burnIn/results$thin,col=4)
      
      varDxE<-scan('resultsETA_8_varB.dat')
      plot(varDxE,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "DxE variance");
      abline(h=results$ETA[[8]]$varB,col=4,lwd=2);
      abline(v=results$burnIn/results$thin,col=4)
      
      ##Convergence diagnostics
      parameters<-data.frame(cbind(varE, varENV, varGEN, varADD, varDOM, varGxE, varAxE, varDxE))
      convergence<-boa.randl(parameters, 0.025, 0.005, 0.95, 0.001)
      cat("-----Raftery and Lewis Convergence diagnostic:", "\n") 
      print(convergence) 
      cat("-----If Depedence factor is lower than 10 the Bayesian analysis converged", "\n")
      cat("-----If Depedence factor is higher than 10 the Bayesian analysis did not converge", "\n")
      cat("-----If Thin calculated is lower than thin used in the Bayesian analysis, it means that the analysis converged", "\n")
      cat("-----If Thin calculated is higher than thin used in the Bayesian analysis, it means that the analysis did not converge", "\n")
      cat("-----If burnIn calculated is lower than burnIn used in the Bayesian analysis, it means that the analysis converged", "\n")
      cat("-----If burnIn calculated is higher than burnIn used in the Bayesian analysis, it means that the analysis did not converge", "\n")
      cat("\n")
      
      ## Estimating the accuracy
      genetic_value_accuracy[j,k]<-cor(genetic_value, matrix(phenotest-fixtest%*%fix_effect[,k]))
      prediction_accuracy[j,k]<-cor(prediction, matrix(phenotest-fixtest%*%fix_effect[,k]))
      genetic_interaction_accuracy[j,k]<-cor(genetic_interaction, matrix(phenotest-fixtest%*%fix_effect[,k]))
    }
    fixed_effects[,j]<-as.matrix(rowMeans(fix_effect))
    environment_effects[,j]<-as.matrix(rowMeans(env_effect))
    environment_SD_effects[,j]<-as.matrix(rowMeans(env_SD_effect))
    genetic_effects[,j]<-as.matrix(rowMeans(gen_effect))
    genetic_SD_effects[,j]<-as.matrix(rowMeans(gen_SD_effect))
    additive_effects[,j]<-as.matrix(rowMeans(add_effect))
    additive_SD_effects[,j]<-as.matrix(rowMeans(add_SD_effect))
    dominance_effects[,j]<-as.matrix(rowMeans(dom_effect))
    dominance_SD_effects[,j]<-as.matrix(rowMeans(dom_SD_effect))
    GeneticxEnvironment_effects[,j]<-as.matrix(rowMeans(GxE_effect))
    GeneticxEnvironment_SD_effects[,j]<-as.matrix(rowMeans(GxE_SD_effect))
    AdditivexEnvironment_effects[,j]<-as.matrix(rowMeans(AxE_effect))
    AdditivexEnvironment_SD_effects[,j]<-as.matrix(rowMeans(AxE_SD_effect))
    DominancexEnvironment_effects[,j]<-as.matrix(rowMeans(DxE_effect))
    DominancexEnvironment_SD_effects[,j]<-as.matrix(rowMeans(DxE_SD_effect))
  }
  
  cat(" -------------------------------------------------------------------", "\n") 
  cat(" Variavel = ",i, "    ",colnames(pheno)[i], "\n") 
  cat(" -------------------------------------------------------------------", "\n") 
  cat("\n") 
  
  cat("-----Genetic accuracy:", "\n") 
  print(genetic_value_accuracy) 
  cat("-----Mean of Genetic accuracy:", "\n") 
  print(mean(genetic_value_accuracy)) 
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(genetic_value_accuracy))) 
  cat("\n") 
  cat("-----Prediction accuracy:", "\n") 
  print(prediction_accuracy)
  cat("-----Mean of Prediction accuracy:", "\n") 
  print(mean(prediction_accuracy))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(prediction_accuracy))) 
  cat("\n") 
  cat("-----Genetic + GxE accuracy:", "\n") 
  print(genetic_interaction_accuracy)
  cat("-----Mean of Genetic + GxE accuracy:", "\n") 
  print(mean(genetic_interaction_accuracy))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(genetic_interaction_accuracy))) 
  cat("\n") 
  cat("-----Residual variance:", "\n") 
  print(residual_variance)
  cat("-----Mean of residual variance:", "\n") 
  print(mean(residual_variance))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(residual_variance))) 
  cat("\n") 
  cat("-----Fixed effects:", "\n") 
  print(fixed_effects)
  cat("-----Mean of fixed effects:", "\n") 
  print(mean(fixed_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(fixed_effects))) 
  cat("\n")
  cat("-----Environmental effects:", "\n") 
  print(environment_effects)
  cat("-----Mean of environmental effects:", "\n") 
  print(mean(environment_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(environment_effects)))
  cat("\n")
  cat("-----Environmental standard deviation:", "\n") 
  print(environment_SD_effects)
  cat("-----Mean of environmental standard deviation:", "\n") 
  print(mean(environment_SD_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(environment_SD_effects)))
  cat("\n")
  plot((rowMeans(environment_effects))^2, ylab='Estimated Squared-Environment Effect',
       type='o',cex=.5,col=4,main='Environment Effects')
  cat("-----Environmental variance:", "\n") 
  print(environment_variance)
  cat("-----Mean of environmental variance:", "\n") 
  print(mean(environment_variance))
  cat("\n")
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(environment_variance)))
  cat("-----Genetic effects:", "\n") 
  print(genetic_effects)
  cat("-----Mean of Genetic effects:", "\n") 
  print(mean(genetic_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(genetic_effects)))
  cat("\n")
  cat("-----Genetic standard deviation:", "\n") 
  print(genetic_SD_effects)
  cat("-----Mean of Genetic standard deviation:", "\n") 
  print(mean(genetic_SD_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(genetic_SD_effects)))
  cat("\n")
  plot((rowMeans(genetic_effects))^2, ylab='Estimated Squared-Genetic Effect',
       type='o',cex=.5,col=4,main='Genetic Effects')
  cat("-----Genetic variance:", "\n") 
  print(genetic_variance)
  cat("-----Mean of Genetic variance:", "\n") 
  print(mean(genetic_variance))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(genetic_variance)))
  cat("\n")
  cat("-----Additive effects:", "\n") 
  print(additive_effects)
  cat("-----Mean of additive effects:", "\n") 
  print(mean(additive_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(additive_effects)))
  cat("\n")
  cat("-----Additive standard deviation:", "\n") 
  print(additive_SD_effects)
  cat("-----Mean of Additive standard deviation:", "\n") 
  print(mean(additive_SD_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(additive_SD_effects)))
  cat("\n")
  plot((rowMeans(additive_effects))^2, ylab='Estimated Squared-Additive Effect',
       type='o',cex=.5,col=4,main='Additive Effects')
  cat("-----Additive variance:", "\n") 
  print(additive_variance)
  cat("-----Mean of Additive variance:", "\n") 
  print(mean(additive_variance))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(additive_variance)))
  cat("\n")
  cat("-----Dominance effects:", "\n") 
  print(dominance_effects)
  cat("-----Mean of dominance effects:", "\n") 
  print(mean(dominance_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(dominance_effects)))
  cat("\n")
  cat("-----Dominance standard deviation:", "\n") 
  print(dominance_SD_effects)
  cat("-----Mean of dominance standard deviation:", "\n") 
  print(mean(dominance_SD_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(dominance_SD_effects)))
  cat("\n")
  plot((rowMeans(dominance_effects))^2, ylab='Estimated Squared-Dominance Effect',
       type='o',cex=.5,col=4,main='Dominance Effects')
  cat("-----Dominance variance:", "\n") 
  print(dominance_variance)
  cat("-----Mean of dominance variance:", "\n") 
  print(mean(dominance_variance))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(dominance_variance)))
  cat("\n")
  cat("-----GeneticxEnvironment effects:", "\n") 
  print(GeneticxEnvironment_effects)
  cat("-----Mean of GeneticxEnvironment effects:", "\n") 
  print(mean(GeneticxEnvironment_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(GeneticxEnvironment_effects)))
  cat("\n")
  cat("-----GeneticxEnvironment standard deviation:", "\n") 
  print(GeneticxEnvironment_SD_effects)
  cat("-----Mean of GeneticxEnvironment standard deviation:", "\n") 
  print(mean(GeneticxEnvironment_SD_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(GeneticxEnvironment_SD_effects)))
  cat("\n")
  plot((rowMeans(GeneticxEnvironment_effects))^2, ylab='Estimated Squared-GxE Effect',
       type='o',cex=.5,col=4,main='GxE Effects')
  cat("-----GeneticxEnvironment variance:", "\n") 
  print(GeneticxEnvironment_variance)
  cat("-----Mean of GeneticxEnvironment variance:", "\n") 
  print(mean(GeneticxEnvironment_variance))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(GeneticxEnvironment_variance)))
  cat("\n")
  cat("-----AdditivexEnvironment effects:", "\n") 
  print(AdditivexEnvironment_effects)
  cat("-----Mean of AdditivexEnvironment effects:", "\n") 
  print(mean(AdditivexEnvironment_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(AdditivexEnvironment_effects)))
  cat("\n")
  cat("-----AdditivexEnvironment standard deviation:", "\n") 
  print(AdditivexEnvironment_SD_effects)
  cat("-----Mean of AdditivexEnvironment standard deviation:", "\n") 
  print(mean(AdditivexEnvironment_SD_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(AdditivexEnvironment_SD_effects)))
  cat("\n")
  plot((rowMeans(AdditivexEnvironment_effects))^2, ylab='Estimated Squared-AxE Effect',
       type='o',cex=.5,col=4,main='AxE Effects')
  cat("-----AdditivexEnvironment variance:", "\n") 
  print(AdditivexEnvironment_variance)
  cat("-----Mean of AdditivexEnvironment variance:", "\n") 
  print(mean(AdditivexEnvironment_variance))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(AdditivexEnvironment_variance)))
  cat("\n")
  cat("-----DominancexEnvironment effects:", "\n") 
  print(DominancexEnvironment_effects)
  cat("-----Mean of DominancexEnvironment effects:", "\n") 
  print(mean(DominancexEnvironment_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(DominancexEnvironment_effects)))
  cat("\n")
  cat("-----DominancexEnvironment standard deviation:", "\n") 
  print(DominancexEnvironment_SD_effects)
  cat("-----Mean of DominancexEnvironment standard deviation:", "\n") 
  print(mean(DominancexEnvironment_SD_effects))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(DominancexEnvironment_SD_effects)))
  cat("\n")
  plot((rowMeans(DominancexEnvironment_effects))^2, ylab='Estimated Squared-DxE Effect',
       type='o',cex=.5,col=4,main='DxE Effects')
  cat("-----DominancexEnvironment variance:", "\n") 
  print(DominancexEnvironment_variance)
  cat("-----Mean of AdditivexEnvironment variance:", "\n") 
  print(mean(DominancexEnvironment_variance))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(DominancexEnvironment_variance)))
  cat("\n")
  cat("-----Deviance information criteria:", "\n") 
  print(DIC)
  cat("-----Mean of Deviance information criteria:", "\n") 
  print(mean(DIC))
  cat("\n")
  
  ##Calculating the geentic value for all full-sib families
  genetic_value<-data.frame(pheno[,1:3],add%*%rowMeans(additive_effects))
  colnames(genetic_value)<-c("Env", "Gen", "Blo", "GV")
  genetic_value1<-genetic_value[order(genetic_value$Env),]
  family_genetic_value<-matrix(nrow=ncol(gen1), ncol=1)
  for (x in 1:ncol(gen1))
  {
    fgv<-matrix(0,nrow=nrow(pheno1), ncol=1)
    for (y in 1:nrow(pheno1))
    {
      if (genetic_value1[y,2]==x){fgv[y,1]=genetic_value1[y,4]}
    }
    family_genetic_value[x,1]<-sum(fgv)/ncol(env1)*ncol(fix1)
  }
  rownames(family_genetic_value)<-seq(1, max(pheno[,2]), 1)
  cat("-----Genetic value for each family:", "\n") 
  print(family_genetic_value)
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(family_genetic_value[,1]))
  cat("\n")
  
  ##Calculating the genetic value for all individuals
  genetic_GxE_value<-data.frame(pheno[,1:3],matrix((gen1%*%rowMeans(genetic_effects)+add%*%rowMeans(additive_effects)+
                                               dom%*%rowMeans(dominance_effects)+GxE%*%rowMeans(GeneticxEnvironment_effects)+
                                               AxE%*%rowMeans(AdditivexEnvironment_effects)+DxE%*%rowMeans(DominancexEnvironment_effects))))
  colnames(genetic_GxE_value)<-c("Env", "Gen", "Blo", "GV")
  genetic_GxE_value1<-genetic_GxE_value[order(genetic_GxE_value$Env),]
  for(l in 1:ncol(env1))
  {
    a<-1+((l-1)*ncol(fix1)*ncol(gen1))
    b<-l*ncol(fix1)*ncol(gen1)
    gv<-data.frame(genetic_GxE_value1[a:b,])
    pheno2<-data.frame(pheno[a:b,])
    cat("-----Individual selection - environment " ,l, ": ", "\n")
    print(gv[order(gv[,4], decreasing = TRUE),])
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(gv[,4]))
    write.table(gv[order(gv[,4], decreasing = TRUE),], paste("genetic_value", i, l), quote = FALSE, row.names = FALSE)
    write.table(emp.hpd(gv[,4]), paste("CI_genetic_value", i), quote = FALSE, row.names = FALSE)
    cat("\n")
    plot(gv[,4]~pheno2[,i+3],xlab='Observed',ylab='Predicted',col=2, main=paste('Environment ', l),
         xlim=c(min(pheno1[i,]), max(pheno1[i,])),ylim=c(min(genetic_interaction),max(genetic_interaction))) 
    abline(a=0,b=1,col=4,lwd=2)
  }
  
  write.table(convergence, paste("convergence", i ,j ,k), quote = FALSE, row.names = FALSE)
  write.table(genetic_value_accuracy, paste("genetic_value_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(genetic_value_accuracy)), paste("CI_genetic_value_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(prediction_accuracy, paste("prediction_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(prediction_accuracy)), paste("CI_prediction_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(genetic_interaction_accuracy, paste("genetic_interaction_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(genetic_interaction_accuracy)), paste("CI_genetic_interaction_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(residual_variance, paste("residual_variance", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(residual_variance)), paste("CI_residual_variance", i), quote = FALSE, row.names = FALSE)
  write.table(fixed_effects, paste("fixed_effectse", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(fixed_effects)), paste("CI_fixed_effects", i), quote = FALSE, row.names = FALSE)
  write.table(environment_effects, paste("environment_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(environment_effects)), paste("CI_environment_effects", i), quote = FALSE, row.names = FALSE)
  write.table(environment_SD_effects, paste("environment_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(environment_SD_effects)), paste("CI_environment_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(environment_variance, paste("environment_variance", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(environment_variance)), paste("CI_environment_variance", i), quote = FALSE, row.names = FALSE)
  write.table(genetic_effects, paste("genetic_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(genetic_effects)), paste("CI_genetic_effects", i), quote = FALSE, row.names = FALSE)
  write.table(genetic_SD_effects, paste("genetic_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(genetic_SD_effects)), paste("CI_genetic_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(genetic_variance, paste("genetic_variance", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(genetic_variance)), paste("CI_genetic_variance", i), quote = FALSE, row.names = FALSE)
  write.table(additive_effects, paste("addtive_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(additive_effects)), paste("CI_addtive_effects", i), quote = FALSE, row.names = FALSE)
  write.table(additive_SD_effects, paste("addtive_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(additive_SD_effects)), paste("CI_addtive_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(additive_variance, paste("addtive_variance", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(additive_variance)), paste("CI_addtive_variance", i), quote = FALSE, row.names = FALSE)
  write.table(dominance_effects, paste("dominance_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(dominance_effects)), paste("CI_dominance_effects", i), quote = FALSE, row.names = FALSE)
  write.table(dominance_SD_effects, paste("dominance_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(dominance_SD_effects)), paste("CI_dominance_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(dominance_variance, paste("dominance_variance", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(dominance_variance)), paste("CI_dominance_variance", i), quote = FALSE, row.names = FALSE)
  write.table(GeneticxEnvironment_effects, paste("GeneticxEnvironment_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(GeneticxEnvironment_effects)), paste("CI_GeneticxEnvironment_effects", i), quote = FALSE, row.names = FALSE)
  write.table(GeneticxEnvironment_SD_effects, paste("GeneticxEnvironment_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(GeneticxEnvironment_SD_effects)), paste("CI_GeneticxEnvironment_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(GeneticxEnvironment_variance, paste("GeneticxEnvironment_variance", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(GeneticxEnvironment_variance)), paste("CI_GeneticxEnvironment_variance", i), quote = FALSE, row.names = FALSE)
  write.table(AdditivexEnvironment_effects, paste("AdditivexEnvironment_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(AdditivexEnvironment_effects)), paste("CI_AdditivexEnvironment_effects", i), quote = FALSE, row.names = FALSE)
  write.table(AdditivexEnvironment_SD_effects, paste("AdditivexEnvironment_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(AdditivexEnvironment_SD_effects)), paste("CI_AdditivexEnvironment_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(AdditivexEnvironment_variance, paste("AdditivexEnvironment_variance", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(AdditivexEnvironment_variance)), paste("CI_AdditivexEnvironment_variance", i), quote = FALSE, row.names = FALSE)
  write.table(DominancexEnvironment_effects, paste("DominancexEnvironment_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(DominancexEnvironment_effects)), paste("CI_DominancexEnvironment_effects", i), quote = FALSE, row.names = FALSE)
  write.table(DominancexEnvironment_SD_effects, paste("DominancexEnvironment_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(DominancexEnvironment_SD_effects)), paste("CI_DominancexEnvironment_SD_effects", i), quote = FALSE, row.names = FALSE)
  write.table(DominancexEnvironment_variance, paste("DominancexEnvironment_variance", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(DominancexEnvironment_variance)), paste("CI_DominancexEnvironment_variance", i), quote = FALSE, row.names = FALSE)
  write.table(DIC, paste("DIC", i), quote = FALSE, row.names = FALSE)
  write.table(family_genetic_value, paste("family_genetic_value", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(family_genetic_value[,1]), paste("CI_family_genetic_value", i), quote = FALSE, row.names = FALSE)
}

#Disconnect from the output file
dev.off()

















