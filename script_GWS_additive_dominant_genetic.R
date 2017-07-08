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
## Model file should have eight rows which each row has one model for the follow effects
## Fixed, environment, genetic (traditional BLUP), additive (markers), dominance (markers),
## GenxEnv, AddxEnv, DomxEnv
## If your model there is no one of these effects type NULL in the model file for the respective line
## Users can choose any model that there is in BGLR package such as FIXED, BRR, BA, BB, BC, and BL

####################################################################
setwd("C:\\Users\\Biometria_pc31\\Google Drive\\Simula??o GWS\\Simula??o RILs")

## Reading genotype file
geno<-read.table("RILgeno.txt", h=F)
geno<-data.frame(geno[,1:500])
## Reading phenotype file
pheno<-as.matrix(read.table("RILvalfen.txt", h=F))
pheno<-data.frame(pheno[1:500,])
## Reading the model file (file which specifies what bayesian model should be used for each factor)
model<-read.table("model.txt", h=F, text=TRUE)
model<-as.matrix(c("FIXED", "BRR", "BRR", "BRR", "BRR", "BRR", "BRR", "BRR"))
library(foreach)
library(argparse)
library(BGLR)
library(rrBLUP)
library(boa)
library(TeachingDemos)
library(HapEstXXR)
library(Matrix)


## Items that should be informed by the User
nfolds=5
niteraction=2
nIter=1000
burmIm=100
thin=1

Fix = TRUE
Env = TRUE
Gen = TRUE
Add = TRUE
Dom = TRUE
Fam = FALSE

## Creating a file called "args" with all files that need to run GWS model
X  = pheno

args <- list('X' = X , 'Fix' = NULL , 'Env' = NULL , 'Gen' = NULL , 'Add' = NULL , 'Dom' = NULL ,
             'GxE' = NULL , 'AxE' = NULL , 'DxE' = NULL , 'fixtest' = NULL, 'envtest' = NULL, 
             'gentest' = NULL, 'addtest' = NULL, 'domtest' = NULL, 'GxEtest' = NULL, 
             'AxEtest' = NULL, 'DxEtest' = NULL, 'fixtrain' = NULL, 'envtrain' = NULL, 
             'gentrain' = NULL, 'addtrain' = NULL, 'domtrain' = NULL, 'GxEtrain' = NULL, 
             'AxEtrain' = NULL, 'DxEtrain' = NULL, 'fixed_effects' = NULL , 'environment_effects' = NULL ,
             'genetic_effects' = NULL , 'additive_effects' = NULL , 'dominance_effects' = NULL , 
             'AdditivexEnvironment_effects' = NULL , 'DominancexEnvironment_effects' = NULL ,
             'GeneticxEnvironment_effects' = NULL , 'fix_effect' = NULL, 'env_effect' = NULL,
             'gen_effect' = NULL, 'add_effect' = NULL, 'dom_effect' = NULL, 'GxE_effect' = NULL, 
             'AxE_effect' = NULL, 'DxE_effect' = NULL, 'train' = NULL, 'fix' = Fix, 'env' = Env, 'gen' = Gen,
             'add' = Add, 'dom' = Dom, 'fam' = Fam )

## Creating a pdf file to save all graphics
pdf.options(family="Helvetica", height=5, width= 10)
pdf("Full model.pdf", pointsize=6)

##calculating alelic frequence 
Z = recode.genotype(M = geno)

newgen=data.frame(maf(Z, marker.label=colnames(Z)))
freq<-cbind(((newgen[,1]+(newgen[,2]/2))/newgen[,4]), ((newgen[,3]+(newgen[,2]/2))/newgen[,4]))
colnames(freq)<-c("p", "q")
rm(Z, newgen)

pheno1<-as.matrix(pheno[,4:ncol(pheno)])
rownames(pheno1)<-seq(1,nrow(pheno1),1)
nvariable=ncol(pheno)-3
Nids = nrow(pheno)

## Creating an incidence matrix for fixed effect
lev  = unique(args$X[,3])
Nfix = length(lev)
if (Fix & Nfix>1) {args$Fix = design.matrix(X = pheno[,3], Nids = Nids)}

## Creating an incidence matrix for environment effect
lev  = unique(args$X[,1])
Nenv = length(lev)
if (Env & Nenv>1) {args$Env = design.matrix(X = pheno[,1], Nids = Nids)}

## Creating an incidence matrix for genetic effect
if (Gen) {args$Gen = design.matrix(X = pheno[,2], Nids = Nids)}

## Creating an incidence matrix for additive genetic effect
if (Add) {args$Add = scale(geno,center=T,scale=F)
rownames(args$Add) <- 1:Nids}

## Creating an incidence matrix for dominance genetic effect
if (Dom) {args$Dom <- scale.dom(X = geno)}

## Creating the incidence matrix of GeneticxEnvironment effect
if (Gen & Env & Nenv>1) {args$GxE <- get.matrix(m1 = args$Gen, m2 = args$Env, nsnp = ncol(args$Gen) , nRow = nrow(args$Gen), nLev = ncol(args$Env))}

## Creating the incidence matrix of AdditivexEnvironment effect
if (Add & Env & Nenv>1) {args$AxE <- get.matrix(m1 = args$Add, m2 = args$Env, nsnp = ncol(args$Add) , nRow = nrow(args$Add), nLev = ncol(args$Env))}

## Creating the incidence matrix of DominancexEnvironment effect
if (Dom & Env & Nenv>1) {args$DxE <- get.matrix(m1 = args$Dom, m2 = args$Env, nsnp = ncol(args$Dom) , nRow = nrow(args$Dom), nLev = ncol(args$Env))}

## Defining how many subsets (folds) will be used to run the analysis
subset<-cut(seq(1,nrow(pheno1)),breaks=nfolds,labels=FALSE)

## Run bayesian analysis
for (i in 1:nvariable)
{
  genetic_value_accuracy<-matrix(nrow=niteraction, ncol=nfolds)
  prediction_accuracy<-matrix(nrow=niteraction, ncol=nfolds)
  genotypic_value_accuracy<-matrix(nrow=niteraction, ncol=nfolds)
  residual_variance<-matrix(nrow=niteraction, ncol=nfolds)
  if (Nenv>1) {environment_variance<-matrix(nrow=niteraction, ncol=nfolds)}
  if (Gen) {genetic_variance<-matrix(nrow=niteraction, ncol=nfolds)}
  if (Add) {additive_variance<-matrix(nrow=niteraction, ncol=nfolds)}
  if (Dom) {dominance_variance<-matrix(nrow=niteraction, ncol=nfolds)}
  if (Nenv>1) {GeneticxEnvironment_variance<-matrix(nrow=niteraction, ncol=nfolds)}
  if (Nenv>1) {AdditivexEnvironment_variance<-matrix(nrow=niteraction, ncol=nfolds)}
  if (Nenv>1) {DominancexEnvironment_variance<-matrix(nrow=niteraction, ncol=nfolds)}
  if (Nfix>1) {args$fixed_effects<-matrix(nrow=ncol(args$Fix), ncol=niteraction)}
  if (Nenv>1) {args$environment_effects<-matrix(nrow=ncol(args$Env), ncol=niteraction)}
  if (Nenv>1) {environment_SD_effects<-matrix(nrow=ncol(args$Env), ncol=niteraction)}
  if (Gen) {args$genetic_effects<-matrix(nrow=ncol(args$Gen), ncol=niteraction)}
  if (Gen) {genetic_SD_effects<-matrix(nrow=ncol(args$Gen), ncol=niteraction)}
  if (Add) {args$additive_effects<-matrix(nrow=ncol(args$Add), ncol=niteraction)}
  if (Add) {additive_SD_effects<-matrix(nrow=ncol(args$Add), ncol=niteraction)}
  if (Dom) {args$dominance_effects<-matrix(nrow=ncol(args$Dom), ncol=niteraction)}
  if (Dom) {dominance_SD_effects<-matrix(nrow=ncol(args$Dom), ncol=niteraction)}
  if (Gen & Nenv>1) {args$GeneticxEnvironment_effects<-matrix(nrow=ncol(args$GxE), ncol=niteraction)}
  if (Gen & Nenv>1) {GeneticxEnvironment_SD_effects<-matrix(nrow=ncol(args$GxE), ncol=niteraction)}
  if (Add & Nenv>1) {args$AdditivexEnvironment_effects<-matrix(nrow=ncol(args$AxE), ncol=niteraction)}
  if (Add & Nenv>1) {AdditivexEnvironment_SD_effects<-matrix(nrow=ncol(args$AxE), ncol=niteraction)}
  if (Dom & Nenv>1) {args$DominancexEnvironment_effects<-matrix(nrow=ncol(args$DxE), ncol=niteraction)}
  if (Dom & Nenv>1) {DominancexEnvironment_SD_effects<-matrix(nrow=ncol(args$DxE), ncol=niteraction)}
  DIC<-matrix(nrow=niteraction, ncol=nfolds)
  for (j in 1:niteraction)
  {
    if (Nfix>1) {args$fix_effect<-matrix(nrow=ncol(args$fix), ncol=nfolds)}
    if (Nenv>1) {args$env_effect<-matrix(nrow=ncol(args$Env), ncol=nfolds)
    env_SD_effect<-matrix(nrow=ncol(args$Env), ncol=nfolds)}
    if (Gen) {args$gen_effect<-matrix(nrow=ncol(args$Gen), ncol=nfolds)
    gen_SD_effect<-matrix(nrow=ncol(args$Gen), ncol=nfolds)}
    if (Add) {args$add_effect<-matrix(nrow=ncol(args$Add), ncol=nfolds)
    add_SD_effect<-matrix(nrow=ncol(args$Add), ncol=nfolds)}
    if (Dom) {args$dom_effect<-matrix(nrow=ncol(args$Dom), ncol=nfolds)
    dom_SD_effect<-matrix(nrow=ncol(args$Dom), ncol=nfolds)}
    if (Gen & Nenv>1) {args$GxE_effect<-matrix(nrow=ncol(args$GxE), ncol=nfolds)
    GxE_SD_effect<-matrix(nrow=ncol(args$GxE), ncol=nfolds)}
    if (Add & Nenv>1) {args$AxE_effect<-matrix(nrow=ncol(args$AxE), ncol=nfolds)
    AxE_SD_effect<-matrix(nrow=ncol(args$AxE), ncol=nfolds)}
    if (Dom & Nenv>1) {args$DxE_effect<-matrix(nrow=ncol(args$DxE), ncol=nfolds)
    DxE_SD_effect<-matrix(nrow=ncol(args$DxE), ncol=nfolds)}
    train<-as.matrix(sample(1:nrow(pheno1), nrow(pheno1), replace=FALSE))
    pheno2<-matrix.creation(X=pheno1[,i], args=args, train=train) 
    if (Nfix>1) {fix2<- matrix.creation(X=args$fix, args=args, train=train)} 
    if (Nenv>1) {env2<-matrix.creation(X=args$Env, args=args, train=train)}
    if (Gen) {gen2<-matrix.creation(X=args$Gen, args=args, train=train)}
    if (Add) {add2<-matrix.creation(X=args$Add, args=args, train=train)}
    if (Dom) {dom2<-matrix.creation(X=args$Dom, args=args, train=train)}
    if (Gen & Nenv>1) {GxE2<-matrix.creation(X=args$GxE, args=args, train=train)}
    if (Add & Nenv>1) {AxE2<-matrix.creation(X=args$AxE, args=args, train=train)}
    if (Dom & Nenv>1) {DxE2<-matrix.creation(X=args$DxE, args=args, train=train)}
    for (k in 1:nfolds)
    {
      timei<-proc.time() #initial time 
      
      testIndexes<-which(subset==i,arr.ind=TRUE) 
      
      phenotrain<-pheno2[-testIndexes, ] ## Phenotype matrix to train 
      phenotest<-as.matrix(pheno2[testIndexes, ]) ## Phenotype effect matrix to validate 
      
      if (Fix & Nfix>1) {args$fixtrain<-as.matrix(fix2[-testIndexes, ]) ## fixed effect matrix to train 
      args$fixtest<-as.matrix(fix2[testIndexes, ])} ## fixed effect matrix to validate 
      
      if (Env & Nenv>1) {args$envtrain<-as.matrix(env2[-testIndexes,]) ## enviroment effect matrix to train 
      args$envtest<-as.matrix(env2[testIndexes,])}  ## enviroment effect matrix to validate 
      
      if (Gen) {args$gentrain<-as.matrix(gen2[-testIndexes,]) ## genetic effect matrix to train 
      args$gentest<-as.matrix(gen2[testIndexes,])} ## genetic effect matrix to validate 
      
      if (Add) {args$addtrain<-as.matrix(add2[-testIndexes,]) ## additive effect matrix to train 
      args$addtest<-as.matrix(add2[testIndexes,])} ## additive effect matrix to validate 
      
      if (Dom) {args$domtrain<-as.matrix(dom2[-testIndexes,]) ## dominance effect matrix to train 
      args$domtest<-as.matrix(dom2[testIndexes,])} ## dominance effect matrix to validate 
      
      if (Gen & Env & Nenv>1) {args$GxEtrain<-as.matrix(GxE2[-testIndexes,]) ## GxE effect matrix to train 
      args$GxEtest<-as.matrix(GxE2[testIndexes,])}  ## GxE effect matrix to validate 
      
      if (Add & Env & Nenv>1) {args$AxEtrain<-as.matrix(AxE2[-testIndexes,]) ## GxE effect matrix to train 
      args$AxEtest<-as.matrix(AxE2[testIndexes,])}  ## GxE effect matrix to validate 
      
      if (Dom & Env & Nenv>1) {args$DxEtrain<-as.matrix(DxE2[-testIndexes,]) ## GxE effect matrix to train 
      args$DxEtest<-as.matrix(DxE2[testIndexes,])}  ## GxE effect matrix to validate 
      
      ETA=ETA.pred(args=args, model=model)
      results<-BGLR(y=phenotrain, response_type= "gaussian",
                    ETA=ETA,
                    nIter=5000, burnIn=1000, thin=1,
                    saveAt = "results")
      
      ## Genetic Estimate Breeding Value
      residual_variance[j,k]= mean(results$varE) #residual variance
      
      if (Fix & Nfix>1) {args$fix_effect[,k]<-as.matrix(results$ETA[[1]]$b)}
      if (Env & Nenv>1) {environment_variance[j,k]= mean(results$ETA[[1]]$varB)
      args$env_effect[,k]<-as.matrix(results$ETA[[1]]$b)
      env_SD_effect[,k]<-as.matrix(results$ETA[[1]]$SD.b)
      varENV<-scan('resultsETA_1_varB.dat')
      }
      if (Gen) {genetic_variance[j,k]= mean(results$ETA[[1]]$varB)
      args$gen_effect[,k]<-as.matrix(results$ETA[[1]]$b)
      gen_SD_effect[,k]<-as.matrix(results$ETA[[1]]$SD.b)
      varGEN<-scan('resultsETA_1_varB.dat')
      }
      if (Add) {additive_variance[j,k]= mean(results$ETA[[1]]$varB)
      args$add_effect[,k]<-as.matrix(results$ETA[[1]]$b)
      add_SD_effect[,k]<-as.matrix(results$ETA[[1]]$SD.b)
      varADD<-scan('resultsETA_1_varB.dat')
      }
      if (Dom) {dominance_variance[j,k]= mean(results$ETA[[1]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[1]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[1]]$SD.b)
      varDOM<-scan('resultsETA_1_varB.dat')
      }
      if (Env & Fix & Nenv>1 & Nfix>1) {environment_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$env_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      env_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      varENV<-scan('resultsETA_2_varB.dat')
      }
      if (Gen & Fix & Nfix>1) {genetic_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$gen_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      gen_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      varGEN<-scan('resultsETA_2_varB.dat')
      }
      if (Gen & Env & Nenv>1) {genetic_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$gen_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      gen_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      GeneticxEnvironment_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$GxE_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      GxE_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      varENV<-scan('resultsETA_1_varB.dat')
      varGEN<-scan('resultsETA_2_varB.dat')
      varGxE<-scan('resultsETA_3_varB.dat')
      }
      if (Add & Fix & Nfix>1) {additive_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$add_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      add_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      varADD<-scan('resultsETA_2_varB.dat')
      }
      if (Add & Env & Nenv>1) {additive_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$add_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      add_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      AdditivexEnvironment_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$AxE_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      AxE_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      varENV<-scan('resultsETA_1_varB.dat')
      varADD<-scan('resultsETA_2_varB.dat')
      varAxE<-scan('resultsETA_3_varB.dat')
      }
      if (Add & Gen) {additive_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$add_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      add_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      varGEN<-scan('resultsETA_1_varB.dat')
      varADD<-scan('resultsETA_2_varB.dat')
      }
      if (Dom & Add) {dominance_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      varADD<-scan('resultsETA_1_varB.dat')
      varDOM<-scan('resultsETA_2_varB.dat')
      }
      if (Dom & Gen) {dominance_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      varGEN<-scan('resultsETA_1_varB.dat')
      varDOM<-scan('resultsETA_2_varB.dat')
      }
      if (Dom & Env & Nenv>1) {dominance_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      DominancexEnvironment_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$DxE_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      DxE_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      varENV<-scan('resultsETA_1_varB.dat')
      varDOM<-scan('resultsETA_2_varB.dat')
      varDxE<-scan('resultsETA_3_varB.dat')
      }
      if (Dom & Fix & Nfix>1) {dominance_variance[j,k]= mean(results$ETA[[2]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[2]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[2]]$SD.b)
      varDOM<-scan('resultsETA_2_varB.dat')
      }
      if (Gen & Env & Fix & Nenv>1 & Nfix>1) {genetic_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$gen_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      gen_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      GeneticxEnvironment_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$GxE_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      GxE_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      varENV<-scan('resultsETA_2_varB.dat')
      varGEN<-scan('resultsETA_3_varB.dat')
      varGxE<-scan('resultsETA_3_varB.dat')
      }
      if (Add & Env & Fix &Nenv>1 & Nfix>1) {additive_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$add_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      add_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      AdditivexEnvironment_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$AxE_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      AxE_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      varENV<-scan('resultsETA_2_varB.dat')
      varADD<-scan('resultsETA_3_varB.dat')
      varAxE<-scan('resultsETA_4_varB.dat')
      }
      if (Add & Gen & Env & Nenv>1) {additive_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$add_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      add_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      GeneticxEnvironment_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$GxE_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      GxE_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      AdditivexEnvironment_variance[j,k]= mean(results$ETA[[5]]$varB)
      args$AxE_effect[,k]<-as.matrix(results$ETA[[5]]$b)
      AxE_SD_effect[,k]<-as.matrix(results$ETA[[5]]$SD.b)
      varENV<-scan('resultsETA_1_varB.dat')
      varGEN<-scan('resultsETA_2_varB.dat')
      varADD<-scan('resultsETA_3_varB.dat')
      varGxE<-scan('resultsETA_4_varB.dat')
      varAxE<-scan('resultsETA_5_varB.dat')
      }
      if (Add & Gen & Fix & Nfix>1) {additive_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$add_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      add_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      varGEN<-scan('resultsETA_2_varB.dat')
      varADD<-scan('resultsETA_3_varB.dat')
      }
      if (Dom & Gen & Fix & Nfix>1) {dominance_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      varGEN<-scan('resultsETA_2_varB.dat')
      varDOM<-scan('resultsETA_3_varB.dat')
      }
      if (Dom & Gen & Env & Nenv>1) {dominance_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      GeneticxEnvironment_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$GxE_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      GxE_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      DominancexEnvironment_variance[j,k]= mean(results$ETA[[5]]$varB)
      args$DxE_effect[,k]<-as.matrix(results$ETA[[5]]$b)
      DxE_SD_effect[,k]<-as.matrix(results$ETA[[5]]$SD.b)
      varENV<-scan('resultsETA_1_varB.dat')
      varGEN<-scan('resultsETA_2_varB.dat')
      varDOM<-scan('resultsETA_3_varB.dat')
      varGxE<-scan('resultsETA_4_varB.dat')
      varDxE<-scan('resultsETA_5_varB.dat')
      }
      if (Dom & Add & Fix & Nfix>1) {dominance_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      varADD<-scan('resultsETA_2_varB.dat')
      varDOM<-scan('resultsETA_3_varB.dat')
      }
      if (Dom & Add & Env & Nenv>1) {dominance_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      AdditivexEnvironment_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$AxE_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      AxE_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      DominancexEnvironment_variance[j,k]= mean(results$ETA[[5]]$varB)
      args$DxE_effect[,k]<-as.matrix(results$ETA[[5]]$b)
      DxE_SD_effect[,k]<-as.matrix(results$ETA[[5]]$SD.b)
      varENV<-scan('resultsETA_1_varB.dat')
      varADD<-scan('resultsETA_2_varB.dat')
      varDOM<-scan('resultsETA_3_varB.dat')
      varAxE<-scan('resultsETA_4_varB.dat')
      varDxE<-scan('resultsETA_5_varB.dat')
      }
      if (Dom & Add & Gen) {dominance_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      varGEN<-scan('resultsETA_1_varB.dat')
      varADD<-scan('resultsETA_2_varB.dat')
      varDOM<-scan('resultsETA_3_varB.dat')
      }
      if (Dom & Env & Fix & Nenv>1 & Nfix>1) {dominance_variance[j,k]= mean(results$ETA[[3]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[3]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[3]]$SD.b)
      DominancexEnvironment_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$DxE_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      DxE_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      varENV<-scan('resultsETA_2_varB.dat')
      varDOM<-scan('resultsETA_3_varB.dat')
      varDxE<-scan('resultsETA_4_varB.dat')
      }
      if (Add & Gen & Env & Fix & Nenv>1 & Nfix>1) {additive_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$add_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      add_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      GeneticxEnvironment_variance[j,k]= mean(results$ETA[[5]]$varB)
      args$GxE_effect[,k]<-as.matrix(results$ETA[[5]]$b)
      GxE_SD_effect[,k]<-as.matrix(results$ETA[[5]]$SD.b)
      AdditivexEnvironment_variance[j,k]= mean(results$ETA[[6]]$varB)
      args$AxE_effect[,k]<-as.matrix(results$ETA[[6]]$b)
      AxE_SD_effect[,k]<-as.matrix(results$ETA[[6]]$SD.b)
      varENV<-scan('resultsETA_2_varB.dat')
      varGEN<-scan('resultsETA_3_varB.dat')
      varADD<-scan('resultsETA_4_varB.dat')
      varGxE<-scan('resultsETA_5_varB.dat')
      varAxE<-scan('resultsETA_6_varB.dat')
      }
      if (Dom & Gen & Env & Fix & Nenv>1 & Nfix>1) {dominance_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      GeneticxEnvironment_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$GxE_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      GxE_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      DominancexEnvironment_variance[j,k]= mean(results$ETA[[6]]$varB)
      args$DxE_effect[,k]<-as.matrix(results$ETA[[6]]$b)
      DxE_SD_effect[,k]<-as.matrix(results$ETA[[6]]$SD.b)
      varENV<-scan('resultsETA_2_varB.dat')
      varGEN<-scan('resultsETA_3_varB.dat')
      varDOM<-scan('resultsETA_4_varB.dat')
      varGxE<-scan('resultsETA_5_varB.dat')
      varDxE<-scan('resultsETA_6_varB.dat')
      }
      if (Dom & Add & Env & Fix & Nenv>1 & Nfix>1) {dominance_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      AdditivexEnvironment_variance[j,k]= mean(results$ETA[[5]]$varB)
      args$AxE_effect[,k]<-as.matrix(results$ETA[[5]]$b)
      AxE_SD_effect[,k]<-as.matrix(results$ETA[[5]]$SD.b)
      DominancexEnvironment_variance[j,k]= mean(results$ETA[[6]]$varB)
      args$DxE_effect[,k]<-as.matrix(results$ETA[[6]]$b)
      DxE_SD_effect[,k]<-as.matrix(results$ETA[[6]]$SD.b)
      varENV<-scan('resultsETA_2_varB.dat')
      varADD<-scan('resultsETA_3_varB.dat')
      varDOM<-scan('resultsETA_4_varB.dat')
      varAxE<-scan('resultsETA_5_varB.dat')
      varDxE<-scan('resultsETA_6_varB.dat')
      }
      if (Dom & Add & Gen & Fix & Nfix>1) {dominance_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      varGEN<-scan('resultsETA_2_varB.dat')
      varADD<-scan('resultsETA_3_varB.dat')
      varDOM<-scan('resultsETA_4_varB.dat')
      }
      if (Dom & Add & Gen & Env & Nenv>1) {dominance_variance[j,k]= mean(results$ETA[[4]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[4]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[4]]$SD.b)
      GeneticxEnvironment_variance[j,k]= mean(results$ETA[[5]]$varB)
      args$GxE_effect[,k]<-as.matrix(results$ETA[[5]]$b)
      GxE_SD_effect[,k]<-as.matrix(results$ETA[[5]]$SD.b)
      AdditivexEnvironment_variance[j,k]= mean(results$ETA[[6]]$varB)
      args$AxE_effect[,k]<-as.matrix(results$ETA[[6]]$b)
      AxE_SD_effect[,k]<-as.matrix(results$ETA[[6]]$SD.b)
      DominancexEnvironment_variance[j,k]= mean(results$ETA[[7]]$varB)
      args$DxE_effect[,k]<-as.matrix(results$ETA[[7]]$b)
      DxE_SD_effect[,k]<-as.matrix(results$ETA[[7]]$SD.b)
      varENV<-scan('resultsETA_1_varB.dat')
      varGEN<-scan('resultsETA_2_varB.dat')
      varADD<-scan('resultsETA_3_varB.dat')
      varDOM<-scan('resultsETA_4_varB.dat')
      varGxE<-scan('resultsETA_5_varB.dat')
      varAxE<-scan('resultsETA_6_varB.dat')
      varDxE<-scan('resultsETA_7_varB.dat')
      }
      if (Dom & Add & Gen & Env & Fix & Nenv>1 & Nfix>1) {dominance_variance[j,k]= mean(results$ETA[[5]]$varB)
      args$dom_effect[,k]<-as.matrix(results$ETA[[5]]$b)
      dom_SD_effect[,k]<-as.matrix(results$ETA[[5]]$SD.b)
      GeneticxEnvironment_variance[j,k]= mean(results$ETA[[6]]$varB)
      args$GxE_effect[,k]<-as.matrix(results$ETA[[6]]$b)
      GxE_SD_effect[,k]<-as.matrix(results$ETA[[6]]$SD.b)
      AdditivexEnvironment_variance[j,k]= mean(results$ETA[[7]]$varB)
      args$AxE_effect[,k]<-as.matrix(results$ETA[[7]]$b)
      AxE_SD_effect[,k]<-as.matrix(results$ETA[[7]]$SD.b)
      DominancexEnvironment_variance[j,k]= mean(results$ETA[[8]]$varB)
      args$DxE_effect[,k]<-as.matrix(results$ETA[[8]]$b)
      DxE_SD_effect[,k]<-as.matrix(results$ETA[[8]]$SD.b)
      varENV<-scan('resultsETA_2_varB.dat')
      varGEN<-scan('resultsETA_3_varB.dat')
      varADD<-scan('resultsETA_4_varB.dat')
      varDOM<-scan('resultsETA_5_varB.dat')
      varGxE<-scan('resultsETA_6_varB.dat')
      varAxE<-scan('resultsETA_7_varB.dat')
      varDxE<-scan('resultsETA_8_varB.dat')
      }
      
      prediction<-pred(args, mod="1", k=k)
      genotypic_value<-pred(args, mod="2", k=k)
      genetic_value<-pred(args, mod="3", k=k)
      statistics<-data.frame(results$fit)
      DIC[j,k]<-statistics[,4]
      
      varE<-scan('resultsvarE.dat')
      plot(varE,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Resildual variance");
      abline(h=results$varE,col=4,lwd=2);
      abline(v=results$burnIn/results$thin,col=4)
      
      if (Env & Nenv>1) {plot(varENV,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Environment variance");
        abline(h=environment_variance[j,k],col=4,lwd=2);
        abline(v=results$burnIn/results$thin,col=4)}
      
      if (Gen) {plot(varGEN,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Genetic variance");
        abline(h=genetic_variance[j,k],col=4,lwd=2);
        abline(v=results$burnIn/results$thin,col=4)}
      
      if (Add) {plot(varADD,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Additive variance");
        abline(h=additive_variance[j,k],col=4,lwd=2);
        abline(v=results$burnIn/results$thin,col=4)}
      
      if (Dom) {plot(varDOM,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "Dominance variance");
        abline(h=dominance_variance[j,k],col=4,lwd=2);
        abline(v=results$burnIn/results$thin,col=4)}
      
      if (Env & Gen & Nenv>1) {plot(varGxE,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "GxE variance");
        abline(h=GeneticxEnvironment_variance[j,k],col=4,lwd=2);
        abline(v=results$burnIn/results$thin,col=4)}
      
      if (Env & Add & Nenv>1) {plot(varAxE,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "AxE variance");
        abline(h=AdditivexEnvironment_variance[j,k],col=4,lwd=2);
        abline(v=results$burnIn/results$thin,col=4)}
      
      if (Env & Dom & Nenv>1) {plot(varDxE,type='o',col=2,cex=.5,ylab=expression(var[e]), main = "DxE variance");
        abline(h=DominancexEnvironment_variance[j,k],col=4,lwd=2);
        abline(v=results$burnIn/results$thin,col=4)}
      
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
      if (Fix & Nfix>1) {genetic_value_accuracy[j,k]<-cor(genetic_value, matrix(phenotest-args$fixtest%*%args$fix_effect[,k]))
      if (env & Nenv>1) {prediction_accuracy[j,k]<-cor(prediction, matrix(phenotest-args$fixtest%*%args$fix_effect[,k]))
      genotypic_value_accuracy[j,k]<-cor(genotypic_value, matrix(phenotest-args$fixtest%*%args$fix_effect[,k]))}
      } else {
        genetic_value_accuracy[j,k]<-cor(genetic_value, matrix(phenotest))
        if (env & Nenv>1) {prediction_accuracy[j,k]<-cor(prediction, matrix(phenotest))
        genotypic_value_accuracy[j,k]<-cor(genotypic_value, matrix(phenotest))}
      }
    }
    args$fixed_effects[,j]<-as.matrix(rowMeans(args$fix_effect, na.rm = TRUE))
    args$environment_effects[,j]<-as.matrix(rowMeans(args$env_effect, na.rm = TRUE))
    environment_SD_effects[,j]<-as.matrix(rowMeans(env_SD_effect, na.rm = TRUE))
    args$genetic_effects[,j]<-as.matrix(rowMeans(args$gen_effect, na.rm = TRUE))
    genetic_SD_effects[,j]<-as.matrix(rowMeans(gen_SD_effect, na.rm = TRUE))
    args$additive_effects[,j]<-as.matrix(rowMeans(args$add_effect, na.rm = TRUE))
    additive_SD_effects[,j]<-as.matrix(rowMeans(add_SD_effect, na.rm = TRUE))
    args$dominance_effects[,j]<-as.matrix(rowMeans(args$dom_effect, na.rm = TRUE))
    dominance_SD_effects[,j]<-as.matrix(rowMeans(dom_SD_effect, na.rm = TRUE))
    args$GeneticxEnvironment_effects[,j]<-as.matrix(rowMeans(args$GxE_effect, na.rm = TRUE))
    GeneticxEnvironment_SD_effects[,j]<-as.matrix(rowMeans(GxE_SD_effect, na.rm = TRUE))
    args$AdditivexEnvironment_effects[,j]<-as.matrix(rowMeans(args$AxE_effect, na.rm = TRUE))
    AdditivexEnvironment_SD_effects[,j]<-as.matrix(rowMeans(AxE_SD_effect, na.rm = TRUE))
    args$DominancexEnvironment_effects[,j]<-as.matrix(rowMeans(args$DxE_effect, na.rm = TRUE))
    DominancexEnvironment_SD_effects[,j]<-as.matrix(rowMeans(DxE_SD_effect, na.rm = TRUE))
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
  cat("-----Genotypic + GxE accuracy:", "\n") 
  print(genotypic_value_accuracy)
  cat("-----Mean of Genotypic + GxE accuracy:", "\n") 
  print(mean(genotypic_value_accuracy))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(genotypic_value_accuracy))) 
  cat("\n") 
  cat("-----Residual variance:", "\n") 
  print(residual_variance)
  cat("-----Mean of residual variance:", "\n") 
  print(mean(residual_variance))
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(rowMeans(residual_variance))) 
  cat("\n") 
  if (Fix & ncol(fixtrain)>1) {cat("-----Fixed effects:", "\n") 
    print(args$fixed_effects)
    cat("-----Mean of fixed effects:", "\n") 
    print(mean(args$fixed_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(args$fixed_effects))) 
    cat("\n")}
  if (Env & ncol(envtrain)>1) {cat("-----Environmental effects:", "\n") 
    print(args$environment_effects)
    cat("-----Mean of environmental effects:", "\n") 
    print(mean(args$environment_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(args$environment_effects)))
    cat("\n")
    cat("-----Environmental standard deviation:", "\n") 
    print(environment_SD_effects)
    cat("-----Mean of environmental standard deviation:", "\n") 
    print(mean(environment_SD_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(environment_SD_effects)))
    cat("\n")
    plot((rowMeans(args$environment_effects))^2, ylab='Estimated Squared-Environment Effect',
         type='o',cex=.5,col=4,main='Environment Effects')
    cat("-----Environmental variance:", "\n") 
    print(environment_variance)
    cat("-----Mean of environmental variance:", "\n") 
    print(mean(environment_variance))
    cat("\n")
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(environment_variance)))}
  if (Gen) {cat("-----Genetic effects:", "\n") 
    print(args$genetic_effects)
    cat("-----Mean of Genetic effects:", "\n") 
    print(mean(args$genetic_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(args$genetic_effects)))
    cat("\n")
    cat("-----Genetic standard deviation:", "\n") 
    print(genetic_SD_effects)
    cat("-----Mean of Genetic standard deviation:", "\n") 
    print(mean(genetic_SD_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(genetic_SD_effects)))
    cat("\n")
    plot((rowMeans(args$genetic_effects))^2, ylab='Estimated Squared-Genetic Effect',
         type='o',cex=.5,col=4,main='Genetic Effects')
    cat("-----Genetic variance:", "\n") 
    print(genetic_variance)
    cat("-----Mean of Genetic variance:", "\n") 
    print(mean(genetic_variance))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(genetic_variance)))
    cat("\n")}
  if (Add) {cat("-----Additive effects:", "\n") 
    print(args$additive_effects)
    cat("-----Mean of additive effects:", "\n") 
    print(mean(args$additive_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(args$additive_effects)))
    cat("\n")
    cat("-----Additive standard deviation:", "\n") 
    print(additive_SD_effects)
    cat("-----Mean of Additive standard deviation:", "\n") 
    print(mean(additive_SD_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(additive_SD_effects)))
    cat("\n")
    plot((rowMeans(args$additive_effects))^2, ylab='Estimated Squared-Additive Effect',
         type='o',cex=.5,col=4,main='Additive Effects')
    cat("-----Additive variance:", "\n") 
    print(additive_variance)
    cat("-----Mean of Additive variance:", "\n") 
    print(mean(additive_variance))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(additive_variance)))
    cat("\n")}
  if (Dom) {cat("-----Dominance effects:", "\n") 
    print(args$dominance_effects)
    cat("-----Mean of dominance effects:", "\n") 
    print(mean(args$dominance_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(args$dominance_effects)))
    cat("\n")
    cat("-----Dominance standard deviation:", "\n") 
    print(dominance_SD_effects)
    cat("-----Mean of dominance standard deviation:", "\n") 
    print(mean(dominance_SD_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(dominance_SD_effects)))
    cat("\n")
    plot((rowMeans(args$dominance_effects))^2, ylab='Estimated Squared-Dominance Effect',
         type='o',cex=.5,col=4,main='Dominance Effects')
    cat("-----Dominance variance:", "\n") 
    print(dominance_variance)
    cat("-----Mean of dominance variance:", "\n") 
    print(mean(dominance_variance))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(dominance_variance)))
    cat("\n")
  }
  if (Env & Gen & ncol(envtrain)>1) {cat("-----GeneticxEnvironment effects:", "\n") 
    print(args$GeneticxEnvironment_effects)
    cat("-----Mean of GeneticxEnvironment effects:", "\n") 
    print(mean(args$GeneticxEnvironment_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(args$GeneticxEnvironment_effects)))
    cat("\n")
    cat("-----GeneticxEnvironment standard deviation:", "\n") 
    print(GeneticxEnvironment_SD_effects)
    cat("-----Mean of GeneticxEnvironment standard deviation:", "\n") 
    print(mean(GeneticxEnvironment_SD_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(GeneticxEnvironment_SD_effects)))
    cat("\n")
    plot((rowMeans(args$GeneticxEnvironment_effects))^2, ylab='Estimated Squared-GxE Effect',
         type='o',cex=.5,col=4,main='GxE Effects')
    cat("-----GeneticxEnvironment variance:", "\n") 
    print(GeneticxEnvironment_variance)
    cat("-----Mean of GeneticxEnvironment variance:", "\n") 
    print(mean(GeneticxEnvironment_variance))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(GeneticxEnvironment_variance)))
    cat("\n")}
  if (Env & Add & ncol(envtrain)>1) {cat("-----AdditivexEnvironment effects:", "\n") 
    print(args$AdditivexEnvironment_effects)
    cat("-----Mean of AdditivexEnvironment effects:", "\n") 
    print(mean(args$AdditivexEnvironment_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(args$AdditivexEnvironment_effects)))
    cat("\n")
    cat("-----AdditivexEnvironment standard deviation:", "\n") 
    print(AdditivexEnvironment_SD_effects)
    cat("-----Mean of AdditivexEnvironment standard deviation:", "\n") 
    print(mean(AdditivexEnvironment_SD_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(AdditivexEnvironment_SD_effects)))
    cat("\n")
    plot((rowMeans(args$AdditivexEnvironment_effects))^2, ylab='Estimated Squared-AxE Effect',
         type='o',cex=.5,col=4,main='AxE Effects')
    cat("-----AdditivexEnvironment variance:", "\n") 
    print(AdditivexEnvironment_variance)
    cat("-----Mean of AdditivexEnvironment variance:", "\n") 
    print(mean(AdditivexEnvironment_variance))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(AdditivexEnvironment_variance)))
    cat("\n")}
  if (Env & Dom & ncol(envtrain)>1) {cat("-----DominancexEnvironment effects:", "\n") 
    print(args$DominancexEnvironment_effects)
    cat("-----Mean of DominancexEnvironment effects:", "\n") 
    print(mean(args$DominancexEnvironment_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(args$DominancexEnvironment_effects)))
    cat("\n")
    cat("-----DominancexEnvironment standard deviation:", "\n") 
    print(DominancexEnvironment_SD_effects)
    cat("-----Mean of DominancexEnvironment standard deviation:", "\n") 
    print(mean(DominancexEnvironment_SD_effects))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(DominancexEnvironment_SD_effects)))
    cat("\n")
    plot((rowMeans(args$DominancexEnvironment_effects))^2, ylab='Estimated Squared-DxE Effect',
         type='o',cex=.5,col=4,main='DxE Effects')
    cat("-----DominancexEnvironment variance:", "\n") 
    print(DominancexEnvironment_variance)
    cat("-----Mean of AdditivexEnvironment variance:", "\n") 
    print(mean(DominancexEnvironment_variance))
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(rowMeans(DominancexEnvironment_variance)))
    cat("\n")}
  cat("-----Deviance information criteria:", "\n") 
  print(DIC)
  cat("-----Mean of Deviance information criteria:", "\n") 
  print(mean(DIC))
  cat("\n")
  
  ##Calculating the genetic value for all families
  if (Fam) {family_genotypic_value= family.pred(args=args)
  cat("-----Genotypic value for each family:", "\n") 
  print(family_genotypic_value)
  cat("-----Highest Posterior Density Intervals:", "\n") 
  print(emp.hpd(family_genotypic_value[,2]))
  cat("\n")}
  
  ##Calculating the genetic value for all individuals
  genotypic_value<-ind.pred(args=args)
  for(l in 1:ncol(args$Env))
  {
    a<-1+((l-1)*ncol(args$fix)*ncol(args$Gen))
    b<-l*ncol(args$fix)*ncol(args$Gen)
    gv<-data.frame(genotypic_value[a:b,])
    pheno2<-data.frame(args$X[a:b,])
    cat("-----Individual selection - environment " ,l, ": ", "\n")
    rank=gv[order(gv[,7], decreasing = TRUE),]
    print(gv[order(gv[,7], decreasing = TRUE),])
    cat("-----Highest Posterior Density Intervals:", "\n") 
    print(emp.hpd(gv[,7]))
    write.table(gv[order(gv[,7], decreasing = TRUE),], paste("genetic_value", i, l), quote = FALSE, row.names = FALSE)
    write.table(emp.hpd(gv[,7]), paste("CI_genetic_value", i), quote = FALSE, row.names = FALSE)
    cat("\n")
    plot(gv[,7]~pheno2[,l+3],xlab='Observed',ylab='Predicted',col=2, main=paste('Environment ', l),
         xlim=c(min(pheno2[,l+3]), max(pheno[,l+3])),ylim=c(min(gv[,7]),max(gv[,7]))) 
    abline(a=0,b=1,col=4,lwd=2)
  }
  
  write.table(convergence, paste("convergence", i ,j ,k), quote = FALSE, row.names = FALSE)
  write.table(genetic_value_accuracy, paste("genetic_value_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(genetic_value_accuracy)), paste("CI_genetic_value_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(prediction_accuracy, paste("prediction_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(prediction_accuracy)), paste("CI_prediction_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(genotypic_value_accuracy, paste("genotypic_value_accuracy", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(rowMeans(genotypic_value_accuracy)), paste("CI_genotypic_value_accuracy", i), quote = FALSE, row.names = FALSE)
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
  write.table(family_genotypic_value, paste("family_genotypic_value", i), quote = FALSE, row.names = FALSE)
  write.table(emp.hpd(family_genotypic_value[,1]), paste("CI_family_genotypic_value", i), quote = FALSE, row.names = FALSE)
}

#Disconnect from the output file
dev.off()