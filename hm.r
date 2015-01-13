# Code for the paper 'Enriching the analysis of genome-wide association studies with hierarchical modeling' AJHG 2007. 

# This function implements a two-stage hierarchical model using a second-stage ML weighted least squares

# firstStageFilename is the name of the file containing three columns where the first column contains the marker IDs, second column contains first stage coefficients, and third contains standard errors from linear or logistic regression.
# covariatesFilename is the name of the file containing the covariates for use in the hierarchical model where the number of rows should match the number of rows in firstStageFilename
# Note that neither file above should contain a header row.
# weightingCol is the index of the column from the covariatesFilename that should be used in conjunction with rho in determining SNP specific shrinkage values.
# tau represents the assumed second-stage residual standard deviations at baseline (i.e. SNPs with no prior evidence). 
# rho represents the assumed residual standard deviation for SNPs with maxium prior evidence.  
# Lower values of tau or rho increase the degree of shrinkage of first stage estiamtes towards second stage estimates.

hm<-function(firstStageFilename,covariatesFilename,weightingCol,tau,rho){
  print("Loading first stage data")
  firstStage<-read.table(firstStageFilename)
  markers<-firstStage[,1]
  betas<-abs(firstStage[,2])
  se<-firstStage[,3]
  print("Loading Z matrix")
  Z<-as.matrix(read.table(covariatesFilename))
  # replace NAs with 0
  Z<-replace(Z,is.na(Z),0)
  rowLength<-length(markers)
  prior_scaling<-(log(tau^2)-log(rho^2))/max(Z[,weightingCol])
  # Compute the weights 
  shrinkage<-1/exp(prior_scaling*Z[,weightingCol])
  # for now defined the weight as the linkage column
  wt<-(se[]^2+tau^2*shrinkage)^-1
  lmsummary<-summary(lm(betas~Z,,,weights=wt))
  outputVector<-c(lmsummary$coefficients[,1],lmsummary$coefficients[,2])
  write.table(matrix(outputVector,length(Z[1,])),file=paste("hm_tau",tau,"_rho",rho,"_coefficients.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE)
  second.betas<-Z %*% lmsummary$coefficients[,1]
  Z_transpose<-t(Z)
  # The following line of code assume that the second stage matrix T is 
  # diagonal in order to speed up computations for vast numbers of SNPs
  zprimewmatrix<-apply(as.matrix(Z_transpose),1,function(x,weights) x*weights,weights=wt)
  secondVarFirstHalf<-Z %*% solve(t(zprimewmatrix) %*% Z) 
  second.Var<-apply(t(secondVarFirstHalf)*Z_transpose,2,sum)
  second.SE<-sqrt(second.Var)
  B<-(se^2)*wt
  hm.betas<-B*second.betas + (1-B)*betas
  H<-second.Var * wt
  C<-se^2 * (1- (1-H) * B )
  hm.SE<-sqrt(C)
  outputVector<-c(markers,hm.betas,hm.SE,second.betas,second.SE)
  write.table(matrix(outputVector,rowLength),file=paste("hm_tau",tau,"_rho",rho,"_hmresults.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE)
}