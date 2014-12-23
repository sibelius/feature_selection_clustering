# tiecorr calculates the correction for tied ranks

tiecorr <- function (rankarray) {
 tie3margsum <- 0
 ranksize <- dim(rankarray)
 for (rowindex in 1:ranksize[1]) {
  tie3rowsum <- 0
  rankindex <- 1
  while (rankindex <= ranksize[2]) {
   tiesum <- sum(rankindex == rankarray[rowindex,])
   if(tiesum > 1) tie3rowsum <- tie3rowsum + (tiesum^3 - tiesum)
   rankindex <- rankindex + 0.5
  }
  tie3margsum <- tie3margsum + tie3rowsum
 }
 return(tie3margsum)
}

# calculates a Z score for a zero-sum contrast on the rank array

BEZ <- function (rankarray,lambda) {
 ranksize <- dim(rankarray)
 L <- 0
 lambda2sum <- 0
 for (col in 1:ranksize[2]) {
  L<-L+lambda[col]*sum(rankarray[,col])
  lambda2sum<-lambda[col]^2 + lambda2sum
 }
 Z<-L/sqrt((ranksize[1]*ranksize[2]*(ranksize[2]+1)*lambda2sum)/12)
 Z<-Z*sqrt(ranksize[1]/tiecorr(rankarray))
 return(Z)
}

# computes Kendall's W from a matrix of either scores or ranks where
# rows are scoring or ranking methods and columns are data objects

kendall.w <- function (x,lambda=NULL,descending=TRUE,ranks=FALSE) {
 if(missing(x))
  stop("Usage: kendall.w(x,lambda=NULL,descending=TRUE,ranks=FALSE)")
 if(!is.data.frame(x) && !is.matrix(x))
  stop("x must be a dataframe or matrix")
 # lookup table for alpha=0.01 critical values for Kendall's W
 # lookup for small N and k starts at [3,3], so use offset of -2 to read
 Wcrit01<-matrix(
  c(NA,NA,NA,NA,NA,.522,.469,.425,.392,.359,.335,.311,.291,.274,.26,.245,.233,.221,
    NA,.768,.644,.553,.491,.429,.39,.351,.328,.306,.284,.262,.24,.227,.216,.204,.193,.182,
    .84,.683,.571,.489,.434,.379,.344,.309,.289,.269,.25,.23,.211,.2,.19,.18,.17,.16,
    .78,.629,.524,.448,.397,.347,.314,.282,.264,.246,.228,.21,.193,.183,.173,.164,.155,.146,
    .737,.592,.491,.419,.371,.324,.293,.263,.246,.229,.212,.195,.179,.169,.16,.152,.144,.136),
    nrow=5,byrow=TRUE)
 # lookup table for alpha=0.05 critical values for Kendall's W
 # lookup for small N and k starts at [3,3], so use offset of -2 to read
 Wcrit05<-matrix(
  c(NA,NA,NA,NA,NA,.376,.333,.3,.275,.25,.232,.214,.2,.187,.176,.166,.158,.15,
    NA,.619,.501,.421,.369,.318,.287,.256,.232,.222,.204,.188,.171,.160,.151,.143,.136,.129,
    .716,.552,.449,.378,.333,.287,.259,.231,.212,.195,.180,.167,.155,.147,.139,.131,.124,.117,
    .66,.512,.417,.351,.305,.267,.24,.215,.201,.186,.172,.158,.145,.137,.13,.123,.116,.109,
    .624,.484,.395,.333,.29,.253,.228,.204,.19,.176,.163,.15,.137,.13,.123,.116,.109,.103),
    nrow=5,byrow=TRUE)
 datadim<-dim(x)
 if(min(datadim) < 3) stop("x must be at least 3x3")
 if(is.null(colnames(x))) cnames<-as.character(1:datadim[2])
 else cnames<-colnames(x)
 if(!is.null(lambda)) max.lambda.len<-max(nchar(unlist(lambda)))
 else max.lambda.len<-4
 col.width<-max(nchar(cnames))
 if(col.width <= max.lambda.len) col.width<-max.lambda.len+1
 cnames<-formatC(cnames,width=col.width)
 if(ranks) {
  # check that all rows sum to the same value
  xsums<-rowSums(x)
  if(!all(xsums == xsums[1]))
   stop("Differing row sums - will not compute properly.")
  rank.mat<-x
 }
 else {
  meanscore<-sapply(x,mean)
  rank.mat <- t(as.matrix(x))
  if(descending) rank.mat <- max(rank.mat) - rank.mat
  exist.tie<-0
  for (i in 1:datadim[1]) rank.mat[,i]<-rank(rank.mat[,i])
  rank.mat <- t(rank.mat)
 }
 exist.tie<-length(unlist(apply(rank.mat,1,unique)))<length(rank.mat)
 meanranks<-apply(rank.mat,2,mean)
 grandmean<-mean(meanranks)
 if(exist.tie) {
  Tj<-tiecorr(rank.mat)
  W<-sum((meanranks-grandmean)^2)/
   ((datadim[2]*(datadim[2]^2-1)-Tj/datadim[1])/12)
 }
 else W<-sum((meanranks-grandmean)^2)/(datadim[2]*(datadim[2]^2-1)/12)
 if(datadim[2] > 7) {
  p.table<-NA
  x2df<-datadim[2]-1
  p.chisq<-pchisq(datadim[1]*(datadim[2]-1)*W,x2df,lower.tail=FALSE)
 }
 else {
  p.table<-ifelse(W > Wcrit01[datadim[2]-2,datadim[1]-2],"<0.01",
   ifelse(W > Wcrit05[datadim[2]-2,datadim[1]-2],"<0.05",">0.05"))
  x2df<-NA
  p.chisq<-NA
 }
 if(!is.null(lambda)) {
  ldim <- dim(lambda)
  if(is.null(ldim)) zstat<-round(BEZ(rank.mat,lambda),3)
  else {
   zstat <- vector("numeric",ldim[1])
   for (i in 1:ldim[1]) {
    zstat[i]<-round(BEZ(rank.mat,lambda[i,]),3)
   }
  }
 }
 else zstat<-NULL
 k.w<-list(W=W,p.table=p.table,p.chisq=p.chisq,x2df=x2df,rank.mat=rank.mat,
  cnames=cnames,meanranks=meanranks,lambda=lambda,zstat=zstat)
 class(k.w)<-"kendall.w"
 return(k.w)
}

print.kendall.w<-function(x,...) {
 cat("\nKendall's W for ordinal data\n")
 cat("W =",x$W)
 if(is.na(x$p.table)) {
  plabel<-paste("  p(X2[",x$x2df,"]) =",sep="",collapse="")
  cat(plabel,x$p.chisq,"\n\n")
 }
 else cat("  p(table) ",x$p.table,"\n\n")
 if(!is.null(x$zstat)) {
  col.width<-ifelse(is.null(x$cnames),8,max(nchar(x$cnames)))
  cat("Contrasts\n")
  cat(x$cnames)
  cat(paste(rep(" ",7),sep="",collapse=""))
  cat("Z\n")
  ldim <- dim(x$lambda)
  if(is.null(ldim))
   cat(formatC(x$lambda,width=col.width),
    rep(" ",8-nchar(x$zstat)),x$zstat,"\n")
  else {
   for(i in 1:length(x$zstat))
    cat(formatC(x$lambda[i,],width=col.width),
     rep(" ",8-nchar(x$zstat[i])),x$zstat[i],"\n")
  }
  cat("\n")
 }
}
