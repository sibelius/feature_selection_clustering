##############################
# X - nxp matrix of data
# y - nx1 vector of classes
#
#############################
library(data.table)
library(ggplot2)
library(LICORS) # kmeanspp -> kmeans++
library(mclust) # adjustedRandIndex, classError, Mclust
library(FSelector) # forward.selection
library(sparcl) # sparse k-means
library(clv) # jaccard index
library(clustvarsel)
library(vscc)

#load("datasets/crabs.Rdata")
#X <- data.table(crabs[,4:8])
#y <- crabs$class

#k <- length(unique(y))

externalValidity <- function(classification, truth) {
    ari <- adjustedRandIndex(classification, truth)
    cer <- classError(classification, truth)$errorRate
    std <- std.ext(classification, truth)

    jaccard <- clv.Jaccard(std)
    return(data.table(ari=ari,cer=cer,jaccard=jaccard))
}

experiment.wrapper <- function(dataset) {
    X <- dataset$X
    y <- dataset$y
    info <- dataset$info
    k <- info$numClust

    t0 <- proc.time()
    tmp <- experiment.simplifiedSilhouette(X,y,k)
    result <- tmp$ev
    subset <- length(tmp$atts)
    t1 <- proc.time()
    print(t1-t0)
    
    results_info <- cbind(dataset$info,result,subset, "simplifiedSilhouette")
    return(results_info)
}

experiment.w.clustvarsel <- function(dataset) {
    X <- dataset$X
    y <- dataset$y
    info <- dataset$info
    k <- info$numClust

    t0 <- proc.time()
    tmp <- experiment.clustvarsel(X,y,k)
    result <- tmp$ev
    subset <- length(tmp$atts)
    t1 <- proc.time()
    print(t1-t0)

    results_info <- cbind(dataset$info,result,subset, "clustvarsel")
    return(results_info)
}

# run all the algorithms in a dataset
# return a data.table with info of dataset and external validation index of algorithms, also the number of selected attributes
experiment.all <- function(dataset) {
    X <- dataset$X
    y <- dataset$y
    info <- dataset$info
    k <- info$numClust

    algorithms <- c(
        "kmeans",
        "modelbased",
        "sparseKmeans",
        "clustvarsel",
        "vscc"
        #"simplifiedSilhouette"
    )
    allVars <- info$numNonNoisy + info$numNoisy
    subset <- c(allVars, allVars)
    runTime <- c()

    results <- list()
    print("kmeans")
    t0 <- proc.time()
    results[[1]] <- experiment.kmeans(X,y,k)
    t1 <- proc.time()
    runTime <- c(runTime, t1[3]-t0[3])

    print("modelbased")
    t0 <- proc.time()
    results[[2]] <- experiment.modelbased(X,y,k)
    t1 <- proc.time()
    runTime <- c(runTime, t1[3]-t0[3])

    print("sparseKmeans")
    t0 <- proc.time()
    tmp <- experiment.sparseKmeans(X,y,k)
    t1 <- proc.time()
    runTime <- c(runTime, t1[3]-t0[3])
    
    results[[3]] <- tmp$ev
    subset <- c(subset, length(tmp$atts))
    
    print("clustvarsel")
    t0 <- proc.time()
    tmp <- experiment.clustvarsel(X,y,k)
    t1 <- proc.time()
    runTime <- c(runTime, t1[3]-t0[3])
   
    results[[4]] <- tmp$ev
    subset <- c(subset, length(tmp$atts))

    print("vscc")
    t0 <- proc.time()
    tmp <- experiment.vscc(X,y,k)
    t1 <- proc.time()
    runTime <- c(runTime, t1[3]-t0[3])
    
    results[[5]] <- tmp$ev
    subset <- c(subset, length(tmp$atts))

    #tmp <- experiment.simplifiedSilhouette(X,y,k)
    #results[[6]] <- tmp$ev
    #subset <- c(subset, length(tmp$atts))
    
    results_all <- rbindlist(results)
    results_info <- cbind(dataset$info,results_all,subset, algorithms, runTime)
    return(results_info)
}

experiment.nClustNotKnown <- function(dataset) {
    X <- dataset$X
    y <- dataset$y
    info <- dataset$info
    k <- info$numClust

    algorithms <- c(
        "clustvarsel",
        "vscc"
        #"simplifiedSilhouette"
    )
    subset <- c()
    clusterSelected <- c()

    t0 <- proc.time()
    results <- list()
    
    print("clustvarsel")
    tmp <- experiment.clustvarsel(X,y)
    results[[1]] <- tmp$ev
    subset <- c(subset, length(tmp$atts))
    clusterSelected <- c(clusterSelected,tmp$numClust)

    print("vscc")
    tmp <- experiment.vscc(X,y)
    results[[2]] <- tmp$ev
    subset <- c(subset, length(tmp$atts))
    clusterSelected <- c(clusterSelected,tmp$numClust)

    #tmp <- experiment.simplifiedSilhouette(X,y,k)
    #results[[6]] <- tmp$ev
    #subset <- c(subset, length(tmp$atts))
    
    t1 <- proc.time()
    print(t1-t0)
    
    results_all <- rbindlist(results)
    results_info <- cbind(dataset$info,results_all,subset, clusterSelected, algorithms)
    return(results_info)
}




################################
# All Features
###############################


# kmeans++
experiment.kmeans <- function(X, y, k) {
    km <- kmeanspp(X, k)
    externalValidity(km$cluster, y)
}

# model based - EM
experiment.modelbased <- function(X, y, k) {
    model <- Mclust(X, k)
    externalValidity(model$classification, y)
}

experiment.modelbasedBIC <- function(X, y) {
    model <- Mclust(X,2:9)
    ev <- externalValidity(model$classification,y)
    return(list(ev=ev,numClust=model$G))
}

###################################
# Feature Selection for Clustering
###################################
source("util.R")
source("simplifiedSilhouette.R")
experiment.simplifiedSilhouette <- function(X, y, k=NULL) {
    if (is.null(k)) {
        minK = 2
        maxK = 9
    } else {
        minK = k
        maxK = k
    }

    bestK <- k
    atts <- forward.search(names(X), function(att) {
        #SSEParallel(X, att, minK, maxK)  
        tmp <- SSEvaluator(X, att, minK, maxK)
        bestK <- tmp$bestK
        return(tmp$result)
    } )
    #atts <- backward.search(names(X), function(att) {
    #    SSEParallel(X, att, minK, maxK)
    #} )

    Xs <- X[,atts,with=FALSE]
    km <- kmeanspp(Xs,k)
    ev <- externalValidity(km$cluster, y)
    return(list(ev=ev,atts=atts))
}

experiment.sparseKmeans <- function(X, y, k) {
    x <- scale(X, TRUE, TRUE)
    # choose tuning parameter
    km.perm <- KMeansSparseCluster.permute(x,K=k,wbounds=seq(3,7,len=15),nperms=5, silent=TRUE)

    km <- KMeansSparseCluster(x, k,wbounds=km.perm$bestw)    

    atts <- names(which(km[[1]]$ws!=0))

    ev <- externalValidity(km[[1]]$Cs, y)
    return(list(ev=ev,atts=atts))
}

experiment.clustvarsel <- function(X, y, k=NULL) {
    if (is.null(k)) {
        k <- 2:9
    }

    c <- clustvarsel(X, k)
    
    if(is.null(c$subset)) {
        Xs <- X
        atts <- names(X)
    }
    else {
        atts <- names(c$subset)
        Xs <- X[,atts,with=FALSE]
    }
    model <- Mclust(Xs, k)

    ev <- externalValidity(model$classification, y)
    return(list(ev=ev,atts=atts,numClust=model$G))
}

experiment.vscc <- function(X, y, k=NULL) {
    if (is.null(k)) {
        k <- 2:9
    }

    vs <- vscc(X, k)

    ev <- externalValidity(vs$bestmodel$classification, y)
    return(list(ev=ev,atts=names(vs$topselected),numClust=vs$bestmodel$G))
}


