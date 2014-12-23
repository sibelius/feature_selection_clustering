library(clusterGeneration)
library(mclust)

#################################################################################################################
# Generate datasets
#
# This script generates artificially datasets to compare different types of feature selection for clustering
#
# The design factors for my experiment are:
# 
# Number of Clusters (6)
#   category1 2,4,6
#   category2 12,14,16
#
# Separation Index (3)
#   0.01 close cluster
#   0.21 separated cluster
#   0.34 well-separated cluster
#
# Number of Non Noisy Variables (6)
#   category1 2,3,5
#   category2 22,23,24
#
# Number of Noisy Variables (6)
#   category1 0,2,4
#   category2 20,30,50
#
# Number of objects within the clusters (3)
#   balanced
#   10% one cluster - rest balanced
#   20% one cluster - rest balanced
#
# Number of replicates (3) - for statistical confidence
#
# Number of data sets generated:
# 6 x 3 x 6 x 6 x 3 x 3 = 2916
#################################################################################################################

# Factors design variables

# Number of clusters
numClust <- c(2,4,6,12,14,16)

# Separation index between clusters and their nearest neighboring clusters
# 0.01 (A=4) close cluster
# 0.21 (A=6) separated cluster
# 0.34 (A=8) well-separated cluster
sepVal <- c(0.01,0.21,0.34)

# Number of non noisy variables
#numNonNoisy <- c(2,3,5,22,24,26)
numNonNoisy <- c(22,24,26)


# Number of noisy variables
numNoisy <- c(0,1,3,20,50)


# Cluster size
# 1 - equal size
# 2 - randomly generated based on rangeN
# 3 - specified by the vector clustSizes
clustszind <- 1

# Number of replications (statistical confidence)
numReplicate <- 3

# Generate the data sets combining the factors
filename <- "dataset"

n <- prod(
        length(numClust),
        length(numNonNoisy),
        length(numNoisy),
        length(sepVal),
        numReplicate)
        
# Generate clusters
system.time(
datasets <- simClustDesign(
    numClust=numClust,
    sepVal=sepVal,
    numNonNoisy=numNonNoisy,
    numNoisy=numNoisy,
    clustszind=clustszind,
    numReplicate=numReplicate,
    fileName=filename)
)

# Generate a cluster size vector based on the percentage of one size of the first cluster, the rest is balanced
calculateClusterSize <- function(perc, numClust, numNonNoisy, numNoisy) {
    p <- numNonNoisy * numNoisy # total number of variables
    minCs <- (100/perc) * p
    restCs <- totalCs - minCs
    balancedCs <- restCs / (numClust-1)

    return(c(minCs, rep(balancedCs, nc-1)))
}

generate_plot <- function(numClust=2, sepVal=0.34, numNonNoisy=2, numNoisy=0) {
    tmp <- genRandomClust(
                numClust=numClust, 
                sepVal=sepVal,  
                numNonNoisy=numNonNoisy, 
                numNoisy=numNoisy, 
                numReplicate=1)
    X <- as.data.table(tmp$datList$test_1)
    Y <- tmp$memList$test_1

    model <- Mclust(X, numClust)
    plot(model, what="classification")

    print(externalValidity(model$classification, Y)) 
}

externalValidity <- function(classification, truth) {
    ARI <- adjustedRandIndex(classification, truth)
    CER <- classError(classification, truth)$errorRate
    return(data.table(ARI=ARI,CER=CER))
}
