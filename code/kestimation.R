######################################
# estimate the number of clusters (k)
######################################

# 1. look at elbow in the sum of squared error (SSE)
wss <- (nrow(X)-1) * sum(apply(X,2,var))
wss[2:15] <- sapply(2:15, function(k) sum(kmeans(X,k)$withinss))
wss <- data.table(wss)
ggplot(wss, aes(x=1:15,y=wss)) + geom_line() + geom_point()

# 2. pamk (medoids estimation) 
library(fpc)
library(cluster)
pamk.best <- pamk(X)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(X, pamk.best$nc))

asw <- numeric(20)
for (k in 2:20)
  asw[[k]] <- pam(X, k)$silinfo$avg.width
k.best <- which.max(asw)
cat("silhouette-optimal number of clusters:", k.best, "\n")

# 3. Calinky criterion
fit <- cascadeKM(scale(X, center = TRUE,  scale = TRUE), 1, 10, iter = 1000)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

# 4 . Bayesian Information Criterion
library(mclust)
# Run the function to see how many clusters
# it finds to be optimal, set it to search for
# at least 1 model and up 20.
d_clust <- Mclust(as.matrix(X), G=1:20)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
plot(d_clust)

# 5. Affine Propagation (AP) clustering
library(apcluster)
d.apclus <- apcluster(negDistMat(r=2), X)
cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")
heatmap(d.apclus)
plot(d.apclus, X)

# 6. Gap Statistic for Estimating the Number of Clusters 
clusGap(X, kmeans, 10, B = 100, verbose = interactive())

# 7. clustergrams

# 8. NbClust
nb <- NbClust(as.matrix(X), min.nc=2, max.nc=15, method = "kmeans", 
        index = "alllong", alphaBeale = 0.1)
nb <- NbClust(as.matrix(X), method="kmeans")
hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

############################
# Dendrogram visualization
############################

# 1. Hierarchical clustering
d_dist <- dist(as.matrix(X))   # find distance matrix 
plot(hclust(d_dist))           # apply hirarchical clustering and plot

# 2. Bayesian clustering method 
library(bclust)
x <- as.matrix(X)
d.bclus <- bclust(x, transformed.par = c(0, -50, log(16), 0, 0, 0))
viplot(imp(d.bclus)$var); plot(d.bclus); ditplot(d.bclus)
dptplot(d.bclus, scale = 20, horizbar.plot = TRUE,varimp = imp(d.bclus)$var, horizbar.distance = 0, dendrogram.lwd = 2)

# p-value hierarchical clustering
library(pvclust)
X.pv <- pvclust(X)
plot(X.pv)

