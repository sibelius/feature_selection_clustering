# My implementation of the simplifiedSilhouette
simplifiedSilhouette <- function(data, centers, cluster, distFun) {
    silh = sapply(1:nrow(data), function(i) {
        a <- distFun(data[i,], centers[cluster[i],])

        b <- min(sapply(1:nrow(centers), function(k) {
            if(k == cluster[i])
                return(Inf)
            else
                return(distFun(data[i,], centers[k,]))
        } ))

        return( (b-a) / max(a,b) )
    } )

    return(mean(silh))
}

euclidian <- function(x1, x2) {
    dist(rbind(as.numeric(x1),as.numeric(x2)))[1]
}

SSEvaluator <- function(data, att, minK=2, maxK=6) {
    data <- data[,att,with=FALSE]
    
    numberRepetitions <- 10
    bestMerit <- -Inf

    pbTotal <- length(minK:maxK ) * numberRepetitions 
    pb <- txtProgressBar(min = 0, max = pbTotal, style = 3)
    pbK <- 0

    for(k in minK:maxK ) {
        for(rpt in 1:numberRepetitions ) {
            km <- kmeans(data,k)

            silh <- simplifiedSilhouette(data, km$centers, km$cluster, euclidian)
            
            if (silh > bestMerit) {
                bestMerit <- silh
                bestPartition <- km
                bestK <- k
            }
            pbK <- pbK + 1
            setTxtProgressBar(pb, pbK)
        }
    }
    #print(att)
    #print(bestK)
    #print(bestMerit)

    return(list(result=bestMerit,bestK=bestK))
}

SSEParallel <- function(data, att, minK = 2, maxK = 6) {
    data <- data[,att, with=FALSE]
    
    numberRepetitions <- 10

    bestMerit <- mclapply(minK:maxK, function(k) {

        bestK <- mclapply(1:numberRepetitions, function(x) {
            km <- kmeans(data, k)

            silh <- simplifiedSilhouette(data, km$centers, km$cluster, euclidian)

            return(list(silh=silh,km=km))
        }, mc.cores=4 )

        silh <- sapply(bestK, function(x) x$silh)
        silhMax <- max(silh)
        return(list(silh=silhMax,km=bestK[[which.max(silh)]]))
    }, mc.cores=4 )

    silh <- sapply(bestMerit, function(x) x$silh)
    silhMax <- max(silh)
    i <- which.max(silh)
    bestK <- minK + i - 1
    bestPartition <- bestMerit[[i]]$km
        
    #print(att)
    #print(bestK)
    #print(silhMax)

    return(silhMax)
}
