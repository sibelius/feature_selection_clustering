library(parallel)
library(data.table)

##------------------------------------------------------------------------------
##' Wrapper around mclapply to track progress
##' 
##' Based on http://stackoverflow.com/questions/10984556
##' 
##' @param X         a vector (atomic or list) or an expressions vector. Other
##'                  objects (including classed objects) will be coerced by
##'                  ‘as.list’
##' @param FUN       the function to be applied to
##' @param ...       optional arguments to ‘FUN’
##' @param mc.preschedule see mclapply
##' @param mc.set.seed see mclapply
##' @param mc.silent see mclapply
##' @param mc.cores see mclapply
##' @param mc.cleanup see mclapply
##' @param mc.allow.recursive see mclapply
##' @param mc.progress track progress?
##' @param mc.style    style of progress bar (see txtProgressBar)
##'
##' @examples
##' x <- mclapply2(1:1000, function(i, y) Sys.sleep(0.01))
##' x <- mclapply2(1:3, function(i, y) Sys.sleep(1), mc.cores=1)
##------------------------------------------------------------------------------
mclapply2 <- function(X, FUN, ..., 
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
    mc.cleanup = TRUE, mc.allow.recursive = TRUE,
    mc.progress=TRUE, mc.style=3) 
{
    if (!is.vector(X) || is.object(X)) X <- as.list(X)

    if (mc.progress) {
        f <- fifo(tempfile(), open="w+b", blocking=T)
        p <- parallel:::mcfork()
        pb <- txtProgressBar(0, length(X), style=mc.style)
        setTxtProgressBar(pb, 0) 
        progress <- 0
        if (inherits(p, "masterProcess")) {
            while (progress < length(X)) {
                readBin(f, "double")
                progress <- progress + 1
                setTxtProgressBar(pb, progress) 
            }
            cat("\n")
            parallel:::mcexit()
        }
    }
    tryCatch({
        result <- mclapply(X, function(...) {
                res <- FUN(...)
                if (mc.progress) writeBin(1, f)
                res
            }, 
            mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
            mc.silent = mc.silent, mc.cores = mc.cores,
            mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
        )

    }, finally = {
        if (mc.progress) close(f)
    })
    result
}

plot.kmeans <- 
function (x, data = NULL, Y = NULL, legend.position = c("right", 
    "bottom", "left", "top", "none"), title = "K-Means Results", 
    xlab = "Principal Component 1", ylab = "Principal Component 2", 
    ...) 
{
    toPlot <- fortify(model = x, data = data)
    legend.position <- match.arg(legend.position)
    
    g <- ggplot(toPlot, aes_string(x = ".x", y = ".y", colour = ".Cluster"))
    if (!is.null(Y)) { 
        toPlot$Y <- factor(Y)
        g <- g + geom_point(aes_string(shape = Y)) + scale_shape_identity()
    } else {
        g <- g + geom_point()
    }

    g <- g + scale_color_discrete("Cluster") + 
        theme(legend.position = legend.position) + labs(title = title, 
        x = xlab, y = ylab)
    print(g)
    return(g)
}

clusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}


# read generated datasets
files <- list.files(pattern=".dat")

getClusterInfo <- function(fileName) {
    info <- gsub("dataset|.dat","",fileName)
    sepVal <- substr(info, 2, 2)

    nc_i <- regexpr("G",info) + 1
    nc_e <- regexpr("v",info) - 1
    numClust <- as.numeric(substr(info,nc_i,nc_e))

    nnn_i <- nc_e + 2
    nnn_e <- regexpr("n",info)-1
    numNonNoisy <- as.numeric(substr(info,nnn_i,nnn_e))

    nn_i <- nnn_e + 3
    nn_e <- regexpr("o",info)-1
    numNoisy <- as.numeric(substr(info,nn_i,nn_e))

    numReplicate <- as.numeric(substr(info, nchar(info), nchar(info)))
    return(data.table(sepVal,numClust,numNonNoisy,numNoisy,numReplicate)) 
}

getInfo <- function(files) {
    info <- mclapply2(files, getClusterInfo, mc.cores=4)
    info <- rbindlist(info)
    return(info)
}

loadDatasets <- function(files) {
    datasets <- mclapply2(files, function(fileName) {
        fileNameMem <- gsub(".dat",".mem",fileName)
        fileNameNoisy <- gsub(".dat",".noisy",fileName)
        
        X <- fread(fileName)
        y <- fread(fileNameMem)$V1  
        noisy <- fread(fileNameNoisy)$V1

        info <- getClusterInfo(fileName)

        return(list(X=X,y=y,noisy=noisy,info=info))
    }, mc.cores=4)
    return(datasets)
}

lapply_pb <- function(X, FUN, ...)
{
 env <- environment()
 pb_Total <- length(X)
 counter <- 0
 pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   

 # wrapper around FUN
 wrapper <- function(...){
   curVal <- get("counter", envir = env)
   assign("counter", curVal +1 ,envir=env)
   setTxtProgressBar(get("pb", envir=env), curVal +1)
   FUN(...)
 }
 res <- lapply(X, wrapper, ...)
 close(pb)
 res
}
