aggExCluster.matrix <- function(s, x, includeSim=FALSE)
{
    noPriorClustering <- (missing(x) || is.null(x))

    if (length(dim(s)) != 2 || (ncol(s) != nrow(s) && noPriorClustering))
        stop("'s' must be a square matrix")

    AggResultObj <- new("AggExResult")

    K <- nrow(s)

    AggResultObj@l <- K

    preserveNames <- (length(rownames(s)) == nrow(s))

    if (noPriorClustering) ## no prior clustering
    {
        AggResultObj@maxNoClusters <- K
        AggResultObj@clusters[[K]] <- as.list(1:K)
        AggResultObj@exemplars[[K]] <- 1:K

        if (preserveNames)
        {
            AggResultObj@labels <- rownames(s)
            names(AggResultObj@exemplars[[K]]) <- rownames(s)

            for (i in 1:K)
                names(AggResultObj@clusters[[K]][[i]]) <- rownames(s)[i]
        }
        else
            AggResultObj@labels <- as.character(1:K)
    }
    else ## prior clustering
    {
        if (x@l != nrow(s))
            stop("data set sizes of 's' and 'x' do not match")

        AggResultObj@sel <- x@sel

        K <- length(x@exemplars)

        if (K < 1)
            stop("'x' empty or corrupted")

        AggResultObj@maxNoClusters <- K
        AggResultObj@clusters[[K]] <- x@clusters
        AggResultObj@exemplars[[K]] <- x@exemplars
        AggResultObj@labels <- paste("Cluster", 1:K)
    }

    if (K < 2)
    {
        warning("there is nothing to cluster")
        return(invisible(AggResultObj))
    }

    objMat <- matrix(NA, K, K) ## matrix of objective values for pairs
    exeMat <- matrix(NA, K, K) ## matrix of joint exemplars
    ## note: only the upper triangle of these matrices is non-NA

    actClust <- AggResultObj@clusters[[K]]
    actExem <- AggResultObj@exemplars[[K]]
    actLabels <- -(1:K)

    AggResultObj@merge <- matrix(NA, K - 1, 2)
    AggResultObj@height <- rep(0, K - 1)

    res <- .Call("aggExClusterC", s, K, actClust, actExem, objMat, exeMat,
        actLabels, AggResultObj@sel, AggResultObj@clusters,
        AggResultObj@exemplars, AggResultObj@merge,
        AggResultObj@height, as.logical(preserveNames)[1])
    
    if (is.element("error", names(res)))
    {
        if (res$error == 1)
        {
            stop("clusters cannot be joined because of missing ",
                "similarity values;\n       maybe increasing the ",
                "cluster size through decreasing\n",
                "       the self similarity 'p' helps.")
        }
        else if (res$error == 2)
        {
            stop("clusters cannot be joined because of missing ",
                "similarity values")
        }
    }
    
    exeMat <- res$exeMat
    objMat <- res$objMat
    
    AggResultObj@merge <- res$merge
    AggResultObj@height <- res$height
    AggResultObj@clusters <- res$clusters
    
    if (length(AggResultObj@sel) > 0)
    {
        colInd <- res$colInd
    }

    ## finally, determine reordering for dendrogram plotting
    AggResultObj@order <- determineOrder(AggResultObj@merge,
        AggResultObj@height, K - 1)

    AggResultObj@call <- deparse(sys.call(-1))

    if (includeSim)
        AggResultObj@sim <- s

    AggResultObj
}

setMethod("aggExCluster", signature("matrix", "missing" ), aggExCluster.matrix)
setMethod("aggExCluster", signature("matrix", "ExClust" ), aggExCluster.matrix)


aggExCluster.Matrix <- function(s, x, includeSim=FALSE)
{
    if (is(s, "sparseMatrix"))
    {
        s <- as.SparseSimilarityMatrix(s)

        rng <- range(s@x)

        fill <- 2 * rng[1] - rng[2]

        s <- as.DenseSimilarityMatrix(s, fill=fill)
    }
    else
        s <- as.DenseSimilarityMatrix(s)

    if (missing(x))
        res <- aggExCluster(s, includeSim=includeSim)
    else
        res <- aggExCluster(s, x, includeSim=includeSim)

    res
}

setMethod("aggExCluster", signature("Matrix", "missing" ), aggExCluster.Matrix)
setMethod("aggExCluster", signature("Matrix", "ExClust" ), aggExCluster.Matrix)


aggExCluster.Clust <- function(s, x, includeSim=TRUE)
{
    if (all(dim(x@sim) <= 1))
        stop("similarity matrix not included in object")

    AggResultObj <- aggExCluster(x@sim, x)

    AggResultObj@call <- deparse(sys.call(-1))

    if (includeSim)
        AggResultObj@sim <- x@sim

    AggResultObj
}

setMethod("aggExCluster", signature("missing" , "ExClust" ), aggExCluster.Clust)


aggExCluster.function <- function(s, x, includeSim=TRUE, ...)
{
    if (is.data.frame(x))
        x <- as.matrix(x[, sapply(x, is.numeric)])

    if (is.matrix(x))
        N <- nrow(x)
    else
        N <- length(x)

    if (N < 2) stop("cannot cluster less than 2 samples")

    if (!is.function(s))
    {
        if (!is.character(s) || !exists(s, mode="function"))
            stop("invalid distance function")

        s <- match.fun(s)
    }

    sim <- s(x=x, ...)

    if (!is.matrix(sim) || (nrow(sim) != N) || ncol(sim) != N)
        stop("computation of similarity matrix failed")

    AggResultObj <- aggExCluster(sim)

    AggResultObj@call <- deparse(sys.call(-1))

    if (includeSim)
        AggResultObj@sim <- sim

    AggResultObj
}

setMethod("aggExCluster", signature("function" , "ANY"), aggExCluster.function)
setMethod("aggExCluster", signature("character", "ANY"), aggExCluster.function)


## auxiliary function for determining the order for dendrogram plotting
## fills up order recursively starting from the last merge
determineOrder <- function(merge, height, k)
{
    I <- merge[k, 1] ## I and J are the clusters merged in the k-th step
    J <- merge[k, 2]

    if (I < 0 && J < 0) ## if both are singletons, list I first
        return(c(-I, -J))
    else if (I < 0) ## if I is a singleton and J is not, list it first
        return(c(-I, determineOrder(merge, height, J)))
    else if (J < 0) ## if J is a singleton and I is not, list it first
        return(c(-J, determineOrder(merge, height, I)))
    else ## if both are non-singleton clusters, list the "tighter" cluster
    {    ## on the left-hand side (see ?hclust)
        if (height[I] > height[J])
            return(c(determineOrder(merge, height, I),
                     determineOrder(merge, height, J)))
        else
            return(c(determineOrder(merge, height, J),
                     determineOrder(merge, height, I)))
    }
}
