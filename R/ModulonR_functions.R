# Define constants for distance and clustering methods
DISTANCE_METHODS <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
CLUSTERING_METHODS <- c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")

# Load required packages quietly
suppressMessages(require(stats))
suppressMessages(require(parallel))
suppressMessages(require(ropls))
suppressMessages(require(operators))

# Function to get hierarchical clusters
#' @title Get Hierarchical Clusters
#' @description
#' This function computes the hierarchical clustering tree of the features in the rows of a matrix 
#' using the `stats::hclust()` method. 
#' @param mat The input matrix for which the features in the rows will be clustered.
#' @param distance.method The distance method to use in `stats::dist()`.
#' @param clustering.method The hierarchical clustering method to use in `stats::hclust()`. Default is "complete" (see ?stats::hclust for more information).
#' @param scale Logical argument indicating whether the rows of the input matrix should be scaled (mean-centered and divided by the standard deviation). Default is TRUE.
#' @param ... Additional parameters for `stats::hclust`.
#' @return An object of class hclust (see ?stats::hclust).
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   # Generate a random matrix with 100 rows and 5 columns
#'   set.seed(123)
#'   random_matrix <- matrix(abs(rnorm(25 * 5)), nrow = 25, ncol = 5)
#'   rownames(random_matrix) <- paste0("F", 1:25)
#'   plot(GetHierarchicalClusters(mat = random_matrix))
#'   }
#' }
#' @rdname GetHierarchicalClusters
#' @export 
GetHierarchicalClusters <- function(mat, 
                                    distance.method = DISTANCE_METHODS, 
                                    clustering.method = CLUSTERING_METHODS,
                                    scale = TRUE,
                                    ...) {
  # Set parameters and sanity checks
  distance.method <- match.arg(distance.method)
  clustering.method <- match.arg(clustering.method)
  
  if (scale) {
    # Scale matrix by rows
    mat.mean <- rowMeans(mat, na.rm = TRUE)
    mat.sd <- apply(mat, 1, sd, na.rm = TRUE)
    mat <- sweep(mat, 1, mat.mean, "-")
    mat <- sweep(mat, 1, mat.sd, "/")
  }
  
  # Get the distance matrix
  mat.dist <- dist(mat, method = distance.method)
  
  # Get the hclust results
  hclust.res <- hclust(mat.dist, method = clustering.method, ...)
  
  return(hclust.res)
}

#' @title Find Feature Clusters
#' @description
#' This function takes a matrix and groups similar features in the rows into clusters through hierarchical clustering. 
#' The total number of clusters can be set using the parameter `cluster.nums`.
#' The function allows the samples (columns) to be aggregated by a group annotation before clustering (see parameter `annotation`).
#' @param mat A matrix with features in rows and samples in columns.
#' @param annotation A character vector containing the annotation of the samples (columns) for data aggregation. 
#' Aggregation is performed using a simple mean. If NULL, no aggregation is performed. Default is NULL.
#' @param cluster.nums A numerical vector specifying the number of clusters. Default is 2:10.
#' @param distance.method The distance method to use in `stats::dist()`.
#' @param clustering.method The hierarchical clustering method to use in `stats::hclust()`. Default is "complete" (see ?stats::hclust for more information).
#' @param ... Additional parameters for `stats::hclust()` (see also GetHierarchicalClusters()).
#' @return A list containing the feature cluster composition. 
#' The names of the list indicate the cluster ID in the following structure:
#' name = <one total number of clusters in cluster.nums>"."<cluster number>
#' The values are lists of genes, i.e., the cluster's signature.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   # Generate a random matrix with 100 rows and 5 columns
#'   set.seed(123)
#'   random_matrix <- matrix(abs(rnorm(25 * 5)), nrow = 25, ncol = 5)
#'   rownames(random_matrix) <- paste0("F", 1:25)
#'   signature.list <- FindFeatureCluster(mat = random_matrix, annotation = c('a','a','a','b','b'))
#'   signature.list[c(2:5)]
#'   }
#' }
#' @rdname FindFeatureCluster
#' @export 
FindFeatureCluster <- function(mat, 
                               annotation = NULL, 
                               cluster.nums = 2:10, 
                               distance.method = DISTANCE_METHODS, 
                               clustering.method = CLUSTERING_METHODS,
                               ...) {
  # Set parameters and sanity checks
  distance.method <- match.arg(distance.method)
  clustering.method <- match.arg(clustering.method)
  cluster.nums <- as.numeric(cluster.nums)
  
  # Aggregate the data if needed
  if (!is.null(annotation)) {
    mat <- aggregate(t(mat),
                     by = list(annotation),
                     FUN = mean)
    rownames(mat) <- mat[, 1]
    mat <- mat[, -1]
    mat <- t(mat)
  }
  
  # Get the hierarchical cluster tree
  hclust.res <- GetHierarchicalClusters(
    mat = mat, 
    distance.method = distance.method, 
    clustering.method = clustering.method, 
    ...)
  
  # Get the cluster(s) identity for every k
  cluster.res <- cutree(hclust.res, k = cluster.nums)
  cluster.res <- data.frame(cluster.res)
  colnames(cluster.res) <- as.character(cluster.nums)
  
  # Convert to a simple 1 level list
  markers.res  <- lapply(cluster.res, function(x) split(rownames(cluster.res), x))
  markers.res <- unlist(markers.res, recursive = FALSE)
  
  return(markers.res)
}

# Wrapper function to calculate signature AUC
#' @title Get Signature AUC
#' @description
#' This function computes the AUC signature score of the signatures in a 
#' signature list using the `AUCell::AUC()` function. 
#' It accepts a list of signatures (see parameter `signature.list`). 
#' Note that this list can be obtained using the function `FindFeatureCluster()`.
#' The input data can be either a matrix with features in rows and samples in columns or the 
#' "rankings" for each sample as defined by `AUCell::AUCell_buildRankings()`.
#' @param signature.list A named list of the signatures.
#' @param mat A matrix with features in rows and samples in columns.
#' @param rankings Cells ranking as defined by `AUCell::AUCell_buildRankings()`. Default is NULL.
#' @param scale A logical value indicating whether to scale the data between 0 and 1. Default is FALSE.
#' @return A dataframe containing the signature score (AUC) for all the signatures in `signature.list`, with signatures in rows and samples in columns.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   # Generate a random matrix with 100 rows and 5 columns
#'   set.seed(123)
#'   random_matrix <- matrix(abs(rnorm(25 * 5)), nrow = 25, ncol = 5)
#'   rownames(random_matrix) <- paste0("F", 1:25)
#'   signature.list <- FindFeatureCluster(mat = random_matrix, annotation = c('a','a','a','b','b'))
#'   GetSignatureAUC(signature.list = signature.list, mat = random_matrix)
#'   }
#' }
#' @seealso 
#'  \code{\link[AUCell]{AUCell_buildRankings}}, \code{\link[AUCell]{AUCell_calcAUC}}, \code{\link[AUCell]{aucellResults-class}}
#' @rdname GetSignatureAUC
#' @export 
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
GetSignatureAUC <- function(signature.list, mat, rankings = NULL, scale = FALSE) {
  # Compute rankings if not provided
  if (is.null(rankings)) {
    if (!is.null(mat)) {
      rankings <- AUCell::AUCell_buildRankings(exprMat = mat, plotStats = FALSE, verbose = FALSE)
    } else {
      stop("Must provide input data matrix or rankings")
    }
  }
  
  # Calculate AUC
  cells.AUC <- AUCell::AUCell_calcAUC(geneSets = signature.list, 
                                      rankings = rankings, 
                                      aucMaxRank = nrow(rankings),
                                      normAUC = FALSE, 
                                      verbose = FALSE)
  
  # Get AUC
  signature.score <- AUCell::getAUC(cells.AUC)
  
  # Create result data.frame
  signature.score <- data.frame(signature.score, check.names = FALSE)
  
  # Scaling if required
  if (scale) {
    signature.score <- t(apply(signature.score, 1, scale_range))
    colnames(signature.score) <- colnames(signature.score)
  }
  
  return(signature.score)
}

#' @title Discriminant Analysis
#' @description This function computes the discriminant analysis (OPLS-DA from the `ropls` library), comparing every class in a categorical variable against the background (i.e., the other classes).
#' @param mat Matrix with features and samples in rows and columns, respectively.
#' @param method Character; the discriminant analysis method. Possible values are: "OPLS" (see `FindMarkersOPLS()`).
#' @param annotation Character vector with the class of the samples.
#' @param BackgroundClasses Classes to compare with the query class. If NULL, all the classes in annotation except the query class will be included as background. Default: NULL
#' @param QueryClasses Character vector with the classes to compare with respect to the background. If NULL, all the classes in annotation will be analyzed (separately) with respect to the background. Default: NULL
#' @return Data frame with the discriminant analysis results for each analyzed class.
#' @rdname DiscriminantAnalysis
#' @export 
DiscriminantAnalysis <- function(mat, method = "OPLS", annotation, BackgroundClasses=NULL, QueryClasses=NULL) {
  if (is.null(BackgroundClasses)){
    BackgroundClasses = unique(annotation)
  }
  if (is.null(QueryClasses)){
    QueryClasses = BackgroundClasses
  }
  
  TRUTH.vector = annotation %in% c(BackgroundClasses,QueryClasses)
  
  mat = mat[,TRUTH.vector]
  annotation = annotation[TRUTH.vector]
  
  
  compare.res.df <- data.frame()
  
  for (population in QueryClasses) {
    print(paste0('Running discriminant analysis for ', population))
    
    # Step 1: Make the background balanced
    ident.df <- data.frame(ident = annotation)
    table_ <- table(ident.df$ident)
    max.n <- max(table_)
    table_[[population]] <- Inf
    min.n <- min(table_)
    
    # Get new ident
    ident.df$new.ident <- "ND"
    ident.df[ident.df$ident == population,]$new.ident <- "ident.1"
    
    for (sub.population in BackgroundClasses[BackgroundClasses != population]) {
      idx.sub.population <- which(ident.df$ident == sub.population)
      idx.random <- sample.int(n = length(idx.sub.population), size = min.n, replace = FALSE)
      idx.sub.population.random <- idx.sub.population[idx.random]
      ident.df[idx.sub.population.random,]$new.ident <- "ident.2"
    }
    
    # Step 2: Discriminant analysis
    marker.res <- RunOPLS(
      mat = mat,
      annotation = ident.df$new.ident,
      ident.1 = "ident.1", 
      ident.2 = "ident.2"
    )
    
    # Add information about cluster ids and idents
    marker.res$feature <- rownames(marker.res)
    marker.res$ident.1 <- population
    marker.res$ident.2 <- NA
    compare.res.df <- rbind(compare.res.df, marker.res)
  }
  
  return(compare.res.df)
}

#' @title Run OPLS
#' @description Wrapper function to run the ropls::opls(). 
#' @param mat Matrix with features and samples in rows and columns, respectively.
#' @param annotation Character vector with the class of the samples.
#' @param ident.1 Character, query class.
#' @param ident.2 Character vector with background classes. If NULL, all the classes in annotation but the query class will be included, Default: NULL.
#' @param scale_weights Boolean; Scale weightMN and weightStarMN.
#' @param scale_vipVn Boolean; Scale vipVn.
#' @return Data frame with a summary of the `ropls::opls()`. The columns are:
#' @return - weightMN: OPLS weights. (see ? ropls::opls)
#' @return - weightStarMN: OPLS projections. (see ? ropls::opls)
#' @seealso 
#'  \code{\link[ropls]{opls}}
#' @rdname RunOPLS
#' @export 
#' @importFrom ropls opls
RunOPLS <- function(mat, annotation, ident.1, ident.2 = NULL, scale_weights = TRUE, scale_vipVn = TRUE, ...) {
  mat <- t(mat)
  
  # Get the classes
  object.classes <- annotation
  
  # Process the information about the ident 1 and 2 (comparison classes)
  if (!is.null(ident.1)) {
    if (is.null(ident.2)) {
      ident.2 <- "Background"
      object.classes[object.classes != ident.1] <- ident.2
    } else {
      object.classes[!(object.classes %in% c(ident.1, ident.2))] <- NA
      mask.2.rm <- is.na(object.classes)
      mat <- mat[!mask.2.rm, ]
      object.classes <- object.classes[!mask.2.rm]
    }
  } else {
    stop("ident.1 is missing, with no default value.")
  }
  
  # The order of the classes matters
  object.classes <- factor(object.classes, levels = c("ident.2", "ident.1"))
  
  # Get the OPLS results
  opls.res <- ropls::opls(
    x = mat, 
    y = object.classes, 
    scaleC = 'standard',  # c("none", "center", "pareto", "standard")
    log10L = FALSE,
    predI = 1,
    orthoI = 1, 
    fig.pdfC = "none",  # no figure is generated
    info.txtC = "none"
  )
  
  # Get the relevant results in a data.frame format
  opls.res.df <- cbind(
    opls.res@weightMN, 
    opls.res@weightStarMN, 
    opls.res@orthoWeightMN,
    opls.res@vipVn, 
    opls.res@orthoVipVn, 
    opls.res@loadingMN, 
    opls.res@orthoLoadingMN, 
    opls.res@coefficientMN
  )
  colnames(opls.res.df) <- c(
    "weightMN", 
    "weightStarMN", 
    "orthoWeightMN",
    "vipVn", 
    "orthoVipVn", 
    "loadingMN", 
    "orthoLoadingMN", 
    "coefficientMN"
  )
  
  # Add variables that were removed due to having variance < 2.2e-16
  ZeroVarVi <- data.frame(row.names = names(opls.res@xZeroVarVi))
  if (nrow(ZeroVarVi) > 0) {
    for (col_ in colnames(opls.res.df)) {
      ZeroVarVi[[col_]] <- NA
    }
    opls.res.df <- rbind(opls.res.df, ZeroVarVi[, colnames(opls.res.df)])
  }
  opls.res.df <- as.data.frame(opls.res.df)
  
  if (scale_weights) {
    sign_w <- sign(opls.res.df$weightMN)
    abs_norm_w <- scale_range(abs(opls.res.df$weightMN))
    norm_w <- sign_w * abs_norm_w
    opls.res.df$weightMN <- norm_w
    
    sign_w <- sign(opls.res.df$weightStarMN)
    abs_norm_w <- scale_range(abs(opls.res.df$weightStarMN))
    norm_w <- sign_w * abs_norm_w
    opls.res.df$weightStarMN <- norm_w
  }
  
  if (scale_vipVn) {
    sign_w <- sign(opls.res.df$vipVn)
    abs_norm_w <- scale_range(abs(opls.res.df$vipVn))
    norm_w <- sign_w * abs_norm_w
    opls.res.df$vipVn <- norm_w
  }
  
  return(opls.res.df)
}



#' @title scale_range
#' @description This function scales a numerical vector from 0 to 1. 
#' @param x Numerical vector.
#' @return Numerical vector. The input numerical vector scaled to range from 0 to 1.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  scale_range(x = c(0.1, 0.2, 0.1, 0.3))
#'  }
#' }
#' @rdname scale_range
#' @export 
scale_range <- function(x) {
  range_ <- range(x, na.rm = TRUE)
  if ((range_[2] - range_[1]) == 0) {
    y <- replicate(n = length(x), expr = 1)
  } else {
    y <- (x - range_[1]) / (range_[2] - range_[1])
  }
  return(y)
}


#' @title Calculate the Global Discriminant Score (GDS)
#' @description This function calculates the average maximum discriminancy (for all classes) at each clustering resolution (k).
#' @param clusters.DA Data frame with the output from the DiscriminantAnalysis() function.
#' @return Data frame with the GDS for each clustering resolution (k).
#' @rdname CalculateGDS
#' @export 
CalculateGDS <- function(clusters.DA = NULL){
  print('Calculating the GDS...')
  if(is.null(clusters.DA)){print('Wrong input');exit}
  # The total number of cluster:
  clusters.DA$NumCluster <- as.numeric(gsub("[.][0-9]+", "", clusters.DA$feature))
  # The cluster id within the total number of cluster
  clusters.DA$ClusterID <- as.numeric(gsub("[0-9]+[.]", "", clusters.DA$feature))
  # Get the test design description
  clusters.DA$TestDesign <- paste0(clusters.DA$ident.1, " VS background")
  clusters.DA$TestDesign <- as.factor(clusters.DA$TestDesign)
  
  
  # Collect the best discriminant score for each cell state and clustering resolution (k)
  clusters.TILs.DA.max.df <- clusters.DA %>% 
    group_by(TestDesign, NumCluster) %>% 
    summarise(
      .groups = 'drop',
      max_weight = max(weightMN),
      max_ClusterID = ClusterID[which(weightMN == max(weightMN))]
    )
  clusters.TILs.DA.max.df <- data.frame(clusters.TILs.DA.max.df)
  clusters.TILs.DA.max.df$group <- paste(clusters.TILs.DA.max.df$NumCluster, clusters.TILs.DA.max.df$TestDesign)
  
  clusters.TILs.DA.max.df <- clusters.DA %>% 
    group_by(TestDesign, NumCluster) %>% 
    summarise(
      .groups ='drop',
      max_weight = max(weightMN),
      max_ClusterID = ClusterID[which(weightMN == max(weightMN))]
    )
  clusters.TILs.DA.max.df <- data.frame(clusters.TILs.DA.max.df)
  
  clusters.TILs.DA.max.df$group <- paste(clusters.TILs.DA.max.df$NumCluster, clusters.TILs.DA.max.df$TestDesign)
  
  # Calculate the Global Discriminant Score (GDS) 
  clusters.TILs.DA.max.avg.df <- clusters.TILs.DA.max.df %>% 
    group_by(NumCluster) %>% 
    summarise(
      .groups = 'drop',
      avg_max_weight = mean(max_weight)
    )
  clusters.TILs.DA.max.avg.df <- data.frame(clusters.TILs.DA.max.avg.df)
  return(clusters.TILs.DA.max.avg.df)
}


#' @title Range [0,1]
#' @description Function to normalize values of a numeric vector as follows: (x-min(x))/(max(x)-min(x)).
#' @param x A numeric vector.
#' @return A list object with the normalized values, and the maximum and minimum values of the input vector.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' range01(c('-50','-5','5','50','100')
#'  }
#' }
#' @rdname range01
#' @export 
range01 = function(x){(x-min(x, na.rm=T))/(max(x,na.rm=T)-min(x, na.rm=T))}

#' @title String split
#' @description Split a vector of composite names into a matrix of simple names.
#' @param x character vector.
#' @param split character to split each element of vector on, see strsplit.
#' @param ... other arguments are passed to strsplit
#' @return Matrix with the elements of the composite names split in columns
#' @details This function is the same as strsplit except that the output value is a matrix instead of a list. 
#' The first column of the matrix contains the first component from each element of x, the second column contains the second components etc.
#' The number of columns is equal to the maximum number of components for any element of x.
#'The motivation for this function in the limma package is handle input columns which are composites of two or more annotation fields.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' x = c("AA_BB_1","AA_BB_2")
#' strsplit2(x,split="_")
#'  }
#' }
#' @rdname strsplit2
#' @export  
strsplit2 = function (x, split, ...) {
  x <- as.character(x)
  n <- length(x)
  s <- strsplit(x, split = split, ...)
  nc <- unlist(lapply(s, length))
  out <- matrix("", n, max(nc))
  for (i in 1:n) {
    if (nc[i]) 
      out[i, 1:nc[i]] <- s[[i]]
  }
  out
}

#' @title Discriminant Analysis Plus
#' @description This function computes the discriminant analysis (OPLS-DA from the `ropls` library), comparing every class in a categorical variable against the background (i.e., the other classes). It provides more information than DiscriminantAnalysis(), as it includes the complete output from the `ropls::opls()` function.
#' @param data Matrix with features and samples in rows and columns, respectively.
#' @param annotation Character vector with the class of the samples.
#' @param BackgroundClasses Classes to compare with the query class. If NULL, all the classes in annotation except the query class will be included as background. Default: NULL
#' @param QueryClasses Character vector with the classes to compare with respect to the background. If NULL, all the classes in annotation will be analyzed (separately) with respect to the background. Default: NULL
#' @param features Character vector with the names of the features to be included in the discriminant analysis. If NULL, all the features in data rows are included. Default: NULL
#' @param reduction.name Name to label the analysis. Default: NULL
#' @param pairwise Logical. If TRUE, a pairwise comparison between classes is performed. If NULL, each query class is compared with respect to the background classes. Default: NULL
#' @return List with complete results from the `ropls::opls()` function for all comparisons.
#' @rdname DiscriminantAnalysis.plus
#' @export 
DiscriminantAnalysis.plus = function(
    data = NULL, 
    annotation = NULL,
    BackgroundClasses = NULL,
    QueryClasses = NULL,
    features = NULL,
    reduction.name = NULL,
    pairwise = NULL
){
  if (is.null(pairwise)){
    pairwise = FALSE
  } 
  
  if (is.null(features)){
    features = rownames(data)
  }
  if (is.null(reduction.name)){
    reduction.name = 'OPLSDA'
  }
  if (is.null(BackgroundClasses)){
    BackgroundClasses = unique(annotation)
  }
  if (is.null(QueryClasses)){
    QueryClasses = BackgroundClasses
  }
  
  TRUTH.vector = annotation %in% c(BackgroundClasses,QueryClasses)
  
  data = data[features,TRUTH.vector]
  annotation = annotation[TRUTH.vector]
  
  #Store original Seurat object idents
  idents.backup = annotation
  
  if(!isTRUE(pairwise)){
    results.list = list()
    for(population in QueryClasses){
      idents.tmp = ifelse(idents.backup == population,"ident.1","ident.2")
      
      opls.res = RunOPLS.plus(data = data,
                                               features = features,
                                               annotation = idents.tmp,
                                               ident.1 = "ident.1", 
                                               ident.2 = "ident.2")
      
      # Complete loadings
      missing = setdiff(features,rownames(opls.res@loadingMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@loadingMN = as.matrix(rbind(opls.res@loadingMN, df.tmp))
      }
      # Complete orthoLoadings
      missing = setdiff(features,rownames(opls.res@orthoLoadingMN))
      if(length(missing)>0){
        df.tmp = data.frame(o1=missing)
        rownames(df.tmp) = df.tmp$o1
        df.tmp$o1 = 0
        opls.res@orthoLoadingMN = as.matrix(rbind(opls.res@orthoLoadingMN, df.tmp))
      }
      # Complete weightStarMN
      missing = setdiff(features,rownames(opls.res@weightStarMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@weightStarMN = as.matrix(rbind(opls.res@weightStarMN, df.tmp))
      }
      # Complete weightMN
      missing = setdiff(features,rownames(opls.res@weightMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@weightMN = as.matrix(rbind(opls.res@weightMN, df.tmp))
      }
      # Complete orthoWeightMN
      missing = setdiff(features,rownames(opls.res@orthoWeightMN))
      if(length(missing)>0){
        df.tmp = data.frame(o1=missing)
        rownames(df.tmp) = df.tmp$o1
        df.tmp$o1 = 0
        opls.res@orthoWeightMN = as.matrix(rbind(opls.res@orthoWeightMN, df.tmp))
      }
      # Complete vipVn
      missing = setdiff(features,names(opls.res@vipVn))
      if(length(missing)>0){
        df.tmp = data.frame(vipVn=missing)
        rownames(df.tmp) = df.tmp$vipVn
        df.tmp$vipVn = 0
        vipVn.tmp = as.table(df.tmp$vipVn)
        names(vipVn.tmp) = rownames(df.tmp)
        opls.res@vipVn = c(opls.res@vipVn, vipVn.tmp)
      }
      # Complete embeddings
      results.list[[population]] = opls.res
    } 
  }
  
  if(isTRUE(pairwise)){
    tests.for = combn(length(QueryClasses), 2, FUN = NULL, simplify = TRUE)
    test.rev = rbind(tests.for[2,],tests.for[1,])
    tests = cbind(tests.for,test.rev)
    
    results.list = list()
    for(i in 1:ncol(tests)){
      population.1 = as.character(QueryClasses[tests[1,i]])
      population.2 = as.character(QueryClasses[tests[2,i]])
      # Step 1: Make the background balanced
      # Get ident.df
      ident.df = data.frame(ident=idents.backup)
      
      ident.df$new.ident = "ident.3"
      ident.df[ident.df$ident == population.1,]$new.ident = "ident.1"
      ident.df[ident.df$ident == population.2,]$new.ident = "ident.2"
      
      opls.res = RunOPLS.plus(data = data,
                                      features = features,
                                      annotation = ident.df$new.ident,
                                      ident.1 = "ident.1", 
                                      ident.2 = "ident.2")
      
      # Complete loadings
      missing = setdiff(features,rownames(opls.res@loadingMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@loadingMN = as.matrix(rbind(opls.res@loadingMN, df.tmp))
      }
      # Complete orthoLoadings
      missing = setdiff(features,rownames(opls.res@orthoLoadingMN))
      if(length(missing)>0){
        df.tmp = data.frame(o1=missing)
        rownames(df.tmp) = df.tmp$o1
        df.tmp$o1 = 0
        opls.res@orthoLoadingMN = as.matrix(rbind(opls.res@orthoLoadingMN, df.tmp))
      }
      # Complete weightStarMN
      missing = setdiff(features,rownames(opls.res@weightStarMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@weightStarMN = as.matrix(rbind(opls.res@weightStarMN, df.tmp))
      }
      # Complete weightMN
      missing = setdiff(features,rownames(opls.res@weightMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@weightMN = as.matrix(rbind(opls.res@weightMN, df.tmp))
      }
      # Complete orthoWeightMN
      missing = setdiff(features,rownames(opls.res@orthoWeightMN))
      if(length(missing)>0){
        df.tmp = data.frame(o1=missing)
        rownames(df.tmp) = df.tmp$o1
        df.tmp$o1 = 0
        opls.res@orthoWeightMN = as.matrix(rbind(opls.res@orthoWeightMN, df.tmp))
      }
      # Complete vipVn
      missing = setdiff(features,names(opls.res@vipVn))
      if(length(missing)>0){
        df.tmp = data.frame(vipVn=missing)
        rownames(df.tmp) = df.tmp$vipVn
        df.tmp$vipVn = 0
        vipVn.tmp = as.table(df.tmp$vipVn)
        names(vipVn.tmp) = rownames(df.tmp)
        opls.res@vipVn = c(opls.res@vipVn,vipVn.tmp)
      }
      # Complete embeddings
      data.matrix = as.matrix(data)
      data.matrix = data.matrix[rownames(opls.res@loadingMN),] 
      data.matrix = t(data.matrix)
      
      mu = colMeans(data.matrix)[rownames(opls.res@loadingMN)]
      sdev = apply(data.matrix,2,sd,na.rm=T)[rownames(opls.res@loadingMN)]
      
      # Center and scale the query data
      data.matrix.centered = t(t(data.matrix[,rownames(opls.res@loadingMN),drop=T]) - mu)
      data.matrix.centered.scaled = t(t(data.matrix.centered) / sdev)
      
      projected.data = data.matrix.centered.scaled %*% opls.res@weightStarMN
      
      E = data.matrix.centered.scaled %*% opls.res@orthoWeightMN
      
      projected.data = projected.data - E
      
      opls.res@scoreMN = projected.data
      
      results.list[[paste0(population.1,'_Vs_',population.2)]] = opls.res
    } 
  }
  
  opls.loadings.list = lapply(results.list,function(x){
    loadings.opls = x@loadingMN
  })
  
  for(name in names(opls.loadings.list)){
    colnames(opls.loadings.list[[name]])=name
  }
  
  opls.embeddings.list = lapply(results.list,function(x){
    embeddings.opls= x@scoreMN
  })
  
  for(name in names(opls.embeddings.list)){
    colnames(opls.embeddings.list[[name]])=name
  }
  
  # Add the missing variables (with 0 variance)
  
  for(name in names(opls.loadings.list)){
    missing = setdiff(features,rownames(opls.loadings.list[[name]]))
    if(length(missing)>0){
      df.tmp = data.frame(p1=missing)
      rownames(df.tmp) = df.tmp$p1
      df.tmp$p1 = 0
      colnames(df.tmp) = name
      opls.loadings.list[[name]] = rbind(opls.loadings.list[[name]], df.tmp)
    }
  }
  
  df.loadings = Reduce(cbind, opls.loadings.list)
  opls.axis.names.loadings = data.frame(Predictive.axis = colnames(df.loadings))
  colnames(df.loadings) = paste0('OPLSDA_',c(1:ncol(df.loadings)))
  opls.axis.names.loadings$Predictive.axis.name = colnames(df.loadings)
  rownames(opls.axis.names.loadings) = opls.axis.names.loadings$Predictive.axis
  
  df.embeddings = Reduce(cbind, opls.embeddings.list)
  opls.axis.names.embeddings = data.frame(Predictive.axis = colnames(df.embeddings))
  colnames(df.embeddings) = paste0('OPLSDA_',c(1:ncol(df.embeddings)))
  opls.axis.names.embeddings$Predictive.axis.name = colnames(df.embeddings)
  rownames(opls.axis.names.embeddings) = opls.axis.names.embeddings$Predictive.axis
  
  results.list[['opls.axis.names.loadings']] = opls.axis.names.loadings
  results.list[['opls.axis.names.embeddings']] = opls.axis.names.embeddings
  
  #saveRDS(file="OPLSDA.results.list.Rds",results.list)
  
  return(results.list)
}

#' @title Discriminant Analysis Plus Seurat
#' @description This function computes the discriminant analysis (OPLS-DA from the `ropls` library), comparing every class in a categorical variable against the background (i.e., the other classes). It provides more information than DiscriminantAnalysis(), as it includes the complete output from the `ropls::opls()` function.
#' @param seurat.object Seurat object
#' @param assay Name of the assay in `seurat.object` to pull the slot data from. If NULL, it takes the name of the active assay in `seurat.object`, default = NULL.
#' @param slot Slot to pull data from, default: 'data'.
#' @param annotation Character vector with the class of the samples.
#' @param BackgroundClasses Classes to compare with the query class. If NULL, all the classes in annotation except the query class will be included as background. Default: NULL
#' @param QueryClasses Character vector with the classes to compare with respect to the background. If NULL, all the classes in annotation will be analyzed (separately) with respect to the background. Default: NULL
#' @param features Character vector with the names of the features to be included in the discriminant analysis. If NULL, all the features in data rows are included. Default: NULL
#' @param reduction.name Name to label the analysis. Default: NULL
#' @param pairwise Logical. If TRUE, a pairwise comparison between classes is performed. If NULL, each query class is compared with respect to the background classes. Default: NULL
#' @return The Seurat object including a list with complete results from the `ropls::opls()`.
#' @seealso 
#'  \code{\link[Seurat]{reexports}}
#' @rdname DiscriminantAnalysis.plus.seurat
#' @export 
#' @importFrom Seurat Idents CreateDimReducObject
DiscriminantAnalysis.plus.seurat = function(seurat.object = NULL, 
                            assay = NULL, 
                            slot = "data",
                            annotation = NULL,
                            features = NULL,
                            BackgroundClasses = NULL,
                            QueryClasses = NULL,
                            reduction.name = NULL,
                            pairwise = NULL) {
  if (is.null(pairwise)){
    pairwise = FALSE
  } 
  if (is.null(assay)){
    assay = seurat.object@active.assay
  } 
  if (is.null(features)){
    DefaultAssay(seurat.object) = assay  
    features = rownames(seurat.object)
  }
  if (is.null(reduction.name)){
    reduction.name = 'OPLSDA'
  }
  
  if (is.null(BackgroundClasses)){
    BackgroundClasses = unique(annotation)
  }
  if (is.null(QueryClasses)){
    QueryClasses = BackgroundClasses
  }
  
  #TRUTH.vector = annotation %in% c(BackgroundClasses,QueryClasses)
  
  #data = data[features,TRUTH.vector]
  #annotation = annotation[TRUTH.vector]
  
  
  
  #Store original Seurat object idents
  idents.backup = Seurat::Idents(object = seurat.object)
  
  # Get the list of classes to compare with the background
  seurat.object@meta.data[, annotation] = as.factor(seurat.object@meta.data[, annotation])
  if (is.null(classes)){
    classes = levels(seurat.object@meta.data[,annotation])
  }
  #Seurat::Idents(seurat.object) = seurat.object[[annotation]]
  Seurat::Idents(seurat.object) = seurat.object@meta.data[,annotation]
  
  
  
  if(!isTRUE(pairwise)){
    results.list = list()
    for(population in QueryClasses){
      # Step 1: Make the background balanced
      # Get ident.df
      ident.df = seurat.object@meta.data[, annotation, drop = F]
      colnames(ident.df) = c("ident")
      
      ident.df$new.ident = "ident.3"
      ident.df[ident.df$ident %in% BackgroundClasses,]$new.ident = "ident.2"
      ident.df[ident.df$ident == population,]$new.ident = "ident.1"
      
      # Set the Idents
      seurat.object@meta.data$new.ident = ident.df[rownames(seurat.object@meta.data),]$new.ident
      Seurat::Idents(seurat.object) = seurat.object$new.ident
      
      opls.res = RunOPLS.plus.seurat(seurat.object = seurat.object,
                                     assay = assay,
                                     features = features,
                                     slot = slot,
                                     ident.1 = "ident.1", 
                                     ident.2 = "ident.2")
      
      # 
      DefaultAssay(seurat.object)=assay
      # Complete loadings
      missing = setdiff(features,rownames(opls.res@loadingMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@loadingMN = as.matrix(rbind(opls.res@loadingMN, df.tmp))
      }
      # Complete orthoLoadings
      missing = setdiff(features,rownames(opls.res@orthoLoadingMN))
      if(length(missing)>0){
        df.tmp = data.frame(o1=missing)
        rownames(df.tmp) = df.tmp$o1
        df.tmp$o1 = 0
        opls.res@orthoLoadingMN = as.matrix(rbind(opls.res@orthoLoadingMN, df.tmp))
      }
      # Complete weightStarMN
      missing = setdiff(features,rownames(opls.res@weightStarMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@weightStarMN = as.matrix(rbind(opls.res@weightStarMN, df.tmp))
      }
      # Complete weightMN
      missing = setdiff(features,rownames(opls.res@weightMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@weightMN = as.matrix(rbind(opls.res@weightMN, df.tmp))
      }
      # Complete orthoWeightMN
      missing = setdiff(features,rownames(opls.res@orthoWeightMN))
      if(length(missing)>0){
        df.tmp = data.frame(o1=missing)
        rownames(df.tmp) = df.tmp$o1
        df.tmp$o1 = 0
        opls.res@orthoWeightMN = as.matrix(rbind(opls.res@orthoWeightMN, df.tmp))
      }
      # Complete vipVn
      missing = setdiff(features,names(opls.res@vipVn))
      if(length(missing)>0){
        df.tmp = data.frame(vipVn=missing)
        rownames(df.tmp) = df.tmp$vipVn
        df.tmp$vipVn = 0
        vipVn.tmp = as.table(df.tmp$vipVn)
        names(vipVn.tmp) = rownames(df.tmp)
        opls.res@vipVn = c(opls.res@vipVn, vipVn.tmp)
      }
      # Complete embeddings
      results.list[[population]] = opls.res
    } 
  }
  
  if(isTRUE(pairwise)){
    tests.for = combn(length(QueryClasses), 2, FUN = NULL, simplify = TRUE)
    test.rev = rbind(tests.for[2,],tests.for[1,])
    tests = cbind(tests.for,test.rev)
    
    results.list = list()
    for(i in 1:ncol(tests)){
      population.1 = as.character(QueryClasses[tests[1,i]])
      population.2 = as.character(QueryClasses[tests[2,i]])
      # Step 1: Make the background balanced
      # Get ident.df
      ident.df = seurat.object@meta.data[, annotation, drop = F]
      colnames(ident.df) = c("ident")
      
      ident.df$new.ident = "ident.3"
      ident.df[ident.df$ident == population.1,]$new.ident = "ident.1"
      ident.df[ident.df$ident == population.2,]$new.ident = "ident.2"
      
      # Set the Idents
      seurat.object@meta.data$new.ident = ident.df[rownames(seurat.object@meta.data),]$new.ident
      Seurat::Idents(seurat.object) = seurat.object$new.ident
      
      opls.res = RunOPLS.plus.seurat(seurat.object = seurat.object,
                                     assay = assay,
                                     features = features,
                                     slot = slot,
                                     ident.1 = "ident.1", 
                                     ident.2 = "ident.2")
      
      # 
      DefaultAssay(seurat.object)=assay
      
      # Complete loadings
      missing = setdiff(features,rownames(opls.res@loadingMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@loadingMN = as.matrix(rbind(opls.res@loadingMN, df.tmp))
      }
      # Complete orthoLoadings
      missing = setdiff(features,rownames(opls.res@orthoLoadingMN))
      if(length(missing)>0){
        df.tmp = data.frame(o1=missing)
        rownames(df.tmp) = df.tmp$o1
        df.tmp$o1 = 0
        opls.res@orthoLoadingMN = as.matrix(rbind(opls.res@orthoLoadingMN, df.tmp))
      }
      # Complete weightStarMN
      missing = setdiff(features,rownames(opls.res@weightStarMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@weightStarMN = as.matrix(rbind(opls.res@weightStarMN, df.tmp))
      }
      # Complete weightMN
      missing = setdiff(features,rownames(opls.res@weightMN))
      if(length(missing)>0){
        df.tmp = data.frame(p1=missing)
        rownames(df.tmp) = df.tmp$p1
        df.tmp$p1 = 0
        opls.res@weightMN = as.matrix(rbind(opls.res@weightMN, df.tmp))
      }
      # Complete orthoWeightMN
      missing = setdiff(features,rownames(opls.res@orthoWeightMN))
      if(length(missing)>0){
        df.tmp = data.frame(o1=missing)
        rownames(df.tmp) = df.tmp$o1
        df.tmp$o1 = 0
        opls.res@orthoWeightMN = as.matrix(rbind(opls.res@orthoWeightMN, df.tmp))
      }
      # Complete vipVn
      missing = setdiff(features,names(opls.res@vipVn))
      if(length(missing)>0){
        df.tmp = data.frame(vipVn=missing)
        rownames(df.tmp) = df.tmp$vipVn
        df.tmp$vipVn = 0
        vipVn.tmp = as.table(df.tmp$vipVn)
        names(vipVn.tmp) = rownames(df.tmp)
        opls.res@vipVn = c(opls.res@vipVn,vipVn.tmp)
      }
      # Complete embeddings
      seurat.object.matrix = as.matrix(GetAssayData(object = seurat.object, assay = assay, slot = slot))
      seurat.object.matrix = seurat.object.matrix[rownames(opls.res@loadingMN),] 
      seurat.object.matrix = t(seurat.object.matrix)
      
      mu = colMeans(seurat.object.matrix)[rownames(opls.res@loadingMN)]
      sdev = apply(seurat.object.matrix,2,sd,na.rm=T)[rownames(opls.res@loadingMN)]
      
      # Center and scale the query data
      seurat.object.matrix.centered = t(t(seurat.object.matrix[,rownames(opls.res@loadingMN),drop=T]) - mu)
      seurat.object.matrix.centered.scaled = t(t(seurat.object.matrix.centered) / sdev)
      
      projected.data = seurat.object.matrix.centered.scaled %*% opls.res@weightStarMN
      
      E = seurat.object.matrix.centered.scaled %*% opls.res@orthoWeightMN
      
      projected.data = projected.data - E
      
      opls.res@scoreMN = projected.data
      
      results.list[[paste0(population.1,'_Vs_',population.2)]] = opls.res
    } 
  }
  
  opls.loadings.list = lapply(results.list,function(x){
    loadings.opls = x@loadingMN
  })
  
  for(name in names(opls.loadings.list)){
    colnames(opls.loadings.list[[name]])=name
  }
  
  opls.embeddings.list = lapply(results.list,function(x){
    embeddings.opls= x@scoreMN
  })
  
  for(name in names(opls.embeddings.list)){
    colnames(opls.embeddings.list[[name]])=name
  }
  
  # Add the missing variables (with 0 variance)
  DefaultAssay(seurat.object) = assay
  for(name in names(opls.loadings.list)){
    missing = setdiff(features,rownames(opls.loadings.list[[name]]))
    if(length(missing)>0){
      df.tmp = data.frame(p1=missing)
      rownames(df.tmp) = df.tmp$p1
      df.tmp$p1 = 0
      colnames(df.tmp) = name
      opls.loadings.list[[name]] = rbind(opls.loadings.list[[name]], df.tmp)
    }
  }
  
  df.loadings = Reduce(cbind, opls.loadings.list)
  opls.axis.names.loadings = data.frame(Predictive.axis = colnames(df.loadings))
  colnames(df.loadings) = paste0('OPLSDA_',c(1:ncol(df.loadings)))
  opls.axis.names.loadings$Predictive.axis.name = colnames(df.loadings)
  rownames(opls.axis.names.loadings) = opls.axis.names.loadings$Predictive.axis
  
  df.embeddings = Reduce(cbind, opls.embeddings.list)
  opls.axis.names.embeddings = data.frame(Predictive.axis = colnames(df.embeddings))
  colnames(df.embeddings) = paste0('OPLSDA_',c(1:ncol(df.embeddings)))
  opls.axis.names.embeddings$Predictive.axis.name = colnames(df.embeddings)
  rownames(opls.axis.names.embeddings) = opls.axis.names.embeddings$Predictive.axis
  
  results.list[['opls.axis.names.loadings']] = opls.axis.names.loadings
  results.list[['opls.axis.names.embeddings']] = opls.axis.names.embeddings
  
  #saveRDS(file="OPLSDA.results.list.Rds",results.list)
  
  DimReducObject = Seurat::CreateDimReducObject(embeddings = as.matrix(df.embeddings),
                                                loadings = as.matrix(df.loadings),
                                                stdev = 0,
                                                key = 'OPLSDA_',
                                                assay = assay)
  
  seurat.object[[reduction.name]] = DimReducObject
  seurat.object@tools[[reduction.name]] = results.list
  
  seurat.object@reductions[[reduction.name]]@cell.embeddings = DimReducObject@cell.embeddings
  seurat.object@reductions[[reduction.name]]@key = DimReducObject@key
  
  Seurat::Idents(object = seurat.object) = idents.backup 
  
  return(seurat.object)
}


#' @title Run OPLS Plus Seurat
#' @description Wrapper function to run the ropls::opls(). 
#' @param seurat.object Seurat object.
#' @param ident.1 Character, query class.
#' @param ident.2 Character vector with background classes. If NULL, all the classes in Seurat::Idents(object=seurat.object) but the query class will be included, Default: NULL.
#' @param features Features to be included in the analysis. If NULL, all the features in in the assay and slot are included , Default: NULL.
#' @param assay Name of the assay in `seurat.object` to pull the slot data from. If NULL, it takes the name of the active assay in `seurat.object`, default = NULL.
#' @param slot Slot to pull data from, default: 'data'.
#' @return Data frame with a summary of the `ropls::opls()`. The columns are:
#' @return - weightMN: OPLS weights. (see ? ropls::opls).
#' @return - weightStarMN: OPLS projections. (see ? ropls::opls).
#' @seealso 
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[ropls]{opls}}
#' @rdname RunOPLS.plus.seurat
#' @export 
#' @importFrom Seurat Idents
#' @importFrom ropls opls
RunOPLS.plus.seurat = function(seurat.object, 
                               ident.1, 
                               ident.2 = NULL,
                               features = NULL,
                               assay = NULL, 
                               slot = "data", ...) {
  # sanity checks + setting parameters
  if (is.null(assay)){
    assay = seurat.object@active.assay
  }
  
  # Get the seurat.object matrix
  seurat.object.matrix = as.matrix(GetAssayData(object = seurat.object, assay=assay, slot=slot))
  seurat.object.matrix = seurat.object.matrix[features,] 
  seurat.object.matrix = t(seurat.object.matrix)
  
  # Get the classes
  seurat.object.classes = Seurat::Idents(object=seurat.object)
  seurat.object.classes = as.character(seurat.object.classes)
  
  # Process the information about the idents 1 and 2 (comparison classes)
  if (!is.null(ident.1)){
    if (is.null(ident.2)){
      ident.2 = "Background"
      seurat.object.classes[seurat.object.classes != ident.1] = ident.2
    } else {
      seurat.object.classes[!(seurat.object.classes %in% c(ident.1, ident.2))] = NA
      
      mask.2.rm = is.na(seurat.object.classes)
      
      seurat.object.matrix = seurat.object.matrix[!mask.2.rm, ]
      seurat.object.classes = seurat.object.classes[!mask.2.rm]
    }
  } else {
    stop("ident.1 is missing, with no default value.")
  }
  
  # ATTENTION. The order of the classes matter. 
  # The background class is the first level and the foreground class is the second level
  seurat.object.classes = factor(seurat.object.classes, levels=c("ident.2", "ident.1"))
  
  # Get the OPLS results
  opls.res = ropls::opls(x = seurat.object.matrix, 
                         y = seurat.object.classes, 
                         scaleC = 'standard', # c("none", "center", "pareto", "standard")
                         log10L = FALSE,
                         predI = 1,
                         orthoI = 1, 
                         fig.pdfC = "none", # no figure is generated
                         info.txtC = "none", # no save file name is provided
                         ...)
  # Return the results
  return(opls.res)
}

#' @title Run OPLS Plus
#' @description Wrapper function to run the ropls::opls(). 
#' @param data Matrix with features and samples in rows and columns, respectively.
#' @param annotation Character vector with the class of the samples.
#' @param ident.1 Character, query class.
#' @param ident.2 Character vector with background classes. If NULL, all the classes in annotation but the query class will be included, Default: NULL.
#' @param features Features to be included in the analysis. If NULL, all the features in in the assay and slot are included , Default: NULL.
#' @return List with complete results from the `ropls::opls()` function for all comparisons.
#' @seealso 
#'  \code{\link[ropls]{opls}}
#' @rdname RunOPLS.plus
#' @export 
#' @importFrom ropls opls
RunOPLS.plus = function(data,annotation,ident.1, ident.2 = NULL, features = NULL, ...) {
  
  if(is.null(features)){
    features = rownames(data)
  }
  
  seurat.object.matrix = data
  seurat.object.matrix = seurat.object.matrix[features,] 
  seurat.object.matrix = t(seurat.object.matrix)
  
  # Get the classes
  seurat.object.classes = annotation
  seurat.object.classes = as.character(seurat.object.classes)
  
  # Process the information about the idents 1 and 2 (comparison classes)
  if (!is.null(ident.1)){
    if (is.null(ident.2)){
      ident.2 = "Background"
      seurat.object.classes[seurat.object.classes != ident.1] = ident.2
    } else {
      seurat.object.classes[!(seurat.object.classes %in% c(ident.1, ident.2))] = NA
      
      mask.2.rm = is.na(seurat.object.classes)
      
      seurat.object.matrix = seurat.object.matrix[!mask.2.rm, ]
      seurat.object.classes = seurat.object.classes[!mask.2.rm]
    }
  } else {
    stop("ident.1 is missing, with no default value.")
  }
  
  # ATTENTION. The order of the classes matters. 
  # The background class is the first level and the foreground class is the second level
  seurat.object.classes = factor(seurat.object.classes, levels=c("ident.2", "ident.1"))
  
  # Get the OPLS results
  opls.res = ropls::opls(x = seurat.object.matrix, 
                         y = seurat.object.classes, 
                         scaleC = 'standard', # c("none", "center", "pareto", "standard")
                         log10L = FALSE,
                         predI = 1,
                         orthoI = 1, 
                         fig.pdfC = "none", # no figure is generated
                         info.txtC = "none", # no save file name is provided
                         ...)
  # Return the results
  return(opls.res)
}



#' @title Modulon Identification
#' @description Wrapper function to perform the first step of a modulon analysis: modulon identification.
#' @param data Matrix of regulon activity with features (regulons) in rows and samples (cells) in columns.
#' @param annotation Character vector with the class (phenotype/state) of the samples (cells).
#' @param BackgroundClasses Classes to compare with the query class. If NULL, all the classes in annotation except the query class will be included as background. Default: NULL
#' @param QueryClasses Character vector with the classes to compare with respect to the background. If NULL, all the classes in annotation will be analyzed (separately) with respect to the background. Default: NULL
#' @param k.range Clustering resolution (k) range to be explored. Default: c(2, 30)
#' @return List with the identified modulons and the global discriminant score (GDS) at every k within the range defined by k.range.
#' @rdname ModulonIdent
#' @export 
ModulonIdent = function(data, annotation,BackgroundClasses=NULL,QueryClasses=NULL, k.range = c(2,30)){
  library(dplyr)
  library(AUCell)
  if (is.null(BackgroundClasses)){
    BackgroundClasses = unique(annotation)
  }
  if (is.null(QueryClasses)){
    QueryClasses = BackgroundClasses
  }
  
  data.backup = data
  annotation.backup = annotation
  
  TRUTH.vector.query = annotation %in% c(QueryClasses)
  TRUTH.vector.background = annotation %in% c(BackgroundClasses)
  TRUTH.vector.query.and.background = annotation %in% c(QueryClasses,BackgroundClasses)
  
  data.query = data[, TRUTH.vector.query]
  annotation.query = annotation[ TRUTH.vector.query]
  
  data.background = data[, TRUTH.vector.background]
  annotation.background = annotation[ TRUTH.vector.background]
  
  data.query.and.background = data[, TRUTH.vector.query.and.background]
  annotation.query.and.background = annotation[ TRUTH.vector.query.and.background]
  
  # Generate TF clusters at different clustering resolutions
  clusters <- FindFeatureCluster(
    distance.method <- "euclidean",
    clustering.method <- "complete",
    mat = data.query, 
    annotation = annotation.query, 
    cluster.nums = k.range)
  # Calculate cluster's signatures scores (AUC)
  clusters.AUC <- GetSignatureAUC(
    mat = data,
    signature.list = clusters)
  # OPLS discriminant analysis
  clusters.DA <- DiscriminantAnalysis(
    mat = clusters.AUC,
    method = "OPLS", 
    annotation = annotation,
    BackgroundClasses = BackgroundClasses,
    QueryClasses = QueryClasses
    )
  # Calculate the global discriminant score (GDS)
  GDS = CalculateGDS(clusters.DA=clusters.DA)
  
  # Identify the best clustering resolution
  Best.k = GDS[which.max(GDS$avg_max_weight),'NumCluster']
  
  # Generate TF clusters at the optimal resolution
  modulons <- FindFeatureCluster(
    distance.method <- "euclidean",
    clustering.method <- "complete",
    mat = data.query, 
    annotation = annotation.query, 
    cluster.nums = Best.k)
  
  return(list(Modulons=modulons,GDS=GDS))
}


#' @title Modulon Selection
#' @description Wrapper function to perform the second step of a modulon analysis: modulon selection.
#' @param data Matrix of regulon activity with features (regulons) in rows and samples (cells) in columns.
#' @param modulons List object with the modulon constituent elements.
#' @param annotation Character vector with the class (phenotype/state) of the samples (cells).
#' @param BackgroundClasses Classes to compare with the query class. If NULL, all the classes in annotation except the query class will be included as background. Default: NULL.
#' @param TargetState Character vector with the names of the classes to be analyzed. If NULL, all the classes in annotation will be analyzed, Default: NULL.
#' @return List object with 3 elements, namely, Modulon_AUC, Modulon_DA, and Selected_Modulon:
#' @return * Modulon_AUC: modulon activity.
#' @return * Modulon_DA: results of the modulon discriminant analysis.
#' @return * Selected_Modulon: top discriminant modulon for each class.
#' @rdname ModulonSelect
#' @export 
ModulonSelect = function(data, modulons, annotation,BackgroundClasses=NULL,TargetState=NULL){
  if(is.null(TargetState)){
    TargetState = unique(annotation)
  }
  if (is.null(BackgroundClasses)){
    BackgroundClasses = unique(annotation)
  }
  
  # Calculate modulon signatures across states
  modulons.AUC <- GetSignatureAUC(
    mat = data,
    signature.list = modulons)
  
  # Modulon signature cells state discriminant analysis
  modulons.DA <- DiscriminantAnalysis(
    mat = modulons.AUC,
    method = "OPLS", 
    annotation = annotation,
    BackgroundClasses = BackgroundClasses,
    QueryClasses = TargetState
  )
  
  modulons.DA.split = split(modulons.DA,modulons.DA$ident.1)
  
  Selected_Modulon = lapply(modulons.DA.split,function(x){
    input.tmp = x
    #input.tmp$label = factor(as.character(strsplit2(input.tmp$feature,'[.]')[,2]),levels = as.character(strsplit2(input.tmp$feature,'[.]')[,2]))
    input.tmp$label = input.tmp$feature
    input.tmp = input.tmp[order(input.tmp$weightStarMN,decreasing = T),]
    input.tmp$order = factor(c(1:nrow(input.tmp)))
    best.modulon = as.character(input.tmp$feature[which.max(input.tmp$weightStarMN)])
    return(best.modulon)
  })
  return(list(Modulon_AUC = modulons.AUC, Modulon_DA = modulons.DA.split, Selected_Modulon = Selected_Modulon))
}


#' @title Modulon Perturbation
#' @description Wrapper function to perform the third step of a modulon analysis: modulon perturbation analysis.
#' @param Regulons List object with the regulon constituent elements.
#' @param Modulons List object with the modulon constituent elements.
#' @param ExpMat Matrix of gene expression with features (genes) in rows and samples (cells) in columns.
#' @param annotation Character vector with the class (phenotype/state) of the samples (cells).
#' @param BackgroundClasses Classes to compare with the query class. If NULL, all the classes in annotation except the query class will be included as background. Default: NULL.
#' @param TargetState Character vector with the names of the classes to be analyzed. If NULL, all the classes in annotation will be analyzed. Default: NULL.
#' @param TargetModulon Name of the target modulon.
#' @param CombSize Number of elements in the combination.
#' @param Weights Data frame with the weights of the regulatory influence between transcription factors and target genes. If NULL, a default weight of 1 is used for all regulations. Default: NULL.
#' @return List of elements, namely, Bipartite_Graph, Combinations, and Modulon_Targets_DA:
#' @return * Bipartite_Graph: Data frame with the bipartite graph between transcription factors and regulated genes.
#' @return * Combinations: Data frame with the ranking of combinations.
#' @return * Modulon_Targets_DA: Gene expression discriminant analysis complete results.
#' @rdname ModulonPert
#' @export 
ModulonPert = function(Regulons, Modulons, ExpMat, annotation,BackgroundClasses = NULL, TargetState, TargetModulon, CombSize, Weights=NULL){
  if (is.null(BackgroundClasses)){
    BackgroundClasses = unique(annotation)
  }
  
  # Collect TF-gene regulation weights
  
  if(!is.null(Weights)){
    Weights$weight.norm = range01(Weights$weight)
    rownames(Weights)= paste0(Weights$TF,' -> ',Weights$Target)
  }  
  
  # Evaluate gene importance for a given cell state
  modulon.query = TargetModulon
  query = Modulons[[modulon.query]]
  targets=unique(unlist(Regulons[query]))
  OPLSDA.results=DiscriminantAnalysis.plus(
    data = ExpMat[targets,], 
    annotation = annotation,
    BackgroundClasses = BackgroundClasses,
    QueryClasses = TargetState,
    features = targets,
    reduction.name = 'OPLSDA_RNA_M_Targets',
    pairwise = F
  )
  
  weightStarMN.tmp = OPLSDA.results[[TargetState]]@weightStarMN
  missing = setdiff(targets,rownames(weightStarMN.tmp))
  if(length(missing)>0){
    df.tmp = data.frame(p1=missing)
    rownames(df.tmp) = df.tmp$p1
    df.tmp$p1 = 0
    weightStarMN.tmp = as.matrix(rbind(weightStarMN.tmp, df.tmp))
  }
  WEIGHTS = weightStarMN.tmp
  WEIGHTS.01 = (range01(abs(WEIGHTS)))*sign(WEIGHTS)
  WEIGHTS.ordered = as.data.frame(WEIGHTS.01[order(WEIGHTS.01, decreasing=T),,drop=F])
  #rownames(WEIGHTS.ordered) = rownames(WEIGHTS.01)[order(WEIGHTS.01, decreasing=T)]
  colnames(WEIGHTS.ordered) = 'Weights'
  target.weights = WEIGHTS.ordered$Weights
  names(target.weights)=rownames(WEIGHTS.ordered)
  
  # Build a TF-gene bipartite graph
  
  if(!is.null(Weights)){
    links.w.weights.norm = Weights
  }
  
  results.tmp = list()
  
  for(regulon.tmp in query){
    links.tmp =  paste0(regulon.tmp,' -> ',unlist(regulons[regulon.tmp]))
    
    if(is.null(Weights)){
      links.w.weights.tmp = data.frame(Interactions = links.tmp, Weights = 1)
    }
    if(!is.null(Weights)){
      links.w.weights.tmp = data.frame(Interactions = links.tmp, Weights = links.w.weights.norm[links.tmp,'weight.norm'])
    }
    
    string.split.tmp = strsplit2(links.tmp,' -> ')
    links.w.weights.tmp$Source = string.split.tmp[,1]
    links.w.weights.tmp$Target = string.split.tmp[,2]
    links.w.weights.tmp[is.na(links.w.weights.tmp)] = 0
    links.w.weights.tmp$Target.Discriminant.Score = as.numeric(target.weights[links.w.weights.tmp$Target])
    links.w.weights.tmp$link.score = links.w.weights.tmp$Weights * links.w.weights.tmp$Target.Discriminant.Score
    results.tmp[[regulon.tmp]] = links.w.weights.tmp
  }
  results.tmp.df = do.call(rbind,results.tmp)
  rownames(results.tmp.df)=results.tmp.df$Interactions
  
  # Rank combinations of KOs by predicted impact
  combination.size=CombSize  
  combinations.df = as.data.frame(combn(query,combination.size))
  colnames(combinations.df)=paste0('Combination_',c(1:ncol(combinations.df)))  
  rownames(combinations.df)=paste0('Element_',c(1:nrow(combinations.df)))  
  
  target.score.weighted.unique = apply(combinations.df,2,function(x){
    targets.tmp=unique(unlist(regulons[x]))
    return(sum(target.weights[targets.tmp]))
  }) 
  
  target.score.weighted = apply(combinations.df,2,function(x){
    targets.tmp=(unlist(regulons[x]))
    return(sum(target.weights[targets.tmp]))
  })
  
  target.score.weighted.w.link.weights = apply(combinations.df,2,function(x){
    combination.tmp = x
    results.tmp.df.subset = results.tmp.df[results.tmp.df$Source %in% combination.tmp,]
    return(sum(results.tmp.df.subset$link.score))
  })
  
  Combinations.Table = data.frame(Combination = colnames(combinations.df),
                                  t(combinations.df),
                                  WCS=target.score.weighted.w.link.weights[colnames(combinations.df)]
  )
  return(list(Bipartite_Graph = results.tmp.df,Combinations = Combinations.Table[order(Combinations.Table$WCS,decreasing = T),],Modulon_Targets_DA=target.weights))
}
