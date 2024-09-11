# Define constants for distance and clustering methods
DISTANCE_METHODS <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
CLUSTERING_METHODS <- c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")

# Load required packages quietly
suppressMessages(require(stats))
suppressMessages(require(parallel))
suppressMessages(require(ropls))

# Function to get hierarchical clusters
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param distance.method PARAM_DESCRIPTION, Default: DISTANCE_METHODS
#' @param clustering.method PARAM_DESCRIPTION, Default: CLUSTERING_METHODS
#' @param scale PARAM_DESCRIPTION, Default: TRUE
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
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

# Function to find feature clusters
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param annotation PARAM_DESCRIPTION, Default: NULL
#' @param cluster.nums PARAM_DESCRIPTION, Default: 2:10
#' @param distance.method PARAM_DESCRIPTION, Default: DISTANCE_METHODS
#' @param clustering.method PARAM_DESCRIPTION, Default: CLUSTERING_METHODS
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
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

# Function to calculate signature AUC
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param signature.list PARAM_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param rankings PARAM_DESCRIPTION, Default: NULL
#' @param scale PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'OPLS'
#' @param annotation PARAM_DESCRIPTION
#' @param classes PARAM_DESCRIPTION, Default: NULL
#' @param TargetState PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname DiscriminantAnalysis
#' @export 
DiscriminantAnalysis <- function(mat, method = "OPLS", annotation,classes=NULL,TargetState=NULL) {
  if (is.null(classes)){
    classes = unique(annotation)
  }
  if (is.null(TargetState)){
    TargetState = classes
  }
  
  compare.res.df <- data.frame()
  
  for (population in TargetState) {
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
    
    for (sub.population in classes[classes != population]) {
      idx.sub.population <- which(ident.df$ident == sub.population)
      idx.random <- sample.int(n = length(idx.sub.population), size = min.n, replace = FALSE)
      idx.sub.population.random <- idx.sub.population[idx.random]
      ident.df[idx.sub.population.random,]$new.ident <- "ident.2"
    }
    
    # Step 2: Discriminant analysis
    marker.res <- FindMarkersOPLS(
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param annotation PARAM_DESCRIPTION
#' @param ident.1 PARAM_DESCRIPTION
#' @param ident.2 PARAM_DESCRIPTION, Default: NULL
#' @param scale_weights PARAM_DESCRIPTION, Default: TRUE
#' @param scale_vipVn PARAM_DESCRIPTION, Default: TRUE
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[ropls]{opls}}
#' @rdname FindMarkersOPLS
#' @export 
#' @importFrom ropls opls
FindMarkersOPLS <- function(mat, annotation, ident.1, ident.2 = NULL, scale_weights = TRUE, scale_vipVn = TRUE, ...) {
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
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


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param clusters.DA PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
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


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname range01
#' @export 
range01 = function(x){(x-min(x, na.rm=T))/(max(x,na.rm=T)-min(x, na.rm=T))}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param split PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION, Default: NULL
#' @param annotation PARAM_DESCRIPTION, Default: NULL
#' @param classes PARAM_DESCRIPTION, Default: NULL
#' @param features PARAM_DESCRIPTION, Default: NULL
#' @param reduction.name PARAM_DESCRIPTION, Default: NULL
#' @param pairwise PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RunOPLSDA
#' @export 
RunOPLSDA = function(
    data = NULL, 
    annotation = NULL,
    classes = NULL,
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
  
  #Store original Seurat object idents
  idents.backup = annotation
  
  # Get the list of classes to compare with the background
  
  if (is.null(classes)){
    classes = unique(annotation)
  }
  
  if(!isTRUE(pairwise)){
    results.list = list()
    for(population in classes){
      idents.tmp = ifelse(idents.backup == population,"ident.1","ident.2")
      
      opls.res = FindMarkersOPLS.ext.wo.seurat(data = data,
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
    tests.for = combn(length(classes), 2, FUN = NULL, simplify = TRUE)
    test.rev = rbind(tests.for[2,],tests.for[1,])
    tests = cbind(tests.for,test.rev)
    
    results.list = list()
    for(i in 1:ncol(tests)){
      population.1 = as.character(classes[tests[1,i]])
      population.2 = as.character(classes[tests[2,i]])
      # Step 1: Make the background balanced
      # Get ident.df
      ident.df = data.frame(ident=idents.backup)
      
      ident.df$new.ident = "ident.3"
      ident.df[ident.df$ident == population.1,]$new.ident = "ident.1"
      ident.df[ident.df$ident == population.2,]$new.ident = "ident.2"
      
      opls.res = FindMarkersOPLS.ext.wo.seurat(data = data,
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
  
  saveRDS(file="OPLSDA.results.list.Rds",results.list)
  
  return(results.list)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seurat.object PARAM_DESCRIPTION, Default: NULL
#' @param assay PARAM_DESCRIPTION, Default: NULL
#' @param slot PARAM_DESCRIPTION, Default: 'data'
#' @param annotation PARAM_DESCRIPTION, Default: NULL
#' @param features PARAM_DESCRIPTION, Default: NULL
#' @param classes PARAM_DESCRIPTION, Default: NULL
#' @param reduction.name PARAM_DESCRIPTION, Default: NULL
#' @param pairwise PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Seurat]{reexports}}
#' @rdname RunOPLSDA_Seurat
#' @export 
#' @importFrom Seurat Idents CreateDimReducObject
RunOPLSDA_Seurat = function(seurat.object = NULL, 
                            assay = NULL, 
                            slot = "data",
                            annotation = NULL,
                            features = NULL,
                            classes = NULL,
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
    for(population in classes){
      # Step 1: Make the background balanced
      # Get ident.df
      ident.df = seurat.object@meta.data[, annotation, drop = F]
      colnames(ident.df) = c("ident")
      
      ident.df$new.ident = "ident.2"
      ident.df[ident.df$ident == population,]$new.ident = "ident.1"
      
      # Set the Idents
      seurat.object@meta.data$new.ident = ident.df[rownames(seurat.object@meta.data),]$new.ident
      Seurat::Idents(seurat.object) = seurat.object$new.ident
      
      opls.res = FindMarkersOPLS.ext(seurat.object = seurat.object,
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
    tests.for = combn(length(classes), 2, FUN = NULL, simplify = TRUE)
    test.rev = rbind(tests.for[2,],tests.for[1,])
    tests = cbind(tests.for,test.rev)
    
    results.list = list()
    for(i in 1:ncol(tests)){
      population.1 = as.character(classes[tests[1,i]])
      population.2 = as.character(classes[tests[2,i]])
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
      
      opls.res = FindMarkersOPLS.ext(seurat.object = seurat.object,
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
  
  saveRDS(file="OPLSDA.results.list.Rds",results.list)
  
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


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seurat.object PARAM_DESCRIPTION
#' @param ident.1 PARAM_DESCRIPTION
#' @param ident.2 PARAM_DESCRIPTION, Default: NULL
#' @param features PARAM_DESCRIPTION, Default: NULL
#' @param assay PARAM_DESCRIPTION, Default: NULL
#' @param slot PARAM_DESCRIPTION, Default: 'data'
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[ropls]{opls}}
#' @rdname FindMarkersOPLS.ext
#' @export 
#' @importFrom Seurat Idents
#' @importFrom ropls opls
FindMarkersOPLS.ext = function(seurat.object, 
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param annotation PARAM_DESCRIPTION
#' @param ident.1 PARAM_DESCRIPTION
#' @param ident.2 PARAM_DESCRIPTION, Default: NULL
#' @param features PARAM_DESCRIPTION, Default: NULL
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[ropls]{opls}}
#' @rdname FindMarkersOPLS.ext.wo.seurat
#' @export 
#' @importFrom ropls opls
FindMarkersOPLS.ext.wo.seurat = function(data,
                                         annotation,       
                                         ident.1, 
                                         ident.2 = NULL,
                                         features = NULL, ...) {
  
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




#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param annotation PARAM_DESCRIPTION
#' @param k.range PARAM_DESCRIPTION, Default: c(2, 30)
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ModulonIdent
#' @export 
ModulonIdent = function(data, annotation, k.range = c(2,30)){
  library(dplyr)
  library(AUCell)
  # Generate TF clusters at different clustering resolutions
  clusters <- FindFeatureCluster(
    distance.method <- "euclidean",
    clustering.method <- "complete",
    mat = data, 
    annotation = annotation, 
    cluster.nums = k.range)
  # Calculate cluster's signatures scores (AUC)
  clusters.AUC <- GetSignatureAUC(
    mat = data,
    signature.list = clusters)
  # OPLS discriminant analysis
  clusters.DA <- DiscriminantAnalysis(
    mat = clusters.AUC,
    method = "OPLS", 
    annotation = annotation)
  # Calculate the global discriminant score (GDS)
  GDS = CalculateGDS(clusters.DA=clusters.DA)
  
  # Identify the best clustering resolution
  Best.k = GDS[which.max(GDS$avg_max_weight),'NumCluster']
  
  # Generate TF clusters at the optimal resolution
  modulons <- FindFeatureCluster(
    distance.method <- "euclidean",
    clustering.method <- "complete",
    mat = data, 
    annotation = annotation, 
    cluster.nums = Best.k)
  
  return(list(Modulons=modulons,GDS=GDS))
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param modulons PARAM_DESCRIPTION
#' @param annotation PARAM_DESCRIPTION
#' @param TargetState PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ModulonSelect
#' @export 
ModulonSelect = function(data, modulons, annotation,TargetState=NULL){
  if(is.null(TargetState)){
    TargetState = unique(annotation)
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
    TargetState = TargetState
  )
  
  modulons.DA.split = split(modulons.DA,modulons.DA$ident.1)
  
  Selected_Modulon = lapply(modulons.DA.split,function(x){
    input.tmp = x
    input.tmp$label = factor(as.character(strsplit2(input.tmp$feature,'[.]')[,2]),levels = as.character(strsplit2(input.tmp$feature,'[.]')[,2]))
    input.tmp = input.tmp[order(input.tmp$weightStarMN,decreasing = T),]
    input.tmp$order = factor(c(1:nrow(input.tmp)))
    best.modulon = as.character(input.tmp$feature[which.max(input.tmp$weightStarMN)])
    return(best.modulon)
  })
  return(list(Modulon_AUC = modulons.AUC, Modulon_DA = modulons.DA.split, Selected_Modulon = Selected_Modulon))
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Regulons PARAM_DESCRIPTION
#' @param Modulons PARAM_DESCRIPTION
#' @param ExpMat PARAM_DESCRIPTION
#' @param annotation PARAM_DESCRIPTION
#' @param TargetState PARAM_DESCRIPTION
#' @param TargetModulon PARAM_DESCRIPTION
#' @param CombSize PARAM_DESCRIPTION
#' @param Weights PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ModulonPert
#' @export 
ModulonPert = function(Regulons, Modulons, ExpMat, annotation, TargetState, TargetModulon, CombSize, Weights=NULL){
  # Collect TF-gene regulation weights
  
  if(!is.null(Weights)){
    Weights$weight.norm = range01(Weights$weight)
    rownames(Weights)= paste0(Weights$TF,' -> ',Weights$Target)
  }  
  
  # Evaluate gene importance for a given cell state
  modulon.query = TargetModulon
  query = Modulons[[modulon.query]]
  targets=unique(unlist(Regulons[query]))
  OPLSDA.results=RunOPLSDA(
    data = ExpMat[targets,], 
    annotation = annotation,
    classes = TargetState,
    features = targets,
    reduction.name = 'OPLSDA_RNA_M_Targets',
    pairwise = F
  )
  
  weightStarMN.tmp = OPLSDA.results[["CD8_Tex"]]@weightStarMN
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
