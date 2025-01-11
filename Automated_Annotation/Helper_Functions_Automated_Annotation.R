############################################
# Load Libraries
############################################

load_libraries <- function(libraries) {
  for (lib in libraries) {
    message("Loading library: ", lib)
    if (!requireNamespace(lib, quietly = TRUE)) {
      message("Installing missing library: ", lib)
      install.packages(lib)
    }
    library(lib, character.only = TRUE)
  }
}

setup_python_environment <- function() {
  library(reticulate)
  
  env_name <- "renv/python/virtualenvs/renv-python-3.11.9"
  required_packages <- c("numpy==1.24.3", "scipy==1.10.1", "magic-impute==3.0.0", "leidenalg", "igraph")
  
  # Check if the Python environment exists
  if (!virtualenv_exists(env_name)) {
    stop("Python virtual environment not found. Please ensure renv.lock is configured correctly and run renv::restore().")
  }
  
  # Activate the environment
  use_virtualenv(env_name, required = TRUE)
  
  # Check which packages are installed
  installed_packages <- py_list_packages()
  missing_packages <- required_packages[!sapply(strsplit(required_packages, "=="), function(x) x[1]) %in% installed_packages$package]
  
  if (length(missing_packages) > 0) {
    cat("Some required packages are missing. Installing them now...\n")
    py_install(missing_packages, pip = TRUE)
  } else {
    cat("All required packages are already installed.\n")
  }
  
  # Import and verify all required modules
  modules_to_import <- c('numpy', 'scipy', 'magic', 'leidenalg', 'igraph')
  imported_modules <- list()
  
  for (module in modules_to_import) {
    tryCatch({
      imported_modules[[module]] <- import(module)
      cat(sprintf("Successfully imported %s\n", module))
    }, error = function(e) {
      cat(sprintf("Failed to import %s: %s\n", module, e$message))
    })
  }
  
  cat("Setup complete.\n")
  
  # Return the imported modules
  return(imported_modules)
}


############################################
# Automated pre-processing functions for mtx
############################################

read_all_mtx_and_create_seurat <- function(directory) {
  # Initialize a vector to store Seurat object names
  pattern_files <- c()
  
  # List all files in the main directory (non-recursive)
  main_files <- list.files(directory, pattern = "(_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz)$", full.names = TRUE, recursive = FALSE)
  
  # Find base names for files in the main directory
  main_base_names <- unique(sapply(strsplit(basename(main_files), "_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz"), `[`, 1))
  
  # Process datasets in the main directory
  for (base in main_base_names) {
    mtx_file <- file.path(directory, paste0(base, "_matrix.mtx.gz"))
    features_file <- file.path(directory, paste0(base, "_features.tsv.gz"))
    cells_file <- file.path(directory, paste0(base, "_barcodes.tsv.gz"))
    
    # Check if all files exist
    if (file.exists(mtx_file) && file.exists(features_file) && file.exists(cells_file)) {
      # Read the matrix
      expression_matrix <- ReadMtx(mtx = mtx_file, features = features_file, cells = cells_file)
      
      # Create a Seurat object
      seurat_object <- CreateSeuratObject(counts = expression_matrix, min.cells = 0, min.features = minFeatures)
      seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
      seurat_object <- subset(seurat_object, subset = nFeature_RNA > minFeatures & percent.mt < Mito_Cutoff & nCount_RNA > minCounts)
      
      # Assign the Seurat object to a new variable with the base name
      assign(base, seurat_object, envir = .GlobalEnv)
      
      # Add the base name to pattern_files vector
      pattern_files <- c(pattern_files, base)
    } else {
      warning(paste("Files for dataset", base, "not found. Skipping."))
    }
  }

  # Now handle subdirectories
  subdirs <- list.dirs(directory, recursive = FALSE, full.names = TRUE)
  
  for (subdir in subdirs) {
    # List files in the subdirectory
    subdir_files <- list.files(subdir, pattern = "(matrix.mtx.gz|features.tsv.gz|barcodes.tsv.gz)$", full.names = TRUE, recursive = FALSE)
    
    # Subfolder name (basename of the subdir)
    subdir_name <- basename(subdir)
    
    # Construct the file paths for matrix, features, and barcodes
    mtx_file <- file.path(subdir, "matrix.mtx.gz")
    features_file <- file.path(subdir, "features.tsv.gz")
    cells_file <- file.path(subdir, "barcodes.tsv.gz")
    
    # Check if all files exist in the subdirectory
    if (file.exists(mtx_file) && file.exists(features_file) && file.exists(cells_file)) {
      # Read the matrix
      expression_matrix <- ReadMtx(mtx = mtx_file, features = features_file, cells = cells_file)
      
      # Create a Seurat object
      seurat_object <- CreateSeuratObject(counts = expression_matrix, min.cells = 0, min.features = minFeatures)
      seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
      seurat_object <- subset(seurat_object, subset = nFeature_RNA > minFeatures & percent.mt < Mito_Cutoff & nCount_RNA > minCounts)
      
      # Assign the Seurat object to a new variable with the subfolder name
      assign(subdir_name, seurat_object, envir = .GlobalEnv)
      
      # Add the subfolder name to pattern_files vector
      pattern_files <- c(pattern_files, subdir_name)
    } else {
      warning(paste("Files for dataset in subdirectory", subdir_name, "not found or incomplete. Skipping."))
    }
  }
  
  # Assign pattern_files to global environment
  assign("pattern_files", pattern_files, envir = .GlobalEnv)
  
  # Return the vector of Seurat object names
  return(pattern_files)
}


Process_Seurat <- function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = FALSE)  # Set verbose to FALSE to reduce console output
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, verbose = FALSE)  # Set verbose to FALSE to reduce console output
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)  # Set verbose to FALSE to reduce console output
  return(obj)
}

#####################################
# Doublet removal functions
#####################################

Assign_ExpectedDoublet_variable <- function(val) {
  
  # Initialize expected_doublet value
  expected_doublet<- NA
  
  # Check the range of value_A and assign the corresponding value to expected_doublet value based on 10x information
  if (val >= 0 & val <= 500) {
    expected_doublet <- 0.004
  } else if (val >= 501 & val <= 1000) {
    expected_doublet <- 0.008
  } else if (val >= 1001 & val <= 2000) {
    expected_doublet <- 0.016
  } else if (val >= 2001 & val <= 3000) {
    expected_doublet <- 0.023
  } else if (val >= 3001 & val <= 4000) {
    expected_doublet <- 0.031
  } else if (val >= 4001 & val <= 5000) {
    expected_doublet <- 0.039
  } else if (val >= 5001 & val <= 6000) {
    expected_doublet <- 0.046
  } else if (val >= 6001 & val <= 7000) {
    expected_doublet <- 0.054
  } else if (val >= 6001 & val <= 8000) {
    expected_doublet <- 0.061
  } else if (val >= 6001 & val <= 9000) {
    expected_doublet <- 0.069
  } else if (val >= 9001) {
    expected_doublet <- 0.076
  }
  
  # Return the expected doublet value
  return(expected_doublet)
}

# Function to process Seurat object
Process_Doublet_Removal <- function(seurat_obj, seurat_names) {
  # Your existing operations
  sweep.res <- paramSweep_v3_alt(seurat_obj, PCs = 1:20, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  val <- length(seurat_obj@assays[["RNA"]]@counts@p)
  expected_doublet <- Assign_ExpectedDoublet_variable(val)
  print(expected_doublet)
  
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(expected_doublet*nrow(seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  print(pK)
  print(nExp_poi)
         
  seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  seurat_obj@meta.data$Singlet <- seurat_obj@meta.data[[8]] == "Singlet"
  seurat_obj <- subset(x = seurat_obj, cells = WhichCells(seurat_obj, expression = Singlet == 'TRUE'))
 
  max_row <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  arrange(desc(pK)) %>%
  head(1)

  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
  x_at_max <- as.numeric(as.character(max_row$pK))

# Create the plot
  
  range_x <- range(bcmvn$pK)
  positions <- seq(range_x[1], range_x[2], length.out = 6)
  
  plot_obj <- ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    ggtitle(seurat_names) +
    geom_point() +
    geom_line() +
    geom_vline(aes(xintercept = x_at_max), linetype="dashed", color = "black") +
    geom_vline(xintercept = positions, color="grey50", linetype="solid") +
    scale_x_continuous(breaks = c(min(bcmvn$pK), x_at_max, max(bcmvn$pK))) +
    theme_minimal() + 
    theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
  
  output_filename <- paste0(output_directory, "Clustering_01/DoubletFinderFunction_", seurat_names, ".png")
  ggsave(filename = output_filename, plot = plot_obj, width = 6, height = 4, bg = 'white')
  
  paste0(output_directory, output_filename)

  return(seurat_obj)
  
}

paramSweep_v3_alt <- function (seu, PCs = 1:10, sct = FALSE, num.cores = 1) 
{
  require(Seurat)
  require(fields)
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
      10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    print(dim(data))
    n.real.cells <- ncol(data)
    print(n.real.cells)
  }
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3, 
      n.real.cells, real.cells, pK, pN, data, orig.commands, 
      PCs, sct, mc.cores = num.cores)
    stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3, 
      n.real.cells, real.cells, pK, pN, data, orig.commands, 
      PCs, sct)
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
      sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

#####################################
# Removal of decontaminated RNA
#####################################

runDecontX <- function(seurat_obj, seed=1){
  counts <- GetAssayData(object = seurat_obj, slot = "counts")
  clusters <- Idents(seurat_obj) %>% as.numeric()

  # Run on only expressed genes
  x <- counts[rowSums(counts)>0,]
  message(sprintf("Running decontX on %s cells with %s non-zero genes...", dim(x)[2], dim(x)[1]))
  decon <- decontX(x, z=clusters, verbose=TRUE, seed=seed)

  # Save desired information back to Seurat Object
  # We will place the estimated 'decontaminated counts' in place of the original counts ('RNA')
  # and keep the original counts as a separate assay called 'origCounts'
  newCounts <- decon$decontXcounts
  # Add back unexpressed genes and sort according to original counts
  newCounts <- rbind(newCounts, counts[rowSums(counts)==0,])[rownames(counts),]
  seurat_obj[["RNA"]]@counts <- as(round(newCounts), "sparseMatrix")
  seurat_obj$estConp <- decon$contamination # Estimated 'contamination proportion, 0 to 1'

  return(seurat_obj)
}

Remove_lowRNA <- function(seurat_obj){
	nPreDecon <- dim(seurat_obj)[2]
	seurat_obj <- subset(seurat_obj, subset = (nFeature_RNA > minFeatures & nCount_RNA > minCounts))
	nPostDecon <- dim(seurat_obj)[2]
	
	message(sprintf("Removed %s cells with too little RNA after decontamination. %s cells remaining.", 
		nPreDecon - nPostDecon, nPostDecon))
	
	return(seurat_obj)
	
}

#####################################
# Cluster scRNA using iterative LSI
#####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(Matrix)
})

################
# Iterative LSI
################

getVarGenes <- function(mat, nvar = 2000, blacklist = NULL){
  # Get the top nvar variable genes present in mat (a gene x sample/cell matrix)
  # If blacklist is present, do not return any genes in the blacklist
  if(!is.null(blacklist)){
    ncount <- nrow(mat)
    mat <- mat[!rownames(mat) %in% blacklist,]
    message(sprintf("Removed %s genes overlapping blacklist prior to selecting variable genes...", ncount - nrow(mat)))
  }
  if(is(mat, "sparseMatrix")){
    varGenes <- rownames(mat)[head(order(sparseMatrixStats::rowVars(mat), decreasing = TRUE), nvar)]
  }else{
    varGenes <- rownames(mat)[head(order(matrixStats::rowVars(mat), decreasing = TRUE), nvar)]
  }
  return(varGenes)
}

sparseLogX <- function(spmat, logtype="log2", scale=FALSE, scaleFactor=10^4){
  stopifnot(any(logtype == c("log", "log2", "log10")))

  if(scale == TRUE){
      spmat <- t(t(spmat)/Matrix::colSums(spmat)) * scaleFactor
  }

  if(is(spmat, "sparseMatrix")){
    matsum <- summary(spmat) # Get the sparse matrix summary
    if(logtype == "log"){
      logx <- log(matsum$x + 1) 
    }else if(logtype == "log2"){
      logx <- log2(matsum$x + 1) 
    }else{
      logx <- log10(matsum$x + 1)
    }
    logmat <- sparseMatrix(i = matsum$i, j = matsum$j, # convert back to sparse matrix
                           x = logx, dims = dim(spmat),
                           dimnames = dimnames(spmat))
  }else{
    if(logtype == "log"){
      logmat <- log(spmat + 1) 
    }else if(logtype == "log2"){
      logmat <- log2(spmat + 1) 
    }else{
      logmat <- log10(spmat + 1) 
    }
  }
  return(logmat)
}


runLSI <- function(mat, nComponents, binarize = FALSE){
  # TF-IDF LSI adapted from Jeff Granja, who adapted from flyATAC (i.e. Cusanovich et al. 2018)
  #
  # Calculate the Term Frequency - Inverse Document Frequency (TF-IDF) for a feature x cell counts
  # matrix, then calculate the Singular Value Decomposition of that matrix, which is then used as
  # input for Surat's SNN clustering
  if(binarize){
    message("Binarizing matrix...")
    # The 'x' slot of the dgCMatrix class contains the non-zero elements of the matrix
    mat@x[mat@x > 0] <- 1
  }
  #Calculate RowSums and ColSums
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)

  # Calculate TF-IDF
  message(sprintf("Calculating TF-IDF with %s features (terms) and %s cells (documents)...", nrow(mat), ncol(mat)))
  start <- Sys.time()
  scaleTo <- 10^4
  tf <- t(t(mat) / colSm)
  idf <- as(ncol(mat) / rowSm, "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% tf
  # Log transform TF-IDF
  tfidf <- sparseLogX(tfidf, logtype="log", scale=TRUE, scaleFactor=scaleTo)
  
  # Clean up
  rm(tf)
  rm(idf)
  invisible(gc())
  
  # Calculate SVD for LSI
  message("Calculating SVD for LSI...")
  svd <- irlba::irlba(tfidf, nv=nComponents, nu=nComponents)
  svdDiag <- Matrix::diag(x=svd$d)
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC", seq_len(ncol(matSVD)))
  
  # Return matSVD and svd
  message(sprintf("LSI complete: %s minutes", round(difftime(Sys.time(), start, units="mins"), 3)))
  if(is.null(rownames(mat))){
    rownames(mat) <- 1:nrow(mat)
  }
  return(
    list(
        matSVD = matSVD, 
        rowSm = rowSm, 
        colSm = colSm, 
        svd = svd, 
        binarize = binarize
        )
    )
}

assign_state <- function(seurat_obj, tags = NULL) {
  # Extract the orig.ident column
  orig_ident <- seurat_obj$orig.ident

  # Function to extract the letter after the underscore
  extract_letter <- function(x) {
    parts <- strsplit(x, "_")[[1]]
    if (length(parts) > 1) {
      return(parts[2])
    } else {
      return(NA)
    }
  }

  # Apply the extraction function to orig.ident
  letters <- sapply(orig_ident, extract_letter)

  # If tags is provided and is a named list, use it for mapping
  if (!is.null(tags) && is.list(tags) && !is.null(names(tags))) {
    disease_state <- sapply(letters, function(l) {
      if (l %in% names(tags)) {
        return(tags[[l]])
      } else {
        return(l)  # Return the letter itself if no matching tag
      }
    })
  } else {
    # If no tags provided or invalid tags, use the letters directly
    disease_state <- letters
  }

  # Add the new column to the Seurat object
  seurat_obj$diseaseStatus <- disease_state

  return(seurat_obj)
}



###########
# Magic Alt
###########

magic_alt <- function(
  data,
  genes = NULL,
  knn = 5,
  knn.max = NULL,
  decay = 1,
  t = 3,
  npca = 100,
  init = NULL,
  t.max = 20,
  knn.dist.method = 'euclidean',
  verbose = 1,
  n.jobs = 1,
  seed = NULL,
  # deprecated args
  k=NULL, alpha=NULL,
  ...
) {
  # check installation
  if (!reticulate::py_module_available(module = "magic") || (is.null(pymagic))) load_pymagic()
  # check for deprecated arguments
  if (!is.null(k)) {
    message("Argument k is deprecated. Using knn instead.")
    knn <- k
  }
  if (!is.null(alpha)) {
    message("Argument alpha is deprecated. Using decay instead.")
    decay <- alpha
  }
  knn <- as.integer(x = knn)
  t.max <- as.integer(x = t.max)
  n.jobs <- as.integer(x = n.jobs)
  if (is.numeric(x = npca)) {
    npca <- as.integer(x = npca)
  } else if (!is.null(x = npca) && is.na(x = npca)) {
    npca <- NULL
  }
  if (is.numeric(x = decay)) {
    decay <- as.double(x = decay)
  } else if (!is.null(x = decay) && is.na(x = decay)) {
    decay <- NULL
  }
  if (is.numeric(x = t)) {
    t <- as.integer(x = t)
  } else if (is.null(x = t) || is.na(x = t)) {
    t <- 'auto'
  }
  if (is.numeric(x = seed)) {
    seed <- as.integer(x = seed)
  } else if (!is.null(x = seed) && is.na(x = seed)) {
    seed <- NULL
  }
  if (is.numeric(x = verbose)) {
    verbose <- as.integer(x = verbose)
  }
  if (!methods::is(object = data, "Matrix")) {
    data <- as.matrix(x = data)
  }
  if (is.null(genes)) {
    gene_names <- colnames(x = data)
} else if (length(genes) == 1 && genes %in% c("all_genes", "pca_only")) {
    # Handle special conditions 
    if(genes == "all_genes") {
        gene_names <- colnames(x = data)
    } else if(genes == "pca_only") {
        gene_names <- paste0("PC", 1:npca)
    }
} else if (is.character(genes)) {
    # If genes is a character vector of gene names
    if (!all(genes %in% colnames(x = data))) {
      warning(paste0("Genes ", genes[!(genes %in% colnames(data))], " not found.", collapse = ", "))
    }
    gene_indices <- which(colnames(data) %in% genes)
    gene_names <- colnames(data)[gene_indices]
    genes <- as.integer(gene_indices - 1)
} else if (is.numeric(genes)) {
    # If genes is numeric, adjust indices if needed
    gene_names <- colnames(x = data)[genes]
    genes <- genes - 1
}
  # store parameters
  params <- list(
    "data" = data,
    "knn" = knn,
    "knn.max" = knn.max,
    "decay" = decay,
    "t" = t,
    "npca" = npca,
    "knn.dist.method" = knn.dist.method
  )
  # use pre-initialized values if given
  operator <- NULL
  if (!is.null(x = init)) {
    if (!methods::is(init, "magic")) {
      warning("object passed to init is not a phate object")
    } else {
      operator <- init$operator
      operator$set_params(
        knn = knn,
        knn_max = knn.max,
        decay = decay,
        t = t,
        n_pca = npca,
        knn_dist = knn.dist.method,
        n_jobs = n.jobs,
        random_state = seed,
        verbose = verbose
      )
    }
  }
  if (is.null(x = operator)) {
    operator <- pymagic$MAGIC(
      knn = knn,
      knn_max = knn.max,
      decay = decay,
      t = t,
      n_pca = npca,
      knn_dist = knn.dist.method,
      n_jobs = n.jobs,
      random_state = seed,
      verbose = verbose
    )
  }
  result <- operator$fit_transform(
    data,
    genes = genes,
    t_max = t.max
  )
  colnames(x = result) <- gene_names
  rownames(x = result) <- rownames(data)
  result <- as.data.frame(x = result)
    result <- list(
      "result" = result,
      "operator" = operator,
      "params" = params
    )
  class(x = result) <- c("magic", "list")
  return(result)
}

groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

#########################
# Matrix Helper Functions
#########################

getFreqs <- function(x){
  # Return a named vector of frequencies of x
  tab <- table(x) %>% as.data.frame.table()
  frqs <- tab$Freq
  names(frqs) <- tab[,1]
  frqs[order(frqs, decreasing=TRUE)]
}

invertList <- function(lst){
  # Swap names and values in list
  split(rep(names(lst), lengths(lst)), unlist(lst))
}

# Helper function for getting the mean of matrix groups
groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}





# Helper function for applying a custom function to matrix groups
groupFun <- function (mat, fun, groups = NULL, ...){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        # Apply custom function to columns matching group
        fun(mat[, which(groups == x), drop = F], ...)
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}


unmelt <- function(long_df, row_col, col_col, val_col){
  # 'unmelt' a long form data.frame
  # Requires columns of long_df to be named
  ##################################
  # long_df = long-format data.frame
  # row_col = name of column with id's that will become rows
  # col_col = name of column with id's that will become columns
  # val_col = name of column with values that will fill in matrix
  wide_df <- pivot_wider(long_df, id_cols=all_of(c(row_col, col_col)), names_from=all_of(col_col), values_from=all_of(val_col)) %>% as.data.frame()
  rownames(wide_df) <- wide_df[,1]
  wide_df <- wide_df[,2:ncol(wide_df)]
  return(wide_df)
}

unmelt_alt <- function(long_df, row_col, col_col, val_col){
  wide_df <- long_df %>%
    pivot_wider(id_cols=row_col, 
                names_from=col_col, 
                values_from=val_col) %>% 
    as.data.frame()

  rownames(wide_df) <- wide_df[[row_col]]
  wide_df <- wide_df[,-1]
  return(wide_df)
}

prettyOrderMat <- function(mat, scale=TRUE, cutOff=1, lmat=NULL, clusterCols=TRUE){
  # Reorder mat in a prettier way for plotting
  # Adapted from Jeff's ArchR .binarySort
  ###################################
  # mat = matrix (like) object to sort
  # scale = should mat be scaled before building logical mat
  # cutOff = cutoff for lmat
  # lmat = logical matrix for ordering rows (binary sorting)
  # clusterCols = should columns be clustered?
  mat <- as.matrix(mat)

  if(is.null(lmat)){
    # Compute row Z-scores
    if(scale){
      lmat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    }else{
      lmat <- mat
    }
    # Logical matrix of values passing cutoff 
    lmat <- lmat >= cutOff
  }

  # Transpose:
  mat <- t(mat)
  lmat <- t(lmat)

  # Identify column ordering:
  if(clusterCols){
    hc <- hclust(dist(mat))
    colIdx <- hc$order
    mat <- t(mat[colIdx,])
    lmat <- t(lmat[colIdx,])
  }else{
    mat <- t(mat)
    lmat <- t(lmat)
    hc <- NULL
  }

  # Identify row ordering:
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  mat <- mat[rowIdx,]

  return(list(mat=mat, hclust=hc))
}


fractionXbyY <- function(x, y, add_total=FALSE, xname="x", yname="y", ylab="proportion"){
  # Create a dataframe and count occurrences of x within y groups
  XbyYdf <- data.frame("xgroup" = x, "ygroup" = y) %>%
    group_by(xgroup, ygroup) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    pivot_wider(names_from=xgroup, values_from=n, values_fill=list(n=0)) %>%
    as.data.frame()
  
  # Set rownames to the first column and remove it from the dataframe
  rownames(XbyYdf) <- XbyYdf[, 1]
  XbyYmat <- as.matrix(XbyYdf[, -1])
  
  # Check if there is only one cluster
  if (ncol(XbyYmat) == 1) {
    # Ensure correct dimensions if only one cluster
    XbyYmat <- matrix(XbyYmat, nrow=1, dimnames=list(NULL, colnames(XbyYdf[, -1])))
    rownames(XbyYmat) <- rownames(XbyYdf)
  }
  
  # If add_total is TRUE, add a row with column sums
  if (add_total) {
    total_row <- colSums(XbyYmat)
    XbyYmat <- rbind(XbyYmat, total_row)
    rownames(XbyYmat)[nrow(XbyYmat)] <- "total"
  }
  
  # Calculate proportions and reshape the matrix
  XbyYmat <- (XbyYmat / rowSums(XbyYmat, na.rm = TRUE)) %>% reshape2::melt()
  colnames(XbyYmat) <- c(xname, yname, ylab)
  
  # Convert xname and yname columns to factors
  XbyYmat[, xname] <- as.factor(XbyYmat[, xname])
  XbyYmat[, yname] <- as.factor(XbyYmat[, yname])
  
  return(XbyYmat)
}




sparseLogX <- function(spmat, logtype="log2", scale=FALSE, scaleFactor=10^4){
  # Adapted from https://rdrr.io/github/YosefLab/FastProjectR/src/R/Utilities.R
  # Takes the log2(x + 1) of a potentially sparse matrix without creating an 
  # intermediate dense matrix
  ###################################
  # spmat = sparse matrix
  # logtype = which log to use (log, log2, log10)
  # scale = should matrix be scaled first?
  # scaleFactor = the amount to depth-normalize to
  stopifnot(any(logtype == c("log", "log2", "log10")))

  if(scale == TRUE){
      spmat <- t(t(spmat)/Matrix::colSums(spmat)) * scaleFactor
  }

  if(is(spmat, "sparseMatrix")){
    matsum <- summary(spmat) # Get the sparse matrix summary
    if(logtype == "log"){
      logx <- log(matsum$x + 1) 
    }else if(logtype == "log2"){
      logx <- log2(matsum$x + 1) 
    }else{
      logx <- log10(matsum$x + 1)
    }
    logmat <- sparseMatrix(i = matsum$i, j = matsum$j, # convert back to sparse matrix
                           x = logx, dims = dim(spmat),
                           dimnames = dimnames(spmat))
  }else{
    if(logtype == "log"){
      logmat <- log(spmat + 1) 
    }else if(logtype == "log2"){
      logmat <- log2(spmat + 1) 
    }else{
      logmat <- log10(spmat + 1) 
    }
  }
  return(logmat)
}

#####################################
# Expression matrix functions
#####################################


averageExpr <- function(counts_mat, group_vec, log=TRUE){
  # Can't use Seurat's average Expression function since it expects certain types of data
  # This function works on the raw counts matrix.
  # If log==TRUE, returns log2(group means + 1)

  # First, get counts per 10k
  cp10k <- sparseLogX(counts_mat, logtype="log2", scale=TRUE, scaleFactor=10000)

  # Get group means
  group_means <- groupMeans(cp10k, group_vec)

  # Log transform, if necessary
  if(log){
    group_means <- log2(group_means + 1)
  }

  return(group_means)
}



pctExpr <- function(mat) {
  # Check if the matrix is empty or has zero columns
  if (ncol(mat) == 0) {
    stop("The input matrix has zero columns.")
  }
  
  if (nrow(mat) == 0) {
    stop("The input matrix has zero rows.")
  }
  
  # Check if the matrix is either numeric or a sparse matrix
  if (!is.numeric(mat) && !inherits(mat, "dgCMatrix")) {
    stop("The input matrix must contain numeric or sparse matrix data.")
  }
  
  # Calculate the fraction of cells with non-zero expression for each gene
  nexpr <- apply(mat, 1, function(x) sum(x > 0))
  pct_expr <- nexpr / ncol(mat)
  
  return(pct_expr)
}


avgAndPctExpressed <- function(count_mat, groups, feature_normalize=FALSE, min_pct=NULL){
  # Create dataframe of average and percent of cells expressing each feature in matrix
  # Takes the raw counts matrix as input (features x cells)

  # First, get average expression in log2CP10k space
  message("Calculating group average expression...")
  avgExpr <- averageExpr(count_mat, groups)

  # If indicated, normalize mean expression to maximum detected
  if(feature_normalize){
    message("Normalizing to maximum...")
    # First, need to drop features with zero expression
    avgExpr <- avgExpr[rowSums(avgExpr) > 0,]
    avgExpr <- avgExpr / apply(avgExpr,1,max)
  }

  # Next, get percent of cells expressing each feature
  message("Calculating percent of cells per group expressing feature...")
  pctExprMat <- groupFun(count_mat, pctExpr, groups=groups) * 100

  # Filter genes that are too lowly expressed, if indicated
  if(!is.null(min_pct)){
    message(sprintf("Filtering features expressed by less than %s%% of cells in any cluster...", min_pct))
    pctExprMat <- pctExprMat[apply(pctExprMat, 1, function(x) max(x) > min_pct),]
  }

  # Melt each matrix and then merge
  message("Reshaping and returning...")
  expMelt <- reshape2::melt(avgExpr)
  colnames(expMelt) <- c("feature", "grp", "avgExpr")
  pctMelt <- reshape2::melt(pctExprMat)
  colnames(pctMelt) <- c("feature", "grp", "pctExpr")

  # Merge dfs:
  mergedDF <- base::merge(expMelt, pctMelt, by=c("feature", "grp"))
  return(mergedDF)
}




#####################################
# Matrix Correlation tools
#####################################

# Rcpp script for performing fast row-wise pearson correlations:

corMatrix <- function(mat, nonredundant=TRUE, subsetx=NULL, subsety=NULL, nThreads=1){
    # Calculate correlation matrix and associated statistics using Rcpp
    ###########################
    # mat: matrix with rows to be correlated
    # nonredunant: if true, will not calculate redundant correlations (requires)
    # subsetx: vector of rownames or indices to use for calculating correlations
    # subsety: vector or rownames or indices to use for calculating correlations

    # First, make sure no rows are zero
    if(any(Matrix::rowSums(mat) == 0)){
      stop("Error: matrix contains rows of all 0's!")
    }
    if(is.null(rownames(mat))){
        rownames(mat) <- 1:nrow(mat)
    }
    mat <- as.matrix(mat)

    if(is.null(subsetx)){
      subsetx <- rownames(mat)
    }
    if(is.null(subsety)){
      subsety <- rownames(mat)
    }
    xmat <- mat[subsetx,]
    ymat <- mat[subsety,]

    # Get indices to correlate
    message("Determining indices to correlate...")
    maxRows <- max(nrow(xmat), nrow(ymat))
    # idx <- combn(1:maxRows, 2) %>% t() # Native combn is actually very slow...
    idx <- RcppAlgos::comboGeneral(1:maxRows, 2, nThreads=nThreads)
    idx <- idx[idx[,1] <= nrow(xmat),]
    idx <- idx[idx[,2] <= nrow(ymat),]
    xidx <- idx[,1]
    yidx <- idx[,2]

    df <- data.frame(
        "x" = rownames(xmat)[xidx],
        "y" = rownames(ymat)[yidx]
        )
    message(sprintf("Calculating %s correlations...", nrow(df)))
    df$Correlation <- rowCorCpp(xidx, yidx, xmat, ymat)
    message("Finished. Calculating statistics...")
    df$TStat <- (df$Correlation / sqrt((pmax(1-df$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(mat)-2))) #T-statistic P-value
    df$Pval <- 2*pt(-abs(df$TStat), ncol(mat) - 2)
    df$FDR <- p.adjust(df$Pval, method = "fdr")
    df <- df[, c("x", "y", "Correlation", "FDR")]
    return(df)
}

#Generates subclusters for each broad cluster in dataset after iterative LSI

makeSubClusts <- function(obj, ident, subgroups, outdir){
  Idents(obj) <- ident
  for(subg in subgroups){
    subsubdir <- paste0(outdir, sprintf("/%s", subg))
    dir.create(subsubdir, showWarnings = FALSE, recursive = TRUE)
    subObj <- subset(obj, idents = c(subg))
    counts <- GetAssayData(object = subObj, slot = "counts")
    newObj <- CreateSeuratObject(counts = counts, project = subg, min.cells = 0, min.features = 200)
    old.meta <- subObj@meta.data
    # Drop selected columns from old metadata
    old.cols <- colnames(old.meta)
    drop.cols <- old.cols[grepl("^RNA_snn", old.cols)]
    newObj@meta.data <- old.meta[,old.cols %ni% drop.cols]
    message(sprintf("Subcluster %s has %s cells", subg, dim(newObj)[2]))
    saveRDS(newObj, file = paste0(subsubdir, "/", subg, ".rds"))
  }
}

cor2Matrices <- function(mat1, mat2, subset1=NULL, subset2=NULL){
    # Calculate row correlations and associated statistics of two distinct matrices using Rcpp
    ###########################
    # mat1: first matrix with rows to be correlated
    # mat2: second matrix with rows to be correlated
    # subset1: vector of rownames or indices to use for calculating correlations
    # subset2: vector or rownames or indices to use for calculating correlations

    # First, make sure no rows are zero
    if(any(Matrix::rowSums(mat1) == 0)){
      stop("Error: matrix 1 contains rows of all 0's!")
    }
    if(any(Matrix::rowSums(mat2) == 0)){
      stop("Error: matrix 2 contains rows of all 0's!")
    }
    if(is.null(rownames(mat1))){
        rownames(mat1) <- 1:nrow(mat1)
    }
    mat1 <- as.matrix(mat1)
    if(is.null(rownames(mat2))){
        rownames(mat2) <- 1:nrow(mat2)
    }
    mat2 <- as.matrix(mat2)

    if(is.null(subset1)){
      subset1 <- rownames(mat1)
    }
    if(is.null(subset2)){
      subset2 <- rownames(mat2)
    }
    mat1 <- mat1[subset1,]
    mat2 <- mat2[subset2,]

    # Get indices to correlate
    message("Determining indices to correlate...")
    idx <- expand.grid(rownames(mat1), rownames(mat2))
    idx1 <- match(idx[,1], rownames(mat1))
    idx2 <- match(idx[,2], rownames(mat2))

    df <- data.frame(
        "x" = rownames(mat1)[idx1],
        "y" = rownames(mat2)[idx2]
        )
    message(sprintf("Calculating %s correlations...", nrow(df)))
    df$Correlation <- rowCorCpp(idx1, idx2, mat1, mat2)
    message("Finished. Calculating statistics...")
    df$TStat <- (df$Correlation / sqrt((pmax(1-df$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(mat1)-2))) #T-statistic P-value
    df$Pval <- 2*pt(-abs(df$TStat), ncol(mat1) - 2)
    df$FDR <- p.adjust(df$Pval, method = "fdr")
    df <- df[, c("x", "y", "Correlation", "FDR")]
    return(df)
}


corPairwise <- function(mat1, mat2){
    # Calculate row correlations and associated statistics of two paired matrices using Rcpp
    # Assumes that rows in mat1 are already paired to those in mat2 (i.e. row 1 of mat 1 will be 
    # correlated to row 1 of mat 2, and so on)
    ###########################
    # mat1: first matrix with rows to be correlated
    # mat2: second matrix with rows to be correlated

    # First, make sure no rows are zero
    if(any(Matrix::rowSums(mat1) == 0)){
      stop("Error: matrix 1 contains rows of all 0's!")
    }
    if(any(Matrix::rowSums(mat2) == 0)){
      stop("Error: matrix 2 contains rows of all 0's!")
    }
    if(is.null(rownames(mat1))){
        rownames(mat1) <- 1:nrow(mat1)
    }
    mat1 <- as.matrix(mat1)
    if(is.null(rownames(mat2))){
        rownames(mat2) <- 1:nrow(mat2)
    }
    mat2 <- as.matrix(mat2)

    df <- data.frame(
        "ix" = rownames(mat1)
        )
    message(sprintf("Calculating %s correlations...", nrow(df)))
    df$Correlation <- rowCorCpp(1:nrow(mat1), 1:nrow(mat2), mat1, mat2)
    message("Finished. Calculating statistics...")
    df$TStat <- (df$Correlation / sqrt((pmax(1-df$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(mat1)-2))) #T-statistic P-value
    df$Pval <- 2*pt(-abs(df$TStat), ncol(mat1) - 2)
    df$FDR <- p.adjust(df$Pval, method = "fdr")
    df <- df[, c("ix","Correlation", "FDR")]
    return(df)
}

####################
# Plotting Functions
####################

plotClusterQC <- function(obj, subgroup, plotDir, pointSize=1.0, barwidth=0.9, sampleCmap=NULL, diseaseCmap=NULL){
  
  # Plot basic clustering plots
  # Set colormap
  qualcmap <- cmaps_BOR$stallion
  quantcmap <- cmaps_BOR$solarExtra
  namedSampCmap <- TRUE
  namedDiseaseCmap <- TRUE

  if(is.null(sampleCmap)){
    sampleCmap <- qualcmap
    namedSampCmap <- FALSE
  }
  if(is.null(diseaseCmap)){
    diseaseCmap <- qualcmap
    namedDiseaseCmap <- FALSE
  }

  ### Bar plot cluster counts ###
  tabDF <- base::table(obj$Clusters) %>% as.data.frame
  colnames(tabDF) <- c("Clusters", "count")

  pdf(paste0(plotDir, sprintf("/clusterBarPlot_%s.pdf", subgroup)))
  print(qcBarPlot(tabDF, cmap=qualcmap, barwidth=barwidth))
  dev.off()

  clustBySamp <- fractionXbyY(obj$Clusters, obj$sample, add_total=TRUE, xname="Cluster", yname="sample")

  pdf(paste0(plotDir, sprintf("/clustBySampleBarPlot_%s.pdf", subgroup)))
  print(stackedBarPlot(clustBySamp, cmap=sampleCmap, namedColors=namedSampCmap, barwidth=barwidth))
  dev.off()

  ### Stacked bar plot fraction disease in clusters ###
  clustByDisease <- fractionXbyY(obj$Clusters, obj$diseaseStatus, add_total=TRUE, xname="Cluster", yname="Disease")

  pdf(paste0(plotDir, sprintf("/clustByDiseaseBarPlot_%s.pdf", subgroup)))
  print(stackedBarPlot(clustByDisease, cmap=diseaseCmap, namedColors=namedDiseaseCmap, barwidth=barwidth))
  dev.off()

  ### Cluster UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$Clusters)
  # Randomize cells before plotting UMAP
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/cluster_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=qualcmap, point_size=pointSize))
  dev.off()

  ### Sample UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$sample)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/sample_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=sampleCmap, namedColors=namedSampCmap, point_size=pointSize))
  dev.off()

  ### Disease UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$diseaseStatus)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/disease_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=diseaseCmap, namedColors=namedDiseaseCmap, point_size=pointSize))
  dev.off()
}

plotClusterQC_subgroup <- function(obj, subgroup, plotDir, pointSize=1.0, barwidth=0.9, sampleCmap=NULL, diseaseCmap=NULL){
  
  # Plot basic clustering plots
  # Set colormap
  qualcmap <- cmaps_BOR$stallion
  quantcmap <- cmaps_BOR$solarExtra
  namedSampCmap <- TRUE
  namedDiseaseCmap <- TRUE

  if(is.null(sampleCmap)){
    sampleCmap <- qualcmap
    namedSampCmap <- FALSE
  }
  if(is.null(diseaseCmap)){
    diseaseCmap <- qualcmap
    namedDiseaseCmap <- FALSE
  }

  ### Bar plot cluster counts ###
  tabDF <- base::table(obj$FineClust) %>% as.data.frame
  colnames(tabDF) <- c("Clusters", "count")

  pdf(paste0(plotDir, sprintf("/clusterBarPlot_%s.pdf", subgroup)))
  print(qcBarPlot(tabDF, cmap=qualcmap, barwidth=barwidth))
  dev.off()

  clustBySamp <- fractionXbyY(obj$FineClust, obj$sample, add_total=TRUE, xname="Cluster", yname="sample")

  pdf(paste0(plotDir, sprintf("/clustBySampleBarPlot_%s.pdf", subgroup)))
  print(stackedBarPlot(clustBySamp, cmap=sampleCmap, namedColors=namedSampCmap, barwidth=barwidth))
  dev.off()

  ### Stacked bar plot fraction disease in clusters ###
  clustByDisease <- fractionXbyY(obj$FineClust, obj$diseaseStatus, add_total=TRUE, xname="Cluster", yname="Disease")

  pdf(paste0(plotDir, sprintf("/clustByDiseaseBarPlot_%s.pdf", subgroup)))
  print(stackedBarPlot(clustByDisease, cmap=diseaseCmap, namedColors=namedDiseaseCmap, barwidth=barwidth))
  dev.off()

  ### Cluster UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$FineClust)
  # Randomize cells before plotting UMAP
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/cluster_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=qualcmap, point_size=pointSize))
  dev.off()

  ### Sample UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$sample)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/sample_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=sampleCmap, namedColors=namedSampCmap, point_size=pointSize))
  dev.off()

  ### Disease UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$diseaseStatus)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/disease_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=diseaseCmap, namedColors=namedDiseaseCmap, point_size=pointSize))
  dev.off()
}

###Plotting configurations###
# Plotting settings

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(magrittr)
  library(ggplot2)
  library(ggrastr)
})


theme_BOR <- function(base_size=14, base_family="Helvetica", border = TRUE) {
  library(grid)
  library(ggthemes)
  # Should plots have a bounding border?
  if(border){
    panel.border <- element_rect(fill = NA, color = "black", size = 0.7)
    axis.line <- element_blank()
  }else{
    panel.border <- element_blank()
    axis.line <- element_line(color = "black", size = 0.5)
  }
  
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = panel.border,
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = axis.line,
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text()
    ))
  
}

scale_fill_BOR <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_BOR <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#---------------------------
# Colormaps
#---------------------------

cmaps_BOR <- list(
  # Many of these adapted from ArchR ColorPalettes.R by Jeff Granja or colors.R from BuenColors
  # https://github.com/GreenleafLab/ArchR/blob/master/R/ColorPalettes.R
  # https://github.com/caleblareau/BuenColors/blob/master/R/colors.R
  
  ## Sequential colormaps:
  solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', 
                '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'), #buencolors
  sunrise = c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
              "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
  horizon = c('#000075', '#2E00FF', '#9408F7', '#C729D6', '#FA4AB5', 
              '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFFF60'),
  horizonExtra =c("#000436", "#021EA9", "#1632FB", "#6E34FC", "#C732D5",
                  "#FD619D", "#FF9965", "#FFD32B", "#FFFC5A"),
  blueYellow = c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                  "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
  sambaNight = c('#1873CC','#1798E5','#00BFFF','#4AC596','#00CC00',
                  '#A2E700','#FFFF00','#FFD200','#FFA500'), #buencolors
  wolfgang_basic = c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", 
                    "#1D91C0", "#225EA8", "#253494", "#081D58"), #buencolors
  wolfgang_extra = c("#FFFFFF", "#FCFED3", "#E3F4B1", "#ABDEB6", "#60C1BF", 
                    "#2A9EC1", "#206AAD", "#243996", "#081D58"), #buencolors
  whitePurple = c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                  '#8c6bb1','#88419d','#810f7c','#4d004b'),
  whiteBlue = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
                '#3690c0','#0570b0','#045a8d','#023858'),
  whiteViolet = c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1', 
                  '#DD3497', '#AE017E', '#7A0177', '#49006A'),
  comet = c("#E6E7E8","#3A97FF","#8816A7","black"),

  flame_flame = c('#000033', '#0000A5', '#1E00FB', '#6F00FD', '#C628D6', 
    '#FE629D', '#FF9B64', '#FFD52C', '#FFFF5F'), # buencolors

  flame_short = c('#000033', '#0000A5', '#1E00FB', '#6F00FD', '#C628D6', 
    '#FE629D', '#FF9B64', '#FFD52C'), # Stop short of yellow (better for tracks, etc.)

  #7-colors
  greenBlue = c('#e0f3db','#ccebc5','#a8ddb5','#4eb3d3','#2b8cbe',
                '#0868ac','#084081'),
  
  #6-colors
  beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30","#F7962E","#FCEE2B"),
  
  #5-colors
  fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A"),
  greyMagma = c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"),
  fireworks2 = c("black", "#2488F0","#7F3F98","#E22929","#FCB31A"),
  purpleOrange = c("#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F"),
  beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30","#F7962E","#FCEE2B"),
  zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"), #wesanderson
  darjeeling = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"), #wesanderson
  rushmore = c("#E1BD6D", "#EABE94", "#0B775E","#35274A" , "#F2300F"), #wesanderson
  FantasticFox1 = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20"), #wesanderson
  BottleRocket2 = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E"), #wesanderson
  Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B"), #wesanderson
  fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A"),

  # Divergent sequential:
  coolwarm = c("#4858A7", "#788FC8", "#D6DAE1", "#F49B7C", "#B51F29"),
  brewer_yes = c("#053061", "#2971B1", "#6AACD0","#C1DDEB", "#F7F7F7", 
                  "#FACDB5", "#E58267", "#BB2933", "#67001F"), #buencolors
  brewer_celsius = c("#313695", "#5083BB", "#8FC3DD", "#D2ECF4", "#FFFFBF", 
                      "#FDD384", "#F88D51", "#DE3F2E", "#A50026"), #buencolors
  flame_blind = c("#0DB2AA", "#0AD7D3", "#00FFFF", "#B1FFFE", "#FFFFFF", 
                  "#FFA3EC", "#FF00D8", "#BD00EC", "#5F00FF"), #buencolors
  solar_flare = c('#3361A5', '#2884E7', '#1BA7FF', '#76CEFF', '#FFFFFF', 
                  '#FFE060', '#FA8E24', '#DA2828', '#A31D1D'), #buencolors
  brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', 
                '#FACDB5', '#E58267', '#BB2933', '#67001F'), #buencolors
  
  ## Qualitative colormaps:
  
  # see: https://carto.com/carto-colors/
  cartoPrism = c('#7F3C8D', '#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310', 
                  '#008695', '#CF1C90', '#F97B72', '#4B4B8F'),
  cartoSafe = c('#88CCEE', '#CC6677', '#DDCC77', '#117733', '#332288', '#AA4499', '#44AA99',
                 '#999933', '#882255', '#661100', '#6699CC'),
  cartoBold = c('#7F3C8D' ,'#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310',
                 '#008695', '#CF1C90', '#f97b72', '#4b4b8f'),
  cartoAntique = c('#855C75', '#D9AF6B', '#AF6458', '#736F4C', '#526A83', '#625377', '#68855C',
                    '#9C9C5E', '#A06177', '#8C785D', '#467378'),
  cartoPastel = c('#66C5CC', '#F6CF71', '#F89C74', '#DCB0F2', '#87C55F', '#9EB9F3', '#FE88B1',
                   '#C9DB74', '#8BE0A4', '#B497E7', '#D3B484'),
  cartoVivid = c('#E58606', '#5D69B1', '#52BCA3', '#99C945', '#CC61B0', '#24796C', '#DAA51B',
                  '#2F8AC4', '#764E9F', '#ED645A', '#CC3A8E'),
  # 15 color
  circus = c("#D52126", "#88CCEE", "#FEE52C", "#117733", "#CC61B0", "#99C945", "#2F8AC4", "#332288", 
              "#E68316", "#661101", "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74"),
  iron_man = c('#371377','#7700FF','#9E0142','#FF0080', '#DC494C',"#F88D51","#FAD510","#FFFF5F",'#88CFA4',
               '#238B45',"#02401B","#0AD7D3","#046C9A", "#A2A475", 'grey35'),
  # The following 3 were designed by Ryan Corces.
  stallion = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B", "#FEE500","#8A9FD1","#C06CAB", "#D8A767",
               "#90D5E4", "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8", "#6E4B9E","#0C727C", "#7E1416", "#E6C2DC"),
  calm = c("#7DD06F", "#844081", "#688EC1", "#C17E73", "#484125", "#6CD3A7", "#597873","#7B6FD0", "#CF4A31", "#D0CD47",
           "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D", "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736"),
  kelly = c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
            "#FF7A5C", "#53377A", "#FF8E00","#B32851", "#F4C800", "#7F180D", "#93AA00", "#593315", "#F13A13")
)

#--------------------------
# Colormap helper functions
#--------------------------

mostDifferentColors <- function(cols, n=20, colorspace="Lab", startingCols=NULL){
  stopifnot(length(cols) > n)
  rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue=255)
  
  # Convert sRGB to another colorspace (more 'perceptually uniform' colorspace, e.g. "Lab")
  rgbCols <- t(col2rgb(cols))
  conv <- grDevices::convertColor(rgbCols, from="sRGB", to=colorspace, scale.in=255)
  
  # Now select n 'furthest neighbors' colors
  # This performs an iterative procedure for picking colors that maximize
  # 'distance' to already selected colors. The first color is picked randomly.
  # If starting cols provided, add these to the list of picked cols
  if(!is.null(startingCols)){
    stConv <- grDevices::convertColor(t(col2rgb(startingCols)), from="sRGB", to=colorspace, scale.in=255)
    pickedColors <- list()
    for(i in seq_len(nrow(stConv))){
      pickedColors[[i]] <- stConv[i,]
    }
    remainingColors <- conv
  }else{
    idx <- sample(1:nrow(conv), 1)
    pickedColors <- list(conv[idx,])
    remainingColors <- conv[-idx,]
  }
  pickedLen <- length(pickedColors)
  
  # Iteratively add the furthest color from the selected colors
  for(i in seq(pickedLen, n - 1)){
    distList <- list()
    for(j in seq_along(pickedColors)){
      colJ <- pickedColors[[j]]
      distMat <- dist(rbind(colJ, remainingColors), method="euclidean") %>% as.matrix
      distList[[j]] <- distMat[2:nrow(distMat),1]
    }
    # Maximize the minimum distance between each color
    distMat <- do.call(cbind, distList)
    distMins <- apply(distMat, 1, FUN = min)
    idx <- which(max(distMins) == distMins)
    pickedColors[[i + 1]] <- remainingColors[idx,]
    remainingColors <- remainingColors[-idx,]
  }
  pickedLab <- do.call(rbind, pickedColors)
  pickedRgb <- round(grDevices::convertColor(pickedLab, from = colorspace, to = "sRGB", scale.out = 255),0)
  hex <- apply(pickedRgb, 1, rgb2hex)
  hex
}


mostSimilarColors <- function(color, colorOptions, n=5, colorspace="Lab"){
  stopifnot(length(colorOptions) > n)
  rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)

  colorOptions <- colorOptions[colorOptions != color]
  
  # Convert sRGB to another colorspace (more 'perceptually uniform' colorspace)
  rgb <- t(col2rgb(color))
  rgbCols <- t(col2rgb(colorOptions))
  fullMat <- rbind(rgb, rgbCols)
  rownames(fullMat) <- 1:nrow(fullMat)
  conv <- grDevices::convertColor(fullMat, from = "sRGB", to = colorspace, scale.in = 255)

  # Calcualte distances and pick n most similar to starting color
  distMat <- dist(conv, method = "euclidean") %>% as.matrix()
  pickedIdx <- distMat[1,2:ncol(distMat)] %>% sort() %>% head(.,n=n) %>% names() %>% as.integer()
  colorOptions[pickedIdx-1]
}


pairwiseColorInterpolations <- function(cols, colorspace = "Lab"){
  # Get all pairwise interpolations between a vector of colors
  rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
  interpolate <- function(c1, c2, colorspace){
    rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
  }
  paired <- sapply(cols, function(x) sapply(cols, function(y) interpolate(x, y, colorspace)))
  unique(as.vector(paired))
}


getColorMap <- function(cmap, n, type='qualitative'){
  stopifnot(n >= 1)
  # Return a character vector of n colors based on
  # the provided colormap. If n > length(cmap), do
  # some smart interpolation to get enough colors
  names(cmap) <- NULL # Having names on colors causes problems for some plotting routines

  if(type == 'qualitative'){
    # If qualitative colormap, do 'mostDifferent' interpolation
    if(length(cmap) < n){
      cmap <- mostDifferentColors(
        pairwiseColorInterpolations(cmap), 
        colorspace = "Apple RGB", n = n, startingCols = cmap
      )
    }

  }else{
    # Otherwise, return sequential colors based on provided palette
    colfunc <- colorRampPalette(cmap)
    cmap <- colfunc(n)
  }
  cmap[1:n]
}

plotColorMap <- function(cols){
  # Plot each of the colors in a colormap
  cols <- base::unname(cols)
  n <- length(cols)
  df <- data.frame(
    x = seq_len(n),
    y = rep(1, n),
    z = factor(seq_len(n))
  )
  p <- (
    ggplot(df, aes(x=x,y=y,color=z)) 
    + geom_tile(aes(fill=z))
    + theme_BOR()
    + scale_color_manual(values = cols)
    + scale_fill_manual(values = cols)
  )
  p
}

#-------------------
# Plotting functions
#-------------------

plotUMAP <- function(df, dataType = "qualitative", cmap = NULL, covarLabel = "", point_size=0.5, 
  namedColors=FALSE, plotTitle=NULL, colorLims=NULL, na.value="grey35", useRaster=TRUE){
  # Given a df containing the UMAP x and y coords and a third column, 
  # plot the UMAP

  if(useRaster){
    p <- (
      ggplot(df, aes(x = df[,1], y = df[,2], color = df[,3]))
      + geom_point_rast(size = point_size)
      + theme_BOR()
      + theme(
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        aspect.ratio=1
        )
      + xlab("UMAP1")
      + ylab("UMAP2")
      )
  }else{
    p <- (
      ggplot(df, aes(x = df[,1], y = df[,2], color = df[,3]))
      + geom_point(size = point_size)
      + theme_BOR()
      + theme(
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        aspect.ratio=1
        )
      + xlab("UMAP1")
      + ylab("UMAP2")
      )
  }

  # Set plot title
  if(!is.null(plotTitle)){
    p <- p + ggtitle(plotTitle)
  }else{
    p <- p + ggtitle(sprintf("n = %s", nrow(df)))
  }
  # If colormap provided, update colors
  if(!is.null(cmap)){
    if(namedColors){
      # Named colormap corresponding to discrete values in third column
      p <- p + scale_color_manual(values=cmap, limits=names(cmap), name=covarLabel, na.value=na.value)
      p <- p + guides(fill = guide_legend(title=covarLabel), 
                      colour = guide_legend(override.aes = list(size=5)))
    }else{
      # Remove names
      names(cmap) <- NULL
      if(dataType == "qualitative"){
        # check to make sure you have enough colors for qualitative mapping
        nvals <- length(unique(df[,3]))
        cmap <- getColorMap(cmap, n=nvals)
        p <- p + scale_color_manual(values=cmap, name=covarLabel, na.value=na.value)
        p <- p + guides(fill = guide_legend(title=covarLabel), 
                          colour = guide_legend(override.aes = list(size=5)))
      }else{
        if(!is.null(colorLims)){
          p <- p + scale_color_gradientn(colors=cmap, name=covarLabel, limits=colorLims, na.value=na.value)
        }else{
          p <- p + scale_color_gradientn(colors=cmap, name=covarLabel, na.value=na.value)
        }
      }
    }
  }
  p
}


plotEachQualitative <- function(umapDF, colors=NULL, defaultColor="red", bgColor="grey", cmap = camps_BOR$stallion, pointSize=0.5){
  # Plot separate UMAPs for each variable in the third column
  # Return plots in a named list
  plotList <- list("all" = plotUMAP(umapDF, dataType = "qualitative", cmap = cmap, point_size=pointSize))
  unq <- unique(umapDF[,3]) %>% as.character()
  if(is.null(colors)){
    colors <- rep(defaultColor, length(unq))
  }
  for(i in seq_along(unq)){
    s <- unq[i]
    sc <- colors[i]
    message(sprintf("Plotting %s...", s))
    pDF <- umapDF
    pDF[,3] <- ifelse(pDF[,3] == s, s, "other") %>% factor(., order = TRUE, levels = c(s, "other"))
    # Sort sample to front
    pDF <- pDF[order(pDF[,3], decreasing=TRUE),]
    nCells <- sum(pDF[,3] == s)
    pTitle <- sprintf("%s; n = %s", s, nCells)
    plotList[[s]] <- plotUMAP(pDF, dataType = "qualitative", cmap = c(sc, bgColor), plotTitle=pTitle, point_size=pointSize)
  }
  plotList
}


qcHistFilter <- function(df, cmap = NULL, bins=100, border_color="black", lower_lim=NULL, upper_lim=NULL){
  # Histogram, with cutoffs if included

  # Fix colormap if provided
  if(is.null(cmap)){
    cmap <- "blue"
  }
  
  p <- (
    ggplot(df, aes(x=df[,2]))
    + geom_histogram(bins=bins, fill=cmap, color=border_color)
    + xlab(colnames(df)[2])
    + ylab("Frequency")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 0.6,
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
  )

  x <- df[,2]

  if(!is.null(lower_lim)){
    if(is.finite(lower_lim)){
      p <- p + geom_vline(aes(xintercept=lower_lim), color="red", linetype="dashed")
      thresh <- round((sum(x < lower_lim) / length(x)) * 100, 2)
      p <- p + annotate("text",  x=-Inf, y = Inf, label = paste0(thresh, "%"), vjust=1, hjust=-1)
    }
  }
  if(!is.null(upper_lim)){
    if(is.finite(upper_lim)){
      p <- p + geom_vline(aes(xintercept=upper_lim), color="red", linetype="dashed")
      thresh <- round((sum(x > upper_lim) / length(x)) * 100, 2)
      p <- p + annotate("text",  x=Inf, y = Inf, label = paste0(thresh, "%"), vjust=1, hjust=1)
    }
  }
  p
}


qcBarPlot <- function(df, cmap = NULL, border_color="black", barwidth=0.5){
  # Plot a bar plot (df is a 2+ column dataframe with column 1 = x and column 2 = y)
  nsamp <- nrow(df)
  # Fix colormap if provided
  if(!is.null(cmap)){
    if(length(cmap) > 1){
      cmap <- getColorMap(cmap, n = nsamp)
    }
  }else{
    cmap <- "blue"
  }
  
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2]))
    + geom_bar(stat = "identity", fill = cmap, width=barwidth, color=border_color)
    + scale_fill_manual(values = cmap)
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
  )
  p
}


qcViolinPlot <- function(df, cmap = NULL, makeLog = FALSE){
  # Plot a violin plot
  nsamp <- length(unique(df[,1]))
  aspectRatio <- 6/nsamp
  # Assume that the first column is the sample and the second column is the variable of interest
  if(makeLog){
    df[,2] <- log10(df[,2])
    colnames(df)[2] <- paste0("log10 ", colnames(df)[2]) 
  }
  
  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2], color = df[,1]))
    + geom_violin(aes(fill = df[,1]))
    + geom_boxplot(width = 0.8, alpha = 0)
    + scale_color_manual(values = cmap)
    + scale_fill_manual(values = alpha(cmap, 0.2))
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  p
  
  # Adjust colors if necessary:
  if(!is.null(cmap)){
    cmap <- getColorMap(cmap, n = nsamp)
  }else{
    cmap <- rep("blue", times = nsamp)
  }
  p <- suppressMessages(p + scale_color_manual(values = cmap))
  p <- suppressMessages(p + scale_fill_manual(values = alpha(cmap, 0.3)))
  p
}


stackedBarPlot <- function(df, xcol = 1, fillcol = 2, ycol = 3, cmap = NULL, border_color="black", covarLabel = "", namedColors=FALSE, barwidth=0.5){
  # Plot a stacked bar plot
  # Expects a 'melted' dataframe/matrix as input
  nsamp <- length(unique((df[,xcol])))
  # Assume that we want to show all xaxis labels
  xID <- unique((df[,xcol]))

  # Fix colormap if provided
  if(!namedColors){
    if(!is.null(cmap)){
      cmap <- getColorMap(cmap, n = length(unique((df[,fillcol]))))
    }else{
      cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique((df[,fillcol]))))
    }
  }
  
  p <- (
    ggplot(df, aes(x=df[,xcol], y=df[,ycol], fill=df[,fillcol]))
    + geom_bar(stat = "identity", position="fill", width=barwidth, color=border_color)
    + xlab(xcol)
    + ylab(ycol)
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
  )

  # If colormap provided, update colors
  if(namedColors){
    # Named colormap corresponding to discrete values in third column
    p <- p + scale_color_manual(values = cmap, limits = names(cmap), name = covarLabel)
    p <- p + scale_fill_manual(values = cmap, limits = names(cmap), name = covarLabel)
  }else{
    p <- p + scale_color_manual(values = cmap, name = covarLabel)
    p <- p + scale_fill_manual(values = cmap, name = covarLabel)
  }
  p 
}


groupedBarPlot <- function(df, xcol=1, ycol=2, fillcol=3, cmap = NULL, border_color="black", barwidth=0.5){
  # Plot a bar plot
  nsamp <- nrow(df)
  ngroups <- length(unique(df[,fillcol]))
  # Fix colormap if provided
  if(!is.null(cmap)){
    cmap <- getColorMap(cmap, n = ngroups, type='qualitative')
  }else{
    cmap <- getColorMap(cmaps_BOR$stallion, n=ngroups, type='qualitative')
  }
  
  p <- (
    ggplot(df, aes(x=df[,xcol], y=df[,ycol], fill=df[,fillcol]))
    + geom_bar(
      stat = "identity", 
      position=position_dodge2(width=barwidth + barwidth/(ngroups*2), preserve="single"), 
      width=barwidth, color=border_color
      )
    + scale_fill_manual(values = cmap)
    + xlab(xcol)
    + ylab(ycol)
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = 8/(nsamp + nsamp/ngroups), # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
    + guides(
      fill = guide_legend(title=fillcol)
      )
  )
  p
}


dotPlot <- function(df, xcol, ycol, color_col, size_col, xorder=NULL, yorder=NULL, cmap=NULL, 
  color_label=NULL, size_label=NULL, aspectRatio=NULL, sizeLims=NULL, colorLims=NULL){
  # Plot rectangular dot plot where color and size map to some values in df
  # (Assumes xcol, ycol, color_col and size_col are named columns)

  # If neither x or y col order is provided, make something up
  # Sort df:
  if(is.null(xorder)){
    xorder <- unique(df[,xcol]) %>% sort()
  }
  if(is.null(yorder)){
    yorder <- unique(df[,ycol]) %>% sort()
  }
  if(is.null(aspectRatio)){
    aspectRatio <- length(yorder)/length(xorder) # What is the best aspect ratio for this chart?
  }
  df[,xcol] <- factor(df[,xcol], levels=xorder)
  df[,ycol] <- factor(df[,ycol], levels=yorder)
  df <- df[order(df[,xcol], df[,ycol]),]

  # Make plot:
  p <- (
    ggplot(df, aes(x=df[,xcol], y=df[,ycol], color=df[,color_col], size=ifelse(df[,size_col] > 0, df[,size_col], NA)))
    + geom_point()
    + xlab(xcol)
    + ylab(ycol)
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,0,0.25,1), "cm"), 
            aspect.ratio = aspectRatio,
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + guides(
      fill = guide_legend(title=""), 
      colour = guide_legend(title=color_label, override.aes = list(size=5)),
      size = guide_legend(title=size_label)
      )
  )
  if(!is.null(cmap)){
    if(!is.null(colorLims)){
      p <- p + scale_color_gradientn(colors=cmap, limits=colorLims, oob=scales::squish, name = "")
    }else{
      p <- p + scale_color_gradientn(colors=cmap, name = "")
    }
  }
  if(!is.null(sizeLims)){
    p <- p + scale_size_continuous(limits=sizeLims)
  }
  p
}


# This is used primarily for making colormaps for ComplexHeatmap
makeColFun <- function(start, end, cmap, midpoint = NULL){
  # Make a color ramp function from provided start and end breaks,
  # and optionally a midpoint
  cmapLen <- length(cmap)
  if(!is.null(midpoint)){
    interpolate <- function(c1, c2, colorspace = "Lab"){
      rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
    }
    if(length(cmap) %% 2 == 0){
      # Interpolate middle colors if necessary to get midpoint
      preMidIdx <- floor(cmapLen / 2)
      midCol <- interpolate(cmap[preMidIdx], cmap[preMidIdx + 1])
      cmap <- c(cmap[1:preMidIdx], midCol, cmap[(preMidIdx + 1):cmapLen])
      cmapLen <- length(cmap)
    }
    midIdx <- ceiling(cmapLen / 2)
    breaks <- c(seq(start, midpoint, length.out = midIdx), seq(midpoint, end, length.out = midIdx)[2:midIdx])
  } else {
    breaks <- seq(start, end, length.out = cmapLen)
  }
  colorRamp2(breaks, cmap)
}


# Heatmap wrapper:
BORHeatmap <- function(
  mat, # Data to plot (matrix or dataframe)
  limits = NULL, # Enforced limits for colormap (2 dimensional array)
  clusterCols = TRUE, # Should columns be clustered
  clusterRows = TRUE, # Should rows be clustered
  labelCols = FALSE, # Should columns be labeled
  labelRows = FALSE, # Should rows be labeled
  dataColors = NULL, # Colormap for plotting data
  dataColorMidPoint = NULL, # The data value to be the middle of the color map
  customRowLabel = NULL,
  customRowLabelIDs = NULL,
  customColLabel = NULL,
  customColLabelIDs = NULL,
  customLabelWidth = 0.15,
  useRaster = TRUE, # Should heatmap be rasterized
  rasterDevice = "CairoPNG",
  rasterQuality = 5, # Raster quality. Higher is {better?}
  fontSize = 6, # Font size for labels
  showColDendrogram = FALSE, # Should the column dendrogram be shown
  showRowDendrogram = FALSE, # Should the row dendrogram be shown
  borderColor = NA, # Color for lines between cells
  mapname = " ", # 'Name' to give heatmap
  legendTitle = " ", # Name of legend
  ...
){
  
  #Packages
  suppressPackageStartupMessages(require(ComplexHeatmap))
  suppressPackageStartupMessages(require(circlize))
  
  # Make sure mat is actually a matrix
  if(!is.matrix(mat)){
    message("'mat' needs to be a matrix. Converting...")
    mat <- as.matrix(mat)
  }
  
  # Prepare color function
  if(!is.null(limits)){
    ll <- limits[1]
    ul <- limits[2]
  }else{
    ll <- min(mat, na.rm=TRUE)
    ul <- max(mat, na.rm=TRUE)
  }
  # If no colormap provided, use solarExtra
  if(is.null(dataColors)){
    dataColors <- c("1"='#3361A5', "2"='#248AF3', "3"='#14B3FF', 
                    "4"='#88CEEF', "5"='#C1D5DC', "6"='#EAD397', 
                    "7"='#FDB31A', "8"= '#E42A2A', "9"='#A31D1D')
  }
  dataColFun <- makeColFun(ll, ul, dataColors, midpoint = dataColorMidPoint)
  
  message("Preparing Heatmap...")
  hm <- Heatmap(
    # Main components:
    matrix = mat,
    name = mapname,
    col = dataColFun,
    
    # Legend options:
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      legend_width = unit(1, "cm"),
      title = legendTitle
    ),
    rect_gp = gpar(col = borderColor), 
    
    # Column options:
    show_column_names = labelCols,
    cluster_columns = clusterCols,
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    #column_names_gp = gpar(fontsize = fontSize), 
    
    # Row options:
    show_row_names = labelRows,
    cluster_rows = clusterRows,
    show_row_dend = showRowDendrogram,
    clustering_method_rows = "ward.D2",
    #row_names_gp = gpar(fontsize = fontSize), 
    
    # Raster info:
    use_raster = useRaster,
    raster_device = rasterDevice,
    raster_quality = rasterQuality,

    # Other
    ...
  )
  
  # Add row labels if provided:
  if(!is.null(customRowLabel)){
    if(is.null(customRowLabelIDs)){
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    hm <- hm + rowAnnotation(
      link = anno_mark(at = customRowLabel, labels = customRowLabelIDs, labels_gp = gpar(fontsize = fontSize)),
      width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs)
    )
  }
  
  return(hm)
}



################################################################################
# Volcano / MA plots
################################################################################

volcanoPlot <- function(df, cmap=NULL, cmap_style='qualitative', title=NULL, covarLabel="", 
  namedColors=FALSE, colorColName="color", minxmax=NULL, minymax=NULL, point_size=1){
  # Plot a volcano plot of differential genes
  # df is a 2+ column df:
  # col 1 = x axis (e.g. fold change)
  # col 2 = significance (e.g. adj log10 pval)
  # col 3 = labels (should points be labeled)
  # col 4 = point color
  # If camp is a named vector, will match values in column 4
  # min(x/y)max indicates the lowest allowable value of (x/y)max
  n <- nrow(df)
  na_col <- "grey88"

  # Convert data.table back to df if necessary (data.tables cause problems)
  df <- as.data.frame(df)

  # Get color col
  color_col <- match(colorColName, colnames(df))

  # Set colors
  if(ncol(df) < 4){
    df[,4] <- NA
  }
  nc <- length(unique(df[,color_col]))
  if(is.null(cmap)){
    cmap <- getColorMap(cmaps_BOR$stallion, n=nc)
  }
  
  # Plot a volcano plot
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2], group=df[,color_col]))
    + geom_point_rast(aes(color=df[,color_col]), size=point_size)
    #+ geom_point(aes(color=df[,color_col]), size=point_size)
    + geom_text_repel(aes(label=df[,3]), size=3, max.overlaps=Inf)
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + ggtitle(title)
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            plot.margin=unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio=1,
            #legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(expand=c(0, 0.75)) # Make bars start at the axis
  )
  if(namedColors){
    # Named colormap corresponding to discrete values in third column
    p <- p + scale_color_manual(values=cmap, limits=names(cmap), 
      name=covarLabel, na.value=na_col, drop=FALSE)
    p <- p + guides(fill = guide_legend(title=covarLabel), 
                    colour = guide_legend(override.aes=list(size=5)))
  }else{
    p <- p + scale_color_manual(values=colors, na.value=na_col)
  }
  # Enforce x and y lims if indicated
  if(!is.null(minxmax)){
    xrng <- layer_scales(p)$x$get_limits()
    xmin <- xrng[1]
    xmax <- xrng[2]
    xmin <- ifelse(abs(xmin) < minxmax, -minxmax, xmin)
    xmax <- ifelse(xmax < minxmax, minxmax, xmax)
    p <- p + xlim(xmin, xmax)
  }
  if(!is.null(minymax)){
    yrng <- layer_scales(p)$y$get_limits()
    ymin <- yrng[1]
    ymax <- yrng[2]
    ymax <- ifelse(ymax < minymax, minymax, ymax)
    suppressMessages(
      p <- p + scale_y_continuous(expand = c(0, 0), limits=c(0, ymax*1.05))
    ) 
  }
  p
}


MAPlot <- function(df, cmap=NULL, cmap_style='qualitative', title=NULL, covarLabel="", 
  namedColors=FALSE, colorColName="color", minxmax=NULL, minymax=NULL, point_size=1,
  set_xmin=NULL, set_xmax=NULL, set_ymin=NULL, set_ymax=NULL, border=FALSE){
  # Plot a MA plot of differential genes
  # df is a 2+ column df:
  # col 1 = x axis (e.g. base mean expression)
  # col 2 = y axis (e.g. fold change)
  # col 3 = labels (should points be labeled)
  # col 4 = point color
  # If camp is a named vector, will match values in column 5
  # min(x/y)max indicates the lowest allowable value of (x/y)max
  # sig_cutoff
  n <- nrow(df)
  na_col <- "grey88"

  # Convert data.table back to df if necessary (data.tables cause problems)
  df <- as.data.frame(df)

  # Get color col
  color_col <- match(colorColName, colnames(df))

  # Set colors
  if(ncol(df) < 4){
    df[,4] <- NA
  }
  nc <- length(unique(df[,color_col]))
  if(is.null(cmap)){
    cmap <- getColorMap(cmaps_BOR$stallion, n=nc)
  }
  
  # Plot a MA plot
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2], group=df[,color_col]))
    + geom_point_rast(aes(color=df[,color_col]), size=point_size)
    #+ geom_point(aes(color=df[,color_col]), size=point_size)
    + geom_text_repel(aes(label=df[,3]), 
      size=3, max.overlaps=Inf,
      box.padding=0.75, force=0.75
      )
    + geom_hline(yintercept=0.0, linetype="dashed")
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + ggtitle(title)
    + theme_BOR(border=border)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            plot.margin=unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio=1,
            #legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    #+ scale_y_continuous(expand=c(0, 0.75)) # Make bars start at the axis
  )
  if(namedColors){
    # Named colormap corresponding to discrete values in third column
    p <- p + scale_color_manual(values=cmap, limits=names(cmap), 
      name=covarLabel, na.value=na_col, drop=FALSE)
    p <- p + guides(fill = guide_legend(title=covarLabel), 
                    colour = guide_legend(override.aes=list(size=5)))
  }else{
    p <- p + scale_color_manual(values=colors, na.value=na_col)
  }
  # Get current x and y limits
  xrng <- layer_scales(p)$x$get_limits()
  xmin <- xrng[1]
  xmax <- xrng[2]
  yrng <- layer_scales(p)$y$get_limits()
  ymin <- yrng[1]
  ymax <- yrng[2]
  if(!is.null(set_xmin)){
    xmin <- set_xmin
  }
  if(!is.null(set_xmax)){
    xmax <- set_xmax
  }
  if(!is.null(set_ymin)){
    ymin <- set_ymin
  }
  if(!is.null(set_ymax)){
    ymax <- set_ymax
  }

  # Enforce x and y lims if indicated
  if(!is.null(minxmax)){
    xmin <- ifelse(xmin > -minxmax, -minxmax, xmin)
    xmax <- ifelse(xmax < minxmax, minxmax, xmax)
  }
  if(!is.null(minymax)){
    ymin <- ifelse(ymin > -minymax, -minymax, ymin)
    ymax <- ifelse(ymax < minymax, minymax, ymax)
  }
  # Reset x and y limits
  p <- p + xlim(xmin, xmax)
  suppressMessages(
    p <- p + scale_y_continuous(expand = c(0, 0), limits=c(ymin*1.05, ymax*1.05))
  ) 
  p
}

MitoCluster_Plot <- function(seurat_obj, mito_pattern = "^MT-") {
  # Ensure the necessary libraries are loaded
  library(Seurat)
  library(ggplot2)
  
  # Calculate percent mitochondrial
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mito_pattern)
  
  # Create a dataframe for plotting
  plot_data <- data.frame(
    Cluster = seurat_obj@meta.data$FineClust,
    PercentMito = seurat_obj@meta.data$percent.mt
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Cluster, y = PercentMito)) +
    geom_violin(trim = FALSE, fill = "lightblue", color = "blue") +
    geom_boxplot(width = 0.1, fill = "white", color = "darkblue", outlier.shape = NA) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cluster", y = "Percent Mitochondrial RNA", 
         title = "Mitochondrial RNA Content by Cluster")
  
  # Return the plot
  return(p)
}

plotUMAP_Labels <- function(df, dataType = "qualitative", cmap = NULL, covarLabel = "", point_size = 0.05, 
                                 namedColors = FALSE, plotTitle = NULL, colorLims = NULL, na.value = "grey35", useRaster = TRUE) {
  
  # Rename columns for clarity
  colnames(df)[1:3] <- c("UMAP1", "UMAP2", "Cluster")
  
  # Calculate centroids of clusters
  centroids <- df %>%
    group_by(Cluster) %>%
    summarise(UMAP1 = mean(UMAP1, na.rm = TRUE), UMAP2 = mean(UMAP2, na.rm = TRUE)) %>%
    ungroup()
  
  # Split centroids into left and right based on UMAP1 value
  median_UMAP1 <- median(centroids$UMAP1)
  left_centroids <- centroids %>% filter(UMAP1 < median_UMAP1) %>% arrange(UMAP2)
  right_centroids <- centroids %>% filter(UMAP1 >= median_UMAP1) %>% arrange(UMAP2)
  
  # Create label positions for left and right columns with increased distance
  label_positions <- data.frame(
    Cluster = c(left_centroids$Cluster, right_centroids$Cluster),
    UMAP1 = c(rep(min(df$UMAP1) - diff(range(df$UMAP1)) * 0.2, nrow(left_centroids)), 
              rep(max(df$UMAP1) + diff(range(df$UMAP1)) * 0.2, nrow(right_centroids))),
    UMAP2 = c(seq(min(df$UMAP2), max(df$UMAP2), length.out = nrow(left_centroids)),
              seq(min(df$UMAP2), max(df$UMAP2), length.out = nrow(right_centroids)))
  )
  
  # Adjust the position of the small colored dot
  dot_positions <- label_positions
  dot_positions$UMAP1 <- dot_positions$UMAP1 + ifelse(dot_positions$UMAP1 < median_UMAP1, 0.05, -0.05) * diff(range(df$UMAP1))
  
  # Merge centroids and dot positions to get correct lines
  segment_data <- merge(centroids, dot_positions, by = "Cluster", suffixes = c("_centroid", "_dot"))
  
  # Base plot with original aesthetics
  if (useRaster) {
    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
      geom_point_rast(size = point_size) +
      theme_minimal() +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 0.7,
        legend.position = "none"
      ) +
      xlab("UMAP1") +
      ylab("UMAP2")
  } else {
    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
      geom_point(size = point_size) +
      theme_minimal() +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 0.7,
        legend.position = "none"
      ) +
      xlab("UMAP1") +
      ylab("UMAP2")
  }
  
  # Add segments from centroids to the dots
  p <- p + geom_segment(data = segment_data, aes(x = UMAP1_centroid, y = UMAP2_centroid, xend = UMAP1_dot, yend = UMAP2_dot), 
                        color = "grey50")
  
  # Add text labels
  p <- p + geom_text(data = label_positions, aes(x = UMAP1, y = UMAP2, label = Cluster), 
                     colour = "black", size = 3, fontface = "bold", hjust = ifelse(label_positions$UMAP1 < median_UMAP1, 1, 0))
  
  # Plot a small dot with the correct color on the inside of each label
  p <- p + geom_point(data = dot_positions, aes(x = UMAP1, y = UMAP2, color = Cluster), 
                      size = 3, show.legend = FALSE)
  
  # Set plot title
  if (!is.null(plotTitle)) {
    p <- p + ggtitle(plotTitle)
  } else {
    p <- p + ggtitle(sprintf("n = %s", nrow(df)))
  }
  
  # If colormap provided, update colors
  if (!is.null(cmap)) {
    if (namedColors) {
      p <- p + scale_color_manual(values = cmap, limits = names(cmap), name = covarLabel, na.value = na.value)
    } else {
      names(cmap) <- NULL
      if (dataType == "qualitative") {
        nvals <- length(unique(df$Cluster))
        cmap <- getColorMap(cmap, n = nvals)  # Assuming getColorMap() is a function that generates a color map
        p <- p + scale_color_manual(values = cmap, name = covarLabel, na.value = na.value)
      } else {
        if (!is.null(colorLims)) {
          p <- p + scale_color_gradientn(colors = cmap, name = covarLabel, limits = colorLims, na.value = na.value)
        } else {
          p <- p + scale_color_gradientn(colors = cmap, name = covarLabel, na.value = na.value)
        }
      }
    }
  }

  # Increase the width of the plotting area
  x_range <- range(df$UMAP1)
  x_range_extended <- c(x_range[1] - diff(x_range) * 0.7, x_range[2] + diff(x_range) * 0.7)
  p <- p + coord_cartesian(xlim = x_range_extended)
  
  return(p)
}

plotUMAP_Labels_columns <- function(df, dataType = "qualitative", cmap = NULL, covarLabel = "", point_size = 0.05, 
                                 namedColors = FALSE, plotTitle = NULL, colorLims = NULL, na.value = "grey35", useRaster = TRUE) {
  
  # Rename columns for clarity
  colnames(df)[1:3] <- c("UMAP1", "UMAP2", "Cluster")
  
  # Calculate centroids of clusters
  centroids <- df %>%
    group_by(Cluster) %>%
    summarise(UMAP1 = mean(UMAP1, na.rm = TRUE), UMAP2 = mean(UMAP2, na.rm = TRUE)) %>%
    ungroup()
  
  # Split centroids into left and right based on UMAP1 value
  median_UMAP1 <- median(centroids$UMAP1)
  left_centroids <- centroids %>% filter(UMAP1 < median_UMAP1) %>% arrange(UMAP2)
  right_centroids <- centroids %>% filter(UMAP1 >= median_UMAP1) %>% arrange(UMAP2)
  
  # Create label positions for left and right columns with increased distance
  label_positions <- data.frame(
    Cluster = c(left_centroids$Cluster, right_centroids$Cluster),
    UMAP1 = c(rep(min(df$UMAP1) - diff(range(df$UMAP1)) * 0.2, nrow(left_centroids)), 
              rep(max(df$UMAP1) + diff(range(df$UMAP1)) * 0.2, nrow(right_centroids))),
    UMAP2 = c(seq(min(df$UMAP2), max(df$UMAP2), length.out = nrow(left_centroids)),
              seq(min(df$UMAP2), max(df$UMAP2), length.out = nrow(right_centroids)))
  )
  
  # Adjust the position of the small colored dot
  dot_positions <- label_positions
  dot_positions$UMAP1 <- dot_positions$UMAP1 + ifelse(dot_positions$UMAP1 < median_UMAP1, 0.05, -0.05) * diff(range(df$UMAP1))
  
  # Merge centroids and dot positions to get correct lines
  segment_data <- merge(centroids, dot_positions, by = "Cluster", suffixes = c("_centroid", "_dot"))
  
  # Base plot with original aesthetics
  if (useRaster) {
    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
      geom_point_rast(size = point_size) +
      theme_minimal() +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 0.7,
        legend.position = "none"
      ) +
      xlab("UMAP1") +
      ylab("UMAP2")
  } else {
    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
      geom_point(size = point_size) +
      theme_minimal() +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 0.7,
        legend.position = "none"
      ) +
      xlab("UMAP1") +
      ylab("UMAP2")
  }
  
  # Add segments from centroids to the dots
  p <- p + geom_segment(data = segment_data, aes(x = UMAP1_centroid, y = UMAP2_centroid, xend = UMAP1_dot, yend = UMAP2_dot), 
                        color = "grey50")
  
  # Add text labels
  p <- p + geom_text(data = label_positions, aes(x = UMAP1, y = UMAP2, label = Cluster), 
                     colour = "black", size = 3, fontface = "bold", hjust = ifelse(label_positions$UMAP1 < median_UMAP1, 1, 0))
  
  # Plot a small dot with the correct color on the inside of each label
  p <- p + geom_point(data = dot_positions, aes(x = UMAP1, y = UMAP2, color = Cluster), 
                      size = 3, show.legend = FALSE)
  
  # Set plot title
  if (!is.null(plotTitle)) {
    p <- p + ggtitle(plotTitle)
  } else {
    p <- p + ggtitle(sprintf("n = %s", nrow(df)))
  }
  
  # If colormap provided, update colors
  if (!is.null(cmap)) {
    if (namedColors) {
      p <- p + scale_color_manual(values = cmap, limits = names(cmap), name = covarLabel, na.value = na.value)
    } else {
      names(cmap) <- NULL
      if (dataType == "qualitative") {
        nvals <- length(unique(df$Cluster))
        cmap <- getColorMap(cmap, n = nvals)  # Assuming getColorMap() is a function that generates a color map
        p <- p + scale_color_manual(values = cmap, name = covarLabel, na.value = na.value)
      } else {
        if (!is.null(colorLims)) {
          p <- p + scale_color_gradientn(colors = cmap, name = covarLabel, limits = colorLims, na.value = na.value)
        } else {
          p <- p + scale_color_gradientn(colors = cmap, name = covarLabel, na.value = na.value)
        }
      }
    }
  }

  # Increase the width of the plotting area
  x_range <- range(df$UMAP1)
  x_range_extended <- c(x_range[1] - diff(x_range) * 0.7, x_range[2] + diff(x_range) * 0.7)
  p <- p + coord_cartesian(xlim = x_range_extended)
  
  return(p)
}


library(ggplot2)
library(dplyr)
library(patchwork)

create_combined_dotplot <- function(subgroups, output_directory, custom_palette, featureSetsList) {
  all_plots <- list()
  plot_info <- list()
  total_groups <- 0

  # Process each subgroup
  for (subgroup in subgroups) {
    file_wd <- file.path(output_directory, "Clustering_04", subgroup, paste0(subgroup, ".rds"))
    if (!file.exists(file_wd)) {
      message(sprintf("File for subgroup '%s' not found. Skipping.", subgroup))
      next
    }
    
    obj <- readRDS(file_wd)
    result <- process_subgroup(subgroup, obj, custom_palette, featureSetsList)
    if (!is.null(result)) {
      all_plots[[subgroup]] <- result$plot
      plot_info[[subgroup]] <- list(n_genes = result$n_genes, n_groups = result$n_groups)
      total_groups <- total_groups + result$n_groups
    }
  }

  # Combine all plots into a single figure
  if (length(all_plots) > 0) {
    # Calculate the relative widths based on the number of groups in each subplot
    relative_widths <- sapply(plot_info, function(x) x$n_groups)
    
    # Add legend to the rightmost plot
    subgroup_names <- names(all_plots)
    all_plots[[subgroup_names[length(subgroup_names)]]] <- all_plots[[subgroup_names[length(subgroup_names)]]] +
      theme(
        legend.position = "right",
        legend.text = element_text(face = "bold", size = 8, family = "Arial"),
        legend.title = element_text(face = "bold", size = 10, family = "Arial")
      )
    
    # Combine plots using patchwork with dynamic widths
    combined_plot <- wrap_plots(all_plots, ncol = length(all_plots), widths = relative_widths) +
      plot_layout(guides = 'collect') &
      theme(plot.margin = margin(5, 5, 5, 5))
    
    # Save the combined plot
    output_filename <- file.path(output_directory, "Combined_Dotplots.png")
    ggsave(filename = output_filename, plot = combined_plot, 
           dpi = 300, height = 10, width = 0.5 * total_groups, limitsize = FALSE)
    
    message("Combined plot saved successfully.")
    
    return(combined_plot)
  } else {
    message("No plots were generated. Check the individual subgroup processing for errors.")
    return(NULL)
  }
}

create_combined_dotplot <- function(subgroups, output_directory, custom_palette, featureSetsList) {
  all_plots <- list()
  plot_info <- list()
  total_groups <- 0

  # Process each subgroup
  for (subgroup in subgroups) {
    file_wd <- file.path(output_directory, "Clustering_04", subgroup, paste0(subgroup, ".rds"))
    if (!file.exists(file_wd)) {
      message(sprintf("File for subgroup '%s' not found. Skipping.", subgroup))
      next
    }
    
    obj <- readRDS(file_wd)
    result <- process_subgroup(subgroup, obj, custom_palette, featureSetsList)
    if (!is.null(result)) {
      all_plots[[subgroup]] <- result$plot
      plot_info[[subgroup]] <- list(n_genes = result$n_genes, n_groups = result$n_groups)
      total_groups <- total_groups + result$n_groups
    }
  }

  # Combine all plots into a single figure
  if (length(all_plots) > 0) {
    # Calculate the relative widths based on the number of groups in each subplot
    relative_widths <- sapply(plot_info, function(x) x$n_groups)
    
    # Add legend to the rightmost plot
    subgroup_names <- names(all_plots)
    all_plots[[subgroup_names[length(subgroup_names)]]] <- all_plots[[subgroup_names[length(subgroup_names)]]] +
      theme(
        legend.position = "right",
        legend.text = element_text(face = "bold", size = 8, family = "Arial"),
        legend.title = element_text(face = "bold", size = 10, family = "Arial")
      )
    
    # Combine plots using patchwork with dynamic widths
    combined_plot <- wrap_plots(all_plots, ncol = length(all_plots), widths = relative_widths) +
      plot_layout(guides = 'collect') &
      theme(plot.margin = margin(5, 5, 5, 5))
    
    # Save the combined plot
    output_filename <- file.path(output_directory, "Combined_Dotplots.png")
    ggsave(filename = output_filename, plot = combined_plot, 
           dpi = 300, height = 10, width = 0.5 * total_groups, limitsize = FALSE)
    
    message("Combined plot saved successfully.")
    
    return(combined_plot)
  } else {
    message("No plots were generated. Check the individual subgroup processing for errors.")
    return(NULL)
  }
}

process_subgroup <- function(subgroup, obj, custom_palette, featureSetsList) {
  message(sprintf("Processing subgroup: %s", subgroup))
  
  tryCatch({
    count_mat <- GetAssayData(object = obj, slot = "counts")
    
    if (!is.numeric(count_mat)) {
      count_mat <- as.matrix(count_mat)
      mode(count_mat) <- "numeric"
      message("Converted count_mat to numeric matrix.")
    }
    
    if (!"cell_type" %in% colnames(obj@meta.data)) {
      stop("'cell_type' column not found in object metadata.")
    }
    
    avgPctMat <- avgAndPctExpressed(count_mat, obj$cell_type, feature_normalize = TRUE, min_pct = 5)
    
    required_cols <- c("feature", "grp", "avgExpr", "pctExpr")
    if (!all(required_cols %in% colnames(avgPctMat))) {
      stop("avgPctMat is missing one or more required columns: ", 
           paste(setdiff(required_cols, colnames(avgPctMat)), collapse = ", "))
    }
    
    list_index <- get_list_index(subgroup)
    subGenes <- featureSetsList[[list_index]] %>% unlist()
    
    avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]
    
    if (nrow(avgPctMat) == 0) {
      message("No matching genes found for this subgroup. Skipping plot creation.")
      return(NULL)
    }
    
    wide_df <- unmelt_alt(avgPctMat, row_col = "feature", col_col = "grp", val_col = "avgExpr")
    wide_df <- prettyOrderMat(wide_df)
    grp_order <- colnames(wide_df$mat)
    gene_order <- rev(rownames(wide_df$mat))
    
    p <- ggplot(avgPctMat, aes(x = grp, y = feature)) +
      geom_point(aes(color = avgExpr, size = pctExpr)) +
      scale_color_gradientn(colors = custom_palette, name = "Avg Expression") +
      scale_size_continuous(range = c(0, 6), limits = c(0, 100), breaks = c(25, 50, 75, 100), name = "% Expressed") +
      scale_y_discrete(limits = gene_order, expand = expansion(add = 0.5)) +
      scale_x_discrete(limits = grp_order, expand = expansion(add = 0.5)) +
      labs(title = subgroup) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text(face = "bold", size = 8, angle = 90, hjust = 1, vjust = 0.5, family = "Arial"),
        axis.text.y = element_text(face = "bold", size = 8, family = "Arial"),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(20, 5, 5, 5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      coord_cartesian(clip = "off")
    
    return(list(plot = p, n_genes = length(gene_order), n_groups = length(grp_order)))
  }, error = function(e) {
    message(sprintf("Error processing subgroup '%s': %s", subgroup, e$message))
    return(NULL)
  })
}
