#' Inserts 0 columns and rows after last row/column to symmetrize matrix.
#'
#' @param mtx matrix in dense format to be symmetrized
#' @param N positive integer; additional argument for symmetrizing matrix to desired N x N dimension; N need not be larger than \code{ncol(mtx)} or \code{nrow(mtx)} in which case submatrix \code{mtx[1:N,1:N]} will be extraced
#'
#' @return N by N matrix which is either submatrix of \code{mtx} or \code{mtx} extended with 0's row and/or columns
#'
#' @examples
#' mtx1 <- matrix(1:24, ncol = 4)
#' mtx2 <- matrix(1:24, nrow = 4)
#' print(mtx1)
#' print(mtx2)
#' balance(mtx1)
#' balance(mtx2)
#' balance(mtx1, N = 8)
#' balance(mtx1, N = 3)
#'
#' @export
balance <- function(mtx, N = NULL){
  if(!is.null(N)){
    stopifnot(N > 0)
  }
  mtx.n.rows <- nrow(mtx)
  mtx.n.cols <- ncol(mtx)
  # balance if nrows != ncols
  if(mtx.n.rows < mtx.n.cols){
    # add 0-s rows at the end
    rep(0, mtx.n.cols) %>%
      replicate(mtx.n.cols - mtx.n.rows,.) %>%
      t() %>% rbind(mtx, .) -> mtx
  } else if(mtx.n.rows > mtx.n.cols){
    # add 0-s columns at the end
    rep(0, mtx.n.rows) %>%
      replicate(mtx.n.rows - mtx.n.cols,.) %>%
      cbind(mtx, .) -> mtx
  }
  # add/remove columns and rows, so dim(mtx) == c(N,N)
  if(!is.null(N)){
    if(N > nrow(mtx)){
      rep(0, nrow(mtx)) %>%
        replicate(N - nrow(mtx), .) %>%
        cbind(mtx, .) -> mtx
      rep(0, ncol(mtx)) %>%
        replicate(N - nrow(mtx), .) %>%
        t() %>% rbind(mtx, .) -> mtx
    } else if(N < nrow(mtx)){
      mtx <- mtx[1:N,1:N]
    }
  }
  return(mtx)
}

#' Converts matrix given in dense format to sparse format data frame.
#'
#' This function only keeps non-zero cells. In case given dense matix is symmetric \code{dense2sparse} will return upper triangular part of the matrix (i.e. where rows <= columns)
#'
#' @param mtx matrix in dense format
#' @param add.diagonal logical, if true an additional column indicating diagonal of each cell will be appended to resulting data frame
#' @param name character, additional argument, if specified column with name will be appended to resulting data frame
#'
#' @return data.frame with columns \code{c("i","j","val")} and optionally \code{c("diagonal","name")} columns; every row of resulting dataframe corresponds to cell in given dense matrix with i-th row, j-th column and value val
#'
#' @examples
#' dense2sparse(matrix(1:24, ncol = 3))
#' dense2sparse(matrix(1:24, ncol = 3), name = "some.matrix")
#' dense2sparse(matrix(1:24, ncol = 3), add.diagonal = FALSE)
#' # symmetric matrix
#' mtx.sym <- matrix(1:25, ncol = 5)
#' mtx.sym <- mtx.sym + t(mtx.sym)
#' dense2sparse(matrix(mtx.sym))
#'
#' @export
dense2sparse <- function(mtx, add.diagonal = TRUE, name = NULL){
  Matrix::Matrix(mtx, sparse = TRUE) %>%
    Matrix::summary() %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("i","j","val")) -> m
  if(add.diagonal){
    m$diagonal <- abs(m$i - m$j)
  }
  if(is.null(name)){
    return(m)
  } else {
    return(cbind(m, name = name))
  }
}

#' Converts sparse matrix to dense matrix.
#'
#' This function handles can distinguish between symmetric and non symmetric matrices when given sparse matrix can be interpreted as square.
#'
#' @param sparse.as.df data.frame with sparse matrix containing 3 columns: i, j, val
#' @param N positive integer, optional desired dimension of dense matrix enforced with \code{balancing} function
#' @param balancing logical, if perform balancing; see \code{balancing} description for information on function behaviour when N is NULL
#' @param missing.cells.na logical, if TRUE then substitute cells with i,j coordinates missing in given sparse matrix with NA, otherwise substitute with 0
#'
#' @seealso \code{\link{balancing}} for balancing
#'
#' @examples
#' # produce dense matrix
#' mtx.dense <- matrix(1:24, ncol = 3)
#' # add some zeros
#' mtx.dense <- rbind(matrix.dense, c(0,100,0))
#' # construct sparse matrix
#' mtx.sparse <- dense2sparse(mtx.dense, add.diagonal = FALSE)
#' # get back dense matrix
#' sparse2dense(mtx.sparse, balancing = FALSE)
#' sparse2dense(mtx.sparse, balancing = FALSE, missing.cells.na = TRUE)
#' # symmetric matrices
#' mtx.sym <- matrix(1:16, ncol = 4)
#' mtx.sym <- mtx.sym + t(mtx.sym)
#' sparse2dense(dense2sparse(mtx.sym, add.diagonal = FALSE))
#' # square non symmetric matrices
#' mtx.sq <- matrix(1:16, ncol = 4)
#' mtx.sq[c(2,3,4,7,12)] <- 0
#' sparse2dense(mtx.sparse)
#' sparse2dense(mtx.sparse, missing.cells.na = TRUE)
#'
#' @export
sparse2dense <- function(sparse.as.df, N = NULL, balancing = TRUE, missing.cells.na = FALSE){
  colnames(sparse.as.df) <- c("i","j","val")
  with(
    sparse.as.df,
    Matrix::sparseMatrix(i=as.numeric(i),
                         j=as.numeric(j),
                         x=as.numeric(val))
  ) %>%
    as.matrix() -> dense.mtx
  if((!is.null(N)) | balancing){
    dense.mtx <- balance(dense.mtx, N = N)
  }
  shape <- dim(dense.mtx)
  if(
    (shape[1] == shape[2]) &
    (all(sparse.as.df$i >= sparse.as.df$j) | all(sparse.as.df$i <= sparse.as.df$j))
  ){
    # if matrix is square symmetric then symmetrize it
    dense.mtx.diag <- diag(dense.mtx)
    dense.mtx <- dense.mtx + t(dense.mtx)
    diag(dense.mtx) <- dense.mtx.diag
  }
  # lastly change 0 to NA if missing.cells.na
  if(missing.cells.na & any(dense.mtx == 0)){
    # convert missing cells to NA instead of 0
    if.missing <- matrix(0, nrow = nrow(dense.mtx), ncol = ncol(dense.mtx))
    if.missing <- if.missing == 0
    # set non missing cells to FALSE
    if.missing[cbind(sparse.as.df$i, sparse.as.df$j)] <- FALSE
    dense.mtx[if.missing] <- NA
  }
  return(dense.mtx)
}

#' Efficiently slices k-th diagonal from matrix A.
#'
#' @param A matrix
#' @param k integer diagonal to be sliced, 1 is main diagonal, positive k will yield diagonals from lower triangular matrix while negative k will yield diagonals from upper triangular
#'
#' @return numeric, vector containing entries on k-th diagonal of matrix A
#'
#' @examples
#' mtx <- matrix(1:25, ncol = 5)
#' superdiag(mtx, k = 1)
#' superdiag(mtx, k = 2)
#' superdiag(mtx, k = -2)
#'
#' @export
superdiag <- function(A, k){
  stopifnot(k != 0)
  if(k < 0){
    A <- t(A)
    k <- abs(k)
  }
  n <- nrow(A)
  len <- n - k + 1
  i <- 1:len
  j <- k:n
  indices <- (j - 1) * n + i
  return(A[indices])
}

#' Removes unmappable regions (all zeros columns and rows) from dense matrix.
#'
#' @param dense.mtx numeric matrix (contact map) in dense format
#'
#' @return list containing 3 elements: indices of removed rows, indices of removed columns and matrix without unmappable regions
#'
#' @examples
#' # construct matrix
#' mtx <- matrix(1:32, ncol = 4)
#' mtx[c(5,7,8),] <- 0
#' mtx[,c(2,4)] <- 0
#' l <- remove_unmappable(mtx)
#' print(l[["indices.rows"]])
#' print(l[["indices.cols"]])
#' print(l[["matrix"]])
#'
#' @export
remove_unmappable <- function(dense.mtx){
  idx.rows <- which(rowSums(dense.mtx) == 0)
  idx.cols <- which(colSums(dense.mtx) == 0)
  list(idx.rows, idx.cols, dense.mtx[-idx.rows,-idx.cols]) %>%
    magrittr::set_names(c("indices.rows", "indices.cols","matrix"))
}

#' Restores deleted regions (cells) in vector, for example unmappable regions in pc vector.
#'
#' @param vec numeric vector
#' @param idx indices of elements in initial vector (before cells removal) - the one which is to be restored
#' @param empty.elem numeric or NA how to fill missing (restored) cells
#'
#' @examples
#' # create vector with zeros
#' v <- c(1,2,3,0,0,0,2,2,0,9,8,0)
#' # get indices of 0 elements
#' idx <- which(v == 0)
#' v.without <- v[-idx]
#' restore_unmappable_vec(v.without, idx)
#' restore_unmappable_vec(v.without, idx, empty.elem = 0)
#'
#' @export
restore_unmappable_vec <- function(vec, idx, empty.elem = NA){
  v <- rep(empty.elem, length(vec) + length(idx))
  v[-idx] <- vec
  return(v)
}

#' Restores unmappable regions (all zeros columns and rows) in dense matrix.
#'
#' @param dense.mtx.mappable numeric matrix
#' @param idx.rows integer vector, row indices of 0's elements in initial matrix (before rows removal) - the one which is to be restored
#' @param idx.cols integer vector, column indices of 0's elements in initial matrix (before columns removal) - the one which is to be restored, by default it is equal to idx.rows
#' @param empty.elem numeric or NA how to fill missing (restored) cells
#'
#' @return numeric matrix in dense format
#'
#' @examples
#' # construct matrix with 0's row 5,7,8 and 0's columns 2,4
#' m <- matrix(1:32, ncol = 4); m[c(5,7,8),] <- 0; m[,c(2,4)] <- 0
#' l <- remove_unmappable(mtx)
#' restore_unmappable_mtx(l[["matrix"]], l[["indices.rows"]], idx.cols = l[["indices.cols"]])
#' restore_unmappable_mtx(l[["matrix"]], l[["indices.rows"]], idx.cols = l[["indices.cols"]], empty.elem = 0)
#'
#' @export
restore_unmappable_mtx <- function(dense.mtx.mappable, idx.rows, idx.cols = NULL, empty.elem = NA){
  if(is.null(idx.cols)){
    idx.cols <- idx.rows
  }
  dense.mtx <- matrix(empty.elem, nrow = nrow(dense.mtx.mappable) + length(idx.rows),
                      ncol = ncol(dense.mtx.mappable) + length(idx.cols))
  dense.mtx[-idx.rows,-idx.cols] <- dense.mtx.mappable
  if(!is.null(colnames(dense.mtx.mappable))){
    colnames(dense.mtx) <- restore_unmappable_vec(
      colnames(dense.mtx.mappable), idx.cols, empty.elem = "unmappable"
    )
  }
  if(!is.null(rownames(dense.mtx.mappable))){
    rownames(dense.mtx) <- restore_unmappable_vec(
      rownames(dense.mtx.mappable), idx.rows, empty.elem = "unmappable"
    )
  }
  return(dense.mtx)
}

#' Reads given npz file and returns list with dense matrices.
#'
#' This function requires python and numpy.
#'
#' @param path character, string specifying path to npz file
#' @param mtx.names character vector specyfing subset of matrices names to read from npz dict like file; by default all matrices are loaded
#'
#' @seealso \url{https://www.python.org/} for python, \url{https://www.numpy.org/} for numpy
#'
#' @examples
#' # get sample npz file name
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
#' mtx.dense.list <- read_dense(mtx.fname) # reads all chromosomes
#' # limiting number of chromosomes
#' mtx.dense.sublist <- read_dense(mtx.fname, mtx.names = c("18","19")) # only read chromosome 18 and 19
#'
#' @export
read_dense <- function(path, mtx.names = "all"){
  np <- reticulate::import("numpy")
  npz <- np$load(path)
  if(all(mtx.names == "all")){
    mtx.names <- magrittr::use_series(npz,"files")
  }
  mtx.names %>%
    lapply(function(name){
      npz %>%
        magrittr::use_series("f") %>%
        magrittr::extract2(name)
    }) %>%
    magrittr::set_names(mtx.names)
}

#' Reads npz file with matrices and converts them sparse matrices.
#'
#' @param path character, string specifying path to npz file
#' @param mtx.names character vector specyfing subset of matrices names to read from npz dict like file; by default all matrices are loaded
#' @param sparse.format logical if FALSE then this function is equivalent to \code{\link{read_dense}} function
#'
#' @return list with matrices in sparse format
#'
#' @seealso \code{\link{dense2sparse}} for conversion of dense matrix to sparse format, \code{\link{read_dense}} on reading npz files
#'
#' @examples
#' # get sample npz file name
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
#' mtx.sparse.list <- read_npz(mtx.fname) # reads all chromosomes
#' mtx.sparse.sublist <- read_npz(mtx.fname, mtx.names = c("18","19")) # only read chromosome 18 and 19
#'
#' @export
read_npz <- function(path, mtx.names = "all", sparse.format = TRUE){
  np <- reticulate::import("numpy")
  npz <- np$load(path)
  if(all(mtx.names == "all")){
    mtx.names <- magrittr::use_series(npz,"files")
  }
  mtx.names %>%
    lapply(function(name){
      npz %>%
        magrittr::use_series("f") %>%
        magrittr::extract2(name) -> dense.mtx
      if(sparse.format){
        return(dense2sparse(dense.mtx, name = name))
      } else {
        return(dense.mtx)
      }
    }) %>%
    magrittr::set_names(mtx.names)
}

#' Saves list with dense matrices to npz compressed file.
#'
#' @param mtx.list list containg dense matrices, list should have names
#' @param path character string specifying path together with filename to save matrices
#'
#' @return None
#'
#' @seealso \url{https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.savez_compressed.html} for python numpy method used to save matrices in npz compressed file
#'
#' @examples
#' # get sample npz file name
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
#' mtx.dense.sublist <- read_dense(mtx.fname, mtx.names = c("18","19"))
#' print(names(mtx.dense.sublist))
#' print(typeof(mtx.dense.sublist))
#' print(typeof(mtx.dense.sublist[["18"]]))
#' print(dim(mtx.dense.sublist[["18"]]))
#' # save
#' save_npz(mtx.dense.sublist, "example.npz")
#'
#' @export
save_npz <- function(mtx.list, path){
  # check if list contains dense matrices
  stopifnot(all(sapply(mtx.list, is.matrix)))
  np <- reticulate::import("numpy")
  c(path, mtx.list) %>%
    do.call(np$savez_compressed, .)
}

#' Reads dimension of every matrix in npz file.
#'
#' @param path character, string specifying path to npz file
#' @param mtx.names character vector specyfing subset of matrices names to read from npz dict like file; by default all matrices are loaded
#'
#' @return matrix where row names are names of matrices from npz file and columns are rows and cols; each cell of the matrix contains number of rows and columns respectively that contact map with given name consits of
#'
#' @examples
#' # get sample npz file name
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
#' sizes <- read_size(mtx.fname)
#' print(sizes)
#'
#' @export
read_size <- function(path, mtx.names = "all"){
  np <- reticulate::import("numpy")
  npz <- np$load(path)
  if(all(mtx.names == "all")){
    mtx.names <- magrittr::use_series(npz,"files")
  }
  mtx.names %>%
    sapply(function(name){
      npz %>%
        magrittr::use_series("f") %>%
        magrittr::extract2(name) %>%
        dim() %>%
        as.numeric()
    }) %>% t() %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("n.rows","n.cols"))
}

#' Calculates ranges of each consecutive compartment along given pc vector.
#'
#' Entries with the same sign (i.e. positive or negative) comprise the same compartment. Positives are assigned to A compartment and negatives to B compartment.
#'
#' @param pc numeric, compartment vector (eigenvector)
#'
#' @return data.frame where each row corrspond to interval of consecutive same sign values of eigenvector; columns are start, end and compartment
#'
#' @examples
#' # make artificial eigenvector
#' ev <- c(-0.3,-0.5,0.2,0.3,0.4,0.4,-0.5,0.2,0.1,0.3,-0.9,-0.7)
#' compartment_ranges(ev)
#'
#' @export
compartment_ranges <- function(pc){
  pc.char <- ifelse(c1 > 0, "A", "B")
  lapply(c("A","B"), function(x){
    v <- which(pc.char == x)
    split(v, cumsum(c(1, diff(v) != 1))) %>%
      sapply(function(y){
        c(min(y),max(y))
      }) %>% t() %>%
      as.data.frame() %>%
      set_colnames(c("start","end")) %>%
      inset("compartment", value = x)
  }) %>%
    do.call("rbind",.) -> df
  df[order(df$start, df$end), ] %>%
    set_rownames(NULL)
}

#' Matches compartment with entries of contact map given in sparse format.
#'
#' Compartment vector indices correspond to bins along chromosome.
#'
#' @param pc numeric, compartment vector or principal component (eigenvector)
#' @param sparse.mtx data.frame, matrix in sparse format containg manadatory fields i and j
#'
#' @return data.frame - contact map in sparse format with additional columns for every row (dense matrix cell): compartment.i, compartment.j, compartment
#'
#' @examples
#' # make artificial eigenvector --> its length is 12 so it represents artificial chromosome with 12 bins
#' pc <- c(-0.3,-0.5,0.2,0.3,0.4,0.4,-0.5,0.2,0.1,0.3,-0.9,-0.7)
#' # make artificial sparse contact map
#' sparse.mtx <- data.frame(i = c(1,3,5,7,8,9,11), j = c(7,2,1,1,3,10,9), val = c(10,8,3,1,1,2,20))
#' print(sparse.mtx)
#' pc2mtx(pc, sparse.mtx)
#'
#' @export
pc2mtx <- function(pc, sparse.mtx){
  stopifnot(length(unique(sparse.mtx$name)) == 1)
  compartments <- data.frame(compartment = ifelse(pc > 0, "A", "B"), stringsAsFactors = FALSE)
  base::merge(sparse.mtx, compartments, by.x = "i", by.y ="row.names") %>%
    base::merge(compartments, by.x = "j", by.y ="row.names", suffixes = c(".i",".j")) %>%
    magrittr::inset(
      "compartment",
      value = paste0(
        magrittr::use_series(.,"compartment.i"),
        magrittr::use_series(.,"compartment.j")
      )
    ) -> df
  df[df$compartment == "BA", "compartment"] <- "AB"
  return(df[c(2,1,3:ncol(df))])
}

#' Constructs Toeplitz matrix with means on diagonals (calculated as arithmetic mean of diagonal).
#'
#' @param dense.mtx matrix in dense format
#'
#' @return numerical dense matrix where all entries on i-th diagonal contain mean of i-th diagonal
#'
#' @seealso \code{\link{toeplitz}} for more information about toeplitz matrices
#'
#' @examples
#' mtx1 <- toeplitz(c(1,2,3,4))
#' mean_expected(mtx1)
#' mtx2 <- matrix(1:16, ncol = 4)
#' mean_expected(mtx2)
#'
#' @export
mean_expected <- function(dense.mtx){
  shape <- dim(dense.mtx)
  stopifnot(shape[1] == shape[2])
  seq(nrow(dense.mtx)) %>%
    sapply(function(x){
      superdiag(dense.mtx, x) %>%
        mean()
    }) %>%
    toeplitz() %>%
    magrittr::divide_by(dense.mtx, .)
}

#' PCA analysis of Hi-C contact maps.
#'
#' Performs PCA on Hi-C contact map as described in Liebermann-Aiden et al 2009. More specifically it runs following routines on dense matrix:
#' \itemize{
#'  \item{}{removes unmappable regions (all zeros rows and columns)}
#'  \item{}{divides each diagonal of every cell by its corresponding mean of cells on diagonal}
#'  \item{}{converts matrix from 2. into PCC matrix}
#'  \item{}{performs PCA on such matrix}
#'  \item{}{fills in unmappable regions into PCA object vector/matrix components}
#' }
#'
#' @param dense.mtx numeric matrix - Hi-C contact map
#' @param ... optional arguments passed to \code{prcomp}
#'
#' @return PCA object returned by \code{prcomp} function.
#'
#' @seealso \code{\link{prcomp}} for how is PCA performed, Lieberman-Aiden E. et al., 2009 "Comprehensive mapping of long-range interactions reveals folding principles of the human genome." for compartment detection in Hi-C contact maps.
#'
#' @examples
#' # load Hi-C contact maps from npz file
#' # get sample npz file name
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
#' mtx.sparse.list <- read_npz(mtx.fname, sparse.format = TRUE)
#' # get matrix for selected chromosome
#' mtx <- mtx.sparse.list[["18"]]
#' # do PCA
#' pca <- do_pca(mtx)
#' print(pca)
#' # it is also possible to visualize results
#' pairs(pca)
#'
#' @export
do_pca <- function(dense.mtx, ...){
  l <- remove_unmappable(dense.mtx)
  l[["matrix"]] %>%
    mean_expected() %>%
    magrittr::inset(is.na(.), value = 0) %>%
    cor() %>%
    prcomp(...) -> pca
  pca$sdev <- restore_unmappable_vec(pca$sdev, l[["indices.rows"]])
  pca$center <- restore_unmappable_vec(pca$center, l[["indices.rows"]])
  pca$rotation <- restore_unmappable_mtx(pca$rotation, l[["indices.rows"]])
  pca$x <- restore_unmappable_mtx(pca$x, l[["indices.rows"]])
  return(pca)
}

#' Convert sparse contact map to atomic interactions.
#'
#' Converts contact map in sparse data frame format into data frame with interaction, where one row is single interaction. This is usefull for bootstraping Hi-C interactions.
#'
#' @param sparse.mtx data.frame Hi-C contact map in sparse format with mandatory columns i, j, val
#'
#' @return data.frame with 2 columns, which are i-th (x) and j-th (y) coordinates of every interaction.
#'
#' @seealso \code{\link{dense2sparse}}
#'
#' @examples
#' # create data frame with artificial interactions, where val
#' # indicates total number of interactions between bins i and j
#' sparse.mtx <- data.frame(i = c(1,3,5,7,8,9,11), j = c(7,2,1,1,3,10,9), val = c(10,8,3,1,1,2,20))
#' print(sparse.mtx)
#' sparse2interactions(sparse.mtx)
#'
#' @export
sparse2interactions <- function(sparse.mtx){
  sparse.mtx %>%
    magrittr::extract(
      rep(rownames(.), magrittr::use_series(., "val")),
      -match("val", colnames(.))
    ) %>%
    magrittr::set_rownames(NULL)
}

#' Hi-C interactions bootstrapping.
#'
#' Randomly samples interactions from data frame with atomic interactions into \code{length(ratio)} data frames with interactions in such a way that i-th dataframe have \code{ratio[i]} fraction of interactions of initial atomic interactions data frame.
#'
#' @param interactions data frame containing row with i and j coordinate for every single interaction
#' @param ratio numeric vector indicating on how many atomic interaction sets should initial atomic interaction set be divided; each entry of ratio vector contains fraction of interaction to be put in corresponding atomic interactions set; ratio vector must sum to 1 and all its entries must be larger than one
#'
#' @return data frame representing sparse Hi-C maps (with i, j, val columns) belonging to corresponding bootstrapped interactions subset - ratio.number column indicates index of fraction from ratio vector
#'
#' @examples
#' # create data frame with artificial interactions, where val
#' # indicates total number of interactions between bins i and j
#' sparse.mtx <- data.frame(i = c(1,3,5,7,8,9,11), j = c(7,2,1,1,3,10,9), val = c(10,8,3,1,1,2,20))
#' atomic.interactions <- sparse2interactions(sparse.mtx)
#' print(head(atomic.interactions))
#' b1 <- bootstrap_interactions(atomic.interactions)
#' print(b1)
#' ratios.desired <- c(0.4,0.3,0.2,0.1)
#' b2 <- bootstrap_interactions(atomic.interactions, ratio = ratios.desired)
#' print(b2)
#' sum.interactions <- nrow(atomic.interactions)
#' ratios.sampled <- sapply(split(b2, b2$ratio.number), function(x){ sum(x$val) / sum.interactions })
#' print(ratios.desired)
#' print(ratios.sampled)
#'
#' @export
bootstrap_interactions <- function(interactions, ratio = c(0.5, 0.5)){
  stopifnot(sum(ratio) == 1)
  stopifnot(all(ratio > 0))
  n <- nrow(interactions)
  floor(ratio[seq(1,length(ratio)-1)] * n) %>%
    c(magrittr::subtract(n,sum(.))) -> groups.sizes
  seq(1,length(groups.sizes)) %>%
    sapply(function(x){
      rep(x, groups.sizes[x])
    }, simplify = FALSE) %>%
    unlist() %>%
    sample() %>%
    split(interactions,.) -> interactions.ratios
  names(interactions.ratios) %>%
    lapply(function(x){
      # convert atomic interactions to sparse matrix
      aggregate(N ~ ., data = transform(interactions.ratios[[x]], N = 1), FUN = length) %>%
        magrittr::inset(, "ratio.number", value = x) -> bootstrapped
      idx <- match("N", colnames(bootstrapped))
      colnames(bootstrapped)[idx] <- "val"
      return(bootstrapped)
    }) %>%
    do.call("rbind",.)
}

#' Bootstraps interactions from contact map given in sparse format.
#'
#' For details of bootstrapping procedure see \code{\link{bootstrap_interactions}}.
#'
#' @param sparse.mtx data.frame Hi-C contact map in sparse format with mandatory columns i, j, val
#' @param ratio numeric vector indicating on how many atomic interaction sets should initial atomic interaction set be divided; each entry of ratio vector contains fraction of interaction to be put in corresponding atomic interactions set; ratio vector must sum to 1 and all its entries must be larger than one
#'
#' @return list with data frames (Hi-C maps in sparse format) containing sampled interactions (according to specified ratio vector)
#'
#' @seealso \code{\link{bootstrap_interactions}} for details of Hi-C interactions bootstrapping procedure
#'
#' @examples
#' sparse.mtx <- data.frame(i = c(1,3,5,7,8,9,11), j = c(7,2,1,1,3,10,9), val = c(10,8,3,1,1,2,20))
#' bootstrapped <- bootstrap_sparse(sparse.mtx)
#' print(bootstrapped)
#'
#' @export
bootstrap_sparse <- function(sparse.mtx, ratio = c(0.5, 0.5)){
  # divides Hi-C contact map on length(ration) maps according to specified ratios,
  # where each map contains fraction of contacts from original contact map
  # initial contact map is given in sparse data.frame format:
  # i j val chromosome
  stopifnot(sum(ratio) == 1)
  stopifnot(all(ratio > 0))
  mtx.by.name <- split(sparse.mtx, sparse.mtx$name)
  names(mtx.by.name) %>%
    lapply(function(name){
      sparse2interactions(mtx.by.name[[name]]) %>%
        bootstrap_interactions(ratio = ratio)
    }) %>%
    do.call("rbind",.) %>%
    split(magrittr::use_series(.,"ratio.number"))
}

#' Maps interactions to TADs.
#'
#' Maps cells of a contact map given in sparse format to TADs. User must provide interactions data frame and TADs for single and the same chromosome, otherwise the function will throw an error or behaviour will undefined.
#'
#' @param mtx.sparse data.frame with Hi-C contact map in sparse format; must have i and j columns, i.e. cell coordinates
#' @param tads data frame containg TADs, for single chromosome; it is assumed that first TAD on every chromosome start at 0 and end of TAD equals start of next TAD (in case if there is no break between consecutive TADs); also user is responsible for converting TAD boundary coordinates to bins
#' @param cols character vector with columns from mtx.sparse or tads data.frames to be included in resulting data frame; by default only val column is included
#'
#' @return data.frame with following columns: i, j, val, tad.id, start, end, name (and additional if cols specified), where each row corresponds to cell in Hi-C contact map with i, j coordinates and cells are assigned to TADs (cells which do not belong to any TAD are discarded)
#'
#' @examples
#' # create artificial interactions set
#' sparse.mtx <- data.frame(i = c(1,3,4,4,8,9,11), j = c(7,2,4,5,3,10,9), val = c(10,8,3,1,1,2,20), compartment = c("AA","AA","AB","AB","BB","BB","AB"))
#' # create artificial TAD set in 20000 bp resolution
#' resolution <- 20000
#' tads <- data.frame(start = c(0,2,10) * resolution, end = c(2,8,13) * resolution, name = as.character(c(1,1,1)))
#' print(tads)
#' # convert basepairs to bins
#' tads$start <- tads$start / resolution
#' tads$end <- tads$end / resolution
#' print(tads)
#' # map interactions to TADs
#' interactions2tads(sparse.mtx, tads)
#' # map interactions to TADs and keep also compartment columns
#' interactions2tads(sparse.mtx, tads, cols = c("val","compartment"))
#'
#' @export
interactions2tads <- function(mtx.sparse, tads, cols = c("val")){
  # map cells to TADs
  tads.parsed <- tads
  tads.parsed$start <- tads.parsed$start + 1
  tads.parsed$tad.id <- as.numeric(rownames(tads.parsed))
  boundaries <- c(tads.parsed$start, tads.parsed$end) %>%
    unique() %>% sort()
  v1 <- findInterval(mtx.sparse$i, boundaries)
  v2 <- findInterval(mtx.sparse$j, boundaries)
  cell.tad <- mtx.sparse[v1 == v2,]
  cell.tad$start <- sapply(v1[v1 == v2], function(x){
    boundaries[x]
  })
  merged <- merge(cell.tad, tads.parsed, by = c("start")) %>%
    extract(c(c("i","j","tad.id","start","end"), cols))
  stopifnot(all((merged$i >= merged$start) & (merged$i <= merged$end)))
  stopifnot(all((merged$j >= merged$start) & (merged$j <= merged$end)))
  return(merged)
}

#' Fits bilinear model to set of x,y points.
#'
#' @param x.vec numeric vector of x coordinates
#' @param y.vec numeric vector of y coordinates, must be the same length as x
#' @param truncate.left positive integer - number of points to exclude from left hand side
#' @param truncate.right positive integer - number of points to exclude from right hand side
#'
#' @return list with two componenets: numeric 2 by 2 matrix of coefficients, where row indicate model (left or right) and columns are intercept and slope; numeric vector intersection.x with: x coordinate of first point in left model closest to y (from right hand side), x coordinate of intersection point between left and right models and x coordinate of first point in right model closest to y (from left hand side); these points may be used as cutoff
#'
#' @seealso The code was taken from \url{https://stackoverflow.com/questions/15874214/piecewise-function-fitting-with-nls-in-r} (with minor modifications).
#'
#' @examples
#' # fix parameters of left linear model
#' a.left <- 0
#' b.left <- 0.8
#' # fix parameters of right linear model
#' a.right <- 15
#' b.right <- 0.1
#' # make models
#' x.left <- 1:20
#' y.left <- a.left + b.left * x.left
#' x.right <- 25:45
#' y.right <- a.right + b.right * x.right
#' # add some noise
#' y.left <- y.left + rnorm(length(y.left))
#' y.right <- y.right + rnorm(length(y.right))
#' # get y vector
#' x <- c(x.left, x.right)
#' y <- c(y.left, y.right)
#' # find best fit bilinear model
#' bf.model <- best_fit_bilinear(x, y)
#' print(bf.model[["coefficients"]])
#' print(bf.model[["intersection.x"]])
#' # plot results: points
#' plot(x, y, cex = 0.1)
#' # plot left model
#' abline(a = a.left, b = b.left, col = "blue")
#' # plot right model
#' abline(a = a.right, b = b.right, col = "green")
#' # plot left model fit
#' abline(a = bf.model[["coefficients"]]["left","intercept"], b = bf.model[["coefficients"]]["left","slope"], col = "blue", lty = 2)
#' # plot left model fit
#' abline(a = bf.model[["coefficients"]]["right","intercept"], b = bf.model[["coefficients"]]["right","slope"], col = "green", lty = 2)
#'
#' @export
best_fit_bilinear <- function(x.vec, y.vec, truncate.left = 0, truncate.right = 0){
  stopifnot(length(x.vec) == length(y.vec))
  idx.inf <- which(y.vec == -Inf | y.vec == Inf)
  x <- x.vec[-idx.inf][(1 + truncate.left):(length(x.vec) - length(idx.inf) - truncate.right)]
  y <- y.vec[-idx.inf][(1 + truncate.left):(length(y.vec) - length(idx.inf) - truncate.right)]
  fit.bilinear <- function(Cx){
    lhs <- function(x) ifelse(x < Cx,Cx-x,0)
    rhs <- function(x) ifelse(x < Cx,0,x-Cx)
    fit <- lm(y ~ lhs(x) + rhs(x))
    fs <- summary(fit)
    c(fs$r.squared, fs$coef[1], fs$coef[2], fs$coef[3])
  }
  # optimize
  score.fun <- function(X) -(fit.bilinear(X)[1])
  opt.x <- optimize(score.fun, interval = c(x[3],x[length(x)-2]))
  opt.params <- c(opt.x$minimum, fit.bilinear(opt.x$minimum))
  # get coefficients for left linear model
  a.left <- opt.params[3] + opt.params[1] * opt.params[4] # intercept
  b.left <- -opt.params[4] # slope
  # get coefficients for right linear model
  a.right <- opt.params[3] - opt.params[1] * opt.params[5] # intercept
  b.right <- opt.params[5] # slope

  # closest point to left linear model from right side
  dev.l <- (a.left + b.left * x) - y
  xl <- tail(which(dev.l <= 0), 1)
  # intersection of both linear models
  xi <- (a.right - a.left) / (b.left - b.right)
  # closest point to right linear model from left side
  dev.r <- (a.right + b.right * x) - y
  xr <- which(dev.r <= 0)[1]

  # format result result
  intersection.x <- set_names(c(xl, xi, xr), c("left", "both", "right"))
  coefs <- matrix(c(a.left, b.left, a.right, b.right), ncol = 2, byrow = TRUE)
  colnames(coefs) <- c("intercept", "slope")
  rownames(coefs) <- c("left", "right")
  # return result
  set_names(list(coefs, intersection.x), c("coefficients", "intersection.x"))
}

#' Finds local maximas indices.
#'
#' Seeks for local maximas in vector v using simple second order difference.
#'
#' @param v numeric vector
#'
#' @return numeric vector with indices of local maxima elements of vector v (i.e. where elements of second order difference equals -2)
#'
#' @examples
#' v <- c(1,2,3,4,5,6,5,3,2,-1,-4,-3,-1,2,4,10,12,11,5,3)
#' idx <- local.max(v)
#' print(idx)
#' print(v[idx])
local.max <- function(v) which(diff(sign(diff(v))) == -2) + 1

#' Finds local minimas indices.
#'
#' Seeks for local minimas in vector v using simple second order difference.
#'
#' @param v numeric vector
#'
#' @return numeric vector with indices of local minima elements of vector v (i.e. where elements of second order difference equals 2)
#'
#' @examples
#' v <- c(1,2,3,4,5,6,5,3,2,-1,-4,-3,-1,2,4,10,12,11,5,3)
#' idx <- local.min(v)
#' print(idx)
#' print(v[idx])
local.min <- function(v) which(diff(sign(diff(v))) == 2) + 1

#' Find first local maximum.
#'
#' Seeks first local maximum starting from the right hand side of given vector and moving towards left hand side.
#'
#' @param v numeric vector
#'
#' @return numeric value of first local maximum in \code{v} starting from the end of v and going left; NA if v is nonincresing starting at last element of v and going left (i.e. v is either constant or there is minmum first)
#'
#' @seealso \code{\link{CopulaHiC::local.min}}, \code{\link{CopulaHiC::local.max}} for finding local minma and maxima in vector v
#'
#' @examples
#' # maximum in 7 (index 2) starting from 1 (index 8) and moving with decreasing indices
#' v1 <- c(5,6,7,6,4,3,2,1)
#' # no maximum or minimum
#' v2 <- c(2,2,2,2,2,2)
#' # minimum at 2 (index 3) and no maximum starting at 6 (index 7)
#' v3 <- c(4,3,2,3,4,5,6)
#' # starting from rightmost element (3 with index 14) and going left, minimum first (in -3, index 9) then maximum (in 6, index 4)
#' v4 <- c(3,4,5,6,4,3,2,-1,-3,-2,0,1,2,3)
#' print(left.max(v1))
#' print(left.max(v2))
#' print(left.max(v3))
#' print(left.max(v4))
left.max <- function(v){
  fmin <- local.min(rev(v))[1]
  fmax <- local.max(rev(v))[1]
  if(is.na(fmax)){
    return(NA)
  } else if(is.na(fmin) | (fmin > fmax)){
    return(v[length(v) - fmax + 1])
  } else {
    return(NA)
  }
}

#' Find first local minimum.
#'
#' Seeks first local minimum starting from the left hand side of given vector and moving towards right hand side.
#'
#' @param v numeric vector
#'
#' @return numeric value of first local minimum in \code{v} starting from the begining of v and going right; NA if v is nondecreasing starting at first element of v and going right (i.e. v is either constant or there is maximum first)
#'
#' @seealso \code{\link{CopulaHiC::local.min}}, \code{\link{CopulaHiC::local.max}} for finding local minma and maxima in vector v
#'
#' @examples
#' # maximum first (in 7, index 3)
#' v1 <- c(5,6,7,6,4,3,2,1)
#' # no maximum or minimum
#' v2 <- c(2,2,2,2,2,2)
#' # minimum at 2 (index 3) and no maximum starting at 6 (index 7)
#' v3 <- c(4,3,2,3,4,5,6)
#' print(left.max(v1))
#' print(left.max(v2))
#' print(left.max(v3))
right.min <- function(v){
  fmin <- local.min(v)[1]
  fmax <- local.max(v)[1]
  if(is.na(fmin)){
    return(NA)
  } else if(is.na(fmax) | (fmax > fmin)){
    return(v[fmin])
  } else {
    return(NA)
  }
}

#' Determines TAD boundaries using Insulation Score.
#'
#' For details on Insulation Score approach of TAD boundaries detection see Crane et al. 2015 "Condensin-driven remodelling of X chromosome topology during dosage compensation", this function is implementation of \url{http://bioinformaticsinstitute.ru/sites/default/files/tad_calling_methods_part_1_-_sidorov_16-dec-2016.pdf} with minor adaptations to handle missing values.
#'
#' @param dense.mtx numeric matrix in dense format representing Hi-C contact map
#' @param resolution numeric resolution of Hi-C contact map in base pairs
#' @param window.bp numeric size of sliding window in base pairs
#' @param delta.bp numeic size of delta window in base pairs
#'
#' @return data frame with TAD positions containing start, end columns
#'
#' @examples
#' # get Hi-C contact map
#' sparse.mtx <- CopulaHiC::sample_hic_maps[["MSC-HindIII-1_40kb-raw"]][["18"]]
#' dense.mtx <- sparse2dense(sparse.mtx[c("i","j","val")], N = 1952)
#' # get tads
#' tads <- map2tads(dense.mtx)
#' # plot results
#' plot_contact_map(dense.mtx)
#' plot_regions(tads)
#'
#' @export
map2tads <- function(dense.mtx, resolution = 40000, window.bp = 1000 * 1000, delta.bp = 200 * 1000){
  w.bins <- window.bp / resolution
  d.bins <- delta.bp / resolution
  unmappable <- remove_unmappable(dense.mtx)
  mtx <- unmappable[["matrix"]]
  N <- nrow(mtx)
  bins <- seq(1,N)
  # 1. get insulations scores, first w.bins will be equal to iscore[w.bins+1] and
  # last w.bins wil be equal to iscore[N-w.bins] - to avoid NA, Nan, Inf values
  iscore <- sapply(bins[(w.bins+1):(N-w.bins)], function(i){
    mean(mtx[(i-w.bins):i-1, (i+1):(i+w.bins)], na.rm = TRUE)
  })
  iscore <- c(rep(iscore[1], w.bins), iscore, rep(iscore[length(iscore)], w.bins))
  mean.is <- mean(iscore, na.rm = TRUE)
  normalized.is <- log2(iscore / mean.is)
  # 2. get deltas
  sapply(bins[2:(length(bins)-1)], function(i){
    l.idx <- seq(i-w.bins, i-1)
    r.idx <- seq(i+1, i+1+w.bins)
    mean(normalized.is[l.idx[l.idx > 0]], na.rm = TRUE) -
      mean(normalized.is[r.idx[r.idx <= N]], na.rm = TRUE)
  }) -> deltas
  deltas <- c(deltas[1], deltas, deltas[length(deltas)])
  # get gaps, i.e. intervals with NA and convert them to what is on the left side
  non.na <- which(!is.na(deltas))
  l <- which(diff(non.na) != 1)
  starts <- non.na[l]
  ends <- non.na[l+1] - 1
  if(length(starts) > 0){
    for(j in seq(1,length(starts))){
      deltas[(starts[j]+1):ends[j]] <- starts[j]
    }
  }
  # 3. get minimum insulation score
  # get indices of delta closest to 0
  idx <- which(diff(sign(deltas)) != 0)
  # for each zero delta find:
  # a) its first left maximum (or NA if it decreases on the left from this zero delta)
  # b) its first left minimum (or NA if it increases on the right from this zero delta)
  sapply(idx, function(i){
    c(i, left.max(deltas[1:i]), right.min(deltas[i:N]))
  }) -> Si
  # now get only these indices for which deltaMax - deltaMin > 0.1
  boundaries <- unique(c(1, Si[1, which((Si[2,] - Si[3,]) > 0.1)], length(deltas)))
  v <- rep(NA, length(deltas))
  v[boundaries] <- 1
  v.restored <- restore_unmappable_vec(v, unmappable[["indices.rows"]])
  br <- which(!is.na(v.restored))
  # convert to data frame
  data.frame(start = c(br[1], br[2:(length(br)-1)] + 1), end = br[2:length(br)])
}

#' Finds intersection between 2 sets of TADs.
#'
#' Intersection is defined for the same chromosomes, for example for 2 TAD sets from chromosome 18 if there is TAD in set1 with coordinates 100, 120 (start, end respectively) and in set2 there is a TAD with coordinates 110, 140 then their intersection is interval with coordinates: 110, 120.
#'
#' @param tads1 data frame of TAD set 1 with columns: start, end
#' @param tads2 data frame of TAD set 2 with columns: start, end
#'
#' @return data frame with start, end columns containing TAD intersection intervals
#'
#' @seealso \code{\link{intervals::Intervals}}, \code{\link{intervals::interval_intersection}} for functions used to find intersection between 2 sets of intervals
#'
#' @examples
#' # simple artificial data set of TADs
#' tads1 <- data.frame(start = c(1,10,30), end = c(5,15,35))
#' tads2 <- data.frame(start = c(3,14,28), end = c(12,18,40))
#' tads.intersection <- intersect_tads(tads1, tads2)
#' print(tads.intersection)
#'
#' @export
intersect_tads <- function(tads1, tads2){
  itads1 <- intervals::Intervals(tads1[c("start","end")])
  itads2 <- intervals::Intervals(tads2[c("start","end")])
  intervals::interval_intersection(itads1, itads2) %>%
    as.data.frame() %>% magrittr::set_colnames(c("start","end"))
}
