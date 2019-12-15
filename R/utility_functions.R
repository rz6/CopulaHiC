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
#' @param return.idx logical whether to return indices of this diagonal
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
superdiag <- function(A, k, return.idx = FALSE){
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
  if(return.idx){
    return(indices)
  }
  return(A[indices])
}

#' Convert interactions into zscores
#'
#' Given single (i.e. one chromosome) contact map in sparse format convert its interactions into zscores diagonal-wise.
#'
#' @param mtx.sparse data frame representing contact map
#' @param substitute.val logical indicating whether zscores should be substituted into val
#'
#' @return data frame with interaction zscores calculated over every diagonal separatelly
#'
#' @examples
#'
#' @export
interactions2zscores <- function(mtx.sparse, substitute.val = TRUE){
  if("name" %in% colnames(mtx.sparse)){
    # check if single matrix
    stopifnot(length(unique(as.character(mtx.sparse$name))) == 1)
  }
  by(mtx.sparse, mtx.sparse$diagonal, function(df){
    if(substitute.val){
      magrittr::inset(df, "val", value = base::scale(df$val)[,1])
    } else {
      magrittr::inset(df, "zscore.val", value = base::scale(df$val)[,1])
    }
  }) %>% do.call("rbind",.)
}

#' Calculate proportion of interactions
#'
#' Appends column with proportion of interactions per diagonal.
#'
#' @param mtx.sparse data frame representing contact map
#'
#' @return data frame with interactions proportion calculated over every diagonal separatelly
#'
#' @examples
#'
#' @export
proportions <- function(mtx.sparse){
  if("name" %in% colnames(mtx.sparse)){
    # check if single matrix
    stopifnot(length(unique(as.character(mtx.sparse$name))) == 1)
  }
  by(mtx.sparse, mtx.sparse$diagonal, function(df){
    magrittr::inset(df, "prop", value = df$val / sum(df$val))
  }) %>% do.call("rbind",.)
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
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
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
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
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
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
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
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
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
      magrittr::set_colnames(c("start","end")) %>%
      magrittr::inset("compartment", value = x)
  }) %>%
    do.call("rbind",.) -> df
  df[order(df$start, df$end), ] %>%
    magrittr::set_rownames(NULL)
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
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
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

#' Move contacts to given pair of interacting regions
#'
#' Randomly selects \code{N} contacts from \code{i} xor \code{j} bins of given \code{mtx.sparse} Hi-C contact map and moves them to cell (\code{i}, \code{j}). This method of producing artificial LRI preserves coverage of initial Hi-C contact map.
#'
#' @param mtx.sparse data frame - Hi-C contact map in sparse format with columns: i, j, val, diagonal, name
#' @param i numeric first region (row of Hi-C contact map)
#' @param j numeric second region (column of Hi-C contact map)
#' @param N numeric number of contacts to be moved from regions \code{i} xor \code{j} to interaction (\code{i},\code{j})
#'
#' @return data frame with the same columns as \code{mtx.sparse}
#'
#' @export
move_contacts <- function(mtx.sparse, i, j, N){
  # select contacts for sampling, i.e.: i-th row xor j-th column
  v <- xor(mtx.sparse$i == i, mtx.sparse$j == j)
  contacts <- sparse2interactions(mtx.sparse[v,])
  n <- nrow(contacts)
  stopifnot(n > N)
  remaining.contacts <- mtx.sparse[!v,]
  sampled <- sample(seq(1, n), replace = FALSE, size = N)
  not.sampled <- contacts[-sampled,]
  sparse.not.sampled <- aggregate(N ~ ., data = transform(not.sampled, N = 1), FUN = length)
  cn <- colnames(sparse.not.sampled)
  colnames(sparse.not.sampled)[which(cn == "N")] <- "val"
  remaining.contacts[(remaining.contacts$i == i) & (remaining.contacts$j == j),"val"] <- remaining.contacts[(remaining.contacts$i == i) & (remaining.contacts$j == j),"val"] + N
  rbind(
    remaining.contacts,
    sparse.not.sampled[c("i","j","val","diagonal","name")]
  )
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
bootstrap_interactions <- function(interactions, ratio = c(0.5, 0.5), with.replacement = FALSE){
  if(with.replacement){
    # for sampling with replacement ratio should be a vector containing number of interactions to draw
    stopifnot(all(ratio %% 1 == 0))
    lapply(ratio, function(int.number){
      interactions[sample(nrow(interactions), int.number, replace = TRUE),]
    }) %>% magrittr::set_names(seq(1,length(ratio))) -> interactions.ratios
  } else {
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
  }
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
#' @param ratio numeric vector indicating on how many atomic interaction sets should initial atomic interaction set be divided; the values inside this vector depend on which sampling is used (with or without replacement); each entry of ratio vector contains fraction of interaction to be put in corresponding atomic interactions set; ratio vector must sum to 1 and all its entries must be larger than one
#' @param with.replacement logical indicating whether to perform sampling with or without (default) replacement
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
bootstrap_sparse <- function(sparse.mtx, ratio = c(0.5, 0.5), with.replacement = FALSE){
  # divides Hi-C contact map on length(ration) maps according to specified ratios,
  # where each map contains fraction of contacts from original contact map
  # initial contact map is given in sparse data.frame format:
  # i j val chromosome
  if(with.replacement){
    stopifnot(all(ratio %% 1 == 0))
  } else {
    stopifnot(sum(ratio) == 1)
  }
  stopifnot(all(ratio > 0))
  sparse2interactions(sparse.mtx) %>%
    bootstrap_interactions(ratio = ratio, with.replacement = with.replacement) %>%
    split(magrittr::use_series(.,"ratio.number"))
}

#' Bootstraps interactions from given Hi-C dataset.
#'
#' Dataset containing number of Hi-C contact maps is used to sample interactions with or without replacement.
#'
#' @param path character path to Hi-C dataset in npz format
#' @param ratio list with vectors - ratios for \code{\link{bootstrap_sparse}}; if contains only one element (one ratio vector) then the same ratio vector is used for all contact maps in given Hi-C dataset; otherwise names of elements (ratio vectors) in list must match those in given Hi-C dataset
#' @param with.replacement logical which type of sampling
#' @param N numeric number of repetitions, i.e. number of bootstraps; each bootstrap will have number of maps equal to length of corresponding entry in ratio list
#' @param mtx.names character vector with subset of Hi-C maps names to be selected for analysis, by default all matrices are used
#' @param n.cores numeric number of cores to be used for parallel processing
#'
#' @return list containing bootstraps of corresponding matrices of Hi-C dataset
#'
#' @examples
#' # say we have 2 Hi-C datasets: IMR90-MboI-1 and MSC-HindIII-1 in 40kb
#' npz1 <- read_npz("IMR90-MboI-1_40kb-raw.npz")
#' npz2 <- read_npz("MSC-HindIII-1_40kb-raw.npz")
#' # we want to produce 2*4 bootstraps of IMR90:
#' # 4 with the same number of interactions as in IMR90-MboI-1 and
#' # 4 with the same number of interactions as in MSC-HindIII-1
#' # first calculate number of interactions in both datasets on all chromosomes
#' nm <- intersect(names(npz1), names(npz2))
#' ratio <- lapply(nm, function(x) c(sum(npz1[[x]]$val), sum(npz2[[x]]$val)))
#' names(ratio) <- nm
#' # now produce bootstraps - pairs of bootstrap maps such that the number of interactions corresponds to first and second datasets can be used to asses technical variability including different sequencing depth
#' bts <- bootstrap_dataset("IMR90-MboI-1_40kb-raw.npz", N = 4, with.replacement = TRUE, ratio = ratio, save.path = "~/bootstrapped")
#'
#' @export
bootstrap_dataset <- function(path, ratio = list(c(0.5, 0.5)), with.replacement = FALSE, N = 3, mtx.names = "all", n.cores = 1, save.path = NULL){
  stopifnot(typeof(ratio) == "list")
  sizes <- read_size(path, mtx.names = mtx.names)
  maps.sparse <- read_npz(path, mtx.names = mtx.names)
  maps.names <- names(maps.sparse)
  if(is.null(names(ratio))){
    ratio <- rep(ratio, length(maps.names))
    names(ratio) <- maps.names
  }
  all.names <- intersect(maps.names, names(ratio))
  lapply(all.names, function(x){
    parallel::mclapply(1:N, function(x, m, r, wr){
      bootstrap_sparse(m, ratio = r, with.replacement = wr) %>%
        lapply(function(df) df[c("i","j","val","diagonal","name")])
    }, maps.sparse[[x]], ratio[[x]], with.replacement, mc.cores = n.cores) %>%
      magrittr::set_names(1:N) %>%
      unlist(recursive = FALSE) %>%
      magrittr::set_names(names(.) %>% gsub("\\.","-",.))
  }) %>% magrittr::set_names(all.names) -> bootstrapped
  # if save path is given then save to file
  if(!is.null(save.path)){
    dir.create(save.path, showWarnings = FALSE)
    name.prefix <- sub('\\..*$','', basename(path))
    bootstrapped.swap <- datasets <- do.call(function(...) Map(list, ...), bootstrapped)
    for(dataset.name in names(bootstrapped.swap)){
      fname <- paste0(name.prefix,"_bootstrap_",dataset.name,".npz")
      dataset <- bootstrapped.swap[[dataset.name]]
      save_npz(lapply(names(dataset), function(x){
        sparse2dense(dataset[[x]], N = sizes[x,"n.rows"])
      }) %>% magrittr::set_names(names(dataset)), paste0(save.path,"/",fname))
    }
    return(NULL)
  }
  return(bootstrapped)
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
    magrittr::extract(c(c("i","j","tad.id","start","end"), cols))
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
  if(length(idx.inf) > 0){
    x <- x.vec[-idx.inf][(1 + truncate.left):(length(x.vec) - length(idx.inf) - truncate.right)]
    y <- y.vec[-idx.inf][(1 + truncate.left):(length(y.vec) - length(idx.inf) - truncate.right)]
  } else {
    x <- x.vec[(1 + truncate.left):(length(x.vec) - length(idx.inf) - truncate.right)]
    y <- y.vec[(1 + truncate.left):(length(y.vec) - length(idx.inf) - truncate.right)]
  }
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
  intersection.x <- magrittr::set_names(c(xl, xi, xr), c("left", "both", "right"))
  coefs <- matrix(c(a.left, b.left, a.right, b.right), ncol = 2, byrow = TRUE)
  colnames(coefs) <- c("intercept", "slope")
  rownames(coefs) <- c("left", "right")
  # return result
  magrittr::set_names(list(coefs, intersection.x), c("coefficients", "intersection.x"))
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
#' @seealso \code{\link[DIADEM]{local.min}}, \code{\link[DIADEM]{local.max}} for finding local minma and maxima in vector v
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
#' @seealso \code{\link[DIADEM]{local.min}}, \code{\link[DIADEM]{local.max}} for finding local minma and maxima in vector v
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
#' For details on Insulation Score approach of TAD boundaries detection see Crane et al. 2015 "Condensin-driven remodelling of X chromosome topology during dosage compensation" methods section, paragraph TAD calling (insulation square analysis).
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
#' sparse.mtx <- DIADEM::sample_hic_maps[["MSC-HindIII-1_40kb-raw"]][["18"]]
#' dense.mtx <- sparse2dense(sparse.mtx[c("i","j","val")], N = 1952)
#' # get tads
#' tads <- map2tads(dense.mtx)
#' # plot results
#' plot_contact_map(dense.mtx)
#' plot_regions(tads)
#'
#' @export
map2tads <- function(dense.mtx, resolution = 40000, window.bp = 500 * 1000, delta.bp = 100 * 1000, without_unmappable = TRUE){
  w.bins <- window.bp / resolution
  d.bins <- delta.bp / resolution
  N <- nrow(dense.mtx)
  bins <- seq(1,N)
  # 1. get insulations scores
  sapply(bins[(w.bins+1):(N-w.bins)], function(i){
    mean(dense.mtx[(i-w.bins):i-1, (i+1):(i+w.bins)], na.rm = TRUE)
  }) -> iscore
  mean.is <- mean(iscore, na.rm = TRUE)
  # below may have -Inf where iscore == 0
  normalized.is <- log2(iscore / mean.is)
  # so replace -Inf with the smallest value from normalized.is (but not -Inf) times 2
  normalized.is[normalized.is == -Inf] <- min(normalized.is[normalized.is != -Inf]) * 2
  # 2. get deltas
  n <- length(normalized.is)
  sapply(1:n, function(i){
    l.idx <- seq(i-w.bins, i-1)
    r.idx <- seq(i+1, i+1+w.bins)
    if(i == 1){
      -mean(normalized.is[r.idx[r.idx <= N]], na.rm = TRUE)
    } else if(i == n){
      mean(normalized.is[l.idx[l.idx > 0]], na.rm = TRUE)
    } else{
      mean(normalized.is[l.idx[l.idx > 0]], na.rm = TRUE) - mean(normalized.is[r.idx[r.idx <= N]], na.rm = TRUE)
    }
  }) -> deltas
  # append NA for first and last window bins
  na.bins <- rep(NA, w.bins)
  deltas <- c(na.bins, deltas, na.bins)
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
  boundaries <- Si[1, which((Si[2,] - Si[3,]) > 0.1)]
  tads <- data.frame(start = boundaries[1:(length(boundaries)-1)], end = boundaries[2:length(boundaries)])
  if(without_unmappable){
    unmappable_idx <- remove_unmappable(dense.mtx)[["indices.rows"]]
    split(unmappable_idx, cumsum(c(1, diff(unmappable_idx) != 1))) %>%
      sapply(FUN = function(v){ c(min(v),max(v)) }) %>% as.vector() -> unmappable_boundaries
    c(1, unmappable_boundaries, nrow(dense.mtx)) %>%
      matrix(ncol = 2, byrow = TRUE) -> mappable
    boundaries <- sort(c(boundaries, unmappable_boundaries))
    lapply(1:nrow(mappable), function(x){
      b <- boundaries[(boundaries >= mappable[x,1]) & (boundaries <= mappable[x,2])]
      data.frame(start = b[1:(length(b)-1)], end = b[2:length(b)])
    }) %>% do.call("rbind",.) -> tads
  }
  return(tads)
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
#' @seealso \code{\link[intervals]{Intervals}}, \code{\link[intervals]{interval_intersection}} for functions used to find intersection between 2 sets of intervals
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

#' Extracts neighborhood around i,j cells of given matrix.
#'
#' Given matrix \code{mtx.dense} and positive integer \code{k} it extracts cross pattern cells:
#' \itemize{
#'  \item{}{i - k, j}
#'  \item{}{i, j - k}
#'  \item{}{i - (k-1), j}
#'  \item{}{i, j - (k-1)}
#'  \item{}{.}
#'  \item{}{.}
#'  \item{}{.}
#'  \item{}{i - 1, j}
#'  \item{}{i, j - 1}
#'  \item{}{i + 1, j}
#'  \item{}{i, j + 1}
#'  \item{}{.}
#'  \item{}{.}
#'  \item{}{.}
#'  \item{}{i + (k - 1), j}
#'  \item{}{i, j + (k - 1)}
#'  \item{}{i + k, j}
#'  \item{}{i, j + k}
#'  \item{}{i, j}
#' }
#' or donut shape. The above cells are columns in resulting data frame, so it has in total 2 + 4 * k + 1 columns: first two are cell coordinates in initial matrix, i.e.: i and j, remaining columns are the above values (in the order specified above). Selection of cross size \code{k} enforces that cells belonging to rows, columns or diagonals (starting at main diagonal) with indices 1 to \code{k} or \code{n-k+1} to \code{n} can't be assigned any value.
#'
#' @param mtx.dense numeric symmetric matrix
#' @param k numeric radius size
#' @param max.d numeric fraction of domains to be taken, i.e.: if \code{mtx.dense} is of size \code{n} (i.e. it have \code{n} diagonals) then \code{floor(max.d * n)} diagonals, starting at main diagonal will be taken and remaining will be discarded.
#' @param m numeric how many diagonals (starting at main diagonal) to discard, by definition of cross this must be fixed to m > k to have effect as m will be finally be selected as: \code{max(k,m)}
#' @param which.ngb character indicating the shape of neighborhood to extract, either "cross" or "donut"
#' @param include.ij logical whether to include the value of cell for which neighborhood is extracted
#'
#' @examples
#'
#' @export
neighborhood <- function(mtx.dense, k = 1, max.d = 0.15, m = 1,
                         which.ngb = c("cross","donut")[1], include.ij = FALSE){
  stopifnot(which.ngb %in% c("cross","donut"))
  n <- nrow(mtx.dense)
  N <- floor(n * max.d)
  # copy matrix and zero m diagonals, k left/right and top/bottom columns and rows respectively
  mtx <- mtx.dense
  i <- row(mtx)
  j <- col(mtx)
  mtx[i < j] <- 0
  idx.diags <- c(seq(1,max(m,k)), seq(N+1,n))
  idx.diags.mtx <- unlist(sapply(seq(1,n), function(x){
    unlist(sapply(idx.diags[idx.diags + x - 1 <= n], function(y){
      n * (x - 1) + (x - 1) + y
    }))
  }))
  mtx[idx.diags.mtx] <- 0
  if(k > 0){
    mtx[seq(1,k),] <- 0
    mtx[,seq(1,k)] <- 0
    mtx[seq(n-k+1,n),] <- 0
    mtx[,seq(n-k+1,n)] <- 0
  }
  # get indices of non zero cells
  idx <- which(mtx != 0, arr.ind = TRUE)
  # for each non zero cell get k cells to the left/right and k cells above/below
  if(k > 0){
    if(which.ngb == "cross"){
      v <- c(seq(-k,-1), seq(1,k))
      ij.idx <- rbind(
        data.frame(Var1 = v, Var2 = rep(0, length(v))),
        data.frame(Var1 = 0, Var2 = 0),
        data.frame(Var1 = rep(0, length(v)), Var2 = v)
      )
    }else{
      ij.idx <- expand.grid(seq(-k,k), seq(-k,k))
    }
    if(!include.ij)
      ij.idx <- ij.idx[-ceiling(nrow(ij.idx) / 2),]
    cnms <- apply(ij.idx, 1, function(x) paste0("i", x[1], ".j", x[2]))
    ###
    apply(idx, 1, function(x){
      apply(ij.idx, 1, function(ij){
        mtx.dense[x["row"] + ij[1], x["col"] + ij[2]]
      })
      # if(which.ngb == "cross"){
      #   if(include.ij)
      #     v <- seq(-k,k)
      #   else
      #     v <- c(seq(-k,-1), seq(1,k))
      #   sapply(v, function(y){
      #     c(mtx.dense[x["row"]+y,x["col"]], mtx.dense[x["row"],x["col"]+y])
      #   }) %>% unlist()
      # }else{
      #   ij.idx <- expand.grid(seq(-k,k), seq(-k,k))
      #   if(!include.ij)
      #     ij.idx <- ij.idx[-ceiling(nrow(ij.idx) / 2),]
      #   apply(ij.idx, 1, function(ij){
      #     mtx.dense[x["row"] + ij[1], x["col"] + ij[2]]
      #   })
      # }
    }) %>% t() %>%
      magrittr::set_colnames(cnms) -> ij.crosses
  } else {
    mtx.dense[idx] -> ij.crosses
  }
  df <- as.data.frame(cbind(idx, ij.crosses))
  #df <- df[c("col","row", colnames(df)[3:ncol(df)])]
  colnames(df)[1:2] <- c("i","j")
  return(df)
}

#' Aggregates diagonals of based on similarity of XY distribution
#'
#' Merges 2 contact maps in sparse format and then aggregates it diagonals. The aggregation process is performed as follows:
#' \itemize{
#'   \item{}{start with k = 1 (first diagonal), group k consists of observations: val.x * val.y (from k-th diagonal)}
#'   \item{}{take merge candidate i.e.: k' = k + 1, , group k' consists of observations: val.x * val.y (from k'-th diagonal)}
#'   \item{}{perform chi-square test between group k and k', H0: observations from both groups come from the same distribution, H1: they come from different distributions}
#'   \item{}{if rejected set k = k' (and pool diagonals k to k'-1 into single group), otherwise set k' = k' + 1}
#'   \item{}{repeat until all diagonals are processed}
#' }
#' After algorithm is finished diagonals are assigned to groups (pools).
#'
#' @param mtx1.sparse data frame, contact map in sparse format
#' @param mtx2.sparse data frame, contact map in sparse format
#' @param agg.diags logical if false leave each diagonal in separate group (no pooling)
#' @param alpha numeric significance threshold for chi square test
#'
#' @examples
#'
#' @export
aggregate_diagonals <- function(mtx1.sparse, mtx2.sparse, agg.diags = TRUE, which.test = c("energy","KS")[1],
                                alpha = 0.05, exclude.outliers = FALSE){
  stopifnot(which.test %in% c("energy","KS"))
  # helper function
  get_next <- function(mtx.by.diag, k = 1, alpha = 0.05, wt = "energy", bn = 1000){
    k.range <- names(mtx.by.diag)
    n <- max(as.numeric(k.range))
    df1 <- mtx.by.diag[[as.character(k)]]
    while (TRUE) {
      k <- k + 1
      if(k > n){
        return(k)
      }
      if(as.character(k) %in% k.range){
        df2 <- mtx.by.diag[[as.character(k)]]
        if(wt == "energy"){
          # Energy test is prefered as it tackles multivariate data, unfortunately for large vectors of observations it takes long time to compute
          #p <- cramer::cramer.test(as.matrix(df1[c("val.x","val.y")]), as.matrix(df2[c("val.x","val.y")]))$p.value
          md1 <- as.matrix(df1[df1$outlier == "no", c("val.x","val.y")])
          md2 <- as.matrix(df2[df2$outlier == "no", c("val.x","val.y")])
          p <- energy::eqdist.etest(rbind(md1, md2), c(nrow(md1), nrow(md2)), R = bn)$p.value
        } else{
          p <- ks.test(df1$val.x * df1$val.y, df2$val.x * df2$val.y)$p.value
        }
        if(p <= alpha){
          return(k)
        }
      }
    }
  }
  ################
  ################
  mg <- base::merge(mtx1.sparse, mtx2.sparse, by = c("i", "j", "diagonal", "name"))
  mg$outlier <- "no"
  # exclude.outliers <- FALSE
  # if(exclude.outliers){
  #   by(mg, mg$diagonal, function(df){
  #     tryCatch({
  #       #fit <- robustreg::robustRegBS(val.y ~ val.x, data = df)
  #       #df[which(fit$weights <= 0), "outlier"] <- "yes"
  #       fit <- robustbase::lmrob(val.y ~ val.x, df, setting = "KS2014")
  #       fit.summary <- summary(fit)
  #       thr <- fit.summary$control$eps.outlier
  #       df[fit.summary$rweights <= thr, "outlier"] <- "yes"
  #     }, error = function(e){
  #       return(df)
  #     })
  #     return(df)
  #   }) %>% do.call("rbind", .) -> mg
  # }
  # split by diagonal
  mgs <- split(mg, mg$diagonal)
  diags <- as.numeric(names(mgs))
  if(agg.diags){
    n <- max(diags)
    k <- 1
    v <- c(k)
    while(k < n){
      k <- get_next(mgs, k = k, alpha = alpha, wt = which.test)
      v <- c(v, k)
    }
    dpool <- data.frame(diagonal = diags, dpool = findInterval(diags, v))
  } else{
    dpool <- data.frame(diagonal = diags, dpool = diags)
  }
  return(dpool)
}

#' Constructs GLM to model per-diagonal Hi-C contact dependancy
#'
#' Models Hi-C contacts using (robust) Negative Binomial (or Poisson when the data is underdispersed) regression. Given the fact that Hi-C data suffers from contact decay bias this method is intendent to model each diagonal separately. By default this function uses robust Negative Binomial regression to model interaction dependencies.
#'
#' @param df data frame with predictor, response, outlier columns
#' @param robust.nb logical whther to use robust fitting procedure (see details)
#' @param overdisp.test.pval numeric significance threshold for testing ovedispersion
#' @param max.nobs numeric maximum number of observations (points), i.e. sample size to be taken for robust NB regression estimation (see details)
#' @param nrep numeric number of repetitions for subsampling (see details)
#'
#' @return object of class glm or MASS::glm.nb
#'
#' @details If \code{robust.nb} is true then this function uses robust Negative Binomial estimation method developed in \insertCite{aeberhard2014robust}{DIADEM}. This function uses the code of glmrob.nb function written by William Aeberhard, which is available at: https://github.com/williamaeberhard/glmrob.nb.
#'
#' At first overdispersion test is performed to decide if Negative Binomial or Poisson regression should be used. If \code{robust.nb} is true the estimation may consume huge amounts of memory for large sample sizes (like for example 400000 points). In order to prevent that whenever the sample size exceeds \code{max.nobs} initial sample is subsampled to \code{max.nobs} size and model is estimated on subsample. This procedure is repeated \code{nrep} times and final parameter estimate equals average over subsampled estimates.
#'
#' @seealso \code{\link{glm}}, \code{\link[MASS]{glm.nb}} to see how GLM are constructed, \code{\link[AER]{dispersiontest}} to see how overdispersion is tested
#'
#' @references{
#'   \insertRef{aeberhard2014robust}{DIADEM}
#' }
#'
#' @export
constructGLM <- function(df, robust.nb = TRUE, overdisp.test.pval = 0.01, max.nobs = 20000, nrep = 10){
  tryCatch({
    # fit robust Poisson regression and test for overdispersion
    model <- robustbase::glmrob(response ~ predictor, data = df[df$outlier == "no",], family = poisson(link = "log"))
    params <- c(0, as.numeric(model$coefficients))
    overdisp.test <- AER::dispersiontest(model)
    if(overdisp.test$p.value <= overdisp.test.pval){
      # overdispersion exists - fit better model - NB regression
      if(robust.nb){
        # check sample size
        df.noout <- df[df$outlier == "no",]
        if(nrow(df.noout) > max.nobs){
          # subsample nrep times and average estimates
          sapply(1:nrep, function(i){
            subdf <- df.noout[sample(1:nrow(df.noout), max.nobs),]
            model <- glmrob.nb(subdf$response, subdf$predictor,
                               bounding.func = 'T/T', weights.on.x='none',
                               c.tukey.beta = 4, c.tukey.sig = 4)
            # dispersion, intercept, beta
            as.numeric(c(model$coef[1], model$coef[2:3]))
          }) %>% apply(1, mean) -> params
        }else{
          model <- glmrob.nb(df.noout$response, df.noout$predictor,
                             bounding.func = 'T/T', weights.on.x='none',
                             c.tukey.beta = 4, c.tukey.sig = 4)
          # dispersion, intercept, beta
          params <- as.numeric(c(model$coef[1], model$coef[2:3]))
        }
      } else{
        model <- MASS::glm.nb(response ~ predictor, data = df[df$outlier == "no",], link = "log")
        if(is.null(model$th.warn)){
          # if model converged - get NB regression params, otherwise keep Poisson regression params
          params <- c(1 / model$theta, as.numeric(model$coefficients))
        }
      }
    }
  }, error = function(e){
    params <- c(NA, NA, NA)
  })
  magrittr::set_names(params, c("dispersion", "intercept", "beta"))
}

#' Computes significance of data given the model
#'
#' Uses either Negative Binomial or Poisson to calculate pvalue of Y | X or X | Y.
#'
#' @param dataset data frame with predictor, response, outlier columns
#' @param model numeric, 3 element vector: dispersion, intercept, beta
#'
#' @return numeric vector with p-values (upper tail)
#'
#' @seealso \code{\link{pnbinom}}, \code{\link{ppois}} to see how to calculate significance
#'
#' @export
significance <- function(dataset, model){
  stopifnot(all(c("predictor","response") %in% colnames(dataset)))
  # get model paarmeters
  dispersion <- model["dispersion"]
  intercept <- model["intercept"]
  beta <- model["beta"]
  if(is.na(dispersion)){
    pval <- rep(NA, nrow(dataset))
  } else if(dispersion){
    theta <- 1 / dispersion
    # if there is overdispersion then NB regression
    sapply(1:nrow(dataset), function(i){
      pnbinom(dataset[i,"response"], size = theta,
              mu = exp(intercept + beta * log(dataset[i,"predictor"])),
              lower.tail = FALSE)
    }) -> pval
  } else{
    # if no overdispersion then poisson regression
    sapply(1:nrow(dataset), function(i){
      ppois(dataset[i,"response"],
            lambda = exp(intercept + beta * log(dataset[i,"predictor"])),
            lower.tail = FALSE)
    }) -> pval
  }
  return(pval)
}

#' Simulate corresponding interactions based on given interactions vector and model.
#'
#' @param interaction.vector integer vector with interactions
#' @param model numeric, 3 element vector: dispersion, intercept, beta
#' @param N numeric number of interactions to sample for every interaction in \code{interaction.vector}
#'
#' @return data frame with columns: predictor (interaction.vector) and response (simulated corresponding interactions)
#'
#' @examples
#'
#' @export
simulate_contacts <- function(interaction.vector, model, N = 1){
  dispersion <- model["dispersion"]
  intercept <- model["intercept"]
  beta <- model["beta"]
  counts <- table(interaction.vector) * N
  lapply(names(counts), function(val){
    n <- as.numeric(counts[val])
    predictor <- as.numeric(val)
    if(is.na(dispersion)){
      response <- rep(NA, n)
    } else if(dispersion){
      # NB regression
      theta <- 1 / dispersion
      response <- rnbinom(n, size = theta,
                          mu = exp(intercept + beta * log(predictor)))
    } else{
      # Poisson regression
      response <- rpois(n, exp(intercept + beta * log(predictor)))
    }
    data.frame(predictor, response) %>%
      magrittr::set_colnames(c("predictor", "response"))
  }) %>% do.call("rbind", .)
}

#' Randomly selects cells for artificial long range interactions.
#'
#' Toghether with \link{\code{DIADEM::simulate_map.HiCglm}} function this can be used to create Hi-C contact maps with artificial long range differential interactions (see examples).
#'
#' @param mtx numeric Hi-C contact map in dense format
#' @param N numeric number of cells to sample
#' @param max.diag numeric specifies diagonal range to sample interactions from
#' @param radius numeric number of cells to include in artificial interaction around (left, right, bottom, top or corners) with sampled interaction
#' @param alri.min.dist numeric minimum distance (in number of cells left, right, bottom, top or diagonally) between any pair of sampled interactions (between their centers)
#'
#' @return matrix with id, i, j columns indicating artificial long range interactions
#'
#' @details Interaction is a square centered on i,j with height i - radius, ..., i, ..., i + radius and width j - radius, ..., j, ..., j + radius
#'
#' The minmum distance between any pair of interactions \code{alri.min.dist} is measured as Euclidean distance between centers of 2 intractions (like i1, j1 and i2, j2)
#'
#' @export
sample_alri <- function(mtx, N = 100, max.diag = NA, radius = 1, alri.min.dist = 2 * radius + 1){
  if(is.na(max.diag)){
    max.diag <- nrow(mtx)
  }
  stopifnot(max.diag %in% 1:nrow(mtx))
  alri <- mtx
  i <- row(alri)
  j <- col(alri)
  # cut off interactions to keep proper range of diagonals and radius
  bmtx <- lower.tri(alri) | (i + (max.diag - radius) <= j) | (i >= j - (radius + 3))
  bmtx[1:radius,] <- TRUE
  bmtx[, seq(ncol(alri) - radius + 1, ncol(alri))] <- TRUE
  alri[bmtx] <- 0
  lapply(1:N, function(x){
    if(sum(alri) == 0){
      return(expand.grid(x, NA, NA))
    }
    alri.sparse <- dense2sparse(alri)
    sampled.ij <- as.numeric(alri.sparse[sample(nrow(alri.sparse), 1),c("i","j")])
    # add extra radius cells
    v <- seq(-alri.min.dist, alri.min.dist)
    ij <- expand.grid(sapply(v, function(y) sampled.ij[1] + y),
                      sapply(v, function(y) sampled.ij[2] + y))
    ij <- ij[(ij[,1] >= 1) & (ij[,1] <= nrow(alri)) & (ij[,2] >= 1) & (ij[,2] <= nrow(alri)),]
    alri[as.matrix(ij)] <<- 0
    expand.grid(x, sapply(seq(-radius, radius), function(y) sampled.ij[1] + y),
                sapply(seq(-radius, radius), function(y) sampled.ij[2] + y))
  }) %>% do.call("rbind",.) %>%
    as.matrix() %>% magrittr::set_colnames(c("id","i","j"))
}
