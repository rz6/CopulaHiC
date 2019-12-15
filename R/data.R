#' Sample data containing Hi-C contact maps as well as HiCcomparator and HiCglm objects.
#'
#' The list contains 5 elements:
#' \describe{
#'  \item{"IMR90-MboI-1"}{List of raw Hi-C contact maps in sparse format (i.e.: data frame with columns: i,j,val,diagonal,name) from GSE63525 study. The list only contains contact map of chromosome 21, first replicate.}
#'  \item{"MSC-HindIII-1"}{List of raw Hi-C contact maps in sparse format (i.e.: data frame with columns: i,j,val,diagonal,name) from GSE52457 study. The list only contains contact map of chromosome 21, first replicate.}
#'  \item{"HiCcomparator"}{HiCcomparator object constructed from the above 2 datasets with parameters as indicated in vignette.}
#'  \item{"HiCglm"}{HiCglm object constructed from the above 2 datasets with parameters as indicated in vignette.}
#'  \item{"chromosome.sizes"}{List with names matching names (chromosomes) from first two datasets and values equal to their sizes (number of bins in contact map).}
#' }
#'
#' @docType data
#'
#' @usage data(sample_hic_data)
"sample_hic_data"
