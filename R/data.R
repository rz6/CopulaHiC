#' Npz file with sample Hi-C contact maps dataset.
#'
#' Npz file containing Hi-C contact maps of human MSC in 40kb resolution. The data comes from Dixon et al., 2015 study and was processed by Imakaev et al., 2012 pipeline without iterative correction step (so it is raw data). It contains contact maps for chromosomes 17, 18, 19.
#'
#' @name MSC-HindIII-1_40kb-raw_maps
#'
#' @seealso \url{https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.savez_compressed.html} for npz format description, Dixon et al., 2015 "Chromatin architecture reorganization during stem cell differentiation" for study where this dataset comes from, Imakaev et al., 2012 "Iterative correction of Hi-C data reveals hallmarks of chromosome organization." for Iterative Correction of Hi-C contact maps and \url{https://mirnylab.bitbucket.io/hiclib/index.html} for its python implementation
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52457}
#'
#' @examples
#' # get file name
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
NULL

#' Npz file with sample Hi-C contact maps dataset.
#'
#' Npz file containing Hi-C contact maps of human IMR90 in 40kb resolution. The data comes from Rao et al., 2014 study and was processed by Imakaev et al., 2012 pipeline without iterative correction step (so it is raw data). It contains contact maps for chromosomes 17, 18, 19.
#'
#' @name IMR90-MboI-1_40kb-raw_maps
#'
#' @seealso \url{https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.savez_compressed.html} for npz format description, Rao et al., 2014 "A three-dimensional map of the human genome at kilobase resolution reveals prinicples of chromatin looping" for study where this dataset comes from, Imakaev et al., 2012 "Iterative correction of Hi-C data reveals hallmarks of chromosome organization." for Iterative Correction of Hi-C contact maps and \url{https://mirnylab.bitbucket.io/hiclib/index.html} for its python implementation
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525}
#'
#' @examples
#' # get file name
#' mtx.fname <- system.file("extdata", "IMR90-MboI-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
NULL

#' Csv file with sample TADs dataset.
#'
#' Contains TAD boundaries of \link{MSC-HindIII-1_40kb-raw_maps} determined using Insulation Score (Crane et al. 2015 "Condensin-driven remodelling of X chromosome topology during dosage compensation") with parameters of window size 1Mbp and delta window size 200 Kbp.
#'
#' @name MSC-HindIII-1_40kb-raw_tads
#'
#' @examples
#' # get file name
#' tads.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.tadIS", package = "CopulaHiC", mustWork = TRUE)
NULL

#' Sample Hi-C contact maps dataset.
#'
#' The list contains 4 entries - Hi-C datasets: IMR90-MboI-1_40kb-raw, IMR90-MboI-2_40kb-raw, MSC-HindIII-1_40kb-raw, MSC-HindIII-2_40kb-raw in 40kb resolution. This data was taken from 2 studies: Rao et al. 2014 "A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping." (IMR90 datasetss) and Dixon 2015 et al. "Chromatin architecture reorganization during stem cell differentiation" (MSC datasets). Data was processed using ICE pipeline \url{https://mirnylab.bitbucket.io/hiclib/index.html} without iterative correction. Each entry in list is a list with Hi-C maps in sparse format (data frames). The data is truncated - it only contains chromosomes 18 and 19.
#'
#' @docType data
#'
#' @usage data(sample_hic_maps)
#'
#' @format list of length 4 with Hi-C contact maps of 2 cell lines, each in 2 replicates. Each entry of list is a 2 element list with data frames. Every data frame is Hi-C contact map of single chromosome in sparse format:
#' \describe{
#'   \item{i}{row (y) coordinate}
#'   \item{j}{column (x) coordinate}
#'   \item{val}{number of contact between i-th and j-th region}
#'   \item{diagonal}{distance between i and j, i.e. abs(i-j) and diagonal of Hi-C contact map which this cell belongs}
#'   \item{name}{Hi-C contact map name, chromosome in this case}
#' }
"sample_hic_maps"

#' Sample TAD dataset.
#'
#' The list contains a set of TADs determined on \code{\link{sampleHiCmaps}} dataset using Insulation Score with window size 1Mbp and delta window size 200Kbp
#'
#' @docType data
#'
#' @usage data(sample_tads)
#'
#' @format list of length 4 with TADs of 2 cell lines, each in 2 replicates. Each entry in list is a data frame with TADs of chromsome 18 and 19:
#' \describe{
#'   \item{start}{start bin of a TAD (first TAD starts at 0)}
#'   \item{end}{end bin of a TAD (for consecutive TADs it is equal to next TAD start bin)}
#'   \item{name}{name of Hi-C contact map on which this TADs were determined - chromosome in this case}
#' }
"sample_tads"
