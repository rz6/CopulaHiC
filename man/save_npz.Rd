% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{save_npz}
\alias{save_npz}
\title{Saves list with dense matrices to npz compressed file.}
\usage{
save_npz(mtx.list, path)
}
\arguments{
\item{mtx.list}{list containg dense matrices, list should have names}

\item{path}{character string specifying path together with filename to save matrices}
}
\value{
None
}
\description{
Saves list with dense matrices to npz compressed file.
}
\examples{
# get sample npz file name
mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
mtx.dense.sublist <- read_dense(mtx.fname, mtx.names = c("18","19"))
print(names(mtx.dense.sublist))
print(typeof(mtx.dense.sublist))
print(typeof(mtx.dense.sublist[["18"]]))
print(dim(mtx.dense.sublist[["18"]]))
# save
save_npz(mtx.dense.sublist, "example.npz")

}
\seealso{
\url{https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.savez_compressed.html} for python numpy method used to save matrices in npz compressed file
}
