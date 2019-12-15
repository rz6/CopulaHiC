# DIADEM #
------------------------------

Authors: Rafal Zaborowski, Bartek Wilczynski

Institution: Faculty of Mathematics, Informatics and Mechanics, University of Warsaw, Poland

License: MIT + file LICENSE

For more information please contact r.zaborowski@mimuw.edu.pl or bartek@mimuw.edu.pl

## Quick summary ##
-------------------

DIADEM is R package for differential analysis of Hi-C data. It takes advantage of significant correlations of main diagonals between different Hi-C data sets (cell lines, experiments, etc.). The number of diagonals (maximum genomic distance between interacting regions) depends on chromosome and data quality but usually will equal to about 5% of total number of bins in given chromosome contact map. DIADEM uses GLM to model relationship between corresponding cells of a pair of Hi-C datasets at given genomic distance and then quantifies deviatons from the model in probabilistic way. The only required input are raw Hi-C contact map files in numpy [npz](https://kite.com/python/docs/numpy.lib.npyio.NpzFile) format.

For more details, examples and quick start refer to vignette (invoke `vignette(package="DIADEM")`). You can also browse documentation of individual functions or objects within the package using standard R syntax (i.e.: `help(foo)` or `?foo`) or have a look at reference manual (invoke `devtools::build_manual()`).

The indepth description of our model together with detailed analysis and motivation is described in manuscript available at: <https://www.biorxiv.org/content/10.1101/654699v1>.

## Prerequisites ##
-------------------

The code is written in R, but data storage is done with numpy, so main requirements are (versions for which tests were performed are given in parenthesis):

*  R (3.4.4)
*  python (2.7.12)
*  numpy (1.12.1)

Additionally following R packages are required:

*  magrittr
*  reticulate
*  Matrix
*  fields
*  parallel
*  igraph
*  raster
*  latex2exp
*  intervals
*  robustreg
*  MASS

Following additional packages are required for examples and plotting:

*  ggplot2
*  reshape2
*  gridExtra

NOTE: Some of the above R packages require GSL (GNU Scientific Library). Before installation make sure that libgsl-dev is installed (`sudo apt-get install libgsl-dev` on Ubuntu).

## Installation ##
-------------------

Two ways of installation are possible (both require R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html) to be installed):

1. from github repository:

    ```r
    devtools::install_github("rz6/DIADEM", build_manual = FALSE, build_vignettes = FALSE)
    ```

2. from source: clone (:warning: NOTE: it must be cloned with --recursive flag, i.e.: `git clone --recursive https://github.com/rz6/DIADEM.git`) repository - by default to directory: diadem, cd to directory containing cloned repo, open R and run:
 
    ```r
    devtools::install("diadem", build_vignettes = TRUE)
    ```
    
:warning: This repository contains submodule, which must be cloned as well for package to compile. Therefore this repository MUST be cloned with --recursive flag.

## Usage ##
-----------

Import DIADEM package and list functions inside it:

```r
library("DIADEM")
getNamespaceExports("DIADEM")
```

A good introduction with some examples and more precise description may be found in vignette. To print it call:
```{r}
vignette(package="DIADEM")
```

#### Sample data ####

DIADEM package contains sample Hi-C contact map as R built-in dataset.

1. It can be accessed as shown below:

    * Hi-C contact maps in sparse format

        ```r
        library("DIADEM")
        # file name of MSC-HindIII-1 (also IMR90-MboI-1 dataset is available)
        data(sample_hic_data, package = "DIADEM")
        msc.df <- sample_hic_data[["MSC-HindIII-1"]]
        # in order to convert contact map to dense format and save in npz file prepare temporary file
        mtx.fname.msc <- file.path(tempdir(), "MSC-HindIII-1_40kb-raw.npz")
        # get chromosome sizes
        chr.sizes <- sample_hic_data[["chromosome.sizes"]]
        # convert to dense matrix format
        l <- lapply(names(msc.df), function(chromosome) sparse2dense(dat2[[chromosome]], N = chr.sizes[[chromosome]]))
        names(l) <- names(msc.df)
        # save to npz file
        save_npz(l, mtx.fname.msc)
        ```

    * Reading Hi-C matrices from npz file

        ```r
        # given file name from previous example one can read matrices in npz format as follows
        sparse.msc <- read_npz(mtx.fname.msc)
        # or in dense format
        dense.msc <- read_npz(mtx.fname.msc, sparse.format = FALSE)
        ```

