# CopulaHiC #
------------------------------

Authors: Rafal Zaborowski, Bartek Wilczynski

Institution: Faculty of Mathematics, Informatics and Mechanics, University of Warsaw, Poland

License: MIT + file LICENSE

For more information please contact r.zaborowski@mimuw.edu.pl or bartek@mimuw.edu.pl

## Quick summary ##
-------------------

CopulaHiC is R package for differential analysis of Hi-C data. It takes advantage of significant correlations of main diagonals between different Hi-C data sets (cell lines, experiments, etc.) - usually first 200 to 300 diagonals from main diagonal are considered. CopulaHiC uses copulas to model these dependancies and then quantifies deviatons from the model in a probabilistic way. The only required input are raw Hi-C contact maps files in numpy [npz](https://kite.com/python/docs/numpy.lib.npyio.NpzFile) format.

For more details, examples and quick start refer to vignette (CopulaHiC_vignette.html file). You can also browse documentation of individual functions or objects within the package using standard R syntax (i.e.: `help(foo)` or `?foo`) or have a look at reference manual (CopulaHiC.pdf file).

The indepth description of our model together with detailed analysis and motivation will be described in manuscript (in preparation - available soon).

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
*  fitdistrplus
*  VineCopula
*  parallel
*  igraph
*  raster
*  latex2exp
*  intervals

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
    devtools::install_github("rz6/CopulaHiC")
    ```

2. from source: clone repository - by default to directory: copulahic, cd to directory containing cloned repo, open R and run:
 
    ```r
    devtools::install("copulahic")
    ```
    
## Usage ##
-----------

Import CopulaHiC package and list functions inside it:

```r
library("CopulaHiC")
getNamespaceExports("CopulaHiC")
```

#### Sample data ####

CopulaHiC package contains sample data of Hi-C contact maps and TADs. Both of them are in 2 formats:  

1. Numpy npz file (Hi-C maps) and csv file (TADs) - to access it run:

    * Hi-C contact maps npz files

        ```r
        # file name of MSC-HindIII-1_40kb-raw_maps
        mtx.fname.msc <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
        # file name of IMR90-MboI-1_40kb-raw_maps
        mtx.fname.imr90 <- system.file("extdata", "IMR90-MboI-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
        # load data
        maps.msc <- read_npz(mtx.fname.msc)
        maps.imr90 <- read_npz(mtx.fname.msc)
        ```

    * TADs csv files

        ```r
        # file name of MSC-HindIII-1_40kb-raw_tads
        tads.fname.msc <- system.file("extdata", "MSC-HindIII-1_40kb-raw.tadIS", package = "CopulaHiC", mustWork = TRUE)
        # load data
        tads.msc <- read.csv(tads.fname.msc)
        ```

2. list with Hi-C contact maps and TADs available in package environment:

    * Hi-C contact maps

        ```r
        # list available Hi-C maps
        print(names(CopulaHiC::sample_hic_maps))
        # get one
        map.imr90 <- CopulaHiC::sample_hic_maps[["IMR90-MboI-1_40kb-raw"]]
        print(names(map.imr90))
        print(head(map.imr90[["18"]]))
        ```

    * TADs

        ```r
        # list available TADs
        print(names(CopulaHiC::sample_tads))
        # get one
        tads.imr90 <- CopulaHiC::sample_tads[["IMR90-MboI-1_40kb-raw"]]
        print(head(tads.imr90))
        ```

