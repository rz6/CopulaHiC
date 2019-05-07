# DIADEM #
------------------------------

Authors: Rafal Zaborowski, Bartek Wilczynski

Institution: Faculty of Mathematics, Informatics and Mechanics, University of Warsaw, Poland

License: MIT + file LICENSE

For more information please contact r.zaborowski@mimuw.edu.pl or bartek@mimuw.edu.pl

## Quick summary ##
-------------------

DIADEM is R package for differential analysis of Hi-C data. It takes advantage of significant correlations of main diagonals between different Hi-C data sets (cell lines, experiments, etc.). The number of diagonals (maximum genomic distance between interacting regions) depends on chromosome and data quality but usually will equal to about 5% of total number of bins in given chromosome contact map. DIADEM uses GLM to model relationship between corresponding cells of a pair of Hi-C datasets at given genomic distance and then quantifies deviatons from the model in probabilistic way. The only required input are raw Hi-C contact map files in numpy [npz](https://kite.com/python/docs/numpy.lib.npyio.NpzFile) format.

For more details, examples and quick start refer to vignette (DIADEM_vignette.html file). You can also browse documentation of individual functions or objects within the package using standard R syntax (i.e.: `help(foo)` or `?foo`) or have a look at reference manual (DIADEM.pdf file).

The indepth description of our model together with detailed analysis and motivation will be described in manuscript (in preparation - available soon).

NOTE: DIADEM_vignette.html file can't be properly rendered from within github. Therefore it is recommended to download the file (or just clone the repo) and then open it.

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
    devtools::install_github("rz6/DIADEM")
    ```

2. from source: clone repository - by default to directory: diadem, cd to directory containing cloned repo, open R and run:
 
    ```r
    devtools::install("diadem")
    ```
    
## Usage ##
-----------

Import DIADEM package and list functions inside it:

```r
library("DIADEM")
getNamespaceExports("DIADEM")
```

#### Sample data ####

DIADEM package contains sample data of Hi-C contact maps and TADs. Both of them are in 2 formats:  

1. Numpy npz file (Hi-C maps) and csv file (TADs) - to access it run:

    * Hi-C contact maps npz files

        ```r
        # file name of MSC-HindIII-1_40kb-raw_maps
        mtx.fname.msc <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
        # file name of IMR90-MboI-1_40kb-raw_maps
        mtx.fname.imr90 <- system.file("extdata", "IMR90-MboI-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
        # load data
        maps.msc <- read_npz(mtx.fname.msc)
        maps.imr90 <- read_npz(mtx.fname.msc)
        ```

    * TADs csv files

        ```r
        # file name of MSC-HindIII-1_40kb-raw_tads
        tads.fname.msc <- system.file("extdata", "MSC-HindIII-1_40kb-raw.tadIS", package = "DIADEM", mustWork = TRUE)
        # load data
        tads.msc <- read.csv(tads.fname.msc)
        ```

2. list with Hi-C contact maps and TADs available in package environment:

    * Hi-C contact maps

        ```r
        # list available Hi-C maps
        print(names(DIADEM::sample_hic_maps))
        # get one
        map.imr90 <- DIADEM::sample_hic_maps[["IMR90-MboI-1_40kb-raw"]]
        print(names(map.imr90))
        print(head(map.imr90[["18"]]))
        ```

    * TADs

        ```r
        # list available TADs
        print(names(DIADEM::sample_tads))
        # get one
        tads.imr90 <- DIADEM::sample_tads[["IMR90-MboI-1_40kb-raw"]]
        print(head(tads.imr90))
        ```

