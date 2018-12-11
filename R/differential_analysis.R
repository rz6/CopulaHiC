#' An S3 class to represent object for Hi-C maps comparisons.
#'
#' HiC comparator object stores Hi-C contact maps from 2 experiments and (optionally) TADs
#' and allows for convenient access to contact matrices, A/B compartments or TADs. HiCcomparator
#' is constructed from npz files containing Hi-C maps in python dict with numpy matrices.
#' Additionally TAD set may be given to HiCcomparator (as list of data frames, where data frames
#' names match those of Hi-C matrices names). One can also choose to determine TADs based on
#' given Hi-C contact maps - only first, only second or determine both and take intersecting
#' intervals between them.
#'
#' @param path1 character - path to npz file containing first set of Hi-C maps
#' @param path2 character - path to npz file containing second set of Hi-C maps
#' @param tads list (optional), set of TADs as named list of data frames, each with at least start, end columns
#' @param mtx.names character vector with subset of Hi-C maps names to be selected for analysis, by default all matrices are used
#' @param which.tads numeric indicating what to do if no TADs are specified: 1 - determine TADs from first set of Hi-C maps, 2 - determine TADs from second set of Hi-C maps, 3 - determine from both sets and then take their intersection, 4 - do not determine TADs
#' @param do.pca logical whether to perform PCA for given maps and determine A/B compartments
#'
#' @return S3 object of class HiCcomparator
#'
#' @seealso \code{\link{read_npz}} for reading npz files, \code{\link{do_pca}} on how A/B compartments are determined, \code{\link{map2tads}} how TADs are determined
#'
#' @examples
#' # get path of first sample maps
#' mtx1.fname <- system.file("extdata", "IMR90-MboI-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
#' # get path of second sample maps
#' mtx2.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "CopulaHiC", mustWork = TRUE)
#' # get sample TADs
#' tads <- CopulaHiC::sample_tads[c("IMR90-MboI-1_40kb-raw", "MSC-HindIII-1_40kb-raw")]
#' # construct HiCcomparator object for chromosomes 18 and 19
#' hic.comparator <- HiCcomparator(mtx1.fname, mtx2.fname, tads, mtx.names = c("18","19"))
#' # plot A/B compartments for first and second map in chromosome 19
#' plot_pc_vector(hic.comparator$pc1.maps1[["19"]]) # first map
#' plot_pc_vector(hic.comparator$pc1.maps2[["19"]]) # second map
#'
#' @export
HiCcomparator <- function(path1, path2, tads = NULL, mtx.names = "all", which.tads = 4, do.pca = FALSE){
  # read HiC contact maps
  maps1 <- read_npz(path1, mtx.names = mtx.names)
  maps2 <- read_npz(path2, mtx.names = mtx.names)
  if(any(sapply(maps1, function(df){ any(df$val %% 1 != 0) }))){
    warning(paste0("HiC map specified as ", path1, " have non integer elements! This may influence comparative analysis!"))
  }
  if(any(sapply(maps2, function(df){ any(df$val %% 1 != 0) }))){
    warning(paste0("HiC map specified as ", path2, " have non integer elements! This may influence comparative analysis!"))
  }
  # get datasets names
  data.names <- c(sub('\\..*$','', basename(path1)), sub('\\..*$','', basename(path2)))
  # get sizes of each chromosome
  sizes1 <- read_size(path1, mtx.names = mtx.names) %>% split(rownames(.))
  sizes2 <- read_size(path2, mtx.names = mtx.names) %>% split(rownames(.))
  names.sizes <- intersect(names(sizes1), names(sizes2))
  sizes <- lapply(names.sizes, function(x){
    rbind(
      sizes1 %>% magrittr::extract2(x) %>% magrittr::set_rownames(data.names[1]),
      sizes2 %>% magrittr::extract2(x) %>% magrittr::set_rownames(data.names[2])
    )
  }) %>% magrittr::set_names(names.sizes)
  # do PCA and get first principal component to annotate cells with A/B compartments
  if(do.pca){
    names(maps1) %>%
      lapply(function(x){
        sparse2dense(maps1[[x]], N = sizes[[x]][data.names[1], "n.rows"]) %>%
          do_pca() %>% magrittr::use_series("x") %>% magrittr::extract(,"PC1")
      }) %>% magrittr::set_names(names(maps1)) -> maps1.pc1
    names(maps2) %>%
      lapply(function(x){
        sparse2dense(maps2[[x]], N = sizes[[x]][data.names[2], "n.rows"]) %>%
          do_pca() %>% magrittr::use_series("x") %>% magrittr::extract(,"PC1")
      }) %>% magrittr::set_names(names(maps2)) -> maps2.pc1
  } else {
    maps1.pc1 <- NULL
    maps2.pc1 <- NULL
  }
  # get TADs
  if(is.null(tads)){
    if(which.tads == 1){
      domains <- lapply(names(maps1), function(x){
        sparse2dense(maps1[[x]], N = sizes[[x]][data.names[1], "n.rows"]) %>%
          map2tads()
      })
    } else if(which.tads == 2){
      domains <- lapply(names(maps2), function(x){
        sparse2dense(maps2[[x]], N = sizes[[x]][data.names[2], "n.rows"]) %>%
          map2tads()
      })
    } else if(which.tads == 3){
      # get intersection and then
      domains <- lapply(intersect(names(maps1), names(maps2)), function(x){
        sparse2dense(maps1[[x]], N = sizes[[x]][data.names[1], "n.rows"]) %>%
          map2tads() -> t1
        sparse2dense(maps2[[x]], N = sizes[[x]][data.names[2], "n.rows"]) %>%
          map2tads() -> t2
        intersect_tads(t1, t2)
      })
    } else {
      # leave TADs as NULL which means use DP in detection of differential regions
      domains <- NULL
    }
  } else {
    if(all(mtx.names == "all")){
      domains <- tads
    } else {
      domains <- tads[mtx.names]
    }
  }

  structure(
    list(maps1, maps2, maps1.pc1, maps2.pc1, domains, data.names, sizes) %>%
      magrittr::set_names(c("maps1", "maps2", "pc1.maps1", "pc1.maps2", "tads", "data.names", "maps.dims")),
    class = "HiCcomparator"
  )
}

#' @export
merge <- function(x, ...) UseMethod("merge", x)

#' Finds regions, which interacts in both experiments (maps).
#'
#' Merges Hi-C contact maps data frames 1 and 2 of HiCcomparator object by i, j, diagonal, name columns.
#'
#' @param x HiCcomparator object
#' @param include.zero.cells logical, whether to include cells, which have non zero number of contacts in one map, but not the other, not recommended
#'
#' @return list of data frames with merged contact maps data in sparse format
#'
#' @seealso \code{\link{HiCcomparator}} on how to construct HiCcomparator object
#'
#' @examples
#' # first create HiCcomparator object - see ?HiCcomparator for examples
#' merged <- merge(hic.comparator)
#'
#' @export
merge.HiCcomparator <- function(x, include.zero.cells = FALSE){
  maps.names <- intersect(names(x$maps1), names(x$maps2))
  lapply(maps.names, function(y){
    map1 <- x$maps1[[y]]
    map2 <- x$maps2[[y]]
    if(include.zero.cells){
      merged <- base::merge(map1, map2, by = c("i","j","diagonal","name"), all = TRUE)
      merged[is.na(merged$val.x), "val.x"] <- 0
      merged[is.na(merged$val.y), "val.y"] <- 0
    } else {
      merged <- base::merge(map1, map2, by = c("i","j","diagonal","name"))
    }
    if(is.null(x$pc1.maps1) | is.null(x$pc1.maps2)){
      return(merged)
    }
    cn <- c("compartment.i","compartment.j","compartment")
    cbind(
      merged,
      pc2mtx(x$pc1.maps1[[y]], merged) %>%
        magrittr::extract(cn) %>% magrittr::set_colnames(paste(cn,"x",sep = ".")),
      pc2mtx(x$pc1.maps2[[y]], merged) %>%
        magrittr::extract(cn) %>% magrittr::set_colnames(paste(cn,"y",sep = "."))
    ) %>% magrittr::inset(c("one.compartment.x", "one.compartment.y"), value = "+") -> merged
    merged[merged$compartment.x == "AB", "one.compartment.x"] <- "-"
    merged[merged$compartment.y == "AB", "one.compartment.y"] <- "-"
    return(merged)
  }) %>% magrittr::set_names(maps.names)
}

#' @export
dominating_signal <- function(x, ...) UseMethod("dominating_signal", x)

#' Calculates coverage or decay of Hi-C maps.
#'
#' Computes coverages or decays of every Hi-C maps in both data sets of given HiCcomparator object. Coverage is defined as sum of contacts on given bin. Decay is sum or mean of contacts for every diagonal.
#'
#' @param hic.comparator object of type HiCcomparator
#'
#' @return dataframe with following columns: (i, sum.contacts, mean.contacts, sd.contacts, name, dataset), which can be used to conveniently visualise coverages or decays (see examples)
#'
#' @seealso \code{\link{HiCcomparator}} on how to construct HiCcomparator object
#'
#' @examples
#' # first create HiCcomparator object - see ?HiCcomparator for examples
#' coverage <- dominating_signal(hic.comparator)
#' # visualise results
#' library("ggplot2")
#' ggplot(coverage, aes(x = i, y = sum.contacts, color = dataset)) +
#' geom_point(size = 0.5) +
#'  geom_smooth(alpha = 0.5) +
#'  facet_wrap(~ name, ncol = 1, scales = "free") +
#'  theme(legend.position = "bottom")
#' # get decay
#' decay <- dominating_signal(hic.comparator, which.signal = "decay")
#' ggplot(decay[decay$diagonal != 0,], aes(x = diagonal, y = mean.contacts, color = dataset)) +
#' geom_point(size = 0.5) +
#'  scale_x_log10() +
#'  scale_y_log10() +
#'  facet_wrap(~ name, ncol = 1, scales = "free") +
#'  theme(legend.position = "bottom")
#'
#' @export
dominating_signal.HiCcomparator <- function(hic.comparator, which.signal = c("coverage","decay")[1]){
  stopifnot(which.signal %in% c("coverage","decay"))
  data.names <- hic.comparator$data.names
  variable <- if(which.signal == "coverage") "i" else "diagonal"
  f <- formula(paste0(". ~ ",variable))
  lapply(names(hic.comparator$maps1), function(n){
    agg <- aggregate(f, data = hic.comparator$maps1[[n]][c(variable,"val")], FUN = function(x){ c(sum.contacts = sum(x), mean.contacts = mean(x), sd.contacts = sd(x)) }) %>%
      magrittr::inset("name", value = n)
  }) %>% do.call("rbind",.) %>% magrittr::inset("dataset", value = data.names[1]) -> signal1
  lapply(names(hic.comparator$maps2), function(n){
    agg <- aggregate(f, data = hic.comparator$maps2[[n]][c(variable,"val")], FUN = function(x){ c(sum.contacts = sum(x), mean.contacts = mean(x), sd.contacts = sd(x)) }) %>%
      magrittr::inset("name", value = n)
  }) %>% do.call("rbind",.) %>% magrittr::inset("dataset", value = data.names[2]) -> signal2
  signal <- rbind(signal1, signal2)
  cbind(signal[,-2], signal[,2])
}

#' @export
decay_correlation <- function(x, ...) UseMethod("decay_correlation", x)

#' Calculates correlations between diagonals.
#'
#' Computes correlations (Pearson, Spearman, Kendall) and significances of corresponding diagonals between 2 Hi-C maps of HiCcomparator object.
#'
#' @param hic.comparator object of type HiCcomparator
#'
#' @return dataframe with following columns: diagonal, pcc, pearson.pval, rho, spearman.pval, tau, kendall.pval, name which can be used to conveniently visualise dependancy between 2 Hi-C maps being compared (see examples)
#'
#' @seealso \code{\link{HiCcomparator}} on how to construct HiCcomparator object
#'
#' @examples
#' first create HiCcomparator object - see ?HiCcomparator for examples
#' library("ggplot2")
#' library("reshape2")
#' decay.cors <- decay_correlation(hic.comparator)
#' # wide to long
#' decay.cors.long <- reshape2::melt(decay.cors[c("name","diagonal","pcc","rho","tau")], id.vars = c("name","diagonal"), variable.name = "correlation", value.name = "coefficient")
#' # remove 0 diagonal (as it is non informative anyways) and illustrate results
#' ggplot(decay.cors.long[decay.cors.long$diagonal != 0,],
#'   aes(x = diagonal, y = coefficient, color = correlation)) +
#'  geom_point(size = 0.3) +
#'  facet_wrap(~ name, ncol = 1, scales = "free") +
#'  theme(legend.position = "bottom")
#'
#' @export
decay_correlation.HiCcomparator <- function(hic.comparator){
  merged <- merge(hic.comparator)
  lapply(names(merged), function(n){
    merged.n <- merged[[n]]
    merged.n.by.diagonal <- split(merged.n, merged.n$diagonal)
    sapply(names(merged.n.by.diagonal), function(d){
      df <- merged.n.by.diagonal[[d]]
      if(nrow(df) < 5){
        # if there are less than 5 observations there is no point in calculating correlation
        cors <- rep(NA, 6)
      } else {
        cors <- sapply(c("pearson","spearman","kendall"), function(x){
          c(
            cor(df$val.x, df$val.y, method = x),
            cor.test(df$val.x, df$val.y, method = x, alternative = "two.sided")$p.value
          )
        })
      }
      c(as.numeric(d), as.vector(cors))
    }) %>% t() %>% as.data.frame() %>%
      magrittr::set_colnames(c("diagonal","pcc","pearson.pval","rho","spearman.pval","tau","kendall.pval")) %>%
      magrittr::inset("name", value = n)
  }) %>% do.call("rbind",.)
}

#' An S3 object to represent Hi-C copula model.
#'
#' Models diagonal-wise dependencies between Hi-C data sets with copulas. Model is constructed as follows:
#' \itemize{
#'   \item{}{merge maps1 with maps2}
#'   \item{}{for each diagonal in diagonals}
#'   \itemize{
#'     \item{}{take all points from this diagonal, such that they are non zero in map1 (X) and non zero in map2 (Y)}
#'     \item{}{model X with gamma distribution}
#'     \item{}{model Y with gamma distribution}
#'     \item{}{transform X and Y to U and V ~ Uniform(0,1) - see \code{\link{VineCopula::pobs}}}
#'     \item{}{model F(U,V) with copula (see \code{\link{VineCopula::BiCopSelect}})}
#'   }
#' }
#' Before fitting the model it's recommended to first inspect correlations between analyzed Hi-C maps before fixing this variable. As the ratio of noise / signal in Hi-C data increases rapidly with decay it's unadvised to use all diagonals for modelling. The number of diagonals to be used will depend on chromosome length, resolution and data quality. As a rule of thumb the number of diagonals should not exceed 0.2 times length of chromosome.
#'
#' @param hic.comparator object of type HiCcomparator
#' @param diagonals fraction or numeric vector or character "all" which diagonals to use to fit models, by default fraction of chromsome length is used to indicate number of diagonals.
#' @param include.zero.cells logical, whether to include cells when merging maps (see \code{\link{merge.HiCcomparator}})
#' @param n.cores numeric number of cores to distribute model computations
#'
#' @return S3 object of class HiCcopula
#'
#' @seealso \code{\link{HiCcomparator}} on how to construct HiCcomparator object, \code{\link{fitdistrplus::fitdist}} on distribution fitting, \code{\link{VineCopula::pobs}} on pseudo observations generation and \code{\link{VineCopula::BiCopSelect}} on finding optimal fit bivariate copula
#'
#' @examples
#' # first create HiCcomparator object - see ?HiCcomparator for examples
#' # construct model
#' hic.copula <- HiCcopula(hic.comparator)
#' # load below packages for visualisation
#' library("fitdistrplus")
#' library("VineCopula")
#' # illustrate gamma fit of X and Y for chromosome 18, diagonal 5
#' plot(copula$model[["18"]][["5"]]$marginal.x)
#' plot(copula$model[["18"]][["5"]]$marginal.y)
#' # illustrate copula fit of U and V for chromosome 18, diagonal 5
#' plot(copula$model[["18"]][["5"]]$bf.copula)
#'
#' @export
HiCcopula <- function(hic.comparator, diagonals = 0.12, include.zero.cells = FALSE, n.cores = 1){
  merged <- merge(hic.comparator, include.zero.cells = include.zero.cells)
  names.merged <- names(merged)
  lapply(names.merged, function(n){
    merged.n <- merged[[n]]
    # split by diagonals
    merged.n.diagonals <- split(merged.n, merged.n$diagonal)
    if(is.character(diagonals)){
      diagonals <- names(merged.n.diagonals)
    } else if((length(diagonals) == 1) & (diagonals < 1)){
      m <- round(diagonals * hic.comparator$maps.dims[[n]][1,1])
      diagonals <- seq(1,m)
    }
    # calculate models for every diagonal
    parallel::mclapply(merged.n.diagonals[as.character(diagonals)], FUN = function(df){
      tryCatch({
        # fit gamma distribution for marginals
        marginal.x <- fitdistrplus::fitdist(df$val.x, distr = "gamma", method = "mle")
        marginal.y <- fitdistrplus::fitdist(df$val.y, distr = "gamma", method = "mle")
        # find best fit copula
        pseudo.observations <- VineCopula::pobs(df[,c("val.x","val.y")])
        u <- pseudo.observations[,1]
        v <- pseudo.observations[,2]
        bf.copula <- VineCopula::BiCopSelect(u,v, familyset = NA)
        # return result as  list with:
        # - marginal distribution of X
        # - marginal distribution of Y
        # - best fit copula
        list(marginal.x, marginal.y, bf.copula) %>%
          magrittr::set_names(c("marginal.x", "marginal.y", "bf.copula"))
      }, error = function(cond){
        # if any error occur return NULL for this diagonal, then one can
        # check in resulting list, for which diagonals error occured:
        # which(sapply(result, is.null))
        return(NULL)
      })
    }, mc.cores = n.cores)
  }) %>% magrittr::set_names(names.merged) -> model
  # return result
  structure(list(hic.comparator$data.names, names(merged), hic.comparator$maps.dims, merged, model) %>%
              magrittr::set_names(c("data.names", "names", "maps.dims", "merged", "model")),
            class = "HiCcopula")
}

#' Compute depletion or enrichment p-value given copula model.
#'
#' Calculates p-value of enrichment or depletion of given point/s with u,v coordinates given copula F(U,V):
#' \itemize{
#'   \item{}{depletion probability is defined as: P(U < u, V > v) = F(U < u, V < 1) - F(U < u, V < v), so its upper left rectangle of copula distribution}
#'   \item{}{enrichment probability is defined as: P(U > u, V < v) = F(U < 1, V < v) - F(U < u, V < v), so its lower right rectangle of copula distribution}
#' }
#' NOTES:
#' \itemize{
#'   \item{}{copula is symmetric with respect to U = V}
#'   \item{}{enrichment is relative to condition u = F(X < x), in other words enriched means chromatin is enriched in interactions in conditon X comparing with Y}
#'   \item{}{P(U,V) distribution can be plotted as U horizontal and increasing from left (start at 0) to right (ends at 1) and V vertical and increasing from bottom (start at 0) to top (ends at 1)}
#'   \item{}{both uv and model is for single (and the same) diagonal}
#'   \item{}{uv must all be points for either top left (depletion) or bottom right (enrichment) part of distribution}
#' }
#' For calculation of copula cdf \code{\link{VineCopula::::BiCopCDF}} is used.
#'
#' @param uv matrix of dimension n x 2 with U,V r.v. ~ Uniform(0,1), see \code{\link{VineCopula::pobs}} for generation of U,V from X,Y; here U,V matrix must be such that either U >= V (lower.right corner) or U <= V (upper.left corner)
#' @param copula.model object of class \code{\link{VineCopula::BiCop}}
#' @param copula.tail character indicating tail (corner) of copula to calculate cdf, either "upper.left" or "lower.right"
#'
#' @return numeric vector of p-value/s of length equal to \code{nrow(uv)}
#'
#' @examples
#' library("copula")
#' library("MASS")
#' # make bivariate standard normal copula of highly correlated variables (0.8)
#' cop <- BiCop(1, 0.8)
#' # convert to package copula class --> for illustration purposes
#' cop.copula.object <- copulaFromFamilyIndex(cop$family, cop$par)
#' # illustrate copula
#' persp(cop.copula.object, dCopula) # 3D
#' contourplot2(cop.copula.object, dCopula, col.regions = terrain.colors) # 2D heatmap
#' # simulate sample of size 10000 and draw copula 2d density plot (heatmap)
#' sample.cop <- data.frame(rCopula(10000, cop.copula.object))
#' plot_copula_density(sample.cop) # 2D heatmap with ggplot
#' # now simulate sample from bivariate standard normal with much lower correlation than that of copula
#' sigma <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
#' xy <- mvrnorm(500, mu = c(0, 0), Sigma = sigma, empirical = T)
#' # convert X,Y to U,V ~ Uniform(0,1) using copula::pobs
#' uv <- pobs(xy)
#' # keep only upper V >= U (\code{copula_pvals} calculates p-values separately for V >= U and U >= V)
#' uv <- uv[uv[,2] >= uv[,1],]
#' # illustrate observations on top of copula density
#' uv.df <- data.frame(uv)
#' colnames(uv.df) <- c("U","V")
#' plot_copula_density(sample.cop) + geom_point(aes(x = U, y = V), data = uv.df, size = 0.5)
#' # calculate p-values of observations given the copula model (use only )
#' uv.df$pval <- copula_pvals(uv, cop)
#' # convert to pvalues to negative log10(pval)
#' uv.df$neg.log.pval <- -log10(uv.df$pval)
#' # illustrate observations significance on top of copula model
#' plot_copula_density(sample.cop) + geom_point(aes(x = U, y = V, color = neg.log.pval), data = uv.df, size = 0.3) + scale_color_gradient(low = "yellow", high = "black", name = "-log10(pval)")
copula_pvals <- function(uv, copula.model, copula.tail = "upper.left"){
  # for some reason if copula family = 214 (Rotated Tawn type 2 180 degrees)
  # the cdf = NaN if u = 1 or v = 1, so sloppy workaround is to set them
  # as close to 1 as possible to not yield NaN, which is 0.9999999999999999
  if(copula.model$family == 214){
    one <- 0.9999999999999999
  } else {
    one <- 1.0
  }
  if(copula.tail == "upper.left"){
    # V >= U
    stopifnot(all(uv[,2] >= uv[,1]))
    u1 <- uv[,1]
    v1 <- rep(one, nrow(uv))
  } else {
    # lower right tail --> U >= V
    stopifnot(all(uv[,1] >= uv[,2]))
    u1 <- rep(one, nrow(uv))
    v1 <- uv[,2]
  }
  cdf1 <- VineCopula::BiCopCDF(u1, v1, family = copula.model$family, par = copula.model$par, par2 = copula.model$par2)
  cdf2 <- VineCopula::BiCopCDF(uv[,1], uv[,2], family = copula.model$family, par = copula.model$par, par2 = copula.model$par2)
  return(cdf1 - cdf2)
}

#' Computes p-values given copula model.
#'
#' Given copula model for diagonal of Hi-C contact map it calculates deviation from this model (p-values) in both directions (i.e. depletion and enrichment). It divides U,V matrix on 3 groups: U > V (enrichment), U == V (no effect), U < V (depletion) and calls \code{\link{copula_pvals}} for each group. Afterwards it merge results.
#'
#' @param uv matrix of dimension n x 2 with U,V r.v. ~ Uniform(0,1), see \code{\link{VineCopula::pobs}} for generation of U,V from X,Y
#' @param copula.model object of class \code{\link{VineCopula::BiCop}}
#' @param mht.correction logical whether to apply Benjamini Hochberg correction for multiple hypothesis testing since we are testing a number of hypothesis (diagonal cells) against null model
#'
#' @return data frame with columns: u, v, effect, p.value, p.value.corrected
#'
#' @seealso \code{\link{copula_pvals}} for p-value calculation
#'
#' @examples
#' # see copula_pvals function
maps_difference_diagonal <- function(uv, copula.model, mht.correction = TRUE){
  list(
    uv[uv[,1] > uv[,2],],
    uv[uv[,1] == uv[,2],],
    uv[uv[,1] < uv[,2],]
  ) %>% magrittr::set_names(c("enrichment", "none", "depletion")) %>%
    magrittr::extract(sapply(., nrow) > 0) -> l
  lapply(names(l), function(x){
    df <- l[[x]]
    if(nrow(df)){
      if(x == "enrichment"){
        copula.tail <- "lower.right"
      } else {
        copula.tail <- "upper.left"
      }
      p.vals <- copula_pvals(df, copula.model, copula.tail = copula.tail)
      cbind(df, effect = x, p.value = p.vals)
    }
  }) %>%
    magrittr::extract(! sapply(., is.null)) %>%
    do.call("rbind",.) -> result
  if(mht.correction){
    result$p.value.corrected <- p.adjust(result$p.value, method = "BH")
  }
  return(result)
}

#' @export
hicdiff <- function(x, ...) UseMethod("hicdiff", x)

#' Computes differential (p-value) map.
#'
#' Given HiCcopula object calculates depletion/enrichment p-values of cell along diagonals w.r.t. background copula models (separate for every diagonal, see \code{\link{HiCcopula} for details}).
#'
#' @param hic.copula object of class HiCcopula
#' @param marginal.distr character wither fit or obs; if fit then fitted gamma distribution for X and Y is used to convert them to U and V respectively
#'
#' @return list of data frames corresponding to Hi-C contact maps; each data frame contain columns: i, j, name, val.x, val.y, u, v, effect, p.value, p.value.corrected
#'
#' @seealso \code{\link{HiCcopula}} on how to construct HiCcopula, \code{\link{maps_difference_diagonal}} and \code{\link{copula_pvals}} on how p-values are calculated
#'
#' @examples
#' # first create HiCcopula object - see ?HiCcopula for examples
#' # calculate p-values and select chromosome of interest (18 in this example)
#' md <- hicdiff(hic.copula)[["18"]]
#' # convert corrected p-values to -log10(pvals)
#' md$neg.log.cor.pvals <- md$p.value.corrected
#' # make neg.log.cor.pvals negative for depleted cells
#' md[md$effect == "depletion", "neg.log.cor.pvals"] <- -md[md$effect == "depletion", "neg.log.cor.pvals"]
#' # convert sparse matrix to dense matrix
#' dense <- sparse2dense(md[c("i","j","neg.log.cor.pvals")], N = hic.copula$maps.dims[["18"]][1,1])
#' # plot results
#' plot_diff_map(dense)
#'
#' @export
hicdiff.HiCcopula <- function(hic.copula, marginal.distr = c("fit","obs")[1]){
  # TODO: add assertion if data (contacts) are reals (i.e.: normalized data) and marginal.distr is obs
  lapply(hic.copula$names, function(n){
    merged <- hic.copula$merged[[n]]
    merged.diagonals <- split(merged, merged$diagonal)
    model <- hic.copula$model[[n]]
    # compute p values for every diagonal in model
    lapply(names(model), function(m){
      model.diagonal <- model[[m]]
      merged.diagonal <- merged.diagonals[[m]]
      if(marginal.distr == "fit"){
        # get gamma parameters
        x.estimate <- model.diagonal$marginal.x[["estimate"]]
        y.estimate <- model.diagonal$marginal.y[["estimate"]]
        # calculate u, v
        u <- pgamma(merged.diagonal$val.x, x.estimate["shape"], rate = x.estimate["rate"])
        v <- pgamma(merged.diagonal$val.y, y.estimate["shape"], rate = y.estimate["rate"])
      } else {
        # calculate observed probability for vector of counts for x
        counts.x <- as.data.frame(table(as.factor(merged.diagonal$val.x)))
        counts.x$Freq %>% divide_by(sum(counts.x$Freq)) %>% magrittr::set_names(counts.x$Var1) -> pobs.x
        u <- pobs.x[as.character(merged.diagonal$val.x)]
        # calculate observed probability for vector of counts for y
        counts.y <- as.data.frame(table(as.factor(merged.diagonal$val.y)))
        counts.y$Freq %>% divide_by(sum(counts.y$Freq)) %>% magrittr::set_names(counts.y$Var1) -> pobs.y
        v <- pobs.y[as.character(merged.diagonal$val.y)]
      }
      uv <- data.frame(u = u, v = v, idx = rownames(merged.diagonal))
      pvals <- maps_difference_diagonal(uv, model.diagonal$bf.copula)
      base::merge(merged.diagonal[c("i","j","name","val.x","val.y")],
                  pvals, by.x = "row.names", by.y = "idx")
    }) %>% do.call("rbind",.)
  }) %>% magrittr::set_names(hic.copula$names)
}

#' Converts significance maps to dense matrices.
#'
#' Given list of significance maps in sparse format produced by \code{\link{hicdiff}} function converts them to dense matrix format.
#'
#' @param hicdiff.list list of data frames, output from \code{\link{hicdiff}} function
#' @param maps.dims list of data frames, usually an attribute of HiCcopula object that was used to produce hicdiff.list
#' @param val.column character indicating which column of sparse significance map to use as cell value for dense matrix (one of p.value or p.value.corrected)
#' @param neg.log logical wheter to apply -log10(.) to val.column
#' @param mark.depleted logical whether to mark depleted cells as negative
#'
#' @return list with dense matrices
#'
#' @seealso \code{\link{hicdiff}} for generation of hicdiff.list
#'
#' @examples
#' # first create HiCcopula object - see ?HiCcopula for examples
#' # then produce hicdiff.list
#' hicdiff.list <- hicdiff(hic.copula)
#' dense.hicdiff.list <- hicdiff2mtx(hicdiff.list)
#' # elements of dense.list can be visualized
#' plot_diff_map(dense.hicdiff.list[["18"]], breaks = 100)
#'
#' @export
hicdiff2mtx <- function(hicdiff.list, maps.dims, val.column = c("p.value.corrected","p.value")[1],
                        neg.log = TRUE, mark.depleted = TRUE){
  names(hicdiff.list) %>%
    lapply(function(n){
      df <- hicdiff.list[[n]]
      md <- maps.dims[[n]]
      subdf <- df[c("i", "j")]
      if(neg.log){
        subdf$val <- -log10(df[,val.column])
      } else {
        subdf$val <- df[,val.column]
      }
      if(mark.depleted){
        subdf[df$effect == "depletion","val"] <- -subdf[df$effect == "depletion","val"]
      }
      sparse2dense(subdf, N = md[1,"n.rows"])
    }) %>% magrittr::set_names(names(hicdiff.list))
}

#' @export
differential_interactions <- function(x, ...) UseMethod("differential_interactions", x)

#' Finds significantly interacting rectangle-like regions.
#'
#' This function works in 3 steps:
#' \itemize{
#'  \item{}{first it calculates differential map,}
#'  \item{}{then it takes negative log10 p-value vector of cells, sorts it, divides on negative and positive parts and lastly fits bilinear model to each of them}
#'  \item{}{finally it retains only those cells that are to the left (right) of intersection point of bilinear model in positive (negative) p-values vector and searches for connected components in both of them separately}
#' }
#' Fitting bilinear model is performed using \code{\link{best_fit_bilinear}} function, while for connected components search \code{\link{raster}} package is used. After detection of significantly interacting regions (connected components) one may further filter list to only retain those with number of non zero cells (n.cells column in interacting.regions data frame) larger than some threshold.
#'
#' @param hic.copula object of class HiCcopula
#' @param plot.models logical if true then plot bilinear model fit for every matrix in hic.copula object; it will plot models for depleted and enriched models separately; if you want to save this results to file open device before calling this function (see for instance \code{\link{pdf}}) and close device after function call (see \code{\link{dev.off}})
#'
#' @return list with number of entries equal to hic.copula$names; each entry is a list with 2 elements: interacting.regions - data frame containing rows with rectangle like regions of significant interactions with coordinates n.cells (number of non zero cells inside rectangle), start.x, end.x, start.y, end.y, effect; connected.components list with cells comprising given connected component; connected components list is named list where each entry name is unique id, which can be mapped to row in interacting.regions (its row names)
#'
#' @seealso \code{\link{best_fit_bilinear}} for fitting bilinear model, \code{\link{raster::raster}} and \code{\link{raster::clump}} for connected components search
#'
#' @examples
#' # first create HiCcopula object - see ?HiCcopula for examples
#' di <- differential_interactions(hic.copula)
#' di18 <- di[["18"]]
#' # if you want to plot results create pvalue map (see hicdiff function) and from that dense map for some chromosome
#' plot_diff_map(dense)
#' # plot regions having at least 10 significant cells in connected component
#' plot_regions(di18[di18$n.cells >= 10,2:6], pal.colors = c("blue","red"))
#'
#' @export
differential_interactions.HiCcopula <- function(hic.copula, plot.models = TRUE){
  maps.difference <- hicdiff(hic.copula)
  lapply(names(maps.difference), function(n){
    df <- maps.difference[[n]]
    df$neg.log.cor.pval <- -log10(df$p.value.corrected)
    df[df$effect == "depletion","neg.log.cor.pval"] <- -df[df$effect == "depletion","neg.log.cor.pval"]
    df <- df[order(df$neg.log.cor.pval),]
    df$X <- seq(1, nrow(df))
    list(df[df$neg.log.cor.pval < 0,], df[df$neg.log.cor.pval > 0,]) %>%
      lapply(function(sgn.df){
        bf <- best_fit_bilinear(sgn.df$X, sgn.df$neg.log.cor.pval)
        x <- bf[["intersection.x"]]["both"]
        # get left and right tails of p-values
        v <- abs(c(sgn.df$neg.log.cor.pval[1], tail(sgn.df$neg.log.cor.pval, 1)))
        if(v[1] > v[2]){
          most.significant <- sgn.df[sgn.df$X < x,]
        } else {
          most.significant <- sgn.df[sgn.df$X > x,]
        }
        if(plot.models){
          plot(sgn.df$X, sgn.df$neg.log.cor.pval, cex = 0.1, xlab = "X", ylab = latex2exp::TeX("-log_{10}(p_{BH})"))
          title(paste0(hic.copula$data.names[1]," vs ",hic.copula$data.names[2],", contact map ",n), cex.main = 1)
          abline(a = bf[["coefficients"]]["left","intercept"], b = bf[["coefficients"]]["left","slope"], col = "green")
          abline(a = bf[["coefficients"]]["right","intercept"], b = bf[["coefficients"]]["right","slope"], col = "blue")
          abline(v = x, col = "red")
        }
        adjacency <- sparse2dense(most.significant[c("i","j","neg.log.cor.pval")],
                                  N = hic.copula$maps.dims[[n]][1,1])
        # convert negative or positive matrix of negative log p-values into adjacency matrix
        adjacency[adjacency != 0] <- 1
        rmat <- raster::raster(adjacency)
        clumps <- matrix(raster::clump(rmat, directions=4),
                         nrow = nrow(adjacency), ncol = ncol(adjacency), byrow = TRUE)
        # get connected components as list
        lapply(seq(1, max(clumps, na.rm=TRUE)), function(m){
          which(clumps == m, arr.ind = TRUE)
        })
      }) -> l
    ids <- seq(1, length(l[[1]]) + length(l[[2]]))
    ids.neg <- ids[1:length(l[[1]])]
    ids.pos <- ids[(length(l[[1]]) + 1):length(ids)]
    # get data frame where each row is interaction with n.cells, start.x, end.x, start.y, end.y, effect and id as row name
    res.df <- rbind(
      sapply(l[[1]], function(cc){
        c(nrow(cc), min(cc[,"col"]), max(cc[,"col"]), min(cc[,"row"]), max(cc[,"row"]))
      }) %>% t() %>% as.data.frame() %>%
        magrittr::set_colnames(c("n.cells","start.x","end.x","start.y","end.y")) %>%
        magrittr::inset("effect", value = "depletion") %>% magrittr::set_rownames(ids.neg),
      sapply(l[[2]], function(cc){
        c(nrow(cc), min(cc[,"col"]), max(cc[,"col"]), min(cc[,"row"]), max(cc[,"row"]))
      }) %>% t() %>% as.data.frame() %>%
        magrittr::set_colnames(c("n.cells","start.x","end.x","start.y","end.y")) %>%
        magrittr::inset("effect", value = "enrichment") %>% magrittr::set_rownames(ids.pos)
    )
    # get list with connected components
    res.cc <- c(magrittr::set_names(l[[1]], ids.neg), magrittr::set_names(l[[2]], ids.pos))
    list(res.df, res.cc) %>%
      magrittr::set_names(c("interacting.regions", "connected.components"))
  }) %>% magrittr::set_names(names(maps.difference))
}
