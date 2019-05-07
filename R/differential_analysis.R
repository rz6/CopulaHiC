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
#' mtx1.fname <- system.file("extdata", "IMR90-MboI-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
#' # get path of second sample maps
#' mtx2.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
#' # get sample TADs
#' tads <- DIADEM::sample_tads[c("IMR90-MboI-1_40kb-raw", "MSC-HindIII-1_40kb-raw")]
#' # construct HiCcomparator object for chromosomes 18 and 19
#' hic.comparator <- HiCcomparator(mtx1.fname, mtx2.fname, tads, mtx.names = c("18","19"))
#' # plot A/B compartments for first and second map in chromosome 19
#' plot_pc_vector(hic.comparator$pc1.maps1[["19"]]) # first map
#' plot_pc_vector(hic.comparator$pc1.maps2[["19"]]) # second map
#'
#' @export
HiCcomparator <- function(path1, path2, tads = NULL, mtx.names = "all", which.tads = 4,
                          do.pca = FALSE, zscore = FALSE, props = FALSE){
  # read HiC contact maps
  maps1 <- read_npz(path1, mtx.names = mtx.names)
  maps2 <- read_npz(path2, mtx.names = mtx.names)
  if(any(sapply(maps1, function(df){ any(df$val %% 1 != 0) }))){
    warning(paste0("HiC map specified as ", path1, " have non integer elements! This may influence comparative analysis!"))
  }
  if(any(sapply(maps2, function(df){ any(df$val %% 1 != 0) }))){
    warning(paste0("HiC map specified as ", path2, " have non integer elements! This may influence comparative analysis!"))
  }
  # calculate number of interactions per diagonal in map 1
  N1 <- lapply(names(maps1), function(x){
    sapply(split(maps1[[x]], maps1[[x]]$diagonal), function(df) sum(df$val))
  }) %>% magrittr::set_names(names(maps1))
  # calculate number of interactions per diagonal in map 2
  N2 <- lapply(names(maps2), function(x){
    sapply(split(maps2[[x]], maps2[[x]]$diagonal), function(df) sum(df$val))
  }) %>% magrittr::set_names(names(maps2))
  if(props){
    lapply(names(maps1), function(x) proportions(maps1[[x]])) %>%
      magrittr::set_names(names(maps1)) -> maps1
    lapply(names(maps2), function(x) proportions(maps2[[x]])) %>%
      magrittr::set_names(names(maps2)) -> maps2
  }
  if(zscore){
    # zscore every diagonal
    lapply(names(maps1), function(x){
      interactions2zscores(maps1[[x]])
    }) %>% magrittr::set_names(names(maps1)) -> maps1
    lapply(names(maps2), function(x){
      interactions2zscores(maps2[[x]])
    }) %>% magrittr::set_names(names(maps2)) -> maps2
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
    list(maps1, maps2, maps1.pc1, maps2.pc1, domains, data.names, sizes, N1, N2) %>%
      magrittr::set_names(c("maps1", "maps2", "pc1.maps1", "pc1.maps2", "tads", "data.names", "maps.dims", "N1", "N2")),
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

#' Constructs GLM to model per-diagonal Hi-C contact dependancy
#'
#' Models Hi-C contacts using Negative Binomial (or Poisson when the data is underdispersed) regression. Given the fact that Hi-C data suffers from contact decay bias this method is intendent to model each diagonal separately.
#'
#' @param df data frame with predictor, response, outlier columns
#' @param link character link to be used in GLM
#' @param try.poisson logical if true then use Poisson regression instead of NB anytime warning is generated when fitting NB model indicating MLE convergence issue
#'
#' @return object of class glm or MASS::glm.nb
#'
#' @seealso \code{\link{glm}}, \code{\link{MASS::glm.nb}} to see how GLM are constructed
#'
#' @export
constructGLM <- function(df, link, try.poisson = TRUE){
  stopifnot(link %in% c("log", "sqrt"))
  tryCatch({
    if(link == "log"){
      model <- MASS::glm.nb(response ~ predictor, data = df[df$outlier == "no",], link = "log")
    } else {
      model <- MASS::glm.nb(response ~ predictor, data = df[df$outlier == "no",], link = "sqrt")
    }
    if(try.poisson & !is.null(model$th.warn)){
      if(link == "log"){
        model <- glm(response ~ predictor, data = df[df$outlier == "no",], family = poisson(link = "log"))
      } else {
        model <- glm(response ~ predictor, data = df[df$outlier == "no",], family = poisson(link = "sqrt"))
      }
    }
  }, error = function(e){
    model <- NULL
  })
  return(model)
}

#' An S3 object to represent differential GLM Hi-C model.
#'
#' Models diagonal-wise dependencies between Hi-C data sets with GLM. Model is constructed as follows:
#' \itemize{
#'   \item{}{merge maps1 with maps2}
#'   \item{}{for each diagonal in diagonals}
#'   \itemize{
#'     \item{}{take all points from this diagonal, such that they are non zero in map1 (X) and non zero in map2 (Y)}
#'     \item{}{remove outliers using robust regression (recommended)}
#'     \item{}{apply transformation to predictor variable (X or Y): log or sqrt depending on link function}
#'     \item{}{model response (Y = f(X) or X = f(Y)) using Negative Binomial distribution with appropriate link function}
#'   }
#' }
#' Before fitting the model it's recommended to first inspect correlations between analyzed Hi-C maps before fixing this variable. As the ratio of noise / signal in Hi-C data increases rapidly with decay it's unadvised to use all diagonals for modelling. The number of diagonals to be used will depend on chromosome length, resolution and data quality. As a rule of thumb the number of diagonals should not exceed 0.1 times length of chromosome.
#'
#' @param hic.comparator object of type HiCcomparator
#' @param diagonals fraction or numeric vector or character "all" which diagonals to use to fit models, by default fraction of chromsome length is used to indicate number of diagonals.
#' @param remove.outliers logical if true try to remove outliers before fitting proper model (NB regression) using robust regression (IRLS) with bisquare weight function
#' @param outlier.weight numeric weight threshold to remove outliers
#' @param link character either log or sqrt - link function for NB regression
#'
#' @return S3 object of class HiCglm
#'
#' @seealso \code{\link{HiCcomparator}} on how to construct HiCcomparator object and \code{\link{robustreg::robustRegBS}} on how robust regression is performed and outliers selection process
#'
#' @examples
#'
#' @export
HiCglm <- function(hic.comparator, diagonals = 0.1, remove.outliers = TRUE, outlier.weight = 0, link = c("log", "sqrt")[1]){
  stopifnot(link %in% c("log", "sqrt"))
  merged <- merge(hic.comparator, include.zero.cells = FALSE)
  nm <- names(merged)
  lapply(nm, function(n){
    merged.n <- merged[[n]]
    # split by diagonals
    md <- split(merged.n, merged.n$diagonal)
    if(diagonals[1] == "all"){
      diagonals <- names(md)
    } else if((length(diagonals) == 1) & (diagonals < 1)){
      m <- round(diagonals * hic.comparator$maps.dims[[n]][1,1])
      diagonals <- as.character(seq(1,m))
    } else {
      diagonals <- as.character(diagonals)
    }
    # build models
    lapply(diagonals, function(k){
      df <- md[[k]]
      # remove outliers if required
      df$outlier <- "no"
      if(remove.outliers){
        tryCatch({
          fit <- robustreg::robustRegBS(val.y ~ val.x, data = df)
          df[which(fit$weights <= outlier.weight), "outlier"] <- "yes"
        }, error = function(e){
          warning(paste0("Failed to perform IRLS for diagonal ",k))
        })
      }
      # build model for enrichment in Y
      df$val.transformed.x <- get(link)(df$val.x) # transform predictor
      model.x <- constructGLM(data.frame(predictor = df$val.transformed.x, response = df$val.y, outlier = df$outlier), link)
      # build model for enrichment in X
      df$val.transformed.y <- get(link)(df$val.y) # transform predictor
      model.y <- constructGLM(data.frame(predictor = df$val.transformed.y, response = df$val.x, outlier = df$outlier), link)
      # results
      list(df, model.x, model.y) %>% magrittr::set_names(c("data","model.x","model.y"))
    }) %>% magrittr::set_names(diagonals)
  }) %>% magrittr::set_names(nm) -> model
  # return result
  structure(list(hic.comparator$data.names, hic.comparator$maps.dims, model) %>%
              magrittr::set_names(c("data.names", "maps.dims", "model")),
            class = "HiCglm")
}

#' Computes significance of data given the model
#'
#' Uses either Negative Binomial or Poisson to calculate pvalue of Y | X or X | Y.
#'
#' @param dataset data frame with predictor, response, outlier columns
#' @param model object of type glm or MASS::glm.nb class
#'
#' @return numeric vector with p-values (upper tail)
#'
#' @seealso \code{\link{pnbinom}}, \code{\link{ppois}} to see how to calculate significance
#'
#' @export
significance <- function(dataset, model){
  stopifnot(model[["family"]]$link %in% c("log", "sqrt"))
  if(model[["family"]]$link == "log"){
    inv.link <- function(v) exp(v)
  }
  if(model[["family"]]$link == "sqrt"){
    inv.link <- function(v) v ** 2
  }
  mu <- predict(model, newdata = dataset)
  mu.muhat <- cbind(dataset$response, mu)
  sapply(1:nrow(mu.muhat), function(i){
    if(model[["family"]]$family == "poisson"){
      ppois(mu.muhat[i,1], inv.link(mu.muhat[i,2]), lower.tail = FALSE)
    } else {
      pnbinom(mu.muhat[i,1], size = model$theta, mu = inv.link(mu.muhat[i,2]), lower.tail = FALSE)
    }
  })
}

#' @export
hicdiff <- function(x, ...) UseMethod("hicdiff", x)

#' Computes differential (p-value/q-value) map.
#'
#' Given HiCglm object calculates enrichment p-values of Y | X or X | Y w.r.t. background models (separate for every diagonal, see \code{\link{HiCglm} for details}).
#'
#' @param hic.glm object of class HiCglm
#' @param marginal.distr character wither fit or obs; if fit then fitted gamma distribution for X and Y is used to convert them to U and V respectively
#'
#' @return list of data frames corresponding to Hi-C contact maps; each data frame contain columns: i, j, name, val.x, val.y, pvalue.x, qvalue.x, pvalue.y, qvalue.y where suffix indicates predictor variable, i.e. pvalue.x indicates E[Y | X] model, so resulting significance means enrichment of Y w.r.t. X
#'
#' @seealso \code{\link{HiCglm}} on how to construct HiCglm, \code{\link{maps_difference_diagonal}} and \code{\link{significance}} on how p-values are calculated
#'
#' @examples
#'
#' @export
hicdiff.HiCglm <- function(hic.glm){
  model <- hic.glm$model
  lapply(names(model), function(n){
    model.n <- model[[n]]
    # compute p values for every diagonal in model
    lapply(names(model.n), function(k){
      l <- model.n[[k]]
      dataset <- l[["data"]]
      model.x <- l[["model.x"]]
      model.y <- l[["model.y"]]
      # enrichment for: Y | X
      px <- significance(data.frame(predictor = dataset$val.transformed.x, response = dataset$val.y),
                         l[["model.x"]])
      # enrichment for: X | Y
      py <- significance(data.frame(predictor = dataset$val.transformed.y, response = dataset$val.x),
                         l[["model.y"]])
      magrittr::inset(dataset, c("pvalue.x", "qvalue.x", "pvalue.y", "qvalue.y"),
                      value = cbind(px, p.adjust(px, method = "BH"), py, p.adjust(py, method = "BH")))
    }) %>% do.call("rbind",.)
  }) %>% magrittr::set_names(names(model))
}

#' Converts significance maps to dense matrices.
#'
#' Given list of significance maps in sparse format produced by \code{\link{hicdiff}} function converts them to dense matrix format.
#'
#' @param hicdiff.list list of data frames, output from \code{\link{hicdiff}} function
#' @param maps.dims list of data frames, usually an attribute of HiCglm object that was used to produce hicdiff.list
#' @param val.column character indicating which column of sparse significance map to use as cell value for dense matrix (either pvalue or qvalue)
#' @param neg.log logical wheter to apply -log transformation to val.column
#' @param which.enrichment character indicating whether to calculate Y enrichment (i.e. E[Y | X] model) - choose y, X enrichment (i.e. E[X | Y] model) - choose x or both
#'
#' @return list with dense matrices containing significance in each cell; when which.enrichment is fixed to both upper triangular map will contain significances of E[Y | X] model and lower triangular will contain significances of E[X | Y]
#'
#' @seealso \code{\link{hicdiff}} for generation of hicdiff.list
#'
#' @examples
#'
#' @export
hicdiff2mtx <- function(hicdiff.list, maps.dims, val.column = c("qvalue","pvalue")[1],
                        neg.log = TRUE, which.enrichment = c("both", "x", "y")[1]){
  stopifnot(val.column %in% c("qvalue","pvalue"))
  stopifnot(which.enrichment %in% c("x", "y", "both"))
  names(hicdiff.list) %>%
    lapply(function(n){
      df <- hicdiff.list[[n]]
      md <- maps.dims[[n]]
      subdf <- df[c("i", "j")]
      vx <- df[,paste0(val.column, ".x")]
      vy <- df[,paste0(val.column, ".y")]
      if(neg.log){
        vx <- -log(vx)
        vy <- -log(vy)
      }
      if(which.enrichment == "y"){
        sparse2dense(magrittr::inset(subdf, "significance", value = vx), N = md[1,"n.rows"])
      } else if(which.enrichment == "x"){
        sparse2dense(magrittr::inset(subdf, "significance", value = vy), N = md[1,"n.rows"])
      } else {
        mtx.x <- sparse2dense(magrittr::inset(subdf, "significance", value = vx), N = md[1,"n.rows"])
        mtx.y <- sparse2dense(magrittr::inset(subdf, "significance", value = vy), N = md[1,"n.rows"])
        mtx <- matrix(0, nrow = md[1,"n.rows"], ncol = md[1,"n.cols"])
        mtx[upper.tri(mtx.x)] <- mtx.x[upper.tri(mtx.x)]
        mtx[lower.tri(mtx.y)] <- mtx.y[lower.tri(mtx.y)]
        return(mtx)
      }
    }) %>% magrittr::set_names(names(hicdiff.list))
}

#' @export
differential_interactions <- function(x, ...) UseMethod("differential_interactions", x)

#' Finds significantly interacting rectangle-like regions.
#'
#' This function works in 3 steps:
#' \itemize{
#'  \item{}{first it calculates differential map,}
#'  \item{}{then it takes negative log significance (qvalue) vector of cells, sorts it and fits bilinear model}
#'  \item{}{finally it retains only those cells that are to the right of intersection point of bilinear model in significance vector (i.e. the most significant cells) and searches for connected components}
#' }
#' Fitting bilinear model is performed using \code{\link{best_fit_bilinear}} function, while for connected components search \code{\link{raster}} package is used. After detection of significantly interacting regions (connected components) one may further filter list to only retain those with number of non zero cells (n.cells column in interacting.regions data frame) larger than some threshold. There are 3 possible ways of selecting significant interactions (cells):
#' \itemize{
#'  \item{}{bilinear model is used to determine significance threshold and then this threshold is compared with \code{pval} parameter - if threshold is less significant than \code{pval} then threshold is substituted with \code{pval} - this is the default behaviour,}
#'  \item{}{only \code{pval} is used as a significance threshold, i.e. hard thresholding,}
#'  \item{}{only bilinear model is used to determine significance threshold (unrecomended, as it may yield non significant interactions).}
#' }
#' When using option 1 and 3 its recommended to plot the fit (enabled by default). An indication of properly determined significance threshold would be when red vertical line (the significance threshold) is located to the right side of grey vertical line.
#'
#' @param hic.glm object of class HiCglm
#' @param plot.models logical if true then plot bilinear model fit for every matrix in hic.glm object; it will plot bilinear fit for E[Y | X] and E[X | Y] models; if you want to save this results to file open device before calling this function (see for instance \code{\link{pdf}}) and close device after function call (see \code{\link{dev.off}})
#' @param pval numeric, pvalue (or qvalue) cutoff to qualify interaction as significant
#' @param sig.thr.selection numeric, if 3 then only use bilinear model fit to establish p-value cutoff for significant interactions, if 2 then select significant interactions using only \code{pval} parameter, if 1 (default) use bilinear model, but if p-value threshold is larger than \code{pval}, use \code{pval} instead
#' @param which.significance character either "qvalue" or "pvalue" indicating, which of the 2 should be used as a measure of interaction significance
#' @param cc.direction specifies criterium for two cells to be considered as neighbors during connected components search, for details see \code{directions} parameter of \code{\link{raster::clump}} function
#'
#' @return list with number of entries equal to hic.glm$names; each entry is a list with 2 elements: interacting.regions - data frame containing rows with rectangle like regions of significant interactions with coordinates n.cells (number of non zero cells inside rectangle), start.x, end.x, start.y, end.y, effect; connected.components list with cells comprising given connected component; connected components list is named list where each entry name is unique id, which can be mapped to row in interacting.regions (its row names); effect column is indicating if interaction refers to Y enrichment (i.e. E[Y | X] model) or X enrichment (i.e. E[X | Y] model)
#'
#' @seealso \code{\link{best_fit_bilinear}} for fitting bilinear model, \code{\link{raster::raster}} and \code{\link{raster::clump}} for connected components search
#'
#' @examples
#'
#' @export
differential_interactions.HiCglm <- function(hic.glm, sig.map = NULL, plot.models = TRUE, pval = 0.05, sig.thr.selection = c(1,2,3)[1], which.significance = c("qvalue","pvalue")[1], cc.direction = c(4,8)[1]){
  if(is.null(sig.map)){
    maps.difference <- hicdiff(hic.glm)
  } else {
    maps.difference <- sig.map
  }
  lapply(names(maps.difference), function(n){
    df <- rbind(
      magrittr::inset(maps.difference[[n]][c("i","j",paste0(which.significance, ".x"))],
                      "effect", value = "Y enrichment") %>% magrittr::set_colnames(c("i","j","significance","effect")),
      magrittr::inset(maps.difference[[n]][c("i","j",paste0(which.significance, ".y"))],
                      "effect", value = "X enrichment") %>% magrittr::set_colnames(c("i","j","significance","effect"))
    )
    df$significance <- -log10(df$significance)
    sig.thr <- -log10(pval)
    # group significantly depleted and enriched cells
    by(df, df$effect, function(df.effect){
      # get appropriate data frame
      dfe <- cbind(df.effect[order(df.effect$significance),], X = seq(1, nrow(df.effect)))
      # find parameter pval based threshold
      x <- dfe[dfe$significance >= sig.thr, "X"][1]
      if(sig.thr.selection != 2){
        # fit bilinear model
        bf <- best_fit_bilinear(dfe$X, dfe$significance)
        x.bm <- bf[["intersection.x"]]["both"]
        if(x.bm < x){
          # bilinear model determined p-value lower than 0.05
          if(sig.thr.selection == 1){
            warning("Bilinear model determined less significant threshold than pval! Switching to pval!")
            x.bm <- x
          } else {
            warning("Bilinear model determined less significant threshold than pval, but sig.thr.selection is set to 3! THIS WILL YIELD NON SIGNIFICANT INTERACTIONS, i.e. p-value > pval!")
          }
        }
        if(plot.models){
          plot(dfe$X, dfe$significance, cex = 0.1, xlab = "X", ylab = latex2exp::TeX("-log_{10}(significance)"))
          title(paste0(hic.glm$data.names[1]," vs ",hic.glm$data.names[2],", ",n,", ",dfe[1,"effect"]), cex.main = 1)
          abline(a = bf[["coefficients"]]["left","intercept"], b = bf[["coefficients"]]["left","slope"], col = "green")
          abline(a = bf[["coefficients"]]["right","intercept"], b = bf[["coefficients"]]["right","slope"], col = "blue")
          legend.labels = c("bilinear model fit 1","bilinear model fit 2")
          colors.labels <- c("green","blue")
          pval.bm <- formatC(10^(-dfe[dfe$X >= bf[["intersection.x"]]["both"],"significance"][1]), format = "e", digits = 2)
          if(x != bf[["intersection.x"]]["both"]){
            legend.labels <- c(legend.labels, paste0("p-value threshold = ",formatC(pval, format = "e", digits = 2)), paste0("intersection fit 1 vs 2, p-value = ", pval.bm))
            colors.labels <- c(colors.labels, "red","gray")
            if(sig.thr.selection == 3 | x < bf[["intersection.x"]]["both"]){
              abline(v = x, col = "gray")
              abline(v = bf[["intersection.x"]]["both"], col = "red")
              legend.labels[3:4] <- legend.labels[4:3]
            } else {
              abline(v = x, col = "red")
              abline(v = bf[["intersection.x"]]["both"], col = "gray")
            }
          } else {
            abline(v = x, col = "red")
            legend.labels <- c(legend.labels, paste0("intersection fit 1 vs 2, p-value = ", pval.bm))
            colors.labels <- c(colors.labels, "red")
          }
          legend("topleft", legend = legend.labels, col = colors.labels, lty = rep(1, length(legend.labels)), bty = "n")
        }
        x <- x.bm
      }
      most.significant <- dfe[dfe$X >= x,]
      # connected components determination
      adjacency <- sparse2dense(most.significant[c("i","j","significance")],
                                N = hic.glm$maps.dims[[n]][1,1])
      # convert negative or positive matrix of negative log p-values into adjacency matrix
      adjacency[adjacency != 0] <- 1
      # take only upper triangular part as the matrix is symmetric
      adjacency[lower.tri(adjacency)] <- 0
      rmat <- raster::raster(adjacency)
      clumps <- matrix(raster::clump(rmat, directions = cc.direction),
                       nrow = nrow(adjacency), ncol = ncol(adjacency), byrow = TRUE)
      # get connected components as list
      lapply(seq(1, max(clumps, na.rm = TRUE)), function(m){
        which(clumps == m, arr.ind = TRUE)
      })
    }, simplify = FALSE) -> l

    ids <- seq(1, length(l[[1]]) + length(l[[2]]))
    ids.neg <- ids[1:length(l[[1]])]
    ids.pos <- ids[(length(l[[1]]) + 1):length(ids)]
    # get data frame where each row is interaction with n.cells, start.x, end.x, start.y, end.y, effect and id as row name
    yenr <- sapply(l[[1]], function(cc){
      c(nrow(cc), min(cc[,"col"]), max(cc[,"col"]), min(cc[,"row"]), max(cc[,"row"]))
    }) %>% t() %>% as.data.frame() %>%
      magrittr::set_colnames(c("n.cells","start.x","end.x","start.y","end.y")) %>%
      magrittr::inset("effect", value = "Y enrichment") %>% magrittr::set_rownames(ids.neg)
    xenr <- sapply(l[[2]], function(cc){
      c(nrow(cc), min(cc[,"col"]), max(cc[,"col"]), min(cc[,"row"]), max(cc[,"row"]))
    }) %>% t() %>% as.data.frame() %>%
      magrittr::set_colnames(c("n.cells","start.x","end.x","start.y","end.y")) %>%
      magrittr::inset("effect", value = "X enrichment") %>% magrittr::set_rownames(ids.pos)
    xenr[c("start.x","end.x","start.y","end.y")] <- xenr[c("start.y","end.y","start.x","end.x")]
    # get list with connected components
    res.cc <- c(magrittr::set_names(l[[1]], ids.neg), magrittr::set_names(l[[2]], ids.pos))
    list(rbind(yenr, xenr), res.cc) %>%
      magrittr::set_names(c("interacting.regions", "connected.components"))
  }) %>% magrittr::set_names(names(maps.difference))
}
