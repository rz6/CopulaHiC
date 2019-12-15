#' Wrapper for fields::image.plot function.
#'
#' Allows to mark -Inf, Inf and NA values in heatmaps with different colors. It will automatically detect such values and assign them color.
#'
#' @param z numeric matrix to be plotted; may contain -Inf, Inf and NA values
#' @param breaks numeric vector of breaks for colorscale
#' @param col character vector of hex color strings; usually generated from some color pallette; it's length must be equal to \code{length(breaks) - 1}
#' @param na.color color for NA values
#' @param neg.inf.color color for -Inf values
#' @param pos.inf.color color for Inf values
#' @param ... additional arguments passed to \code{\link[fields]{image.plot}}
#'
#' @return NULL
#'
#' @seealso \code{\link[fields]{image.plot}} for function which finally handles heatmap plotting
#'
#' @examples
#' # matrix of data
#' mtx <- toeplitz(c(5:1))
#' # make lower triangle part negative
#' mtx[lower.tri(mtx)] <- -mtx[lower.tri(mtx)]
#' # make some cells -Inf
#' mtx[matrix(c(2,1,2,2,3,2), ncol = 2, byrow = TRUE)] <- -Inf
#' # make some cells Inf
#' mtx[matrix(c(3,5,4,5), ncol = 2, byrow = TRUE)] <- Inf
#' # make some cells NA
#' mtx[matrix(c(1,3,2,3,5,3), ncol = 2, byrow = TRUE)] <- NA
#' print(mtx)
#' # prepare breaks --> symmetric, from -5 to 5, spaced by 2, with 0 in the middle
#' # one can also introduce here log, sqrt or other scales by applying proper transformations
#' breaks <- sort(c(0,seq(-5,5,2)))
#' print(breaks)
#' # prepare symmetric color pallette indicating negtive values with blue, middle with white and positive with red
#' colors.pal = c("blue","white","red")
#' pal = colorRampPalette(colors.pal)
#' colors <- pal(length(breaks) - 1L)
#' # finally plot matrix
#' image_plot_na(mtx, breaks, colors)
image_plot_na <- function(z,  breaks, col, na.color = 'black',
                          neg.inf.color = "gold", pos.inf.color = "darkgreen", ...){
  zlim <- c(breaks[1], breaks[length(breaks)])
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  range.fraction <- round((zlim[2] - zlim[1]) * 0.05)
  cbar.coords <- round(c(seq(zlim[1] + range.fraction, 0, length.out = 5),
                         seq(0, zlim[2] - range.fraction, length.out = 5)[2:5]))
  cbar.labels <- cbar.coords
  # handle NA
  newz.na <- zlim[2] + zstep # new z for NA
  breaks <- c(breaks, newz.na)
  z[which(is.na(z > zlim[2]))] <- newz.na # same for newz.na
  zlim[2] <- zlim[2] + zstep # extend top limit to include na
  col <- c(col, na.color) # we construct the new color range by including: na.color
  cbar.coords <- c(cbar.coords, newz.na)
  cbar.labels <- c(cbar.labels, "NA")
  # handle -Inf
  if(any(z == -Inf)){
    # prepend extra color for negative infinity to color palette
    newz.neg.inf <- zlim[1] - zstep
    breaks <- c(newz.neg.inf, breaks)
    z[which(z == -Inf)] <- newz.neg.inf
    zlim[1] <- zlim[1] - zstep
    col <- c(neg.inf.color, col)
    cbar.coords <- c(newz.neg.inf, cbar.coords)
    cbar.labels <- c("-Inf", cbar.labels)
  }
  # handle Inf
  if(any(z == Inf)){
    # append extra color for positive infinity to color palette
    newz.pos.inf <- zlim[2] + zstep
    breaks <- c(breaks, newz.pos.inf)
    z[which(z == Inf)] <- newz.pos.inf
    zlim[2] <- zlim[2] + zstep
    col <- c(col, pos.inf.color)
    cbar.coords <- c(cbar.coords, newz.pos.inf)
    cbar.labels <- c(cbar.labels, "Inf")
  }

  params <- list(...)
  if(!("xlab" %in% names(params))){
    params["xlab"] <- ""
  }
  if(!("ylab" %in% names(params))){
    params["ylab"] <- ""
  }
  params[c("x", "y", "z", "zlim", "col", "breaks")] <- list(seq(1,ncol(z)), seq(1,nrow(z)), z,  zlim, col, breaks)
  colorbar <- TRUE
  if("colorbar" %in% names(params)){
    colorbar <- params[["colorbar"]]
    params[["colorbar"]] <- NULL
  }
  if(colorbar){
    params[["axis.args"]] <- list(at = cbar.coords, labels = cbar.labels)
    do.call(fields::image.plot, params)
  } else {
    do.call(image, params)
  }
}

#' Plots simple contact map or fold-change contact map with log2 scale.
#'
#' @param contact.map numeric matrix (with non-negative entries) to be plotted
#' @param fc logical whether given matrix is fold change
#' @param breaks numeric number of breaks on color scale
#' @param colors.pal colors for pallette
#' @param title character plot title
#' @param useRaster logical if TRUE (default) lower quality image is produced, but the process is much faster and resulting file uses little space
#'
#' @return NULL
#'
#' @seealso \code{\link[DIADEM]{image_plot_na}} for function handling heatmap plotting
#'
#' @examples
#' # get sample contact map (MSC replicate 1) for chromosome 18
#' mtx1.sparse <- DIADEM::sample_hic_maps[["MSC-HindIII-1_40kb-raw"]][["18"]]
#' # convert to dense
#' mtx1 <- sparse2dense(mtx1.sparse[c("i","j","val")], N = 1952)
#' plot_contact_map(mtx1)
#' # now make fold-change matrix
#' # get another sample contact map (MSC replicate 2) for chromosome 18
#' mtx2.sparse <- DIADEM::sample_hic_maps[["MSC-HindIII-2_40kb-raw"]][["18"]]
#' merged <- base::merge(mtx1.sparse, mtx2.sparse, by = c("i", "j"))
#' merged$fc <- merged$val.x / merged$val.y
#' fc.mtx <- sparse2dense(merged[c("i","j","fc")], N = 1952)
#' plot_contact_map(fc.mtx, fc = TRUE, colors.pal = c("blue","white","red"))
#'
#' @export
plot_contact_map <- function(contact.map, fc = FALSE, breaks = 100, zeros.na = TRUE,
                             colors.pal = c("white","red"), useRaster = TRUE, ...){
  pal = colorRampPalette(colors.pal)
  if(zeros.na){
    contact.map[contact.map == 0] <- NA
  }
  if(fc){
    cm <- contact.map
  } else {
    cm <- contact.map + 1
  }
  mtx.dim <- dim(cm)
  # log transformation of color range
  mtx.color.range <- log2(cm[(!is.na(cm)) & (cm > 0) & (abs(cm) != Inf)])
  color.range.min <- min(mtx.color.range)
  color.range.max <- max(mtx.color.range)
  if(fc){
    abs.max <- max(abs(c(color.range.min,color.range.max)))
    color.range.min <- -abs.max
    color.range.max <- abs.max
  }
  colorBreaks <- 2 ^ seq(color.range.min, color.range.max, by = (color.range.max - color.range.min) / breaks)
  colors <- pal(length(colorBreaks) - 1L)
  # if(any(is.na(contact.map))){
  #   na.col.val <- 2 ^ (color.range.max + (color.range.max - color.range.min) / breaks)
  #   cm[is.na(cm)] <- na.col.val
  #   colorBreaks <- c(colorBreaks, na.col.val)
  #   colors <- c(colors, na.color)
  # }

  image_plot_na(cm, colorBreaks, colors, ...)
}

#' Plots differential map.
#'
#' Draws differential map (i.e. one with positive and negative entries).
#'
#' @param mtx.dense numeric matrix with positive and negative entries; if matrix does not contain any negative values use \code{\link{plot_contact_map}} function
#' @param zeros.na logical if TRUE convert zero cells to NA
#' @param breaks numeric number of breaks on color scale
#' @param colors.pal colors for pallette
#' @param color.range numeric vector of length 2 or NULL; if specified gives minimum and maximum values for color scale; this is manual adjustment of scale
#' @param sqrt.transform logical if TRUE apply sqrt transformation: sqrt(pos) to positive elements of matrix and -sqrt(abs(neg)) to negative elements of matrix
#' @param useRaster logical if TRUE (default) lower quality image is produced, but the process is much faster and resulting file uses little space
#'
#' @return NULL
#'
#' @seealso \code{\link[DIADEM]{image_plot_na}} for function handling heatmap plotting
#'
#' @examples
#' # get sample contact map (MSC replicate 1) for chromosome 18
#' mtx1.sparse <- DIADEM::sample_hic_maps[["MSC-HindIII-1_40kb-raw"]][["18"]]
#' # get another sample contact map (MSC replicate 2) for chromosome 18
#' mtx2.sparse <- DIADEM::sample_hic_maps[["MSC-HindIII-2_40kb-raw"]][["18"]]
#' # make differential map
#' merged <- base::merge(mtx1.sparse, mtx2.sparse, by = c("i", "j"))
#' merged$difference <- merged$val.x - merged$val.y
#' dense <- sparse2dense(merged[c("i","j","difference")], N = 1952)
#' # plot
#' plot_diff_map(dense)
#' # plot with sqrt transformation of data
#' plot_diff_map(dense, sqrt.transform = TRUE)
#'
#' @export
plot_diff_map <- function(mtx.dense, zeros.na = TRUE, breaks = 10,
                          colors.pal = c("blue","white","red"),
                          color.range = NULL, sqrt.transform = FALSE,
                          na.color = 'black', neg.inf.color = "gold",
                          pos.inf.color = "darkgreen", useRaster = TRUE, ...){
  mtx <- mtx.dense
  if(zeros.na){
    mtx[mtx == 0] <- NA
  }
  mtx.no.inf <- mtx[abs(mtx) != Inf]
  if(sqrt.transform){
    abs.max <- max(c(
      sqrt(max(abs(mtx.no.inf[mtx.no.inf <= 0]), na.rm = TRUE)),
      sqrt(max(mtx.no.inf[mtx.no.inf >= 0], na.rm = TRUE))
    ))
  } else {
    abs.max <- max(abs(c(min(mtx.no.inf, na.rm = TRUE), max(mtx.no.inf, na.rm = TRUE))))
  }
  if(is.null(color.range)){
    color.range.min <- -abs.max
    color.range.max <- abs.max
  } else {
    color.range.min <- color.range[1]
    color.range.max <- color.range[2]
  }
  if(sqrt.transform){
    n <- (color.range.max - color.range.min) / (2 * breaks)
    colorBreaks <- unique(c(
      seq(color.range.min, 0, by = n),
      seq(0, color.range.max, by = n)
    ))
  } else {
    colorBreaks <- seq(color.range.min, color.range.max, by = (color.range.max - color.range.min) / breaks)
  }
  pal = colorRampPalette(colors.pal)
  colors <- pal(length(colorBreaks) - 1L)

  image_plot_na(mtx, colorBreaks, colors, na.color = na.color,
                neg.inf.color = neg.inf.color, pos.inf.color = pos.inf.color, ...)
}

#' Plots A/B compartment vector.
#'
#' @param pc numeric vector
#'
#' @return NULL
#'
#' @seealso \code{\link{do_pca}} on A/B compartments and how to determine them, \code{\link{HiCcomparator}} object on real data example of compartments and plotting them (examples section)
#'
#' @examples
#' # for real data example check ?HiCcomparator
#' # below is with simulated data
#' v <- sin(seq(1, 10, length.out = 100))
#' plot_pc_vector(v)
#'
#' @export
plot_pc_vector <- function(pc, colors=c("blue","red"), ...){
  x <- seq(length(pc))
  n <- length(x)
  pc[is.na(pc)] <- 0
  y.pos <- pc
  y.neg <- pc
  y.pos[pc<0] <- 0
  y.neg[pc>0] <- 0
  plot(pc~x, type="n", xaxs = "i", ...)
  grid()
  polygon(c(x[1],x,x[n]), c(0,y.pos,0), col=colors[1], border=NA)
  polygon(c(x[1],x,x[n]), c(0,y.neg,0), col=colors[2], border=NA)
}

#' Plots regions (TADs or Long Range interactions) on contact map.
#'
#' @param regions data frame or matrix with 2 (+1) or 4 (+1) columns; 2 columns format is for tads - first column is start bin, second column is end bin; 4 columns format is for LR interactions (rectangle like regions) - columns 1,2 are start and end bin of first interacting region, columns 3,4 are start and end bin of second interacting region; Last column (optional) is category column - like effect column with depletion, no.change, enrichment values to color TADs or LR interactions accordingly
#' @param pal.colors character vector of colors if additional category column is specified - how to color categories
#'
#' @return NULL
#'
#' @seealso \code{\link{plot_contact_map}}, \code{\link{plot_diff_map}} for plotting contact maps or \code{\link{differential_interactions}} for determining and plotting p-value map with differential interactions
#'
#' @examples
#' # plot contact map or differential map - see ?plot_diff_map
#' plot_diff_map(dense, sqrt.transform = TRUE)
#' # then get tads and plot them
#' tads <- DIADEM::sample_tads[["MSC-HindIII-1_40kb-raw"]]
#' tads18 <- tads[tads$name == "18",]
#' plot_regions(tads18[c("start","end")])
#' # for plotting differential interactions see ?differential_interactions
#'
#' @export
plot_regions <- function(regions, pal.colors = NULL, lty = 1, lwd = 0.5){
  # regions - data drame must have 2 or 4 columns and optionally one additional column
  # - 2 columns is for TAD plotting --> in which case first column is start and the other is end
  # - 4 columns is for differential regions plotting --> x.start, x.end, y.start, y.end
  # - optional column if given must be convertible to factor and is category by which regions should be colored
  n.cols <- ncol(regions)
  if((n.cols == 3 | n.cols == 5)){
    cols <- colors()
    n <- length(base::levels(factor(regions[,n.cols])))
    if(is.null(pal.colors) | n != length(pal.colors)){
      # select colors excluding white
      pal.colors <- cols[round(seq(2,length(cols), length.out = 3))]
    }
  } else {
    if(is.null(pal.colors)){
      pal.colors <- c("black")
    }
  }
  if(n.cols == 2){
    regions <- cbind(regions,regions)
  } else if(n.cols == 3) {
    regions <- cbind(regions[,1:2], regions[,1:2], regions[,3])
  }
  if(ncol(regions) == 4){
    regions <- cbind(regions, rep("significant regions", nrow(regions)))
  }
  regions <- magrittr::set_colnames(as.data.frame(regions), c("xs","xe","ys","ye","category"))
  l <- split(regions, regions$category)
  for(i in seq(1,length(l))){
    df <- l[[i]]
    rect(df$xs, df$ys, df$xe, df$ye, col = NA, border = pal.colors[i], lty = lty, lwd = lwd)
  }
}

#' Plots contact map or diff map with inset.
#'
#' Additionally it can plot regions on both orignal image and inset.
#'
#' @param args.map named list of arguments for map plotting function for type of args see \code{\link{plot_contact_map}} and \code{\link{plot_diff_map}}
#' @param xlim numeric 2-element vector o x limits
#' @param ylim numeric 2-element vector o y limits
#' @param which.map character string indicating if contact map or diff map should be plotted
#' @param mar numeric vector of length 4 specifying margins
#' @param args.regions named list of arguments passed to \code{\link{plot_regions}} function, if NULL then don not plot regions
#'
#' @return NULL
#'
#' @seealso \code{\link{plot_contact_map}}, \code{\link{plot_diff_map}} for plotting contact maps and difference maps and \code{\link{plot_regions}} for plotting regions and its arguments
#'
#' @examples
#' # get Hi-C map file
#' mtx.fname <- system.file("extdata", "MSC-HindIII-1_40kb-raw.npz", package = "DIADEM", mustWork = TRUE)
#' # read it and take chromosome 18
#' m.sparse18 <- read_npz(mtx.fname, mtx.names = c("18"))[["18"]]
#' dense <- sparse2dense(m.sparse18[c("i","j","val")], N = 1952)
#' plot_with_inset(list(dense), c(500,800), c(500,800))
#' # get TADs
#' tads <- map2tads(dense)
#' # plot with TADs
#' plot_with_inset(list(dense), c(500,800), c(500,800), args.regions = list(tads))
#'
#' @export
plot_with_inset <- function(args.map, xlim, ylim, which.map = c("contact.map","diff.map")[1],
                            mar = c(1, 2, 1, 0.1), args.regions = NULL){
  # adjust margins for inset
  par(mar = mar)
  # create layout: inset on the left hand side, full size map on the right hand side
  layout(matrix(c(1,2,3,3), ncol = 2, byrow = FALSE),
         widths = c(1,2), heights = c(2,2))
  # plot inset
  if(which.map == "contact.map"){
    do.call(plot_contact_map, c(args.map, list(xlim = xlim, ylim = ylim, colorbar = FALSE)))
  } else {
    do.call(plot_diff_map, c(args.map, list(xlim = xlim, ylim = ylim, colorbar = FALSE)))
  }
  if(!is.null(args.regions)){
    # plot regions
    do.call(plot_regions, args.regions)
  }
  # empty space
  plot.new()
  # plot full size map
  # adjust margins for full size plot
  par(mar = c(1, 0.5, 0.1, 0.1))
  if(which.map == "contact.map"){
    do.call(plot_contact_map, args.map)
  } else {
    do.call(plot_diff_map, args.map)
  }
  # TODO: TADs are shifted on full size map - must fix it
  #if(!is.null(args.regions)){
  #  # plot regions
  #  do.call(plot_regions, args.regions)
  #}
}
