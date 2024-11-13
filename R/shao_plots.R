#' Boxplots from matrix with NON-unique gene symbols
#'
#' Draw boxplots for matrix with unique gene symbols, such as microaray data or matrix use ensembl ids as rownames.
#' Use boxplt_uniq_mult internally.
#'
#' @param x, a gene symbol or a probe/ensembl id.
#' @param is_unique, is \code{"x"} representing a unique id, e.g., proble id.
#' @param gene_anno, match between ids and symbols, MUST have "uniq_id", and "symbol" columns.
#' @param mtx, gene matrix, ids by symbols.
#' @param ..., parameters to
#' @return Save the figures in the sepecified folder.
#' @export
boxplt_mult <- function(x, is_unique = FALSE, gene_anno, mtx, ...) {
  symbol <- uniq_id <- NULL
  if (is_unique) {
    all_id <- x
  } else {
    all_id <- gene_anno[symbol == x, uniq_id] %>%
      intersect(rownames(mtx))
  }

  lapply(all_id, boxplt_uniq_mult, gene_anno = gene_anno, mtx = mtx, ...)
}


#' Boxplots from matrix with unique ids
#'
#' Draw boxplots for matrix with unique ids but not symbols, such as microaray data or matrix use ensembl ids as rownames.
#'
#' @param one_id, a single unique id which is rownames of the input.
#' @param gene_anno, match between ids and symbols, MUST have "uniq_id", and "symbol" columns.
#' @param mtx, input matrix, row names contains \code{"one_id"}.
#' @param s_meta, meta info of samples with grouping info for plotting. There MUST be a "match_id" column where are colnames of mtx (orders don't matter).
#' @param cond, a column name to be used as grouping (e.g., x-axis).
#' @param plt_palette, color palette for plotting.
#' @param test_comp, list for comparisons among groups. E.g., list(c("A", "B")): A vs. B. Student's t-test is used.
#' @param test_method, t-test is used as default.
#' @param paired_test, whether paired t-test is used.
#' @param title, title info.
#' @param subtitle, subtitle info.
#' @param xlab, xlab info.
#' @param ylab, ylab info.
#' @param width, width of the figure.
#' @param height, height of the figure.
#' @param plt_pref, prefix for the figure names.
#' @param plt_dir, folder to save the figure.
#' @return Save the figures in the sepecified folder.
#' @export
boxplt_uniq_mult <- function(one_id,
    gene_anno,
    mtx,
    s_meta,
    cond,
    plt_palette,
    test_comp = list(),
    test_method = "t.test",
    paired_test = FALSE,
    title,
    subtitle = NULL,
    xlab = "",
    ylab = "",
    width = 6,
    height = 7,
    plt_pref = "",
    plt_dir = ".") {
  uniq_id <- value <- NULL
  plt_dta <- data.table(value = mtx[one_id, ], match_id = colnames(mtx)) %>%
    merge(s_meta, by.x = "match_id", by.y = "match_id")

  #- Stop if multiple exist.
  g <- gene_anno[uniq_id == one_id]$symbol
  if (length(g) > 1) browser()

  if (!dir.exists(plt_dir)) dir.create(plt_dir, showWarnings = FALSE, recursive = TRUE)

  p01 <- ggplot(as.data.frame(plt_dta), aes(.data[[cond]], value, fill = .data[[cond]])) +
    geom_boxplot() +
    labs(x = xlab, y = paste(g, one_id, ylab), title = stringr::str_wrap(title, 30), subtitle = stringr::str_wrap(subtitle, 50)) +
    scale_fill_manual(values = plt_palette, guide = "none") +
    geom_signif(comparisons   = test_comp,
                test          = test_method,
                test.args     = list(paired = paired_test),
                step_increase = 0.1,
                color         = "black") +
    theme_bw(16)

  ggsave(paste0(plt_pref, "_", g, "_", one_id, ".png") %>% file.path(plt_dir, .),
         p01,
         width = width,
         height = height)
}


#' Boxplots or barplots from matrix with unique gene symbols
#'
#' Draw barplots or boxplots for matrix with unique gene symbols, such as matrix use symbols as rownames.
#'
#' @param g, a gene to be plotted
#' @param mtx, genes by samples matrix, \code{"g"} used rownames.
#' @param s_meta, meta info of samples with grouping info for ploting. There MUST be a "match_id" column where are colnames of mtx (orders don't matter).
#' @param cond, a column name to be used as grouping (e.g., x-axis).
#' @param plot_type, "boxplot" or "barplot".
#' @param plt_palette, color palette for plotting.
#' @param with_beeswarm, whether add beeswarm points to the boxplots.
#' @param test_comp, list for comparisons among groups. E.g., list(c("A", "B")): A vs. B. Student's t-test is used.
#' @param pre_calulated_test, a data.table MUST with columns "symbol" (unique) and "pvalue".
#' @param test_method, t-test is used as default.
#' @param title, title info.
#' @param subtitle, subtitle info.
#' @param title_wrap, number of characters to wrap the title.
#' @param xlab, xlab info.
#' @param ylab, ylab info.
#' @param ylim, a vector of y limits. c(NA, 1): NA means automatic scaling.
#' @param ylab_rotate, whether rotate the ylab 90 degree.
#' @param paired_test, whether paired t-test is used.
#' @param return_plt, whether return the ggplot object.
#' @param width, width of the figure.
#' @param height, height of the figure.
#' @param plt_pref, prefix for the figure names.
#' @param plt_dir, folder to save the figure.
#' @return Figures in the sepecified folder.
#' @export
bbplt <- function(g,
    mtx,
    s_meta,
    cond,
    plot_type = c("boxplot", "barplot"),
    plt_palette,
    with_beeswarm = FALSE,
    test_comp = list(),
    pre_calulated_test = NULL,
    test_method = "t.test",
    title,
    subtitle = NULL,
    title_wrap = 50,
    xlab = "",
    ylab = "",
    ylim = NULL,
    ylab_rotate = FALSE,
    paired_test = FALSE,
    return_plt = FALSE,
    width = 6,
    height = 7,
    plt_pref = "",
    plt_dir = ".") {
  if (!("ggsignif" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package ggsignif.")
  if (!("ggbeeswarm" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package ggbeeswarm.")

  value <- symbol <- pvalue <- NULL

  if (g %in% rownames(mtx)) {
    plt_dta <- data.table(value = mtx[g, ], match_id = colnames(mtx)) %>%
      merge(s_meta, by.x = "match_id", by.y = "match_id")

    if (!dir.exists(plt_dir)) dir.create(plt_dir, showWarnings = FALSE, recursive = TRUE)

    if (plot_type == "boxplot") {
      p00 <- ggplot(as.data.frame(plt_dta), aes(.data[[cond]], value, fill = .data[[cond]])) +
        geom_boxplot()

      if (with_beeswarm) p00 <- p00 + geom_beeswarm(cex = 3)
    } else {
      p00 <- ggplot(as.data.frame(plt_dta), aes(.data[[cond]], value, fill = .data[[cond]])) +
        stat_summary(geom = "bar", fun = mean, position = "dodge", width = 0.6) +
        stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2) +
        geom_point(position = position_jitter(width = 0.05))
    }

    if (! is.null(ylim)) p00 <- p00 + ylim(ylim)

    if (is.null(pre_calulated_test)) {
      p01 <- p00 +
        labs(x = xlab, y = paste(g, ylab), title = stringr::str_wrap(title, title_wrap), subtitle = stringr::str_wrap(subtitle, 50)) +
        scale_fill_manual(values = plt_palette, guide = "none") +
        geom_signif(comparisons   = test_comp,
                    test          = test_method,
                    test.args     = list(paired = paired_test),
                    step_increase = 0.1,
                    map_signif_level = \(p) paste("p =", signif(p, 4)),
                    color         = "black") +
        theme_bw(16)
    } else {
      pvalue <- signif(pre_calulated_test[symbol == g]$pvalue, 4)
      p01 <- p00 +
        labs(x = xlab, y = paste(g, ylab), title = stringr::str_wrap(title, title_wrap), subtitle = stringr::str_wrap(subtitle, 50)) +
        scale_fill_manual(values = plt_palette, guide = "none") +
        geom_signif(comparisons   = test_comp,
                    annotations = paste("p =", pvalue),
                    step_increase = 0.1,
                    color         = "black") +
        theme_bw(16)
    }


    if (ylab_rotate) p01 <- p01 + theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggsave(paste0(plt_pref, "_", g, ".png") %>% file.path(plt_dir, .),
           p01,
           width  = width,
           height = height)

    if (return_plt) return(p01)
  } else {
    message("(II) ", g, " NOT found!")
  }
}


#' Makes a stream plot
#'
#' Make a stream plot with value on a matrix.
#'
#' @param x a vector of value to be used in x axis.
#' @param y a matrix of m * n, while m equals to the length of x. The columns are the value to be plot.
#' @param order.method a value of c('as.is', 'max', 'first').
#' @param center if TRUE, the stacked polygons will be centered so that the middle, i.e. baseline ('g0'), of the stream is approximately equal to zero.
#' @param frac.rand fraction of the overall data 'stream' range used to define the range of random wiggle (uniform distrubution) to be added to the baseline g0.
#' @param spar setting for smooth.spline function to make a smoothed version of baseline g0.
#' @param xlab label for x axis.
#' @param ylab label for y axis.
#' @param col fill colors for polygons corresponding to y columns (will recycle).
#' @param ylim limts for y axis.
#' @param border border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details).
#' @param lwd border line width for polygons corresponding to y columns (will recycle).
#' @param ... other plot arguments.
#' @return A plot.
#' @details
#'   In order.method, as.is, plot in order of y column; max, plot in order of when each y series reaches maximum value;
#'   first, plot in order of when each y series first value > 0.
#'
#'   Centering is done before the addition of random wiggle to the baseline.
#' @references
#'   Stackover thread, https://stackoverflow.com/a/13087137/349054
#' @examples
#' \dontrun{
#' set.seed(1)
#' m <- 500
#' n <- 50
#' x <- seq(m)
#' y <- matrix(0, nrow=m, ncol=n)
#' colnames(y) <- seq(n)
#' for(i in seq(ncol(y))) {
#'     mu <- runif(1, min=0.25*m, max=0.75*m)
#'     SD <- runif(1, min=5, max=30)
#'     TMP <- rnorm(1000, mean=mu, sd=SD)
#'     HIST <- hist(TMP, breaks=c(0,x), plot=FALSE)
#'     fit <- smooth.spline(HIST$counts ~ HIST$mids)
#'     y[,i] <- fit$y
#' }
#'
#' y <- replace(y, y<0.01, 0)
#'
#' COLS <- rainbow(ncol(y))
#'
#' plt.stream(x,y, axes=FALSE, xlim=c(100, 400), xaxs="i", center=TRUE,
#'            spar=0.2, frac.rand=0.1, col=COLS, border=1, lwd=0.1)
#' }
#' @export
plt_stream <- function(x,
    y,
    order.method = "as.is",
    center       = TRUE,
    frac.rand    = 0.1,
    spar         = 0.2,
    xlab         = "",
    ylab         = "",
    col          = rainbow(ncol(y)),
    ylim         = NULL,
    border       = NULL,
    lwd          = 1,
  ...) {

  if (sum(y < 0) > 0) stop("y cannot contain negative numbers")

  if (is.null(border)) border <- par("fg")

  border <- as.vector(matrix(border, nrow = ncol(y), ncol = 1))
  col    <- as.vector(matrix(col, nrow = ncol(y), ncol = 1))
  lwd    <- as.vector(matrix(lwd, nrow = ncol(y), ncol = 1))

  if (order.method == "max") {
    ord    <- order(apply(y, 2, which.max))
    y      <- y[, ord]
    col    <- col[ord]
    border <- border[ord]
  }

  if (order.method == "first") {
    ord    <- order(apply(y, 2, function(x) min(which(x > 0))))
    y      <- y[, ord]
    col    <- col[ord]
    border <- border[ord]
  }

  bottom.old <- x * 0
  top.old    <- x * 0
  polys      <- vector(mode = "list", ncol(y))

  for (i in seq(polys)) {
    if (i %% 2 == 1) {
      # if odd
      top.new    <- top.old + y[, i]
      polys[[i]] <- list(x = c(x, rev(x)), y = c(top.old, rev(top.new)))
      top.old    <- top.new
    }
    if (i %% 2 == 0) {
      # if even
      bottom.new <- bottom.old - y[, i]
      polys[[i]] <- list(x = c(x, rev(x)), y = c(bottom.old, rev(bottom.new)))
      bottom.old <- bottom.new
    }
  }

  ylim.tmp   <- range(sapply(polys, function(x) range(x$y, na.rm = TRUE)), na.rm = TRUE)
  outer.lims <- sapply(polys, function(r) rev(r$y[(length(r$y) / 2 + 1):length(r$y)]))
  mid        <- apply(outer.lims, 1, function(r) mean(c(max(r, na.rm = TRUE), min(r, na.rm = TRUE)), na.rm = TRUE))

  # center and wiggle
  if (center) {
    g0 <- -mid + runif(length(x), min = frac.rand * ylim.tmp[1], max = frac.rand *
      ylim.tmp[2])
  } else {
    g0 <- runif(length(x), min = frac.rand * ylim.tmp[1], max = frac.rand * ylim.tmp[2])
  }

  fit <- smooth.spline(g0 ~ x, spar = spar)

  for (i in seq(polys)) {
    polys[[i]]$y <- polys[[i]]$y + c(fit$y, rev(fit$y))
  }

  if (is.null(ylim))
    ylim <- range(sapply(polys, function(x) range(x$y, na.rm = TRUE)), na.rm = TRUE)

  plot(x, y[, 1], ylab = ylab, xlab = xlab, ylim = ylim, t = "n", ...)

  for (i in seq(polys)) {
    polygon(polys[[i]], border = border[i], col = col[i], lwd = lwd[i])
  }
}
