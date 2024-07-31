#' @import data.table
#' @import ggplot2
#' @importFrom magrittr %>% %<>% set_rownames set_colnames set_names extract inset
#' @importFrom stats na.omit quantile kmeans sd as.dist cor relevel p.adjust runif smooth.spline lm prcomp setNames var
#' @importFrom utils object.size write.table head
#' @importFrom graphics abline plot par polygon
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom grDevices dev.off pdf rainbow palette.colors
#' @importFrom grid gpar
#' @importFrom ggsignif geom_signif
NULL

#- deal with . in magrittr
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
