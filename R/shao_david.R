#' Flatten DAVID GO annotation
#'
#' Conver the summary DAVID results to one gene per line, filtered by Pvalue or FDR.
#'
#' @param GOres a data.table with DAVID GO annotation.
#' @param p_thres threshold for pvalue, only used if fdr_thres is NULL.
#' @param fdr_thres threshold for FDR.
#' @param anno a data.table with two columns: "enid" and "symbol".
#' @return A data.table with flatted DAVID results, one gene per line.
#' @export
flat_david <- function(GOres, p_thres = 0.001, fdr_thres = NULL, anno = NULL) {
  FDR <- PValue <- GO <- NULL

  GOres[, `:=` (fdr_thres = as.numeric(FDR), PValue = as.numeric(PValue))]

  if (is.null(fdr_thres)) {
      GOres %<>% .[PValue <= p_thres]
  } else {
      GOres %<>% .[FDR <= fdr_thres]
  }

  if (nrow(GOres) == 0) stop("(EE) No GO terms kept under the threshold.")
  res <- apply(GOres, 1, function(x) {
                 cc1 <- x["Genes"]  %>%
                   gsub(" ", "", .) %>%
                   stringi::stri_split_fixed(., pattern = ",", simplify = TRUE) %>%
                   .[1, ]
                 return(data.table(GO = x["Term"], genes = cc1, PValue = as.numeric(x["PValue"]), FDR = as.numeric(x["FDR"])))}) %>%
    rbindlist

  if (!is.null(anno)) res <- merge(anno, res, by.x = "geneid", by.y = "genes") %>% .[order(PValue, GO)]
  return(res)
}


#' GO plots for DAVID GO annotation
#'
#' Put the GO terms in lollipop plots or barplots.
#'
#' @param in_file path to the david GO results.
#' @param height pass to \code{"ggsave"} height.
#' @param width  pass to \code{"ggsave"} width.
#' @param thresh_p P_value threshold.
#' @param thresh_fdr FDR threshold (default NULL and using the P value threshold).
#' @param lollipop whether using lollipop plots.
#' @param cex_size size for the row names.
#' @param col_palette color to be used for barplot.
#' @param output_name name of the output file; by default using the base name of in_file.
#' @param save_as_pdf whether to save as pdf.
#' @param return_plot whether to reture the ggplot object.
#' @return Save a png/pdf file.
#' @export
GO_plot_david <- function(in_file,
    height      = NULL,
    width       = 7,
    thresh_p    = 0.001,
    thresh_fdr  = NULL,
    lollipop    = TRUE,
    cex_size    = 6,
    col_palette = "#0073c2ff",
    output_name = NULL,
    save_as_pdf = FALSE,
    return_plot = FALSE) {
  shortterm <- log10pvalue <- Term <- PValue <- shortterm_break <- FDR <- log10fdr <- NULL
  message("(==) Processing: ", in_file)

  david <- fread(in_file)
  if (nrow(david) != 0) {
    if (is.null(thresh_fdr)) {
    david <- david[, `:=`(shortterm = stringr::str_split(Term, "~", simplify = TRUE)[, 2], log10pvalue = -log10(PValue))] %>%
      .[PValue <= thresh_p] %>%
      .[, shortterm_break := stringr::str_wrap(as.character(shortterm), 50)] %>%
      .[, shortterm_break := factor(shortterm_break, levels = rev(shortterm_break))]
    } else {
    david <- david[, `:=`(shortterm = stringr::str_split(Term, "~", simplify = TRUE)[, 2], log10fdr = -log10(FDR))] %>%
      .[FDR <= thresh_fdr] %>%
      .[, shortterm_break := stringr::str_wrap(as.character(shortterm), 50)] %>%
      .[, shortterm_break := factor(shortterm_break, levels = rev(shortterm_break))]
    }

    f_base_name <- rsplit(basename(in_file), ".", 2)[1]

    if (is.null(output_name)) output_name <- f_base_name

    if (nrow(david) > 0) {
      est_height <- round(nrow(david) / 3)
      if (is.null(height)) height <- ifelse(est_height < 3, 3, est_height)

      if (lollipop) {
        if (is.null(thresh_fdr)) {
          p01 <- ggplot(as.data.frame(david), aes(shortterm_break, log10pvalue)) +
            geom_segment(aes(x = shortterm_break, xend = shortterm_break, y = 0, yend = log10pvalue)) +
            geom_point(color = col_palette, size = 4) +
            coord_flip() +
            theme_bw(cex_size) +
            labs(y = "P value (-log10)", x = "", title = stringr::str_wrap(as.character(output_name), 10))
        } else {
          p01 <- ggplot(as.data.frame(david), aes(shortterm_break, log10fdr)) +
            geom_segment(aes(x = shortterm_break, xend = shortterm_break, y = 0, yend = log10fdr)) +
            geom_point(color = col_palette, size = 4) +
            coord_flip() +
            theme_bw(cex_size) +
            labs(y = "FDR (-log10)", x = "", title = stringr::str_wrap(as.character(output_name), 10))
        }
      } else {
        if (is.null(thresh_fdr)) {
          p01 <- ggplot(as.data.frame(david), aes(shortterm, log10pvalue)) +
            geom_col(fill = alpha(col_palette, 0.5), width = 0.7) +
            coord_flip() +
            theme_bw(cex_size) +
            geom_text(aes(label = shortterm, y = 1),  hjust = -.005, size = 4) +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
            labs(y = "P value (-log10)", x = "", title = stringr::str_wrap(as.character(output_name), 10))
        } else{
          p01 <- ggplot(as.data.frame(david), aes(shortterm, log10fdr)) +
            geom_col(fill = alpha(col_palette, 0.5), width = 0.7) +
            coord_flip() +
            theme_bw(cex_size) +
            geom_text(aes(label = shortterm, y = 1),  hjust = -.005, size = 4) +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
            labs(y = "FDR (-log10)", x = "", title = stringr::str_wrap(as.character(output_name), 10))
        }
      }

      if (is.null(thresh_fdr)) {
        if (save_as_pdf) {
          file_name <- paste0("thresh_p", "_", thresh_p, "_", output_name, ".pdf")
        } else {
          file_name <- paste0("thresh_p", "_", thresh_p, "_", output_name, ".png")
        }
      } else {
        if (save_as_pdf) {
          file_name <- paste0("thresh_fdr", "_", thresh_fdr, "_", output_name, ".pdf")
        } else {
          file_name <- paste0("thresh_fdr", "_", thresh_fdr, "_", output_name, ".png")
        }
      }

      ggsave(file_name, p01, height = height, width = width)

      if (return_plot) return(p01)

    } else {
      message("(II) No enriched GO under the threshold.")
    }
  } else {
    message("(II) Empty file.")
  }
}
