#' @importFrom knitr kable
#' @export
knitr::kable


#' Load a set of libraries
#' @export
loadlib <- function() lapply(c("anatools", "colorout", "magrittr", "data.table"), library, character.only = TRUE) %>% invisible()


#' Get the ids of down/up-regulated genes from DESeq2.
#' @param x path to a file, which must have two columns of "log2FoldChange" and "geneid"
#' @param to_file whether save up-/down-regulated files.
#' @return A list gene ids of up-/down-regulated genes. Related files are writen to the same folder.
#' @export
split_updown <- function(x, to_file = TRUE) {
  log2FoldChange <- geneid <- NULL
  dta    <- fread(x)
  f_name <- basename(x) %>% tools::file_path_sans_ext()

  if (to_file) {
    write.table(dta[log2FoldChange > 0], paste0(f_name, "_up.csv"),   quote = FALSE, row.names = FALSE, sep = ",")
    write.table(dta[log2FoldChange < 0], paste0(f_name, "_down.csv"), quote = FALSE, row.names = FALSE, sep = ",")
  }

  return(list(dta[log2FoldChange > 0, geneid], dta[log2FoldChange < 0, geneid]) %>%
         set_names(paste0(f_name, c("_up", "_down"))))
}


#' Display the corner of data
#' @param x data.
#' @param r number of rows to show.
#' @param c number of columns to show.
#' @export
cn <- function(x, r = 5L, c = 5L) {
  r <- if (nrow(x) < r) nrow(x) else r
  c <- if (ncol(x) < c) ncol(x) else c
  x[1:r, 1:c, drop = FALSE]
}


#' Head
#' @param ... further parameters to \code{"head"}.
#' @export
h <- function(...) head(...)

#' Show the data in a html page
#' @param ... further parameters to \code{"DT::datatable"}.
#' @export
webtable <- function(...) {
  if (!("DT" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package DT.")
  DT::datatable(class = "cell-border display", filter = "top", options = list(DT.TOJSON_ARGS = list(na = "string")), ...)
}


#' Read RSEM result files
#'
#' Read the RSEM output files, either gene results or transcript results, to a
#' a list with count, TPM, fpkm.
#'
#' @param inputfiles path to input files.
#' @param is_gene whethe the input are genes (default) or transcript.
#' @return A list of three type of data, each is a data.table of all samples.
#' @examples
#' \dontrun{
#' inputfiles <- dir("/results", pattern=".*.genes.results$", full.names = TRUE, recursive = TRUE)
#' res <- readrsem(inputfiles)
#' countdata <- res$count_data
#' tpmdata <- res$tpm_data
#' fpkmdata <- res$fpkm_data
#' }
#' @export
readrsem <- function(inputfiles, is_gene = TRUE) {
  id     <- NULL
  id2use <- ifelse(is_gene, "gene_id", "transcript_id")
  cc1    <- lapply(inputfiles, fn_sub_readresm, is_gene)

  message("(II) ", length(cc1), " files read")

  #- check the rownames are identical for records
  cc2 <- lapply(cc1, function(x) x[[id2use]]) %>% setDT %>% setnames(paste0("N", seq_len(ncol(.))))
  if (all(apply(cc2, 2, identical, cc2$N1))) {
    message("(II) All ids matched to each other.")
  } else {
    stop("(EE) Not all gene ids are identical.")
  }

  dt_count <- lapply(cc1, function(x) x[, 2]) %>% unlist(recursive = FALSE) %>% setDT %>% .[, id := cc2$N1] %>% setcolorder("id")
  dt_tpm   <- lapply(cc1, function(x) x[, 3]) %>% unlist(recursive = FALSE) %>% setDT %>% .[, id := cc2$N1] %>% setcolorder("id")
  dt_fpkm  <- lapply(cc1, function(x) x[, 4]) %>% unlist(recursive = FALSE) %>% setDT %>% .[, id := cc2$N1] %>% setcolorder("id")

  return(list(count_data = dt_count, tpm_data = dt_tpm, fpkm_data = dt_fpkm))
}


fn_sub_readresm <- function(x, is_gene = TRUE) {
  message("(==) Processing ", basename(x))
  cc1     <- fread(x)
  experi  <- basename(unlist(rsplit(x, "\\.", 2))[1])
  id2use  <- ifelse(is_gene, "gene_id", "transcript_id")
  coltype <- c("expected_count", "TPM", "FPKM")
  cc2     <- cc1[, c(id2use, coltype), with = FALSE] %>%
    setnames(2:4, paste(experi, coltype, sep = "."))
  return(cc2)
}


#' Split a string from right end
#'
#' Split a string from the righ end and return the result vector. The number of splits could
#' be specified.
#'
#' @param x a string to be splitted.
#' @param s the separator, "fixed".
#' @param n number of desired pieces.
#' @return A matrix with splitted parts.
#' @examples
#' \dontrun{
#' rsplit("A_B_C_D_E", "_", 3)
#' }
#' @export
rsplit <- function(x, s, n = Inf) {
  res <- unlist(stringi::stri_split_fixed(stringi::stri_reverse(x), s, n)) %>%
    vapply(., function(x) stringi::stri_reverse(x), character(1)) %>%
    rev %>%
    set_names(NULL) %>%
    matrix(ncol = n, byrow = T) %>%
    .[nrow(.):1, ]
  return(res)
}


#' Split a string into two part from right end
#'
#' Split a string from the righ end and return the result vector.
#'
#' @param x a string to be splitted.
#' @param s the separator
#' @param n the right location of s.
#' @return A vector with splitted parts.
#' @examples
#' \dontrun{
#' rs2("A_B_C_D_E", "_", 3)
#' # "A_B", "C_D_E"
#' }
#' @export
rs2 <- function(x, s, n) {
  p   <- paste0("[^", s, "]*")
  rx  <- paste0(s, "(?=", paste(rep(paste0(p, s), n - 1), collapse = ""), p, "$)")
  res <- unlist(strsplit(x, rx, perl = TRUE)) %>%
    matrix(ncol = 2, byrow = TRUE) %>%
    .[1:nrow(.), ]
  return(res)
}


#' Extract the differentially expressed genes
#'
#' Get the DE genes from DESeq2 results, removing the NA rows.
#'
#' @param dds a DESeq2 returned results. The input matrix is expected to have the format of
#'        genes by samples, while genes are ensembl gene ids.
#' @param contrast contrast to be used.
#' @param ganno a data.table containing genes annotation, must have two columns with
#'        "geneid", and "symbol".
#' @return A data.table with all results.
#' @examples
#' \dontrun{
#' ll::gde(dds, contrast = c("cond", "CHARGE", "WT"), ganno = unique(ll::mm10enid[, -"transid"]))
#' }
#' @export
# gde <- function(dds, contrast, fc_thresh = 0, count_thresh = 0, fdr_thresh = 0.001, p_thresh = 0.001, ganno) {
gde <- function(dds, contrast, ganno) {
  if (!("DESeq2" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package DESeq2.")
  # padj <- pvalue <- log2FoldChange <- baseMean <- NULL

  res <- DESeq2::results(dds, cooksCutoff = FALSE, contrast = contrast) %>%
    as.data.frame %>%
    as.data.table(keep.rownames = TRUE) %>%
    na.omit %>%
    merge(ganno, ., by.x = "geneid", by.y = "rn", all.y = TRUE)

  col_meta <- SummarizedExperiment::colData(dds) %>% as.data.table(keep.rownames = TRUE)
  col_lite <- col_meta[get(contrast[1]) %in% contrast[2:3]]

  n_count <- DESeq2::counts(dds, normalized = TRUE)
  mtx     <- n_count[, colnames(n_count) %in% col_lite$rn]

  gr_avg <- future.apply::future_apply(mtx, 1, \(x) tapply(x, droplevels(col_lite[[contrast[1]]]), \(y) data.table(avg = mean(y), min_v = min(y), max_v = max(y))) %>% unlist) %>%
    t %>%
    as.data.table(keep.rownames = TRUE)

  # gr_avg <- future.apply::future_apply(mtx[1:5, ], 1, \(x) tapply(x, droplevels(col_lite[[contrast[1]]]), \(y) data.table(avg = mean(y), min_v = min(y), max_v = max(y)))) %>%
  #   t %>%
  #   as.data.table(keep.rownames = TRUE)

  res <- merge(gr_avg, res, by.x = "rn", by.y = "geneid") %>%
    setnames(1, "geneid") %>%
    setcolorder(c("geneid", "symbol"))

  #- Feel not flexibale.
  # if (!is.null(fdr_thresh)) {
  #   filtered_res <- res[padj <= fdr_thresh & abs(log2FoldChange) >= fc_thresh & baseMean >= count_thresh]
  # } else {
  #   filtered_res <- res[pvalue <= p_thresh & abs(log2FoldChange) >= fc_thresh & baseMean >= count_thresh]
  # }

  # ures <- filtered_res[log2FoldChange > 0]
  # dres <- filtered_res[log2FoldChange < 0]
  return(res)
}


#' Object size
#'
#' @param units the units to be used.
#' @return the object size (MB).
#'
#' @export
ob <- function(units = "MB") {
  ob_size <- NULL

  if (length(ls(envir = .GlobalEnv)) != 0) {
    res <- sapply(ls(envir = .GlobalEnv), function(x) format(object.size(get(x)), units = units)) %>%
      {data.table(name = names(.), stringi::stri_split_fixed(., " ", simplify = TRUE))} %>%
      setnames(2:3, c("ob_size", "unit")) %>%
      .[, ob_size := as.numeric(ob_size)] %>%
      .[order(-ob_size)]
    return(res)
  } else {
    message("(II) No objects in the environments.")
  }
}


#' Pattern to grep
#'
#' Build a string from vectors to grep for exact match.
#'
#' @param x a vector of characters to build the string.
#' @param exact exact match of vectors (default TRUE).
#' @return A string to be used in grep or similar regex setting.
#' @export
matchstring <- function(x, exact = TRUE) {
  if (exact) {
    #- with \\b the - is known as separator.
    vapply(unique(x), function(x) paste0("^", x, "$"), character(1)) %>% paste(collapse = "|")
  } else {
    paste(unique(x), collapse = "|")
 }
}


#' Cluster the time-course data by K-means/PAM
#'
#' Perform and plot the K-means or PAM clustering on time course data.
#'
#' @param indata a matrix with gene in rows and time points in columns. Colnames have the pattern of "T_1", "T_2", etc.
#' @param ncl number of k-means clusters.
#' @param timepoint a data table with two column t1 and t2, t1 is the character name of time, and t2 is the numberic time.
#' @param cluster_method K-means or PAM.
#' @param dist_method method to caluclate the metric used in PAM.
#' @param xlab_name name for the xlab.
#' @param ylab_name name for the ylab.
#' @param save_fig whether saving the plot.
#' @param width width of the figure
#' @param height height of the figure
#' @param clustercolour color themes for the smooth line of kmeans clusters.
#' @param figure_prefix prefix for the figure
#' @param show_se whether show the se of the gam smooth line.
#' @return A ggplot object and a csv file with gene clustering information saved in the working folder.
#' @export
tc_cluster <- function(indata,
    ncl,
    timepoint,
    cluster_method = c("kmeans", "pam"),
    dist_method    = function(x) {as.dist((1 - cor(t(x))) / 2)},
    xlab_name,
    ylab_name,
    save_fig       = TRUE,
    width          = 12,
    height         = 10,
    clustercolour  = rainbow(20),
    figure_prefix  = "",
    show_se        = FALSE) {
  if (!("directlabels" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package directlabels")
  if (!("cluster" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package cluster.")
  if (!("ggplot2" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package ggplot2.")

  if (cluster_method == "kmeans") {
    res_kmeans <- kmeans(indata, ncl)
    temdata    <- as.data.frame(indata) %>% inset("genecl", value = paste0("kmeans_", res_kmeans$cluster))
  } else if (cluster_method == "pam") {
    distmat    <- dist_method(indata)
    res_pam    <- cluster::pam(distmat, ncl)
    temdata    <- as.data.frame(indata) %>% inset("genecl", value = paste0("pam_", res_pam$clustering))
  } else {
    stop("(EE) cluster_method should be either kmeans or pam.")
  }

  p01 <- fn_tc_plot(temdata, timepoint, clustercolour, xlab_name, ylab_name, show_se)

  if (figure_prefix != "") figure_prefix <- paste0(figure_prefix, "_")

  if (cluster_method == "kmeans") {
    if (save_fig) ggsave(paste0(figure_prefix, "kmeans_", ncl, ".png"), p01, width = width, height = height)
    write.table(temdata, paste0(figure_prefix, "kmeans_", ncl, ".csv"), quote = FALSE, sep = "\t")
  } else {
    if (save_fig) ggsave(paste0(figure_prefix, "pam_", ncl, ".png"), p01, width = width, height = height)
    write.table(temdata, paste0(figure_prefix, "pam_", ncl, ".csv"), quote = FALSE, sep = "\t")
  }

  return(p01)
}

#- otherwise there will be note for patterns
patterns <- function(...) NULL

fn_tc_plot <- function(indata, timepoint, clustercolour, xlab_name, ylab_name, show_se) {
  N <- value <- genecl <- variable <- sev <- sdv <- t2 <- avgv <- x <- NULL

  indata <- as.data.table(indata, keep.rownames = TRUE)

  data_summary <- melt(indata, id.vars = "genecl", measure.vars = patterns("^T_")) %>%
    .[, .(avgv = mean(value), sdv = sd(value), .N), by = .(genecl, variable)] %>%
    .[, sev := sdv / sqrt(N)]

  data_summary2 <- merge(data_summary, timepoint, by.x = "variable", by.y = "t1")

  p01 <- ggplot(as.data.frame(data_summary2), aes(t2, avgv, group = genecl, color = genecl)) +
    # geom_point() +
    geom_smooth(aes(color = genecl), method = "gam", se = show_se, formula = y ~ s(x, bs = "cr")) +
    theme_classic() +
    scale_color_manual(values = clustercolour) +
    theme(legend.position = "top", legend.title = element_blank()) +
    directlabels::geom_dl(aes(label = N), method = list(directlabels::dl.trans(x = x + .4), "last.points")) +
    labs(x = xlab_name, y = ylab_name)

  return(p01)
}


#' Scale data to arbitrary range
#'
#' Scale vector to sepecified range, by default [0, 1]
#'
#' @param x numeric vector.
#' @param lb low boundary, 0 by default
#' @param ub up boundary, 1 by default.
#' @return A numeric vector.
#' @export
scale_data <- function(x, lb = 0, ub = 1) {
  #- a is the scale factor
  a <- (ub - lb) / (max(x) - min(x))
  #- after scale, how much it needs to move to reach the up boundary.
  b <- ub - a * max(x)
  return(a * x + b)
}


#' Stack the results from FactoMineR::dimdesc
#'
#' Combined the dimdesc results as a single data.table.
#'
#' @param dimdesc_res results from FactoMineR::dimdesc.
#' @return A data.table.
#' @examples
#' \dontrun{
#' #- Input data are samples by features.
#' res_pca     <- FactoMineR::PCA(decathlon[, 1:10], ncp = 5, graph = FALSE)
#' dimdesc_res <- FactoMineR::dimdesc(res_pca, axes = 1:3)
#' }
#' @export
corr_pcvar <- function(dimdesc_res) {
  if (!("FactoMineR" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package FactoMineR")
  p.value <- NULL

  res <- lapply(grep("^Dim", names(dimdesc_res), value = TRUE),
                function(x) as.data.table(dimdesc_res[[x]]$quanti, keep.rownames = TRUE)[, `:=` (pc = x, padj = p.adjust(p.value, "BH"))]) %>%
    rbindlist(fill = TRUE) %>%
    na.omit %>%
    setnames(1, "variable")
  return(res)
}


#' Make a vector unique.
#'
#' Add the sep delim to duplicated elements of vector based on their location in the vector, starting form 1.
#'
#' @param x input vector.
#' @param sep_char suffix.
#' @return A unique vector.
#' @export
make_unique <- function(x, sep_char = "_") {
  i    <- NULL
  cc_1 <- rep("", length(x))
  cc_2 <- names(table(x)[table(x) > 1])

  foreach(i = cc_2) %do% {
    cc_3       <- grep(matchstring(i, exact = TRUE), x)
    cc_1[cc_3] <- paste0(sep_char, as.character(seq_along(cc_3)))
  }

  return(paste0(x, cc_1))
}


#' Discrete colour palettes from the pals package
#'
#' Taken from Seurat, which taken from pals.
#'
#' @param n Number of colours to be generated.
#' @param palette Options are
#'   "alphabet" (26), "alphabet2" (26), "glasbey" (32), "polychrome" (36), and "stepped" (24).
#'   Can be omitted and the function will use the one based on the requested n.
#'
#' @return A vector of colors
#'
#' @export
DiscretePalette <- function(n, palette = NULL) {
  palettes <- list(
    alphabet = c(
      "#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919", "#005C31",
      "#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00",
      "#C20088", "#003380", "#FFA405", "#FFA8BB", "#426600", "#FF0010",
      "#5EF1F2", "#00998F", "#E0FF66", "#740AFF", "#990000", "#FFFF80",
      "#FFE100", "#FF5005"
    ),
    alphabet2 = c(
      "#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356",
      "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F", "#C4451C", "#DEA0FD",
      "#FE00FA", "#325A9B", "#FEAF16", "#F8A19F", "#90AD1C", "#F6222E",
      "#1CFFCE", "#2ED9FF", "#B10DA1", "#C075A6", "#FC1CBF", "#B00068",
      "#FBE426", "#FA0087"
    ),
    glasbey = c(
      "#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300",
      "#FFD300", "#009FFF", "#9A4D42", "#00FFBE", "#783FC1", "#1F9698",
      "#FFACFD", "#B1CC71", "#F1085C", "#FE8F42", "#DD00FF", "#201A01",
      "#720055", "#766C95", "#02AD24", "#C8FF00", "#886C00", "#FFB79F",
      "#858567", "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF",
      "#004CFF", "#F2F318"
    ),
    polychrome = c(
      "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE",
      "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD",
      "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
      "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0",
      "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5",
      "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB"
    ),
    stepped = c(
      "#990F26", "#B33E52", "#CC7A88", "#E6B8BF", "#99600F", "#B3823E",
      "#CCAA7A", "#E6D2B8", "#54990F", "#78B33E", "#A3CC7A", "#CFE6B8",
      "#0F8299", "#3E9FB3", "#7ABECC", "#B8DEE6", "#3D0F99", "#653EB3",
      "#967ACC", "#C7B8E6", "#333333", "#666666", "#999999", "#CCCCCC"
    )
  )
  if (is.null(x = palette)) {
    if (n <= 26) {
      palette <- "alphabet"
    } else if (n <= 32) {
      palette <- "glasbey"
    } else {
      palette <- "polychrome"
    }
  }
  palette.vec <- palettes[[palette]]
  if (n > length(x = palette.vec)) {
    warning("Not enough colours in specified palette")
  }
  palette.vec[seq_len(length.out = n)]
}


#' PC1 location of samples
#'
#' Mainly used for Nanostring data. Get the location of samples on PC1 calculated with signature genes in predefined gene sets.
#'
#' @param mtx samples by genes matrix.
#' @param prb.set.df a data.frame with predefined gene sets, columns are gene sets.
#' @param prunedprobes probes in mtx that are not used.
#' @param min.sd used for scaling.
#' @param mingenes.pc1 min number of genes in a gene set for PCA.
#' @param adjust.for a covariates to be adjusted. Not used by me.
#' @details Extracted and modified from Samira's codes.
#' @return A list of signs (location of sampmles in PCA) and signaturematrix (loadings of signature genes).
#' @export
calc.PC1.scores <- function(mtx,
    prb.set.df,
    prunedprobes = NULL,
    min.sd       = 1,
    mingenes.pc1 = 5,
    adjust.for   = NULL) {
  # prune if needed:
  mtx <- mtx[, !is.element(colnames(mtx), prunedprobes)]

  gene.sets.list <- list()

  for (i in seq_len(ncol(prb.set.df))) {
    if (length(rownames(prb.set.df)[(prb.set.df[, i] != 0) & !is.na((prb.set.df[, i]))]) < 1) {
      gene.sets.list[[i]] <- NA
      names(gene.sets.list)[i] <- colnames(prb.set.df)[i]
    } else {
      gene.sets.list[[i]] <- rownames(prb.set.df)[(prb.set.df[, i] != 0) & !is.na((prb.set.df[, i]))]
      names(gene.sets.list)[i] <- colnames(prb.set.df)[i]
    }
  }

  gene.sets.list        <- lapply(gene.sets.list, function(x) x[!is.na(x)])
  gene.sets.list.final  <- gene.sets.list[lapply(gene.sets.list, length) > 0]


  #>>#########################################################
  #- Not sure whether I need this.
  # adjust for covariates if indicated:
  if (length(adjust.for) > 0) {
    # Check the degrees of freedom:
    temp <- cbind(mtx[, 1], adjust.for)

    if (summary(lm(temp))$df[2] < 1) {
      warning("Warning: Too many parameters or levels to adjust for given the number of observations - Pathway scoring will not be run.")
      return(scores = NULL)
    }

    adjusted <- setNames(mtx, colnames(mtx))  # Need this weird construct because the column names of mtx have a colon in them which gets replaced by a . otherwise
    adjusted[TRUE] <- NA

    for (i in seq_len(dim(adjusted)[2])) {
      # combine with adjust.for to create a mtx frame on which to run a lm:
      temp <- cbind(mtx[, i], adjust.for)
      mod  <- lm(temp)
      # fill in the 'adjusted' mtx frame with the residuals of the lm:
      adjusted[names(mod$resid), i] <- mod$resid
    }

    # remove NA rows:
    adjusted <- adjusted[rowSums(is.na(adjusted)) == 0, ]
    # replace the 'mtx' object:
    mtx <- adjusted
  }
  #<<#########################################################

  # for each gene set, calculate PC1:
  sigs <- matrix(NA, nrow(mtx), length(gene.sets.list))
  dimnames(sigs) <- list(rownames(mtx), names(gene.sets.list))

  signaturematrix <- matrix(0, ncol(mtx), length(gene.sets.list))
  dimnames(signaturematrix) <- list(colnames(mtx), names(gene.sets.list))

  #- PCA for each gene set.
  for (k in seq_along(gene.sets.list)) {
    sharedgenes <- intersect(colnames(mtx), gene.sets.list[[k]])

    if (length(sharedgenes) >= mingenes.pc1) {
      temp <- mtx[, sharedgenes]

      # scale the mtx appropriately:
      scaling.factors <- sqrt(apply(temp, 2, var) + min.sd^2)
      temp2 <- t(t(temp) / scaling.factors)
      # run PCA:
      pct <- prcomp(temp2, scale. = FALSE)
      # reorient the signature if the average gene's weight is negative:
      sigsign <- sign(mean(pct$rotation[, 1]))

      if (sigsign == 0) sigsign <- 1
      # save output:
      sigs[, k] <- pct$x[, 1] * sigsign
      signaturematrix[sharedgenes, k] <- pct$rotation[, 1] * sigsign  ## <---------------- need to undo the scaling!!!!
    }
  }

  sigs <- sigs[, colSums(is.na(sigs)) < nrow(sigs)]
  # signaturematrix <- signaturematrix[,colSums(is.na(sigs)) < nrow(sigs)]
  signaturematrix <- signaturematrix[, which(colSums(signaturematrix) != 0)]
  out <- list(sigs = sigs, signaturematrix = signaturematrix)
  return(out)
}


#' Download and annotate GEO microarray data
#'
#' Download data, add gene symbols, entrez ids if applicable.
#'
#' @param x GEO accession
#' @param GPL GPL id.
#' @param gene_anno gene annotation with three columns: probe_id, symbol, gene_id.
#' @return A list of five elements if gene_anno is provided: ex_f, matrix with filted probes (20% quantile);
#'   ex_anno, annotated matrix; ex_g, matrix with unique symbols;
#'   gset, downloaded GEO data; ex_all, unfiltered matrix.
#'   Or a list of three elements if gene_anno is NOT provided: ex_f, matrix with filted probes (20% quantile);
#'   gset, downloaded GEO data; ex_all, unfiltered matrix.
#' @export
get_geo <- function(x, GPL, gene_anno = NULL) {
  ravg <- NULL

  #- Codes from GEO2R.
  if (!("GEOquery" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package GEOquery.")
  if (!("Biobase" %in% utils::installed.packages()[, "Package"])) stop("(EE) Require the package Biobase.")
  gset <- GEOquery::getGEO(x, GSEMatrix = TRUE, getGPL = FALSE)
  if (length(gset) > 1) idx <- grep(GPL, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]

  ex <- Biobase::exprs(gset)
  qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
  LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

  if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
  }

  #- Filter the data by 20% quantile.
  ex_f <- ex[rowMeans(ex) >= quantile(ex, 0.2), ]

  if (is.null(gene_anno)) {
    return(list(ex_f = ex_f, gset = gset, ex_all = ex))
  } else {
    ex_anno <- as.data.table(ex_f, keep.rownames = TRUE) %>%
      merge(gene_anno, ., by.x = "probe_id", by.y = "rn")

    cc_1 <- copy(ex_anno) %>%
    .[, ravg := rowMeans(.SD), .SDcols = -c("probe_id", "symbol", "gene_id")]
    ex_g <- cc_1[cc_1[, .I[(which.max(ravg))], by = "symbol"]$V1] %>%
      .[, ravg := NULL]
    return(list(ex_f = ex_f, ex_anno = ex_anno, ex_g = ex_g, gset = gset, ex_all = ex))
  }
}


#' Collect and merge non-synonymous mutations and CNA from Depmap and cosmic for cell lines.
#'
#' Get the annoations of cell lines relying on the Depmap IDs.
#'
#' @param g, a vector of interested genes. !! First check if any symbols missed in the Depmap expression data.
#' @param sample_ano, a data.table MUST with columns of cell_line and ModelID. If ModelID is not available, then it is "".
#' @param dp_cm, cosmic id <--> Depmap id.
#' @param dp_expr, Depmap expression data.
#' @param MUT_CNA, whether extracting MUT and CNA data from Depmap and Cosmic.
#' @param dp_mut, Depmap mutations data.
#' @param dp_cna, Depmap CNA data.
#' @param cm_mut, Cosmic mutations data.
#' @param cm_cna, Cosmic CNA data.
#' @return A list of annotations.
#' @examples
#' \dontrun{
#' 'load("Depmap/24Q2/data.DepMap.24Q2.rda")
#' '#- Get the symbols for CNA and expression.
#' 'colnames(cna_tab)  <- gsub("\\s.*$", "", colnames(cna_tab), perl = TRUE)
#' 'colnames(expr_tab) <- gsub("\\s.*$", "", colnames(expr_tab), perl = TRUE)
#' '#- COSMIC data.
#' 'cm_mut <- qs::qread("CellLinesProject_GenomeScreensMutant_v100_GRCh38.qs", nthreads = 4)
#' 'cm_cna <- readRDS("CellLinesProject_CompleteCNA_v100_GRCh38.rds")
#' 'dp_cm  <- readRDS("CellLinesProject_Sample_v100_GRCh38_depmapID.rds")
#' '
#' '#- Mutations and CNA
#' 'g <- c("KRAS", "TP53", "ALK", "BRAF", "CDKN2A", )
#' 'res_1 <- fn_cosmic_mut_cna(g, ano, dp_cm, expr_tab, MUT_CNA = TRUE,
#'   dp_mut_tab, cna_tab, cm_mut, cm_cna)
#'
#' '#- Expression for the CBM genes.
#' 'g <- c("KRAS", "EGFR", "CARD9", "CARD10")
#' 'res_2 <- fn_cosmic_mut_cna(g, ano, dp_cm, expr_tab, MUT_CNA = FALSE)
#' }
#' @export
depmap_cosmic_mut_cna <- function(g,
    sample_ano,
    dp_cm,
    dp_expr,
    MUT_CNA = FALSE,
    dp_mut,
    dp_cna,
    cm_mut,
    cm_cna) {
  #>>#########################################################
  #- !! Expression data.
  #- Remove cells without ModelID.
  cell_line <- ModelID <- Hugo_Symbol <- Chromosome <- Protein_Change <- Variant_Classification <- gr <- value <- NULL
  MUTATION_DESCRIPTION <- MUTATION_AA <- GENE_SYMBOL <- COSMIC_SAMPLE_ID <- depmap_ModelID <- MUT_TYPE <- ii <- jj <- NULL

  exp_tab_2 <- dp_expr[, intersect(g, colnames(dp_expr)), drop = FALSE] %>%
    as.data.table(keep.rownames = TRUE) %>%
    merge(sample_ano[, .(cell_line, ModelID)], ., by.x = "ModelID", by.y = "rn") %>%
    .[sample_ano$cell_line, on = "cell_line"] %>%
    .[, ModelID := NULL] %>%
    setcolorder("cell_line") %>%
    setnames(1, "cell") %>%
    setnames(names(.)[-1], paste0(names(.)[-1], "_log2TPM"))

  #- Round the data.
  sdcol <- grep("_log2TPM", names(exp_tab_2), value = TRUE)
  exp_tab_2[, (sdcol) := lapply(.SD, \(x) round(x, 2)), .SDcols = sdcol]

  #- TPM value.
  mtx <- exp_tab_2[, -"cell"] %>%
    as.matrix() %>%
    {2**. - 1} %>%
    round(2) %>%
    set_colnames(gsub("_log2", "_", colnames(.)))

  #- Combine the log2TPM and TPM
  exp_tab_2 <- cbind(exp_tab_2, as.data.table(mtx))

  #<<#########################################################


  #>>#########################################################
  #- !! Mutations and CNA in Depmap.
  if (MUT_CNA) {
    #- Find ModelID without Mut data.
    dp_mut_no_cell <- dp_mut[!(dp_mut$Variant_Classification %in% c("3'Flank", "5'Flank", "5'UTR", "RNA")), ] %>%
      as.data.table %>%
      setdiff(sample_ano$ModelID[sample_ano$ModelID != ""], .$ModelID)

    if (length(dp_mut_no_cell) > 0) message("ModelID without Mut data:", paste(dp_mut_no_cell, callapse = ";"))

    #- !! Mutations Depmap
    #- Remove the synonymous mutations
    dp_mut_2 <- dp_mut[!(dp_mut$Variant_Classification %in% c("3'Flank", "5'Flank", "5'UTR", "RNA")), ] %>%
      as.data.table %>%
      .[Hugo_Symbol %in% g] %>%
      merge(sample_ano[, .(cell_line, ModelID)], by = "ModelID") %>% #- Keep wanted cells.
      .[, .(ModelID, Chromosome, Hugo_Symbol, Protein_Change, Variant_Classification, cell_line)] %>%
      unique %>%
      .[, Protein_Change := gsub("Ter", "*", Protein_Change)] #- Align with COSMIC for termination, which uses Ter for termination.

    #- wide.
    dp_mut_2_w <- dcast(dp_mut_2, cell_line ~ Hugo_Symbol, value.var = "Protein_Change", fun.aggregate = \(x) paste(x, collapse = ";")) %>%
      setnames(1, "cell") %>%
      setnames(names(.)[-1], paste0(names(.)[-1], "_Mut")) %>%
      .[sample_ano$cell_line, on = "cell"] #- To include all cells

    #- Now treat them as WT, later if ModelID is na, replace with NA.
    #- If those genes are all WT in a cell line, the above "on = "cell"" make NA as well.
    dp_mut_2_w[is.na(dp_mut_2_w)] <- ""

    #- !! CNA Depmap
    #- Mina methods for CNA categories.
    dp_cna_no_cell <- setdiff(sample_ano$ModelID[sample_ano$ModelID != ""], rownames(dp_cna))

    if (length(dp_cna_no_cell) > 0) message("ModelID without CNA data:", paste(dp_cna_no_cell, callapse = ";"))

    dp_cna_2 <- dp_cna[, intersect(g, colnames(dp_cna)), drop = FALSE] %>%
      as.data.table(keep.rownames = TRUE) %>%
      melt(id.vars = "rn") %>%
      merge(sample_ano[, .(cell_line, ModelID)], by.x = "rn", by.y = "ModelID") %>%
      .[, gr := fcase(value < 0.87/2, "DeepDel",
                      value >= 0.87/2 & value < 1.32/2, "ShallowDel",
                      value > 2.64/2 & value <= 3.36/2, "Gain",
                      value > 3.36/2, "AMP",
                      default = "")]

    dp_cna_2_w <- dcast(dp_cna_2, cell_line ~ variable, value.var = "gr") %>%
      setnames(1, "cell") %>%
      setnames(names(.)[-1], paste0(names(.)[-1], "_CNA")) %>%
      .[sample_ano$cell_line, on = "cell"]

    dp_cna_2_w[is.na(dp_cna_2_w)] <- ""

    stopifnot(identical(dp_mut_2_w$cell, dp_cna_2_w$cell))
    stopifnot(identical(exp_tab_2$cell, dp_cna_2_w$cell))

    # dp_ano <- cbind(dp_mut_2_w, dp_cna_2_w[, -"cell"]) %>%
    dp_ano <- Reduce(cbind, list(dp_mut_2_w, dp_cna_2_w[, -"cell"], exp_tab_2[, -"cell"])) %>%
      cbind(ModelID = sample_ano[, ModelID], .)
    #<<#########################################################

    #>>#########################################################
    #- !! COSMIC mutations and CNA.
    #- !! remove synonymous_variant, "p.?"
    # 187324927     c.1_*del         p.0
    #- About 10 mutations in intron_variant leads to protein seq changes.
    #- 5_prime_UTR_variant has protein change muations as well.
    cm_mut_2 <- cm_mut[MUTATION_DESCRIPTION != "synonymous_variant"] %>%
      .[MUTATION_AA != "p.?"] %>%
      .[GENE_SYMBOL %in% g] %>%
      merge(dp_cm[, .(COSMIC_SAMPLE_ID, depmap_ModelID)], by.x = "COSMIC_SAMPLE_ID", by.y = "COSMIC_SAMPLE_ID") %>% #- Add depmap ids.
      merge(sample_ano[, .(cell_line, ModelID)], by.x = "depmap_ModelID", by.y = "ModelID") %>%
      .[, .(depmap_ModelID, COSMIC_SAMPLE_ID, GENE_SYMBOL, MUTATION_AA, cell_line)] %>%
      unique

    cm_mut_2_w <- dcast(cm_mut_2[, .(cell_line, GENE_SYMBOL, MUTATION_AA)] %>% unique, cell_line ~ GENE_SYMBOL, value.var = "MUTATION_AA", fun.aggregate = \(x) paste(x, collapse = ";")) %>%
      setnames(1, "cell") %>%
      setnames(names(.)[-1], paste0(names(.)[-1], "_Mut")) %>%
      .[sample_ano$cell_line, on = "cell"] #- To include all cells

    cm_mut_2_w[is.na(cm_mut_2_w)] <- ""

    cm_cna_2 <- cm_cna[GENE_SYMBOL %in% g] %>%
      .[, gr := MUT_TYPE] %>%
      .[MUT_TYPE == "gain", gr := "AMP"] %>%
      .[MUT_TYPE == "loss", gr := "DeepDel"] %>%
      merge(dp_cm[, .(COSMIC_SAMPLE_ID, depmap_ModelID)], by.x = "COSMIC_SAMPLE_ID", by.y = "COSMIC_SAMPLE_ID") %>%
      merge(sample_ano[, .(cell_line, ModelID)], by.x = "depmap_ModelID", by.y = "ModelID")

    #- !! In rare case there are multiple CNA events
    #    depmap_ModelID COSMIC_SAMPLE_ID COSMIC_CNV_ID COSMIC_GENE_ID GENE_SYMBOL
    #            <char>           <char>        <char>         <char>      <char>
    # 1:     ACH-000496       COSS724868 COSCNV6310492      COSG68225        XPO1
    # 2:     ACH-000496       COSS724868 COSCNV6310746      COSG68225        XPO1
    #    SAMPLE_NAME COSMIC_PHENOTYPE_ID TOTAL_CN MINOR_ALLELE MUT_TYPE
    #         <char>              <char>    <int>        <int>   <char>
    # 1:   NCI-H1792        COSO29915674        5            1     gain
    # 2:   NCI-H1792        COSO29915674        6            1     gain
    #    COSMIC_STUDY_ID CHROMOSOME GENOME_START GENOME_STOP     gr cell_line
    #             <char>     <char>        <int>       <int> <char>    <char>
    # 1:         COSU619          2     61494856    65564634    AMP  NCIH1792
    # 2:         COSU619          2     57777333    61494483    AMP  NCIH1792

    #- NO CNA from COSMIC
    if (nrow(cm_cna_2) != 0) {
      cm_cna_2_w <- dcast(cm_cna_2[, .(cell_line, GENE_SYMBOL, gr)] %>% unique, cell_line ~ GENE_SYMBOL, value.var = "gr") %>%
        setnames(1, "cell") %>%
        setnames(names(.)[-1], paste0(names(.)[-1], "_CNA")) %>%
        .[sample_ano$cell_line, on = "cell"]

      cm_cna_2_w[is.na(cm_cna_2_w)] <- ""
    } else {
      cc_1 <- paste0(g, "_CNA")
      cm_cna_2_w <- data.table(cell = sample_ano$cell_line) %>%
        .[, (cc_1) := ""]
    }

    #<<#########################################################

    # identical(dp$cell, cm$cell)
    stopifnot(identical(dp_ano$cell, cm_mut_2_w$cell))
    stopifnot(identical(dp_ano$cell, cm_cna_2_w$cell))

    #>>#########################################################
    #- Merge mutation data.
    #- overlapping mutation genes
    o_mut <- intersect(grep("_Mut", names(dp_ano), value = TRUE), grep("_Mut", names(cm_mut_2_w), value = TRUE))
    comb_mut <- sapply(o_mut, \(x) {
                    cc_dp <- dp_ano[[x]]
                    cc_cm <- cm_mut_2_w[[x]]

                    foreach(ii = cc_dp, jj = cc_cm, .combine = c) %do% {
                      if (ii == "" & jj == "") {
                        ""
                      } else if (ii != "" & jj == "") {
                        ii
                      } else if (ii == "" & jj != "") {
                        jj
                      } else if (ii != "" & jj != "") {
                        paste(unique(c(unlist(strsplit(ii, ";")), unlist(strsplit(jj, ";")))), collapse = ";")
                      }
                    }
    }) %>% as.data.table
    #<<#########################################################

    #>>#########################################################
    #- Merge CNA data.
    o_cna <- intersect(grep("_CNA", names(dp_ano), value = TRUE), grep("_CNA", names(cm_cna_2_w), value = TRUE))
    comb_cna <- sapply(o_cna, \(x) {
                    cc_dp <- dp_ano[[x]]
                    cc_cm <- cm_cna_2_w[[x]]

                    foreach(ii = cc_dp, jj = cc_cm, .combine = c) %do% {
                      if (ii == "" & jj == "") {
                        ""
                      } else if (ii != "" & jj == "") {
                        ii
                      } else if (ii == "" & jj != "") {
                        jj
                      } else if (ii != "" & jj != "") {
                        paste(unique(c(unlist(strsplit(ii, ";")), unlist(strsplit(jj, ";")))), collapse = ";")
                      }
                    }
    }) %>% as.data.table
    #<<#########################################################

  #- cell, mut from dp, mut from cm, cna from dp
    comb_res <- cbind(dp_ano[, c("ModelID", "cell"), with = FALSE],
                      comb_mut, #- Combined mutations.
                      dp_ano[, setdiff(grep("_Mut", names(dp_ano), value = TRUE), o_mut), with = FALSE], # mutatons in depmap
                      cm_mut_2_w[, setdiff(grep("_Mut", names(cm_mut_2_w), value = TRUE), o_mut), with = FALSE], # mutations in COSMIC
                      comb_cna, #- combined CNA
                      dp_ano[, setdiff(grep("_CNA", names(dp_ano), value = TRUE), o_cna), with = FALSE], # CNA in depmap
                      cm_cna_2_w[, setdiff(grep("_CNA", names(cm_cna_2_w), value = TRUE), o_cna), with = FALSE], # CNA in COSMIC
                      dp_ano[, grep("TPM", names(dp_ano), value = TRUE), with = FALSE])

    #- After combined with COSMIC data, if ModelID is empty, change to N/A
    sdcol_1 <- grep("_Mut|_CNA", names(comb_res), value = TRUE)
    comb_res[ModelID == "", (sdcol_1) := lapply(.SD, \(x) x <- "N/A"), .SDcols = sdcol_1]

    #- ModelID is available but no Mut or CNA, e.g., only screen data in Depmap.
    sdcol_2 <- grep("_Mut", names(comb_res), value = TRUE)
    #- If the annotation is empty, change to N/A. If annotated by COSMIC, wont be updated.
    comb_res[ModelID %in% dp_mut_no_cell, (sdcol_2) := lapply(.SD, \(x) x[x == ""] <- "N/A"), .SDcols = sdcol_2]

    sdcol_3 <- grep("_CNA", names(comb_res), value = TRUE)
    #- If the annotation is empty, change to N/A. If annotated by COSMIC, wont be updated.
    comb_res[ModelID %in% dp_cna_no_cell, (sdcol_3) := lapply(.SD, \(x) x[x == ""] <- "N/A"), .SDcols = sdcol_3]

    return(list(final_results = comb_res,
                depmap_ano    = dp_ano,
                depmap_epxr   = exp_tab_2,
                depmap_mut    = dp_mut_2,
                depmap_mut_w  = dp_mut_2_w,
                depmap_cna    = dp_cna_2,
                depmap_cna_w  = dp_cna_2_w,
                cosmic_mut    = cm_cna_2,
                cosmic_mut_w  = cm_cna_2_w,
                cosmic_cna    = cm_cna_2,
                cosmic_cna_w  = cm_cna_2_w))
  } else {
    return(exp_tab_2)
  }
}
