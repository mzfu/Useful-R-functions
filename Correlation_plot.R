#=======================================================================

# File: Function to Make a Correlation Matrix
# Author: Joy Fu
# Update Date: 07/24/2019

#=======================================================================

#=======================================================================
#   Quick Reference
#   Correlation Plot -- make_heatmap() to visualize correlation
#=======================================================================


#+++++++++++++++++++++++
# Helper Functions
#+++++++++++++++++++++++

# Get lower triangle of the correlation matrix
.get_lower_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[upper.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# Get upper triangle of the correlation matrix
.get_upper_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[lower.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# hc.order correlation matrix
.hc_cormat_order <- function(cormat, hc.method = "complete") {
  dd <- stats::as.dist((1 - cormat) / 2)
  hc <- stats::hclust(dd, method = hc.method)
  hc$order
}

.no_panel <- function() {
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )
}

# Convert a tbl to matrix
.tibble_to_matrix <- function(x){
  x <-  as.data.frame(x)
  rownames(x) <- x[, 1]
  x <- x[, -1]
  as.matrix(x)
}

# Start self-defined plotting function
mycorrplot = function (corr, method = c("square", "circle"), 
                       type = c("full", "lower", "upper"), ggtheme = ggplot2::theme_minimal, title = "", 
                        show.legend = TRUE, legend.title = "Corr", show.diag = FALSE, 
                        colors = c("blue", "white", "red"), outline.color = "gray", 
                        hc.order = FALSE, hc.method = "complete", lab = FALSE, lab_col = "black", 
                        lab_size = 4, p.mat = NULL, sig.level = 0.05, 
                        significant = c("pch", "blank"), pch = 1, pch.col = "black", pch.cex = 10, tl.cex = 12, 
                        tl.col = "black", tl.srt = 45, digits = 2) 
{
  type <- match.arg(type)
  method <- match.arg(method)
  significant <- match.arg(significant)
  if (!is.matrix(corr) & !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  corr <- as.matrix(corr)
  corr <- base::round(x = corr, digits = digits)
  if (hc.order) {
    ord <- .hc_cormat_order(corr)
    corr <- corr[ord, ord]
    if (!is.null(p.mat)) {
      p.mat <- p.mat[ord, ord]
      p.mat <- base::round(x = p.mat, digits = digits)
    }
  }
  if (type == "lower") {
    corr <- .get_lower_tri(corr, show.diag)
    p.mat <- .get_lower_tri(p.mat, show.diag)
  }
  else if (type == "upper") {
    corr <- .get_upper_tri(corr, show.diag)
    p.mat <- .get_upper_tri(p.mat, show.diag)
  }
  corr <- reshape2::melt(corr, na.rm = TRUE)
  colnames(corr) <- c("Var1", "Var2", "value")
  corr$pvalue <- rep(NA, nrow(corr))
  corr$signif <- rep(NA, nrow(corr))
  if (!is.null(p.mat)) {
    p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
    corr$coef <- corr$value
    corr$pvalue <- p.mat$value
    corr$signif <- as.numeric(p.mat$value <= sig.level)
    p.mat <- subset(p.mat, p.mat$value <= sig.level)
    if (significant == "blank") {
      corr$value <- corr$value * corr$signif
    }
  }
  corr$abs_corr <- abs(corr$value) * 10
  p <- ggplot2::ggplot(corr, ggplot2::aes_string("Var1", "Var2", 
                                                 fill = "value"))
  if (method == "square") {
    p <- p + ggplot2::geom_tile(color = outline.color)
  }
  else if (method == "circle") {
    p <- p + ggplot2::geom_point(color = outline.color, shape = 21, 
                                 ggplot2::aes_string(size = "abs_corr")) + 
                                 ggplot2::scale_size(range = c(4, 10)) + ggplot2::guides(size = FALSE)
  }
  p <- p + ggplot2::scale_fill_gradient2(low = colors[1], high = colors[3], 
                                         mid = colors[2], midpoint = 0, limit = c(-1, 1), space = "Lab", 
                                         name = legend.title)
  if (class(ggtheme)[[1]] == "function") {
    p <- p + ggtheme()
  }
  else if (class(ggtheme)[[1]] == "theme") {
    p <- p + ggtheme
  }
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = tl.srt, 
                                                              vjust = 1, size = tl.cex, hjust = 1), axis.text.y = ggplot2::element_text(size = tl.cex)) + 
    ggplot2::coord_fixed()
  label <- round(corr[, "value"], 2)
  if (lab) {
    p <- p + ggplot2::geom_text(ggplot2::aes_string("Var1", 
                                                    "Var2"), label = label, color = lab_col, size = lab_size, family = "serif")
  }
  if (!is.null(p.mat) & significant == "pch") {
    p <- p + ggplot2::geom_point(data = p.mat, ggplot2::aes_string("Var1", 
                                                                   "Var2"), shape = pch, size = pch.cex, color = pch.col)
  }
  if (title != "") {
    p <- p + ggplot2::ggtitle(title)
  }
  if (!show.legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p <- p + .no_panel()
  p
}

# Use this function for making the heatmap
make_heatmap = function(dataset, lab_size = 3, fontsize = 12, save_name, output_path) {
  
  cor_data = round(cor(dataset, use = "complete.obs"), 2)
  p.mat = cor_pmat(dataset)
  
  cor_plot = mycorrplot(cor_data, type = "lower", p.mat = p.mat, lab_size = lab_size, lab = T, significant  = "pch", 
                        legend.title = 'Pearson\nCorrelation', outline.color = "white", tl.cex = 11,
                        colors = c("#6D9EC1", "white", "red")) + 
    theme(
      axis.text.x = element_text(size = fontsize, family = "serif"),
      axis.text.y = element_text(size = fontsize, family = "serif"),
      legend.text = element_text(size = fontsize, family = "serif"),
      legend.title = element_text(size = fontsize, family = "serif"),
      text = element_text(size = fontsize, family = "serif")
    ) 
  ggsave(save_name, plot = cor_plot, path = output_path, width = 8, height = 8)
  
  return(cor_plot)
}
