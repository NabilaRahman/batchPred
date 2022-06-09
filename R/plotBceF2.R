#' Plot Receiver-Operator curves for evaluating the performance of batch correction
#' 
#' \code{refNet} Plots ROC curve using results from \code{computeBceF2} 
#' 
#' @param compute_bcef_object A list of all variables needed to plot ROC curve
#' @param line_color Character vector of colors (Optional). Must match number of curves plotted. 
#' @param plot_title Title of the plot (Optional).
#' @param showLegendBorder If true, draw border around legend. Default=FALSE
#' @param title_size Numerical value specifying size of plot title.
#' @return \code{plotBCeF2} Plots ROC curves of raw and batch corrected datasets.
#'  Allows you to compare the performance of different sets of covariates. 
#' @export

plotBceF2 <- function (compute_bcef_object
                       , line_color = compute_bcef_object[["line.color"]]
                       , plot_title = compute_bcef_object[["plot.title"]]
                       , showLegendBorder = F
                       , title_size=1.5
                       )
{
  if (length(compute_bcef_object[["line.color"]]) != length(line_color)) {
    "Too many/Too few colors specified for ROC curve"
    return(0)
  }
  roc_obj_linear <- compute_bcef_object[["roc_obj_linear"]]
  legend.to.use <- compute_bcef_object[["legend.to.use"]]
  l.auc <- compute_bcef_object[["l.auc"]]
  pvals.df.log <- compute_bcef_object[["pvals.df.log"]]
  gold.standard.bool <- compute_bcef_object[["gold.standard.bool"]]
  plot(roc_obj_linear, col = line_color[1], cex.main = 0.8,
       main = "")

  final.legend.to.use <- paste0(legend.to.use[1], " (auc=",
                                round(l.auc, 3), ")")
  for (i in 2:ncol(pvals.df.log)) {
    roc_obj_tmp = pROC::roc(gold.standard.bool, pvals.df.log[,
                                                             i], quiet = T)
    tmp.auc = pROC::auc(roc_obj_tmp)
    lines(roc_obj_tmp, col = line_color[i])
    final.legend.to.use = c(final.legend.to.use, paste0(legend.to.use[i],
                                                        " (auc=", round(tmp.auc, 3), ")"))
  }
  if (showLegendBorder) {
    legend("bottomright", cex = 1, legend = final.legend.to.use,
           col = line_color, lwd = 0.8, lty = c(1, 1, 1, 1))
  }
  else {
    legend("bottomright", cex = 1, legend = final.legend.to.use,
           col = line_color, lwd = 0.8, lty = c(1, 1, 1, 1),
           bty = "n")
  }
  title(plot_title, adj=0, cex.main=title_size, line=2.5)
}
