#' Generates values for roc curve for evaluating the performance of batch correction
#' 
#' This package corrects for confounders in gene expression datasets using 
#'  multiple linear regression model and then evaluates the improvement in gene 
#'  coexpression of using a gold standard co-expression network. 
#' @param input.edata Matrix of raw expression dataset to be adjusted. Rows 
#'  represent the genes/probes and columns represent samples.
#' @param covariates.df.list List of covariate dataframes. Each dataframe in the 
#'  list consists of columns representing known covariates e.g. Age, Sex, or  
#'  unknown covariates such as principle components. Rows represent samples.
#' @param input.gold.standard Dataframe. A gold standard that includes a gene 
#'  coexpression confidence. First column represents the first gene, the second 
#'  column the second gene and the third columns should be a binary vector where 
#'  1 indicates true associations and 0 indicates false associations. The gene 
#'  IDs in the first two column of \code{input.gold.standard} and the rownames of
#'  \code{input.edata} must be of the same type.
#' @param plot.title Title of the plot.
#' @param roc.curve.legend Character vector of descriptions for each set of 
#'  covariates to plot. Must be in the same order as the dataframes in 
#'  \code{covariates.df.list}. This modifies the roc curve legend.
#' @param line.color Character vector specifying color of roc curves. By default 
#'  sets color for upto 6 roc curves with raw displayed as black. Manually 
#'  specify line.color if plotting more than 6 different sets of covariates, 
#'  i.e. over 6 dataframes in \code{covariates.df.list}.
#' @return \code{computeBCeF2} Generates a list of all variables needed for plotting roc curves of raw and batch corrected datasets.
#'  One batch corrected roc curve is computed for each dataframe specified in 
#'  \code{covariates.df.list}, allowing you to compare the performance of 
#'  different sets of covariates. Use \code{plotBCeF2} to plot these variables
#' @export

computeBCeF2<- function (input.edata, covariates.df.list, input.gold.standard,
          plot.title = "Batch Correction Evaluation", roc.curve.legend = NULL,
          line.color=NULL)
{
  if (!require("pROC")) {
    print("error: pROC package was not installed")
    return(0)
  }
  if (ncol(input.edata) != nrow(input.covariates.df)) {
    print("error: covariates and input.edata has different number of samples")
    return(0)
  }
  
  # Input Data correction
  somecolors<- c("black", "red", "blue", "green","pink","cyan", "magenta")
  if (is.null(line.color) ) {
    line.color <- somecolors[1:(1+length(covariates.df.list))] 
  } else if ( length(covariates.df.list) > 6 )  {
      message("ERROR: Not enough colors specified. Please specify colors using the 'line.color' parameter" )
      stop()
    }
  
  rownames(input.edata) <- as.character(rownames(input.edata))
  colnames(input.gold.standard) = c("Gene1", "Gene2", "Confidence")
  input.gold.standard$Gene1 <- as.character(input.gold.standard$Gene1)
  input.gold.standard$Gene2 <- as.character(input.gold.standard$Gene2)
  input.gold.standard <- as.data.frame(input.gold.standard)
  
  raw.edata = as.matrix(input.edata)
  t.raw.edata = t(raw.edata)
  if (is.null(names(covariates.df.list)) ) { names(covariates.df.list) <- c( 1:length(covariates.df.list) )}
  if (class(covariates.df.list)!="list") { covariates.df.list <- list("Adjusted"=covariates.df.list) }
  
  bins.cors.df<-data.frame()
  all.mod.name<-vector()
  for (i in 0:length(names(covariates.df.list)) ) {
    if ( i > 0 ) {
      subset.covariates.df <- covariates.df.list[[i]]
      factors = sapply(1:ncol(subset.covariates.df), function(x) paste0("subset.covariates.df[,", x, "]"))
      my.formula = reformulate(termlabels = factors, response = "t.raw.edata")
      lm.fitted.edata = stats::lm(my.formula)
      adjusted.reads = t(lm.fitted.edata$residuals)
      mod.name<-names(covariates.df.list)[i]
    } else {
      adjusted.reads = raw.edata
      mod.name <-"Raw"
    }
    
    gold.standard.to.use = input.gold.standard
    rownames(gold.standard.to.use) = c(1:nrow(gold.standard.to.use))
    bins.cors.tmp = gold.standard.to.use
    bins.cors.tmp$binAll.pval = NA
    for (j in 1:nrow(gold.standard.to.use)) {
      geneA.geneID = gold.standard.to.use$Gene1[j]
      geneB.geneID = gold.standard.to.use$Gene2[j]
      if ((geneA.geneID %in% rownames(raw.edata)) & (geneB.geneID %in%
                                                     rownames(raw.edata))) {
        #Calculate Corr of Gene i j of Adjusted Dataset
        geneA.bin.all.vec = geneB.bin.all.vec = 0
        geneA.bin.all.vec = adjusted.reads[geneA.geneID,]
        geneB.bin.all.vec = adjusted.reads[geneB.geneID,]
        result.all = stats::cor.test(geneA.bin.all.vec, geneB.bin.all.vec,
                                     method = "spearman")
        bins.cors.tmp[j, "binAll.pval"] = result.all$p.value
      }
      else {
        print(paste0(j, " is index of not found gene in input.edata"))
        j = j + 1
      }
    }
    colnames(bins.cors.tmp)[4] <- mod.name
    if (i > 0) {
      bins.cors.tmp[,1:3] <- NULL
      bins.cors.df <- cbind(bins.cors.df, bins.cors.tmp)
      all.mod.name<-c(all.mod.name, mod.name)
    } else { 
      bins.cors.df <- bins.cors.tmp  
      }
  }
  
  if (is.null(roc.curve.legend)) {roc.curve.legend <-all.mod.name}
  
  pvals.df = bins.cors.df[!is.na(bins.cors.df$Raw),  ]
  rownames(pvals.df) = c(1:nrow(pvals.df))
  gold.standard.bool = pvals.df[, "Confidence"]
  pvals.df = pvals.df[,c(4:ncol(pvals.df))]
  pvals.df = apply(pvals.df, 2, function(x) p.adjust(x, method = "BH"))
  dataset.col.names <- c("Raw",roc.curve.legend)
  colnames(pvals.df) = dataset.col.names
  pvals.df.log = 0
  pvals.df.log = -log(pvals.df, 10)
  tmp = pvals.df.log
  tmp[is.infinite(tmp)] = 0
  max.val = max(tmp)
  pvals.df.log[is.infinite(pvals.df.log)] <- (max.val + 10)
  legend.to.use = dataset.col.names
  cases.gold = length(which(gold.standard.bool == 1)) 
  control.gold = length(which(gold.standard.bool == 0)) 
  roc_obj_linear = pROC::roc(gold.standard.bool, pvals.df.log[,"Raw"])
  l.auc = pROC::auc(roc_obj_linear)
  auc.vec <- l.auc
  
  # Export relevant important factors
  out <- list()
  out[["roc_obj_linear"]] <- roc_obj_linear
  out[["line.color"]] <- line.color
  out[["plot.title"]] <- plot.title
  out[["legend.to.use"]] <- legend.to.use
  out[["l.auc"]] <- l.auc
  out[["pvals.df.log"]] <- pvals.df.log
  out[["gold.standard.bool"]] <- gold.standard.bool
  return(out)
                   
  # #Plot ROC curves
  # plot(roc_obj_linear
  #      , col = line.color[1]
  #      , cex.main = 0.8
  #      , main=paste0(plot.title)
  #      )
  # final.legend.to.use <- paste0(legend.to.use[1]
  #                              , " (auc=", round(l.auc, 3), ")")
  # for (i in 2:ncol(pvals.df.log)) {
  #   roc_obj_tmp = pROC::roc(gold.standard.bool, pvals.df.log[,i], quiet = T)
  #   tmp.auc = pROC::auc(roc_obj_tmp)
  #   lines(roc_obj_tmp, col = line.color[i])
  #   final.legend.to.use = c(final.legend.to.use, paste0(legend.to.use[i],
  #                                                       " (auc=", round(tmp.auc, 3), ")"))
  # }
  # legend("bottomright", cex = 1, legend = final.legend.to.use,
  #        col = line.color, lwd = 0.8, lty = c(1, 1, 1, 1))
}
