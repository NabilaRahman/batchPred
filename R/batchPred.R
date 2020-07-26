#' Predict relevant covariates for batch effect adjustment
#' 
#' \code{batchPred} predicts relevant batch variables for improving the 
#'  normalisation and the interpretation of heterogenous datasets.
#' 
#' @param input.edata matrix of numeric expression data, where rows are genes
#'  and columns are samples.
#' @param input.covariates.df the covariate dataframe. Each covariate is a 
#'  dataframe column that represents known covariate such as batch number, 
#'  age or hidden covariates such as a principle components.
#' @param input.gold.standard reference data table consisting of known gene 
#'  associations. The number of true positives and false positives should be 
#'  approximately equal. gene IDs in \code{input.gold.standard} must match 
#'  gene IDs in \code{input.edata}.
#' @param threshold numeric. Minimum AUC score improvement required for 
#'  addition to the linear design.  
#' @param cores integer. Number of cores / threads. Number of threads is the 
#' primary bottleneck for computing the initial covariate order.  
#' @return \code{batchPred} returns a table consisting of with the following 
#' columns: 
#' \itemize{
#'  \item Covariates: covariate, whose adjustment show greatest improvement in 
#'  AUC.  
#'  \item LinearModel: Combination of covariates tested. 
#'  \item AUC: Area under roc curve (AUC) value after batch adjustment with 
#'  \code{LinearModel}. 
#'  \item AUCvRaw: difference in AUC value of batch corrected and raw dataset. 
#'  \item covEffectOnAUC: improvement in AUC value compared to previously 
#'  tested \code{LinearModel}.
#'   }
#' @export

batchPred <- function (input.edata
                       , input.covariates.df
                       , input.gold.standard
                       , threshold=0.0001
                       , cores=2
){
  if (!require("pROC")) {
    MESSAGE("ERROR: pROC package was not installed")
    return(0)
  }
  if (ncol(input.edata) != nrow(input.covariates.df)) {
    message("ERROR: covariates and input.edata has different number of samples")
    return(0)
  }
  
  rownames(input.edata) <- as.character(rownames(input.edata))
  colnames(input.gold.standard) = c("Gene1", "Gene2", "Confidence")
  input.gold.standard$Gene1 <- as.character(input.gold.standard$Gene1)
  input.gold.standard$Gene2 <- as.character(input.gold.standard$Gene2)
  input.gold.standard <- as.data.frame(input.gold.standard)
  
  
  #Credit Judith Somekh (https://github.com/jsomekh/BCeF_)
  #Get AUC
  pBCeF <-function( subset.covariates.df, input.edata, input.gold.standard, raw=F ) { 
    raw.edata = as.matrix(input.edata)
    if ( raw == F ) {
      t.raw.edata = t(raw.edata)
      factors = sapply(1:ncol(subset.covariates.df), function(x) paste0("subset.covariates.df[,",
                                                                        x, "]"))
      my.formula = reformulate(termlabels = factors, response = "t.raw.edata")
      lm.fitted.edata = stats::lm(my.formula)
      adjusted.reads = t(lm.fitted.edata$residuals)
    } else {
      adjusted.reads = raw.edata
    }
    
    gold.standard.to.use = input.gold.standard
    rownames(gold.standard.to.use) = c(1:nrow(gold.standard.to.use))
    bins.cors.df = gold.standard.to.use
    bins.cors.df$binAll.pval = NA
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
        bins.cors.df[j, "binAll.pval"] = result.all$p.value
      }
      else {
        print(paste0(j, " is index of not found gene in input.edata"))
        j = j + 1
      }
    }
    pvals.df = bins.cors.df[!is.na(bins.cors.df$binAll.pval), ]  #Rm NA
    rownames(pvals.df) = c(1:nrow(pvals.df))       # Reset df's rowname index
    gold.standard.bool = pvals.df[, "Confidence"]  # Retain postProb in diff var
    pvals.df = pvals.df[, "binAll.pval"]  #Only keep pval columns
    pvals.df = p.adjust(pvals.df, method = "BH") # multiple correction for each column
    names(pvals.df) = "Adjusted"
    pvals.df.log = 0
    pvals.df.log = -log(pvals.df, 10)  #Take log 10 of p-values
    tmp = pvals.df.log  
    tmp[is.infinite(tmp)] = 0   #Set infinite values to 0
    max.val = max(tmp)    #Take most significant p-value
    pvals.df.log[is.infinite(pvals.df.log)] <- (max.val + 10)  # Set non finite val to maxV+10
    roc_obj_tmp = pROC::roc(gold.standard.bool, pvals.df.log, quiet = T) #Plot adjusted Curve
    roc_auc = pROC::auc(roc_obj_tmp)
    
    return(roc_auc)
  }
  
  func.start.time<-Sys.time()
  #Calculate Raw AUC
  raw_auc <-pBCeF( subset.covariates.df=NULL, input.edata, input.gold.standard, raw=T ) 
  
  #Initialise empty variables 
  auc_table = data.frame()
  
  
  
  # Sequentially add covariates and keep formula that increases the AUC
  input.covariates <- colnames(input.covariates.df)
  
  #setup parallel backend to use many processors
  doParallel::registerDoParallel(cores)
  
  tempMatrix <- data.frame()
  var.table <- data.frame( Covariate="Raw"
                           , LinearModel=NA
                           , AUC = raw_auc
                           , AUCvRaw = 0
                           , covEffectOnAUC = 0 )  
  variable_confirm <- vector()
  variable_test <- vector()
  score_confirm <- vector()
  nomatch = F
  
  # Loops until no more variables improve AUC 
  library(foreach)
  while (nomatch == F) { 
    if (length(input.covariates) > 0) {
      start.time <- Sys.time()
      
      finalMatrix <- data.frame()
      # Shuffling of variables
      shuffledIndex <- sample(length(input.covariates), replace=F)
      input.covariates <- input.covariates[shuffledIndex]
      
      #Divide tasks upto cores to not overload memory
      quot <- length(input.covariates) %/% cores
      remainder <- length(input.covariates) %% cores
      cat("Computing max of", length(input.covariates), "AUC scores...")
      if (quot > 0 ) {
        for (j in (1:quot)) {
          mini=((cores*j)-cores)+1
          maxi=cores*j
          tempMatrix <-  foreach::foreach(i=mini:maxi, .combine=rbind) %dopar% {
            variable_test <- c(variable_confirm, input.covariates[i])
            subset.covariates.df<-as.data.frame(input.covariates.df[,variable_test])
            auc_test <- pBCeF( subset.covariates.df, input.edata, input.gold.standard, raw=F )
            tempMat <- data.frame(Variable=input.covariates[i], Score = auc_test, stringsAsFactors = F)
            tempMat
            
          }
          finalMatrix <- rbind(finalMatrix, tempMatrix)
        }
      }
      if (remainder > 0) {
        mini=(cores*quot)+1
        maxi = length(input.covariates)
        tempMatrix <-  foreach::foreach(i=mini:maxi, .combine=rbind) %dopar% {
          variable_test <- c(variable_confirm, input.covariates[i])
          subset.covariates.df<-as.data.frame(input.covariates.df[,variable_test])
          auc_test <- pBCeF( subset.covariates.df, input.edata, input.gold.standard, raw=F )
          tempMat <- data.frame(Variable=input.covariates[i], Score = auc_test, stringsAsFactors = F)
          tempMat
        }
        
        finalMatrix <- rbind(finalMatrix, tempMatrix)
      }
      rm(tempMatrix)
      cat("done\n") 
      maxScore=max(finalMatrix$Score)
      chosenVar = finalMatrix$Variable[finalMatrix$Score == maxScore]
      # Check for Maximum score for initial variable
      if ( length(variable_confirm) == 0  & maxScore > raw_auc & maxScore-raw_auc > threshold  ) {
        variable_confirm <- c( variable_confirm, chosenVar )
        input.covariates <- input.covariates[!input.covariates %in% chosenVar]
        cat(chosenVar, "| AUC improvement:", round(maxScore-raw_auc, 7), "\n" )
        
        var.table <- rbind(var.table, data.frame(Covariate=chosenVar
                                                 , LinearModel=paste(variable_confirm, collapse=",")
                                                 , AUC = maxScore
                                                 , AUCvRaw = maxScore-raw_auc
                                                 , covEffectOnAUC = NA    )  )
        LastScore <- maxScore
        # Check for Maximum score for 2nd variable onwards
      } else if ( length(variable_confirm) > 0  & maxScore > LastScore & maxScore-LastScore > threshold ) { 
        variable_confirm <- c( variable_confirm, chosenVar ) 
        input.covariates <- input.covariates[!input.covariates %in% chosenVar]
        cat(chosenVar, "| AUC improvement:", round(maxScore-LastScore, 7), "\n")
        
        var.table <- rbind(var.table, data.frame(Covariate=chosenVar
                                                 , LinearModel=paste(variable_confirm, collapse=",")
                                                 , AUC = maxScore
                                                 , AUCvRaw = maxScore-raw_auc
                                                 , covEffectOnAUC = maxScore-LastScore)  )
        LastScore <- maxScore
      } else if (length(variable_confirm) == 0  &  maxScore < raw_auc ) {
        message("No Covariate improved AUC of Raw data. Try other covariates.","Exiting...")
        nomatch = T
      } else  {
        message("No other Covariate improves current AUC of: ", round(LastScore,6), ". Exiting.")
        nomatch = T
      }
      
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      cat( "Elapsed Time:", round(time.taken,2), "minutes.\n" )
    } else { 
      break 
    }
  }
  #Calc Total Elapsed Time
  end.time <- Sys.time()
  time.taken <- end.time - func.start.time
  
  message( "Total Elapsed Time:", round(time.taken,2), " minutes" )
  return(var.table)
}
