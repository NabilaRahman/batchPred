#' Generates non-overlapping subsets from a reference gene co-expression network
#' 
#' \code{refNetSubsets} generates non-overlapping subsets from a 
#' reference gene association network.
#' 
#' @param gold_standard_filt Data table of gene-gene network
#'  associations. \code{gold_standard_filt} is basically the output of 
#'  \code{createRefGeneNet}. Column 1 and 2 are gene IDs. Column 3 is 
#'  posterior probability with known edges set to 1. While column 4 defines 
#'  whether a gene relation is a true positive (1) or false positive (0).
#' @param seed Integer. Seed allows reproduction of a sequence of random 
#'  numbers. Important for reproducibly subsetting the reference network.
#' @param number_of_subsets Integer. Specify the number of reference subsets 
#'  to create. It is recommended to run this function once to generate all 
#'  your subsets at once to avoid gene-gene relation overlap between subsets.
#' @param sample_size Integer. Specifies the total number of gene relations 
#'  to extract from \code{gold_standard_filt}. Final number of gene relations
#'  will typically be smaller than this number to ensure a similar
#'  number of true and false positives are generated. Default: 400000
#' @param false_positive cutoff value for false positive gene associations.
#'  Should be the same cutoff used in \code{createRefGeneNet} or lower. Used 
#'  as a starting point for adjusting the number of false positives and 
#'  negatives. Default: 0.025
#' @param true_positive cutoff value for true positive gene associations.
#'  Should be the same cutoff used in \code{createRefGeneNet} or higher. Used
#'  as a starting point for adjusting the number of false positives and 
#'  negatives. Default: 0.5
#' @param threshold Maximum difference allowed between number of true and false 
#' positive. You should aim to have similar number of false and true positive 
#' gene-gene interactions. Set this to NULL calculate an appropriate number 
#' automatically.
#' 
#'  Should be the same cutoff used in \code{createRefGeneNet} or higher. Used
#'  as a starting point for adjusting the number of false positives and 
#'  negatives. Default: 0.5
#' @return A data table or a list of data tables. Generates non-overlapping 
#'  subsets from a reference gene association network and calibrates the 
#'  number of true positives and false positives to be approximately equal.
#' @examples
#' gold_standard_filt <- data.table::data.table(geneID1 = c(6776, 5359,
#'     9039, 10651, 6627, 3950, 1998, 222950, 54732, 92745, 
#'     1974, 30814, 29982, 4967, 10234, 9742, 57484, 5818, 
#'     5590 , 92211)
#'  , geneID2 = c(7318, 3601, 1027, 4190, 6427, 22797, 9215, 
#'     83642, 10388, 1288, 474, 5973, 1735, 374454, 2886, 
#'     1933, 79029, 6666, 7307, 51285)
#'  , postProbability = c(0.62357, 0.784499, 0.535758, 0.824298, 0.773425, 
#'     0.0131272, 0.00365631, 0.0205155, 0.0201044, 0.0212581, 0.0045336, 
#'     0.0115648, 0.00168017, 0.0164424, 0.00369498, 0.00608776, 0.0177325, 
#'     0.0152875, 0.0196584, 0.0222284)
#'  , Confidence = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#'  )
#'	
#' subset_list <- refNetSubsets( gold_standard_filt
#'  , seed = 50
#'  , number_of_subsets = 3
#'  , sample_size = 5
#'  , false_positive = 0.019
#'  , true_positive = 0.5
#'  ) 
#' subset_list[2]
#' @export

refNetSubsets <- function(
  gold_standard_filt
  , seed = 50
  , number_of_subsets = 3
  , sample_size = 400000
  , false_positive = 0.025
  , true_positive = 0.5
  , threshold = NULL
) {
  set.seed(seed)
  colnames(gold_standard_filt) <-c("geneID1", "geneID2", "postProbability", "Confidence")
  gold_standard_filt$geneID1<- as.character(gold_standard_filt$geneID1)
  gold_standard_filt$geneID2 <- as.character(gold_standard_filt$geneID2)
  if (number_of_subsets > 1) {calibratedSamples <- list() } 
  
  for (s.num in 1:number_of_subsets) {
    message(paste0("Sample ", s.num,": "))
    #Subset Rows
    rows <- sample(nrow(gold_standard_filt), size = sample_size, replace=F)
    sample_network_gold_standard <- gold_standard_filt[rows, ]
    sample_network_gold_standard<-sample_network_gold_standard[order(sample_network_gold_standard$postProbability, decreasing = T),]
    #Check that there is at least 1 true and 1 false postive
    if ( nrow(sample_network_gold_standard[sample_network_gold_standard$Confidence==0,]) >= 1 &  nrow(sample_network_gold_standard[sample_network_gold_standard$Confidence==1,]) >= 1) {
      
      #Balance Number of True and False Positives
      TF_summary  <- table(sample_network_gold_standard$Confidence)
      if (is.null(threshold)) {
      TF_opt <- ceiling(min(TF_summary)*0.1) 
      } else {
        TF_opt <- threshold
      }
      
      
      if ( abs(TF_summary[1] - TF_summary[2]) > TF_opt ) {
        TF_diff <- abs(TF_summary[1] - TF_summary[2])
        if ( TF_summary[1] >TF_summary[2] ) {
          type_calibrated = "false_positive"
          temp_summary  <- TF_summary
          temp_diff <- TF_diff
          calibration = false_positive-0.000005
          while ( temp_diff > TF_opt ) {
            sample_network_filt = sample_network_gold_standard[ sample_network_gold_standard$postProbability > calibration, ];
            temp_summary = table(sample_network_filt$Confidence)
            temp_diff = abs(temp_summary[1] - temp_summary[2])
            calibration = calibration - 0.000005
            if(is.na(temp_diff)) { 
              message("ERROR: Calibration generated NA values. Try a lower 'false_positive' threshold or larger sample_size")
              sample_network_filt <- NULL
              break }
          }
        } else if ( TF_summary[1] < TF_summary[2]) {
          type_calibrated = "false_positive"
          temp_summary  <- TF_summary
          temp_diff <- TF_diff
          calibration = true_positive + 0.000005
          while ( temp_diff > TF_opt ) {
            sample_network_filt = sample_network_gold_standard[ sample_network_gold_standard$postProbability > calibration, ];
            temp_summary = table(sample_network_filt$Confidence )
            temp_diff = abs(temp_summary[1] - temp_summary[2])
            calibration = calibration + 0.000005
            if(is.na(temp_diff)) { 
              message("ERROR: Calibration generated NA values. Try a higher 'true_positive' threshold or larger sample_size.")
              sample_network_filt <- NULL
              break }
          }
        }
        message(paste(type_calibrated, "calibration threshold set to:", calibration ))
      } else {
        sample_network_filt = sample_network_gold_standard
      }
      final_summary <- table(sample_network_filt$Confidence)
      message(paste("true_positives_relations:", final_summary[2], "\nfalse_positives_relations:", final_summary[1], "\n" ))
      gold_standard_filt <- gold_standard_filt[-rows, ]
      
      if (number_of_subsets > 1) {calibratedSamples[[s.num]] <- as.data.frame(sample_network_filt)[,c(1,2,4)]
      } else { calibratedSamples <- as.data.frame(sample_network_filt)[,c(1,2,4)] } 
    } else { 
      message( paste0("ERROR: Seed, ", seed, ", returned 0 False or True Positive relations. Try a different seed.\n" ) ) }
    seed = seed+10
  }    
  
  return(calibratedSamples)
}
