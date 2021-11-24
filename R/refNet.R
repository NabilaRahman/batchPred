#' Create reference gene association network
#' 
#' \code{refNet} creates reference gene association network for batch effect 
#'  evaluation. The output is required for \code{create_ref_sample} and 
#'  \code{batchPred} functions.
#' 
#' @param network_gold_standard Datatable of gene-gene network gold standard.
#'  Column 1 and 2 are gene IDs. Column 3 is posterior probability with
#'  known edges set to 1. \code{refNet} is optimised for "full 
#'  networks", which can be be downloaded from 
#'  http://giant.princeton.edu/download/.
#' @param input_edata Character vector of gene ids that 
#'  exist within your expression dataset.
#' @param gold_standard_path (Optional). Path to gene-gene network gold
#'  standard file. Can be specified in place of \code{network_gold_standard}.
#' @param false_positive cutoff value for false positive gene-gene 
#'  interactions. Only false positives with lower post probability are 
#'  listed. 
#' @param true_positive cutoff value for true positive gene-gene 
#'  interactions. Only true positives with higher post probability are 
#'  listed.
#'
#' @return A reference gene association network for genes present in 
#'	\code{input_edata}. This gold standard only retains the highest 
#'	and lowest confidence gene interactions.Output is a data.table with 
#'	columns 1 and 2 being gene ID, column 2 is posterior probability. 
#'	While column 3 is defines whether gene relation is a true positive or 
#'	false positive (0) based on user defined thresholds.
#' @examples
#' #Create Mock input gold standard
#' network_gold_standard<- data.table(
#'     V1 = c(3091,4763,4204,10628,4321,6359,5630,5630,100131077) 
#'     , V2 = c(7352,8698,6855,2475,5155,23513,100131077,751816,751816)
#'     , V3 =  c(1,1,1,1,1,1,0.0384244,0.00357639,0.0097727)
#'     , V4 = c(0.0193132,0.0299962,0.0244092,0.0437484,0.107417,0.0145827,NA
#'       ,NA,NA)
#'     , stringsAsFactors = F )
#'     
#' #Create Mock Gene IDs    
#' input_edata <- c('3091','4763','4204','7352','8698','6855' )
#'	
#'	#Generate Reference Network
#' refNet( network_gold_standard, input_edata, 
#'	gold_standard_path = NULL, false_positive = 0.025, true_positive = 0.5) 
#' @export

refNet <- function(network_gold_standard=NULL, input_edata=input_edata, gold_standard_path = NULL, false_positive = 0.025, true_positive = 0.5) 
{
  start.time <- Sys.time()
  if(!exists("network_gold_standard") & is.null(gold_standard_path) & is.null(ncol(network_gold_standard)) ) {
    message("Warning: 'network_gold_standard' and 'gold_standard_path' missing. Please specify either 'network_gold_standard' or 'gold_standard_path'.")
  } else if ( !is.null(gold_standard_path) && file.exists(gold_standard_path) ) {
    network_gold_standard<-data.table::fread(gold_standard_path, fill=T)
  } else if ( !is.null(gold_standard_path) && !file.exists(gold_standard_path )) {
    message("File specified at 'gold_standard_path' does not exist.")
  }
  
  if(!exists("input_edata") | class(input_edata) != "character" ) {
    message("Warning: 'input_edata' missing. This must be a vector of gene_ids from your dataset of interest")
  } 
  
  # Clean up Reference Network
  network_gold_standard <- network_gold_standard[,c(1:3)]
  colnames(network_gold_standard) <- c("geneID1", "geneID2", "postProbability")
  network_gold_standard$geneID1 <- as.character(network_gold_standard$geneID1)
  network_gold_standard$geneID2 <- as.character(network_gold_standard$geneID2)
  input_edata <- as.character(input_edata)
  
  # Only retains reference gene-gene relation postProbability is > true_positive or < false_positive
  network_gold_standard_filt_prob <- network_gold_standard
  network_gold_standard_filt_prob$Confidence[ network_gold_standard_filt_prob$postProbability > true_positive ] <- 1
  network_gold_standard_filt_prob$Confidence[ network_gold_standard_filt_prob$postProbability < false_positive ] <- 0
  network_gold_standard_filt_prob<-network_gold_standard_filt_prob[!is.na(network_gold_standard_filt_prob$Confidence),]
  message(paste(round(((nrow(network_gold_standard)-nrow(network_gold_standard_filt_prob))/nrow(network_gold_standard))*100,2), "% associations removed due to false/true positive thresholds"))
  
  # Only retains reference gene-gene relation if it exists in input_edata
  gold_standard_filt <- network_gold_standard_filt_prob[network_gold_standard_filt_prob$geneID1 %in% input_edata,]
  gold_standard_filt <- gold_standard_filt[gold_standard_filt$geneID2 %in% input_edata,]
  message(paste( round(((nrow(network_gold_standard_filt_prob)-nrow(gold_standard_filt))/nrow(network_gold_standard))*100,3), "% associations removed because genes do not exist in input_edata") )
  
  #Time Taken
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message( paste("Elapsed Time:", round(time.taken,2), "minutes" ))
  
  # Output Reference Gene Network
  return(gold_standard_filt)
}
