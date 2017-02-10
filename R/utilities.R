#' @importFrom methods new
#' @importFrom stats p.adjust
#' @importFrom utils write.table
#' @import biclust


#' @title Heat Stress Data
#'
#' @description Gene expression data related with the yeast response to heat stress (see Reference).
#'
#' @format A matrix with 5811 genes and 5 conditions
#' @references A. P. Gasch, P. T. Spellman, C. M. Kao, O. Carmel-Harel, M. B. Eisen, G. Storz, D. Votstein, and P. O. Brown. Genomic expression programs in the response of yeast cells to environmental changes. Molecular Biology of the Cell, 11:4241-4257,2000.
#' @name heatdata
NULL

# A S4 Class Object completely similar to a Biclust S4 object apart from the extra information in the \code{info} slot.
#
# @slot Parameters Saves input Parameters in a list.
# @slot RowxNumber Logical Matrix which contains 1 in [i,j] if Row i is in Bicluster j
# @slot NumberxCol Logical Matrix which contains 1 in [i,j] if Col j is in Bicluster i
# @slot Number Number of Biclusters
# @slot info List containing additional info such as p-values, patterns and overlap filter of biclusters
#
# @export
# setClass('CCC',
#          representation = representation(
#            Parameters = 'list',
#            RowxNumber = 'matrix',
#            NumberxCol = 'matrix',
#            Number = 'numeric',
#            info = 'list')
# )


#' @title Discretization used for CCC algorithm.
#' @description Discretization procedure proposed by Ji and Tan (see References). A \code{M x N} input matrix results in a \code{M x (N-1)} discretized matrix.
#' @param matrix Gene Expression Matrix
#' @param threshold Threshold for second step in discretization procedure.
#' @references L. Ji and K. Tan. Mining gene expression data for positive and negative co-regulated gene clusters. Bioinformatics, 20(16):2711-2718, 2004.
#' @references L. Ji and K. Tan. Identifying time-lagged gene clusters using gene expression data. Bioinformatics, 21(4):509-516, 2005.
#' @author Ewoud De Troyer
#' @return Returns a discretized matrix (see Description)
#' @export
CCC_disc <- function(matrix,threshold=1){
  
  if(class(matrix)!="matrix"){stop("matrix should be a matrix")}
  if(is.null(rownames(matrix))){rownames(matrix) <- paste0("Row",c(1:nrow(matrix)))}
  if(is.null(colnames(matrix))){colnames(matrix) <- paste0("Col",c(1:ncol(matrix)))}
  
  # No duplicate row names allowed!
  if(sum(table(rownames(matrix))>1)){stop("No duplicate row names allowed!")}
  
  if(threshold<=0){"threshold should be greater than 0"}
  

  rowA1 <- function(row,threshold){
    out <- rep("",length(row)-1)
    for(i in 1:(length(row)-1)){
      
      if(is.na(row[i])|is.na(row[i+1])){
        out[i] <- NA
      }else{
        if(row[i]==0){
          if(row[i+1]<0){
            out[i] <- "D"
          }else if(row[i+1]>0){
            out[i] <- "U"
          }else{
            out[i] <- "N"
          }
        }else{
          A2 <- (row[i+1]-row[i])/abs(row[i])
          if(A2<=(-threshold)){
            out[i] <- "D"
          }else if(A2>=threshold){
            out[i] <- "U"
          }else{
            out[i] <- "N"
          }
        }
      }
    }
    return(out)
  }
  
  out1 <- t(apply(matrix,MARGIN=1,FUN=rowA1,threshold=threshold))
  
  colnames_temp <- rep("",ncol(matrix)-1)
  for(i in 1:(ncol(matrix)-1)){
    colnames_temp[i] <- paste0(colnames(matrix)[i],"-",colnames(matrix)[i+1])
  }
  
  colnames(out1) <- colnames_temp
  
  return(out1)
}


#' @title Filter CCCinfo in info slot of CCC Result
#' @description Compute adjusted p-values and order biclusters based on them. 
#' @param resbic Result from CCC.
#' @param method Adjust p-values for multiplicity. Can be one of the following: \code{c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")}.
#' @param alpha Significance level.
#' @param filter_overlap Filter Bicluster results based on overlap set in CCC algorithm.
#' @author Ewoud De Troyer
#' @return Returns a filtered and ordered CCCinfo (info slot of CCC result)
#' @examples 
#' \dontrun{
#' data(heatdata)
#' out <- CCC(heatdata)
#' CCCinfo(out,method="bonferroni")
#' }
#' @export
CCCinfo <- function(resbic,method="none",alpha=0.01 ,filter_overlap=FALSE){
  
  if(class(resbic)!="Biclust"){stop("resbic is not a Biclust object")}
  if(resbic@Parameters$Method!="CCC"){stop("resbic is not a result from the CCC algorithm")}
  if(!("CCCinfo"%in%names(resbic@info))){stop("CCCinfo not available in info slot")}
  if(!(method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))){stop("Incorrect method. Should be one of the following : \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\", \"none\"")}

  if(filter_overlap){
    out <- resbic@info$CCCinfo[resbic@info$CCCinfo$filter_overlap,c(1:5)]
  }else{
    out <- resbic@info$CCCinfo[,c(1:5)]
  }
  
  out$adjusted_pvalues <- p.adjust(out$pvalues,method)
  out <- out[which(out$adjusted_pvalues<=alpha),]
  out <- out[order(out$pvalues),]
  
  out <- out[,c("patterns","pvalues","adjusted_pvalues","nGenes","nConditions","from_to_cond")]
  
  if(method=="none"){out <- out[,-3]}
  
  cat(nrow(out),"significant biclusters were discovered.\n\n")
  return(out)
}


#' @title Filter eCCCinfo in info slot of eCCC Result
#' @description Reorders eCCCinfo slot based on the p-values order and applies one of the 2 available filters.
#' @param resbic Result from eCCC.
#' @param filter Choose one of the following filters: \code{"none", "Bonf0.01","Bonf0.01_overlap"}. Respectively no biclusters are filtered, significant biclusters with under 1 percent after Bonferroni correction are filtered and biclusters after Bonferroni (1 percent) with taking the maximum overlap into account are filtered.
#' @author Ewoud De Troyer
#' @return Returns a filtered and ordered eCCCinfo (info slot of eCCC result)
#' @examples 
#' \dontrun{
#' data(heatdata)
#' out <- eCCC_ext(heatdata,minr=3,minc=2)
#' eCCCinfo(out,filter="Bonf0.01_overlap")
#' }
#' @export
eCCCinfo <- function(resbic,filter="none"){
  if(class(resbic)!="Biclust"){stop("resbic is not a Biclust object")}
  if((resbic@Parameters$Method!="eCCC")&(resbic@Parameters$Method!="eCCC_ext")){stop("resbic is not a result from the e-CCC algorithm")}
  if(!("eCCCinfo"%in%names(resbic@info))){stop("eCCCinfo not available in info slot")}
  if(!(filter %in% c("none","Bonf0.01","Bonf0.01_overlap"))){stop("filter needs to be \"none\", \"Bonf0.01\" or \"Bonf0.01_overlap\"")}
  
  
  out <- resbic@info$eCCCinfo
  out <- out[order(out$pvalues_order),]
  if(filter!="none"){
    if(filter=="Bonf0.01"){
      out <- out[out$filter_Bonf0.01,]
      cat(nrow(out),"significant biclusters were discovered.\n\n")
    }else if(filter=="Bonf0.01_overlap"){
      out <- out[out$filter_Bonf0.01_overlap,]
      cat(nrow(out),"significant biclusters were discovered.\n\n")
    }
  }
  
  out <- out[,c(1,3,4,5),drop=FALSE]
  
  return(out)
}


## TO DO:
# - give info on 1 or multiple BC
# - order BC by pvalues + apply pvalue adjustment (check if pvalue is already bonferoni)