

#' @title The CCC Algorithm
#' 
#' @description A R-wrapper which directly calls the original Java code for the CCC algorithm (\url{http://kdbio.inesc-id.pt/software/ccc-biclustering/}) and transforms it to the output format of the \code{Biclust} R package.
#' 
#' @details PLACEHOLDER
#' @author Ewoud De Troyer
#' 
#' @references Sara C. Madeira, Miguel C. Teixeira, Isabel Sa Correia and Arlindo L. Oliveira, "Identification of Regulatory Modules in Time Series Gene Expression Data using a Linear Time Biclustering Algorithms", \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics}
#' 
#' @export
#' @param matrix The input matrix in which the columns are ordered by time.
#' @param overlap Numeric in ]0,1[] containing the maximum percentage of overlapping allowed. The \code{info} slot will contain a vector which filters the discovered Biclusters based on this value.
#' @return A Biclust S4 Class object containing extra information of the CCC algorithm result in the \code{info} slot.
#' 
#' @examples 
#' \dontrun{
#' data(heatdata)
#' out <- CCC(heatdata)
#' info(out,method="bonferroni")
#' }
CCC <- function(matrix,overlap=0.1){
  call <- match.call()
   
  # Check overlap
  if(!(overlap>0 & overlap<1)){stop("overlap needs to be in ]0,1[ interval")}
  
  # Check if matrix is matrix
  if(class(matrix)!="matrix"){stop("matrix parameter should contain a matrix object",call.=FALSE)}

  # Add row- and column names if not available
  if(is.null(rownames(matrix))){rownames(matrix) <- paste0("Row",c(1:nrow(matrix)))}
  if(is.null(colnames(matrix))){colnames(matrix) <- paste0("Col",c(1:ncol(matrix)))}
  
  # Transform matrix to data frame
  matrixdf <- as.data.frame(matrix)
  matrixdf <- cbind(data.frame(GENES=rownames(matrixdf)),matrixdf)
  rownames(matrixdf) <- NULL
  
  # Make a txt file containing the matrix
  matrixloc <- tempfile("matrixloc",fileext=".txt") 
  write.table(matrixdf,file=matrixloc,sep="\t",quote=FALSE,row.names=FALSE,na="")
  # write.table(format(matrixdf, digits=2),file=matrixloc,sep="\t",quote=FALSE,row.names=FALSE)
  
  
  
  # Result location
  resultloc <- dirname(matrixloc)
  
  # Java file location
  javaloc <- paste0("\"",find.package("CCCBiclust")[1],"/java/Test_TCBB_CCC_Biclustering.jar\"")
  

  # Java command and execution
  current_wd <- getwd()
  # javaloc <- paste0("\"",current_wd,"/inst/java/Test_TCBB_CCC_Biclustering.jar\"")
  # javaloc <- gsub("/","\\\\",javaloc)
  
  
  setwd(resultloc)
  command <- paste("java -jar -Xss50M -Xms1024M -Xmx1024M",javaloc,matrixloc,overlap)
  
  out_command <- tryCatch({
    system(command)
  },warning=function(x){return("WARNING")},error=function(x){return("WARNING")})
  
  setwd(current_wd)
  
  if(out_command!="WARNING"){
    out <- CCC2biclust(loc=resultloc,nRows=nrow(matrix),row_names=rownames(matrix),col_names=colnames(matrix),nCols=ncol(matrix),call=call)

  }else{
    
    if(file.exists(paste0(resultloc,"/CCC_BICLUSTERS_SORTED_PVALUE.txt"))){
      warning("CCC Algorithm may not have succeeded. Check if output is not the same as a previous run!")
      out <- CCC2biclust(loc=resultloc,nRows=nrow(matrix),row_names=rownames(matrix),col_names=colnames(matrix),nCols=ncol(matrix),call=call)
    }else{
      warning("CCC Algorithm did not succeed!")
      
      out <- new("Biclust",Parameters=list(Call=call,Method="CCC"),
                 RowxNumber=matrix(FALSE,nrow=nrow(matrix),ncol=1,dimnames=list(rownames(matrix),NULL)),
                 NumberxCol=matrix(FALSE,nrow=1,ncol=ncol(matrix),dimnames=list(NULL,colnames(matrix))),
                 Number=0,
                 info=list())
    }
  }
 
  return(out)
}


CCC2biclust <- function(loc="",nRows,nCols,row_names,col_names,call){
  
  # READ IN `CCC_BICLUSTERS_SORTED_PVALUE`
  CCC_file <- file(paste0(loc,"/CCC_BICLUSTERS_SORTED_PVALUE.txt"))
  CCC_lines <- readLines(CCC_file)
  close(CCC_file)
  
  if(length(CCC_lines)>0){
    # `#`-info lines
    info_index <- which(sapply(CCC_lines,FUN=function(x){substring(x,1,1)})=="#")
    
    
    # Prepare output
    nBC <- length(info_index)
    pvalues <- 1:nBC
    patterns <- rep("",nBC)
    filter <- rep(FALSE,nBC)
    names(pvalues) <- names(patterns) <- names(filter) <- paste0("BC",1:nBC)
    
    RowxNumber <- matrix(FALSE,nrow=nRows,ncol=nBC,dimnames=list(row_names,NULL))
    NumberxCol <- matrix(FALSE,nrow=nBC,ncol=nCols,dimnames=list(NULL,col_names))
    
    
    for(i in 1:length(info_index)){
      info <- strsplit(names(info_index)[i],"\t")[[1]]
      BC <- as.numeric(strsplit(info[1],"_")[[1]][2])
      ngenes <- as.numeric(strsplit(info[4],"=")[[1]][2])
      pvalues[BC] <- as.numeric(strsplit(info[5],"=")[[1]][2])
      
      genenames <- sapply(CCC_lines[(info_index[i]+1):(info_index[i]+ngenes)],FUN=function(x){strsplit(x,"\t")[[1]][1]})
      names(genenames) <- NULL
      RowxNumber[unlist(sapply(genenames,FUN=function(x){return(which(x==row_names))})),BC] <- TRUE
      
      cond_temp <- as.numeric(strsplit(strsplit(info[3],"=")[[1]][2],"-")[[1]])
      NumberxCol[BC,(cond_temp[1]:cond_temp[2])] <- TRUE
      
      patterns[BC] <- strsplit(CCC_lines[info_index[i]+1],"\t")[[1]][2]
      
    }
    
    
    # READ IN `CCC_BICLUSTERS_SORTED_PVALUE_FILTERED_OVERLAPPING`
    if(file.exists(paste0(loc,"/CCC_BICLUSTERS_SORTED_PVALUE_FILTERED_OVERLAPPING.txt"))){
      CCC_file <- file(paste0(loc,"/CCC_BICLUSTERS_SORTED_PVALUE_FILTERED_OVERLAPPING.txt"))
      CCC_lines <- readLines(CCC_file)
      close(CCC_file)
      
      if(length(CCC_lines)>0){
        info_index2 <- which(sapply(CCC_lines,FUN=function(x){substring(x,1,1)})=="#")
        
        for(i in 1:length(info_index2)){
          info <- strsplit(names(info_index2)[i],"\t")[[1]]
          BC <- as.numeric(strsplit(info[1],"_")[[1]][2])
          filter[BC] <- TRUE
        }
      }
    }
    
    result <- new("Biclust",Parameters=list(Call=call,Method="CCC"),
                  RowxNumber=RowxNumber,
                  NumberxCol=NumberxCol,
                  Number=nBC,
                  info=list(
                    CCCinfo=data.frame(patterns=patterns,
                                       pvalues=pvalues,
                                       filter=filter)
                  ))
    
    return(result)
  }else{
    out <- new("Biclust",Parameters=list(Call=call,Method="CCC"),
               RowxNumber=matrix(FALSE,nrow=nRows,ncol=1,dimnames=list(row_names,NULL)),
               NumberxCol=matrix(FALSE,nrow=1,ncol=nCols,dimnames=list(NULL,col_names)),
               Number=0,
               info=list())
  }
  

}

