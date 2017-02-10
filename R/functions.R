### TO DO: TRY LOWER MINR/MINC for eCCC!! + check if switching rows gives different results




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
#' @param overlap Numeric in ]0,1[ containing the maximum percentage of overlapping allowed. The \code{info} slot will contain a logical vector which can filter the biclusters which overlap more than \code{overlap*100} percent. 
#' @return A Biclust S4 Class object containing extra information of the CCC algorithm result in the \code{info} slot.
#' 
#' @examples 
#' \dontrun{
#' data(heatdata)
#' out <- CCC(heatdata)
#' CCCinfo(out,method="bonferroni")
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
  
  # No duplicate row names allowed!
  if(sum(table(rownames(matrix))>1)){stop("No duplicate row names allowed!")}
  
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
  cat("\nTransforming to Biclust Output... ")
  
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
    nGenes <- nConditions <- 1:nBC
    from_to_cond <- rep("",nBC)
    names(pvalues) <- names(patterns) <- names(filter) <- names(nGenes) <- names(nConditions) <- names(from_to_cond) <- paste0("BC",1:nBC)
    
    RowxNumber <- matrix(FALSE,nrow=nRows,ncol=nBC,dimnames=list(row_names,NULL))
    NumberxCol <- matrix(FALSE,nrow=nBC,ncol=nCols,dimnames=list(NULL,col_names))
    
    
    for(i in 1:length(info_index)){
      info <- strsplit(names(info_index)[i],"\t")[[1]]
      BC <- as.numeric(strsplit(info[1],"_")[[1]][2])
      ngenes <- as.numeric(strsplit(info[4],"=")[[1]][2])
      nGenes[BC] <- ngenes 
      pvalues[BC] <- as.numeric(strsplit(info[5],"=")[[1]][2])
      
      genenames <- sapply(CCC_lines[(info_index[i]+1):(info_index[i]+ngenes)],FUN=function(x){strsplit(x,"\t")[[1]][1]})
      names(genenames) <- NULL
      RowxNumber[unlist(sapply(genenames,FUN=function(x){return(which(x==row_names))})),BC] <- TRUE
      
      from_to_cond[BC] <- strsplit(info[3],"= ")[[1]][2]
      cond_temp <- as.numeric(strsplit(strsplit(info[3],"=")[[1]][2],"-")[[1]])
      NumberxCol[BC,(cond_temp[1]:cond_temp[2])] <- TRUE
      
      nConditions[BC] <- length(cond_temp[1]:cond_temp[2])
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
    
    CCCinfo_temp <- data.frame(patterns=patterns,
                               pvalues=pvalues,
                               nGenes=nGenes,
                               nConditions=nConditions,
                               from_to_cond=from_to_cond,
                               filter_overlap=filter)
    CCCinfo_temp$patterns <- as.character(CCCinfo_temp$patterns)
    CCCinfo_temp$from_to_cond <- as.character(CCCinfo_temp$from_to_cond)
    
    result <- new("Biclust",Parameters=list(Call=call,Method="CCC"),
                  RowxNumber=RowxNumber,
                  NumberxCol=NumberxCol,
                  Number=nBC,
                  info=list(
                    CCCinfo=CCCinfo_temp
                  ))
    cat("DONE\n\n")
    return(result)
  }else{
    out <- new("Biclust",Parameters=list(Call=call,Method="CCC"),
               RowxNumber=matrix(FALSE,nrow=nRows,ncol=1,dimnames=list(row_names,NULL)),
               NumberxCol=matrix(FALSE,nrow=1,ncol=nCols,dimnames=list(NULL,col_names)),
               Number=0,
               info=list())
    cat("DONE\n\n")
    return(out)
  }
  

}



#' @title The Extended e-CCC Algorithm
#' 
#' @description A R-wrapper which directly calls the original Java code for the Extended e-CCC algorithm (\url{http://kdbio.inesc-id.pt/software/e-ccc-biclustering/}) and transforms it to the output format of the \code{Biclust} R package.
#' 
#' @details PLACEHOLDER
#' @author Ewoud De Troyer
#' 
#' @references Sara C. Madeira and Arlindo L. Oliveira, "A polynomial time biclustering algorithm for finding genes with approximate expression patterns in gene expression time series", Algorithms for Molecular Biology 2009, 4:8 (4 June 2009)
#' 
#' @export
#' @param matrix The input matrix in which the columns are ordered by time.
#' @param minr Integer containing the row quorum (= minimum number of genes allowed in e-CCC-Biclusters).
#' @param minc Integer containing the column quorum (= minimum number of contiguous time points allowed in e-CCC-Biclusters).
#' @param maxErrors Integer containing the amount of errors allowed, per gene, in the e-CCC-Biclustering algorithm (value of e).
#' @param overlap Numeric in ]0,1[ containing the maximum percentage of overlapping allowed. The \code{info} slot will contain a logical vector which can filter the biclusters which overlap more than \code{overlap*100} percent. 
#' @param missings How to handle missing values:
#' \itemize{
#' \item \code{"remove"}: remove genes with missing values
#' \item \code{"allow"}: allow missing values as valid errors
#' \item \code{"jump"}: 'jump over' missing values
#' }
#' @param anticorrelation Logical value to allow or not allow anticorrelation. If allowed (\code{TRUE}), the algorithm will look for e-CCC-Biclusters with Sign-Changes.
#' @param restrictedErrors Logical value restricting or not restricting errors. 
#' If restricting (\code{TRUE}), errors are restricted to the symbols in the 1-neighbourhood of the symbols in 
#' the alphabet. Since the alphabet {D,N,U} is used in the predefined discretization step provided in this version 
#' of the prototype, the number of neighbours used in the restricted errors extension can only be equal to 1.
#' 
#' @return A Biclust S4 Class object containing extra information of the e-CCC algorithm result in the \code{info} slot.
#' 
#' @examples 
#' \dontrun{
#' data(heatdata)
#' out <- eCCC_ext(heatdata,minr=3,minc=2,missings="jump",anticorrelation = TRUE)
#' eCCCinfo(out,filter="Bonf0.01")
#' out@info$eCCCpatterns$BC53[1:10,]
#' }
eCCC_ext <- function(matrix,minr=1,minc=1,maxErrors=1,overlap=0.1,missings="allow",anticorrelation=FALSE,restrictedErrors=FALSE){
  call <- match.call()
  
  # Check parameters
  if(minr<1){stop("minr should be larger than 1")}
  minr <- as.integer(minr)
  if(minc<1){stop("minc should be larger than 1")}
  minc <- as.integer(minc)
  if(maxErrors<0){stop("maxErrors should be larger than 0")}
  maxErrors <- as.integer(maxErrors)
  if(!(overlap>0 & overlap<1)){stop("overlap needs to be in ]0,1[ interval")}
  
  if(length(missings)!=1){stop("missings should have a single input")}
  if(!(missings%in%c("remove","allow","jump"))){stop("missings should be one of the following: \"remove\",\"allow\",\"jump\"")}
  missings <- switch(missings,remove="R",allow="A",jump="J")
  
  anticorrelation <- ifelse(anticorrelation,"Y","N")
  restrictedErrors <- ifelse(restrictedErrors,"Y","N")
  
  
  # Check if matrix is matrix
  if(class(matrix)!="matrix"){stop("matrix parameter should contain a matrix object",call.=FALSE)}
  
  # Add row- and column names if not available
  if(is.null(rownames(matrix))){rownames(matrix) <- paste0("Row",c(1:nrow(matrix)))}
  if(is.null(colnames(matrix))){colnames(matrix) <- paste0("Col",c(1:ncol(matrix)))}
  
  # No duplicate row names allowed!
  if(sum(table(rownames(matrix))>1)){stop("No duplicate row names allowed!")}
  
  # Transform matrix to data frame
  matrixdf <- as.data.frame(matrix)
  matrixdf <- cbind(data.frame(GENES=rownames(matrixdf)),matrixdf)
  rownames(matrixdf) <- NULL
  
  # Make a txt file containing the matrix
  matrixloc <- tempfile("matrixloc",fileext=".txt") 
  write.table(matrixdf,file=matrixloc,sep="\t",quote=FALSE,row.names=FALSE,na="")

  
  # Result location
  resultloc <- dirname(matrixloc)
  
  # Java file location
  javaloc <- paste0("\"",find.package("CCCBiclust")[1],"/java/Test_AMB_E_CCC_Biclustering_Extended.jar\"")
  
  
  # Java command and execution
  current_wd <- getwd()
  # javaloc <- paste0("\"",current_wd,"/inst/java/Test_AMB_E_CCC_Biclustering_Extended.jar\"")
  # javaloc <- gsub("/","\\\\",javaloc)
  
  
  setwd(resultloc)
  command <- paste("java -jar -Xss50M -Xms1024M -Xmx1024M",javaloc,matrixloc,maxErrors,minr,minc,overlap,missings,anticorrelation,restrictedErrors)
  # java -jar -Xss50M -Xms1024M -Xmx1024M Test_AMB_E_CCC_Biclustering_Extensions.jar yourExpressionMatrix.txt maxErrors rowQuorum columnQuorum overlapping missings anticorrelation restrictedErrors
  
  out_command <- tryCatch({
    system(command)
  },warning=function(x){return("WARNING")},error=function(x){return("WARNING")})
  
  setwd(current_wd)
  
  if(out_command!="WARNING"){
    if( file.exists(paste0(resultloc,"/",gsub(".txt","",basename(matrixloc)),"_1_CCC_BICLUSTERS_SORTED_PVALUE.txt"))){
      out <- eCCCext2biclust(loc=resultloc,matrixloc=matrixloc,nRows=nrow(matrix),row_names=rownames(matrix),col_names=colnames(matrix),nCols=ncol(matrix),call=call)
    }else{
      warning("e-CCC Extended Algorithm did not succeed!")
      
      out <- new("Biclust",Parameters=list(Call=call,Method="eCCC_ext"),
                 RowxNumber=matrix(FALSE,nrow=nrow(matrix),ncol=1,dimnames=list(rownames(matrix),NULL)),
                 NumberxCol=matrix(FALSE,nrow=1,ncol=ncol(matrix),dimnames=list(NULL,colnames(matrix))),
                 Number=0,
                 info=list())
    }
    
  }else{
    if( file.exists(paste0(resultloc,"/",gsub(".txt","",basename(matrixloc)),"_1_CCC_BICLUSTERS_SORTED_PVALUE.txt"))){
      warning("e-CCC Extended Algorithm may not have succeeded. Check if output is not the same as a previous run!")
      out <- eCCCext2biclust(loc=resultloc,matrixloc=matrixloc,nRows=nrow(matrix),row_names=rownames(matrix),col_names=colnames(matrix),nCols=ncol(matrix),call=call)
    }else{
      warning("e-CCC Extended Algorithm did not succeed!")
      out <- new("Biclust",Parameters=list(Call=call,Method="eCCC_ext"),
                 RowxNumber=matrix(FALSE,nrow=nrow(matrix),ncol=1,dimnames=list(rownames(matrix),NULL)),
                 NumberxCol=matrix(FALSE,nrow=1,ncol=ncol(matrix),dimnames=list(NULL,colnames(matrix))),
                 Number=0,
                 info=list())

    }
  }
  
  return(out)
  
}



eCCCext2biclust <- function(loc="",matrixloc,nRows,nCols,row_names,col_names,call){
  cat("\nTransforming to Biclust Output... ")
  
  # READ IN `CCC_BICLUSTERS_SORTED_PVALUE`
  CCC_file <- file(  paste0(loc,"/",gsub(".txt","",basename(matrixloc)),"_1_CCC_BICLUSTERS_SORTED_PVALUE.txt"))
  CCC_lines <- readLines(CCC_file)
  close(CCC_file)
  
  
  if(length(CCC_lines)>0){
    # `#`-info lines
    info_index <- which(sapply(CCC_lines,FUN=function(x){substring(x,1,1)})=="#")
    
    
    # Prepare output
    nBC <- length(info_index)
    pvalues_order <- rep(0,nBC)
    patterns <- rep("",nBC)
    filter_Bonf0.01 <- filter_Bonf0.01_overlap <- rep(FALSE,nBC)
    eCCCpatterns <- vector("list",nBC)
    nGenes <- nConditions <- 1:nBC
    from_to_cond <- rep("",nBC)
    names(nGenes) <- names(nConditions) <- names(from_to_cond) <- names(pvalues_order) <- names(patterns) <- names(filter_Bonf0.01) <- names(filter_Bonf0.01_overlap) <- names(eCCCpatterns) <- paste0("BC",1:nBC)
    
    
    RowxNumber <- matrix(FALSE,nrow=nRows,ncol=nBC,dimnames=list(row_names,NULL))
    NumberxCol <- matrix(FALSE,nrow=nBC,ncol=nCols,dimnames=list(NULL,col_names))
    
  
    for(i in 1:length(info_index)){
      info <- strsplit(names(info_index)[i],"\t")[[1]]
      BC <- as.numeric(strsplit(info[1],"_")[[1]][2])
      ngenes <- as.numeric(strsplit(info[4],"=")[[1]][2])
      nGenes[BC] <- ngenes
      pvalues_order[BC] <- i
      patterns[BC] <- (strsplit(info[5],"=")[[1]][2])
      
      
      genes <- t(sapply(CCC_lines[(info_index[i]+1):(info_index[i]+ngenes)],FUN=function(x){strsplit(x,"\t")[[1]]}))
      rownames(genes) <- NULL
      
      n_error <- sapply(genes[,2],FUN=function(x){
        sum(strsplit(x,"")[[1]]!=strsplit(patterns[BC],"")[[1]])
      })
      names(n_error) <- NULL
      
      df_temp <- data.frame(Index=unlist(sapply(genes[,1],FUN=function(x){return(which(x==row_names))})),Pattern=genes[,2],nErrors=n_error,Conditions=rep(strsplit(info[3],"= ")[[1]][2],ngenes))
      df_temp$Pattern <- as.character(df_temp$Pattern)
      df_temp$Conditions <- as.character(df_temp$Conditions)
      rownames(df_temp) <- genes[,1]
      
      RowxNumber[df_temp$Index,BC] <- TRUE
      
      from_to_cond[BC] <- strsplit(info[3],"= ")[[1]][2]
      cond_temp <- as.numeric(strsplit(strsplit(info[3],"=")[[1]][2],"-")[[1]])
      nConditions[BC] <- length(cond_temp[1]:cond_temp[2])
      NumberxCol[BC,(cond_temp[1]:cond_temp[2])] <- TRUE
      
      eCCCpatterns[[BC]] <- df_temp
    }
    
  
    # Read in extra 2 files for filter info
    
    file1 <-  paste0(loc,"/",gsub(".txt","",basename(matrixloc)),"_CCC_BICLUSTERS_SORTED_PVALUE_FILTERED_PVALUE.txt")
    file2 <-  paste0(loc,"/",gsub(".txt","",basename(matrixloc)),"_CCC_BICLUSTERS_SORTED_PVALUE_FILTERED_PVALUE_OVERLAPPING.txt")
    
    if(file.exists(file1)){
      CCC_file <- file(file1)
      CCC_lines <- readLines(CCC_file)
      close(CCC_file)
      
      if(length(CCC_lines)>0){
        info_index <- which(sapply(CCC_lines,FUN=function(x){substring(x,1,1)})=="#")
        
        BC_filter <- sapply(info_index,FUN=function(x){
          return(as.numeric(strsplit(strsplit(CCC_lines[x],"\t")[[1]][1],"_")[[1]][2]))
        })
        names(BC_filter) <- NULL
        
        filter_Bonf0.01[BC_filter] <- TRUE
        
      }else{
        warning("FILTERED_PVALUE empty!")
      }
    }else{
      warning("FILTERED_PVALUE not available!")
    }
    
    if(file.exists(file2)){
      CCC_file <- file(file2)
      CCC_lines <- readLines(CCC_file)
      close(CCC_file)
      
      if(length(CCC_lines)>0){
        info_index <- which(sapply(CCC_lines,FUN=function(x){substring(x,1,1)})=="#")
        
        BC_filter <- sapply(info_index,FUN=function(x){
          return(as.numeric(strsplit(strsplit(CCC_lines[x],"\t")[[1]][1],"_")[[1]][2]))
        })
        names(BC_filter) <- NULL
        
        filter_Bonf0.01_overlap[BC_filter] <- TRUE
        
      }else{
        warning("FILTERED_PVALUE empty!")
      }
    }else{
      warning("FILTERED_PVALUE not available!")
    }
    
    
    eCCCinfo_temp <- data.frame(patterns=patterns,
                               pvalues_order=pvalues_order,
                               nGenes=nGenes,
                               nConditions=nConditions,
                               from_to_cond=from_to_cond,
                               filter_Bonf0.01=filter_Bonf0.01,
                               filter_Bonf0.01_overlap=filter_Bonf0.01_overlap)
    eCCCinfo_temp$patterns <- as.character(eCCCinfo_temp$patterns)
    eCCCinfo_temp$from_to_cond <- as.character(eCCCinfo_temp$from_to_cond)
    
    result <- new("Biclust",Parameters=list(Call=call,Method="eCCC_ext"),
                  RowxNumber=RowxNumber,
                  NumberxCol=NumberxCol,
                  Number=nBC,
                  info=list(
                    eCCCinfo=eCCCinfo_temp,
                    eCCCpatterns=eCCCpatterns
                  ))

    cat("DONE\n\n")
    return(result)
    
  }else{
    cat("DONE\n\n")
    out <- new("Biclust",Parameters=list(Call=call,Method="eCCC_ext"),
               RowxNumber=matrix(FALSE,nrow=nrow(matrix),ncol=1,dimnames=list(rownames(matrix),NULL)),
               NumberxCol=matrix(FALSE,nrow=1,ncol=ncol(matrix),dimnames=list(NULL,colnames(matrix))),
               Number=0,
               info=list())
    return(out)
  }
  
}






#' @title The e-CCC Algorithm
#' 
#' @description A R-wrapper which directly calls the original Java code for the e-CCC algorithm (\url{http://kdbio.inesc-id.pt/software/e-ccc-biclustering/}) and transforms it to the output format of the \code{Biclust} R package.
#' 
#' @details PLACEHOLDER
#' @author Ewoud De Troyer
#' 
#' @references Sara C. Madeira and Arlindo L. Oliveira, "A polynomial time biclustering algorithm for finding genes with approximate expression patterns in gene expression time series", Algorithms for Molecular Biology 2009, 4:8 (4 June 2009)
#' 
#' @export
#' @param matrix The input matrix in which the columns are ordered by time.
#' @param minr Integer containing the row quorum (= minimum number of genes allowed in e-CCC-Biclusters).
#' @param minc Integer containing the column quorum (= minimum number of contiguous time points allowed in e-CCC-Biclusters).
#' @param maxErrors Integer containing the amount of errors allowed, per gene, in the e-CCC-Biclustering algorithm (value of e).
#' @param overlap Numeric in ]0,1[ containing the maximum percentage of overlapping allowed. The \code{info} slot will contain a logical vector which can filter the biclusters which overlap more than \code{overlap*100} percent. 

#' @return A Biclust S4 Class object containing extra information of the e-CCC algorithm result in the \code{info} slot.
#' 
#' @examples 
#' \dontrun{
#' data(heatdata)
#' out <- eCCC(heatdata,minr=3,minc=2) 
#' eCCCinfo(out,filter="Bonf0.01")
#' out@info$eCCCpatterns$BC53[1:10,]
#' }
eCCC <- function(matrix,minr=1,minc=1,maxErrors=1,overlap=0.1){
  call <- match.call()
  
  # Check parameters
  if(minr<1){stop("minr should be larger than 1")}
  minr <- as.integer(minr)
  if(minc<1){stop("minc should be larger than 1")}
  minc <- as.integer(minc)
  if(maxErrors<0){stop("maxErrors should be larger than 0")}
  maxErrors <- as.integer(maxErrors)
  if(!(overlap>0 & overlap<1)){stop("overlap needs to be in ]0,1[ interval")}
  
  # Check if matrix is matrix
  if(class(matrix)!="matrix"){stop("matrix parameter should contain a matrix object",call.=FALSE)}
  
  # Add row- and column names if not available
  if(is.null(rownames(matrix))){rownames(matrix) <- paste0("Row",c(1:nrow(matrix)))}
  if(is.null(colnames(matrix))){colnames(matrix) <- paste0("Col",c(1:ncol(matrix)))}
  
  # No duplicate row names allowed!
  if(sum(table(rownames(matrix))>1)){stop("No duplicate row names allowed!")}
  
  # Transform matrix to data frame
  matrixdf <- as.data.frame(matrix)
  matrixdf <- cbind(data.frame(GENES=rownames(matrixdf)),matrixdf)
  rownames(matrixdf) <- NULL
  
  # Make a txt file containing the matrix
  matrixloc <- tempfile("matrixloc",fileext=".txt") 
  write.table(matrixdf,file=matrixloc,sep="\t",quote=FALSE,row.names=FALSE,na="")
  
  
  # Result location
  resultloc <- dirname(matrixloc)
  
  # Java file location
  javaloc <- paste0("\"",find.package("CCCBiclust")[1],"/java/Test_AMB_E_CCC_Biclustering.jar\"")
  
  
  # Java command and execution
  current_wd <- getwd()
  # javaloc <- paste0("\"",current_wd,"/inst/java/Test_AMB_E_CCC_Biclustering_Extended.jar\"")
  # javaloc <- gsub("/","\\\\",javaloc)
  
  
  setwd(resultloc)
  command <- paste("java -jar -Xss50M -Xms1024M -Xmx1024M",javaloc,matrixloc,maxErrors,minr,minc,overlap)

  out_command <- tryCatch({
    system(command)
  },warning=function(x){return("WARNING")},error=function(x){return("WARNING")})
  
  setwd(current_wd)
  
  if(out_command!="WARNING"){
    if( file.exists(paste0(resultloc,"/",gsub(".txt","",basename(matrixloc)),"_1_CCC_BICLUSTERS_SORTED_PVALUE.txt"))){
      out <- eCCCext2biclust(loc=resultloc,matrixloc=matrixloc,nRows=nrow(matrix),row_names=rownames(matrix),col_names=colnames(matrix),nCols=ncol(matrix),call=call)
    }else{
      warning("e-CCC Extended Algorithm did not succeed!")
      
      out <- new("Biclust",Parameters=list(Call=call,Method="eCCC_ext"),
                 RowxNumber=matrix(FALSE,nrow=nrow(matrix),ncol=1,dimnames=list(rownames(matrix),NULL)),
                 NumberxCol=matrix(FALSE,nrow=1,ncol=ncol(matrix),dimnames=list(NULL,colnames(matrix))),
                 Number=0,
                 info=list())
    }
    
  }else{
    if( file.exists(paste0(resultloc,"/",gsub(".txt","",basename(matrixloc)),"_1_CCC_BICLUSTERS_SORTED_PVALUE.txt"))){
      warning("e-CCC Extended Algorithm may not have succeeded. Check if output is not the same as a previous run!")
      out <- eCCCext2biclust(loc=resultloc,matrixloc=matrixloc,nRows=nrow(matrix),row_names=rownames(matrix),col_names=colnames(matrix),nCols=ncol(matrix),call=call)
    }else{
      warning("e-CCC Extended Algorithm did not succeed!")
      out <- new("Biclust",Parameters=list(Call=call,Method="eCCC_ext"),
                 RowxNumber=matrix(FALSE,nrow=nrow(matrix),ncol=1,dimnames=list(rownames(matrix),NULL)),
                 NumberxCol=matrix(FALSE,nrow=1,ncol=ncol(matrix),dimnames=list(NULL,colnames(matrix))),
                 Number=0,
                 info=list())
      
    }
  }
  
  return(out)
}


