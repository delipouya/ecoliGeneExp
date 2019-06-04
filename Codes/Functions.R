
# function to check to see if packages are installed. 
#   Install them if they are not, then load them into the R session.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    BiocManager::install(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}



## loading required packages 
Initialize <- function(){
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  options(stringsAsFactors = FALSE)
  listOfPackages <- c('Biobase',  'limma',  'pheatmap', 'gplots', 'reshape2', 
                      'plyr', 'stringr', 'ggplot2', 'WGCNA',  'fastLiquidAssociation', 'parallel')
  
  ipak(listOfPackages)
  
}


## reshaping data for visualization by ggplot 
Reshape <- function( data ){
  data.re <- as.data.frame(unlist(as.data.frame(t( data )) ) )
  colnames(data.re) <- 'geneExp'
  return(data.re)
}




