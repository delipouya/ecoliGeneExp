# ecoliGeneExp
a project on analyzing gene expression of E. coli microarray data

Functions.R >> contains a few simple functions which are needed in other scripts       
Preprocess_Data.R >>  microarray data is preprocessed and cleaned for further analysis          
QualityControl.R >> some QC for double-checking quality of the data. data is already normalized       
difExpAnalysis.R >> performing differential-expression analysis using the limma package;  
   differentially-expressed genes will be used as controllers(X3) during the liquid association analysis          
fastLA.R >> applying fastLA function to the data         
