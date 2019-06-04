## in this sript we'll perform differential expression analysis 
##    using the limma package

#  input: gene-expression matrices of aerobic and Anaerobic data
#  output:  highly up and down regulated genes(probes) to be used as controllers(X3) in FastLA


source('Codes/Functions.R')
Initialize()

##  loading input data
listOfGeneExp <- readRDS('Data/listOfGeneExp.rds')
rawData_aerobic <- listOfGeneExp[['rawData_aerobic']]
rawData_Anaerobic <- listOfGeneExp[['rawData_Anaerobic']]

GeneExp_merged <- cbind(rawData_aerobic,  rawData_Anaerobic)
condition_Names = gsub('.CEL','',gsub('ec_', '', colnames(GeneExp_merged)))
respirations <- c(rep('aerobic',ncol(rawData_aerobic)), rep('Anaerobic', ncol(rawData_Anaerobic)))


### Differential expression analysis 

# make a design matrix
groups <- factor(respirations)
GeneExp_merged.t <- data.frame(t(GeneExp_merged))
GeneExp_merged.t$description <- groups
design <- model.matrix(~description + 0, GeneExp_merged.t)
colnames(design) <- levels(groups)


## fit a linear model 
fit <- lmFit(GeneExp_merged, design)
cont.matrix <- makeContrasts(aerobic - Anaerobic, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
topDiffExp <- topTable(fit2, adjust='fdr', sort.by = 'B', number = 1000) 



#### visualize results
hist(topDiffExp$logFC, main = 'logFC distribution', xlab='logFC', col="#009999") ## logFC distribution
ggplot(topDiffExp, aes(x=logFC,y=-log10(adj.P.Val) ))+
  geom_point(color='dark blue')+theme_bw()+ggtitle('Volcano Plot') ## volcano plot


### subseting highly up and down regulated genes(probes) 
##  to be used as controller(X3) in Liquid association analysis
topDiffExp_upRegulated  <- subset(topDiffExp, logFC>2 & adj.P.Val < 0.05)
topDiffExp_downRegulated  <- subset(topDiffExp, logFC<(-2) & adj.P.Val < 0.05)


saveRDS(topDiffExp_upRegulated, 'Data/topDiffExp_upRegulated.rds')
saveRDS(topDiffExp_downRegulated, 'Data/topDiffExp_downRegulated.rds')


