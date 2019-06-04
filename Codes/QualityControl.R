## in this script, we will reshape the data for visualization 
##    and perform quality control on it 


source('Codes/Functions.R')
Initialize()


listOfGeneExp <- readRDS('Data/listOfGeneExp.rds')
rawData_aerobic <- listOfGeneExp[['rawData_aerobic']]
rawData_Anaerobic <- listOfGeneExp[['rawData_Anaerobic']]


## reshape data for visualization
listOfGeneExp_reshaped <- lapply(listOfGeneExp, function(x) Reshape(x))
sapply(1:length(listOfGeneExp_reshaped), 
       function(i) listOfGeneExp_reshaped[[i]]$respiration <<-c('aerobic','anaerobic','nitrate')[i], simplify = F)

##  head(melt(listOfGeneExp.t[[3]]))
GeneExps_reshaped <- do.call(rbind, listOfGeneExp_reshaped)
rownames(GeneExps_reshaped) <- NULL
ggplot2::ggplot(GeneExps_reshaped, aes(y=geneExp, x=respiration))+geom_violin(aes(fill=respiration))+theme_bw()



GeneExp_merged <- cbind(rawData_aerobic,  rawData_Anaerobic)
condition_Names = gsub('.CEL','',gsub('ec_', '', colnames(GeneExp_merged)))
respirations <- c(rep('aerobic',ncol(rawData_aerobic)), rep('Anaerobic', ncol(rawData_Anaerobic)))
saveRDS(GeneExp_merged, 'Data/GeneExp_merged.rds')

### checking if data needs normalization
pdf('QC.pdf',width = 15, height = 12)

par(mfrow=c(2,1))
boxplot(rawData_aerobic, col='pink1',horizontal=F,names=NULL, outline=F, 
        main = "aerobic samples GE distribution",xlab = "#samples",ylab = "Gene expression value")
boxplot(rawData_Anaerobic, col='slateblue1',horizontal=F,names=NULL, outline=F,
        main = "Anaerobic samples GE distribution",xlab = "#samples",ylab = "Gene expression value")
par(mfrow=c(1,1))


### Principal component analysis

## without scaling
pc <- prcomp(GeneExp_merged)
plot(pc, col ='goldenrod3',main='PCA distribution(before scaling)', xlab='PCA index')
plot(pc$x[,1:2], col='dark blue',main='gene PCA(before scaling)')
pcr <- data.frame(pc$rotation[,1:3], Group = respirations)
ggplot(pcr,aes(x=PC1, y=PC2, color= Group))+geom_point(size=2)+theme_bw()+ggtitle('Sample PCA(before scaling)')


## with scaling
m.scale <- t(scale(t(GeneExp_merged), scale = F))
pc.scaled <- prcomp(m.scale)
plot(pc.scaled, col='cyan4', main='PCA distribution(after scaling)', xlab='PCA index')
plot(pc.scaled$x[,1:2],  main='gene PCA(after scaling)')
pcr <- data.frame(pc.scaled$rotation[,1:3], Group = respirations)
ggplot(pcr,aes(x=PC1, y=PC2, color= Group))+geom_point(size=2)+theme_bw()+ggtitle('Sample PCA(after scaling)')

dev.off()



