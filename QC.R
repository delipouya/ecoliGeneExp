#setwd("C:/Users/Delaram/Desktop/Palsson_microarray_data")
library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(reshape2)
library(plyr)
library(stringr)
library(ggplot2)
library(WGCNA)
library(fastLiquidAssociation)
library(parallel)

rawData = read.csv("GeneExp_data.csv", header = T)
arrayInfo = read.csv("array_info.csv", header = T)


#### rawData preprocessing >> row: gene , columns: conditions 
row.names(rawData) = rawData[,1]
rawData = rawData[,-1]

#### arrayInfo preprocessing  >> rows: conditions , columns: experiment info.
arrayInfo = arrayInfo[-c(214,215),]

#### finding indexs based on Aerobic status 
length(arrayInfo$Aerobic.Anaerobic)
aerobic_index = which(arrayInfo$Aerobic.Anaerobic %in% c("Aerobic", "aerobic"))
anaerobic_index = which(arrayInfo$Aerobic.Anaerobic %in% c("anaerobic" ,"Anaerobic"))
Anaerobic_Nitrate_index = which(arrayInfo$Aerobic.Anaerobic == "Anaerobic/Nitrate")

###### splitting rawData based on aerobic/anaerobic/Anaerobic-Nitrate
rawData_aerobic = rawData[,aerobic_index]
rawData_Anaerobic = rawData[,anaerobic_index]
rawData_Anaerobic_Nitrate = rawData[,Anaerobic_Nitrate_index]

################# QC ##

l <- list(rawData_aerobic, rawData_Anaerobic, rawData_Anaerobic_Nitrate)
l.t <- sapply(1:length(l), function(i) as.data.frame(unlist(as.data.frame(t(l[[i]])) ) ) , simplify = F)
sapply(1:length(l.t), function(i) colnames(l.t[[i]])<<- 'geneExp', simplify = F)
sapply(1:length(l.t), function(i) l.t[[i]]$resp<<-c('aerobic','anaerobic','nitrate')[i])
sum.df <- do.call(rbind, l.t)
merged<- cbind(rawData_aerobic, rawData_Anaerobic)
gr = gsub('.CEL','',gsub('ec_', '', colnames(merged)))
grm <- c(rep('aerobic',ncol(rawData_aerobic)), rep('Anaerobic', ncol(rawData_Anaerobic)))

pdf('heatmap.pdf', width = 30, height=30)
pheatmap(cor(merged), labels_row  = gr, labels_col = gr)
dev.off()

pdf('QC.pdf',width = 15, height = 12)
par(mfrow=c(2,1))
boxplot(rawData_aerobic, col='pink1',horizontal=F,names=NULL, outline=F, 
        main = "aerobic samples GE distribution",xlab = "#samples",ylab = "Gene expression value")
boxplot(rawData_Anaerobic, col='slateblue1',horizontal=F,names=NULL, outline=F,
        main = "Anaerobic samples GE distribution",xlab = "#samples",ylab = "Gene expression value")
par(mfrow=c(1,1))
#ggplot(sum.df, aes(x=resp, y=geneExp, fill = resp)) + 
#  geom_boxplot() + scale_fill_brewer(palette="Dark2")

pc <- prcomp(merged)
plot(pc, col ='goldenrod3',main='PCA distribution(before scaling)', xlab='PCA index')
plot(pc$x[,1:2], col='dark blue',main='gene PCA(before scaling)')
pcr <- data.frame(pc$rotation[,1:3], Group = grm)
ggplot(pcr,aes(x=PC1, y=PC2, color= Group))+geom_point(size=2)+theme_bw()+ggtitle('Sample PCA(before scaling)')

m.scale <- t(scale(t(merged), scale = F))
pc.scaled <- prcomp(m.scale)
plot(pc.scaled, col='cyan4', main='PCA distribution(after scaling)', xlab='PCA index')
plot(pc.scaled$x[,1:2],  main='gene PCA(after scaling)')
pcr <- data.frame(pc.scaled$rotation[,1:3], Group = grm)
ggplot(pcr,aes(x=PC1, y=PC2, color= Group))+geom_point(size=2)+theme_bw()+ggtitle('Sample PCA(after scaling)')

dev.off()
