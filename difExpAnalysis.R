setwd("C:/Users/Delaram/Desktop/Palsson_microarray_data")
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

################# QC

l <- list(rawData_aerobic, rawData_Anaerobic, rawData_Anaerobic_Nitrate)
l.t <- sapply(1:length(l), function(i) as.data.frame(unlist(as.data.frame(t(l[[i]])) ) ) , simplify = F)
sapply(1:length(l.t), function(i) colnames(l.t[[i]])<<- 'geneExp', simplify = F)
sapply(1:length(l.t), function(i) l.t[[i]]$resp<<-c('aerobic','anaerobic','nitrate')[i])
sum.df <- do.call(rbind, l.t)

ggplot(sum.df, aes(x=resp, y=geneExp, fill = resp)) + 
  geom_boxplot() + scale_fill_brewer(palette="Dark2")

merged<- cbind(rawData_aerobic, rawData_Anaerobic)
grm <- c(rep('aerobic',ncol(rawData_aerobic)), rep('Anaerobic', ncol(rawData_Anaerobic)))

boxplot(merged)
pheatmap(cor(merged), labels_row  = grm, labels_col = grm)

pc <- prcomp(merged)
plot(pc, col ='blue')
plot(pc$x[,1:2], col='dark blue')
pcr <- data.frame(pc$rotation[,1:3], Group = grm)
ggplot(pcr,aes(x=PC1, y=PC2, color= Group))+geom_point(size=2)+theme_bw()

m.scale <- t(scale(t(merged), scale = F))
pc.scaled <- prcomp(m.scale)
plot(pc.scaled)
plot(pc.scaled$x[,1:2])
pcr <- data.frame(pc.scaled$rotation[,1:3], Group = grm)
ggplot(pcr,aes(x=PC1, y=PC2, color= Group))+geom_point(size=2)+theme_bw()

#### Differential expression analysis
gr <- factor(grm)
t.merged <- data.frame(t(merged))
t.merged$description <- gr
design <- model.matrix(~description + 0, t.merged)
colnames(design) <- levels(gr)
## fit a linear model 
fit <- lmFit(merged, design)
cont.matrix <- makeContrasts(aerob-Anaerob, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
ex <- topTable(fit2, adjust='fdr', sort.by = 'B', number = Inf) 


## cleaning the results 
b.inc <- grep("_b", row.names(ex), perl=TRUE, value=TRUE)
ex.s <- subset(ex,row.names(ex) %in% b.inc )
ex.s2 <- subset(ex, !row.names(ex) %in% b.inc )
ex.s2 <- subset(ex.s2,row.names(ex.s2) %in% grep("b+", row.names(ex.s2), perl=TRUE, value=TRUE))
ex.s$genes <- str_split_fixed(row.names(ex.s), "_", 3)[,2]
ex.s2$genes <- str_split_fixed(row.names(ex.s2), "_", 2)[,1]
ex <- rbind(ex.s, ex.s2)
length(unique(ex$genes))

####
hist(ex$logFC, main = 'logFC distribution', xlab='logFC', col="#009999")
ex.up  <- subset(ex, logFC>3 & adj.P.Val < 0.05)
ex.up.Genes <- unique(ex.up$genes)
dim(ex.up)

ex.down  <- subset(ex, logFC<(-3) & adj.P.Val < 0.05)
ex.down.Genes <- unique(ex.down$genes)
dim(ex.down)


#### fastLA Analysis 
# input > col:gene , row:observations 

allowWGCNAThreads()
data <- t(merged)
selected <- c(rownames(ex.up),rownames(ex.down))
indices <- which(colnames(data)%in%selected)
result <- fastMLA(data=data,topn=200, nvec=indices, rvalue=1.0, threads=detectCores()) # cut?
CNMcalc <- mass.CNM(data=data, GLA.mat=result, nback= 200)
finalres<-CNMcalc$"top p-values"
write.csv(finalres,file="resfastLA_ecoli.csv")
boot <- CNMcalc$"bootstrap triplets"
write.csv(boot,file="boot.csv")
clust <- makeCluster(32)
GLAeasy <- fastboots.GLA(boot, data=data, clust=clust, boots=30, perm=100, cut=3)
stopCluster(clust)
closeAllConnections()
write.csv(GLAeasy,file=" resfasLAboot.csv")


