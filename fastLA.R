source('QC.R') 
source('difExpAnalysis.R')

#### fastLA Analysis 
# input > col:gene , row:observations #
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


