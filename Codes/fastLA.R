source('Codes/Functions.R')
Initialize()


### importing and combining up and down regulated genes to use as controllers
topDiffExp_upRegulated <- readRDS('Data/topDiffExp_upRegulated.rds')
topDiffExp_downRegulated <- readRDS('Data/topDiffExp_downRegulated.rds')
topDiffExp_merged <- rbind(topDiffExp_upRegulated, topDiffExp_downRegulated)
topDiffExp_genes <- row.names(topDiffExp_merged)

## importing gene-expression matrix 
GeneExp_merged <- readRDS('Data/GeneExp_merged.rds')



#### fastLA Analysis 
# fastMLA takes a gene expression matrix as input in which,
#   columns are gene & rows are observations(conditions)

allowWGCNAThreads()
GeneExp_merged.t <- t(GeneExp_merged)
indices <- which( colnames(GeneExp_merged.t) %in% topDiffExp_genes )
result <- fastMLA(data=GeneExp_merged.t, topn=200, nvec=indices, rvalue=1.0, threads=detectCores()) 


CNMcalc <- mass.CNM(data=GeneExp_merged.t, GLA.mat=result, nback= 200)
finalres<-CNMcalc$"top p-values"
write.csv(finalres,file="resfastLA_ecoli.csv")


boot <- CNMcalc$"bootstrap triplets"
write.csv(boot,file="boot.csv")

clust <- makeCluster(32)
GLAeasy <- fastboots.GLA(boot, data=GeneExp_merged.t, clust=clust, boots=30, perm=100, cut=3)
stopCluster(clust)
closeAllConnections()
write.csv(GLAeasy,file=" resfasLAboot.csv")


