source('QC.R')
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

