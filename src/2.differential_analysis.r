setwd("C:/Users/86188/Desktop/work/brainnetome/project/help/01_Insular_Proteomics/1.analysis_init/src")

library(readxl)
library(DEqMS)
library(ggplot2)
library(clusterProfiler)
library(org.Mmu.eg.db)

#######################################
# 0. Load raw proteomics data (after PD processed)
#######################################
rawdata <- read.csv("../data/0.rawdata/wosp22003/2-Input/protein_matrix.csv")
rawinfo <- read_excel("../data/0.rawdata/wosp22003/2-Input/样本实验信息登记表.xlsx")[1:15, ]
expr.log <- read.table("../data/2.prodata/protomics_expr_processed.tsv")


################################################
## differential analysis (Part1)
#  1. Pairwise model   (one vs. one)
#  2. Enrichment model (one vs. all)
################################################

subregion <- as.factor(as.numeric(rawinfo$实验组别))
design <- model.matrix(~0+subregion)
fit <- lmFit(expr.log, design = design)
eb <- eBayes(fit)

## 1. pairwise model
region_combs <- combn(colnames(design), 2)
region_contrasts <- apply(region_combs, 2, function(x) {
  z <- paste(x, collapse = '-')
  makeContrasts(contrasts = z, levels = design)
})
rownames(region_contrasts) <- colnames(design)
colnames(region_contrasts) <- apply(region_combs, 2, paste, collapse="-")
eb_contrasts <- eBayes(contrasts.fit(fit, region_contrasts))
pvals_contrasts <- eb_contrasts$p.value
fdr_contrasts <- apply(pvals_contrasts, 2, p.adjust, 'BH')
pairwise_model_conuting <- data.frame("FDR" = colSums(fdr_contrasts < 0.05),
                                      "Pval-2" = colSums(pvals_contrasts < 0.01),
                                      "Pval-6" = colSums(pvals_contrasts < 1e-6)
)

## 2. enrichment model
eb0_list <- lapply(unique(rawinfo$实验组别), function(x){
  res <- as.factor(as.numeric(rawinfo$实验组别 == x))
  m <- model.matrix(~res)
  tmplm <- lmFit(expr.log, design = m)
  eBayes(tmplm)
})
pvals0_contrasts <- sapply(eb0_list, function(x){ x$p.value[, 2, drop=FALSE] })
fdr0_contrasts <- apply(pvals0_contrasts, 2, p.adjust, 'BH')
tstats0_contrasts <- sapply(eb0_list, function(x) {x$t[, 2, drop = FALSE]})
logFC0_contrasts <- 
  
  enrichment_model_counting <- data.frame("FDR" = colSums(fdr0_contrasts < 0.05),
                                          "Pval-2" = colSums(pvals0_contrasts < 0.01),
                                          "Pval-6" = colSums(pvals0_contrasts < 1e-6)
  )
rownames(enrichment_model_counting) <- colnames(design)

## DE output
de_counting <- rbind(pairwise_model_conuting, enrichment_model_counting)
de_summary <- list()

pairwise_genelist <- c()
for (idx in 1:6){
  tmp <- topTable(eb_contrasts, coef = idx, number = Inf)[c(1, 4, 5)]
  fdrsig <- tmp$adj.P.Val < 0.05
  tmp$reg <- "No"
  tmp$reg[fdrsig & tmp$logFC>0] = "Up"
  tmp$reg[fdrsig & tmp$logFC<0] = "Down"
  de_summary[[idx]] <- tmp
  pairwise_genelist <- append(pairwise_genelist, rownames(tmp)[tmp$reg != "No"], length(pairwise_genelist))
}
pairwise_genelist <- unique(pairwise_genelist)

enrichment_genelist <- c()
for (idx in 1:4){
  tmp <- topTable(eb0_list[[idx]], number = Inf)[c(1, 4, 5)]
  fdrsig <- tmp$adj.P.Val < 0.05
  tmp$reg <- "No"
  tmp$reg[fdrsig & tmp$logFC>0] = "Up"
  de_summary[[idx+6]] <- tmp
  enrichment_genelist <- append(enrichment_genelist, rownames(tmp)[tmp$reg != "No"], length(enrichment_genelist))
  
  de_counting[idx+6, 1] <- sum(tmp$reg == "Up")
  de_counting[idx+6, 2] <- sum(tmp$P.Value<0.01 & tmp$logFC>0)
  de_counting[idx+6, 3] <- sum(tmp$P.Value<1e-6 & tmp$logFC>0)
}
enrichment_genelist <- unique(enrichment_genelist)

combined_genelist <- unique(append(pairwise_genelist, enrichment_genelist, length(pairwise_genelist)))
write.table(combined_genelist, file="../data/2.prodata/combined_genelist.txt", quote = F, row.names = F, col.names = F)

##############################################
# heatmap
##############################################
expr.log <- read.table("../data/2.prodata/protomics_expr_processed.tsv")
diff_pro_list <- read.table("../data/2.prodata/combined_genelist.txt")$V1
diff_pro_expr <- expr.log[diff_pro_list, ]
colnames(diff_pro_expr) <- gsub(".", "-", gsub("X", "INS ", colnames(diff_pro_expr)), fixed = T)

anno_col = data.frame(row.names = colnames(diff_pro_expr), 
                      "Subregion"=paste("INS ", rawinfo$实验组别))
anno_col_color = list(Subregion=c(`INS 1`='"#C9A92E', `INS 2`='#7BBE31', 
                                  `INS 3`='#2972A9', `INS 4`='#973144'))
anno_col_color = list(Subregion=c("#C9A92E", "#7BBE31", "#2972A9", "#973144"))
# "#C9A92E", "#7BBE31", "#2972A9", "#973144"

library(pheatmap)
tmp <- pheatmap(as.matrix(diff_pro_expr), 
         cluster_rows = T, cluster_cols = F, scale = "row",
         clustering_distance_rows = "correlation",
         show_rownames = F, show_colnames = T,
         #clustering_method = "mcquitty",
         annotation_col = anno_col,
         annotation_colors = anno_col_color[0]
         )

cut_tree = cutree(tmp$tree_row, k=2)
cluster1_expr = diff_pro_expr[cut_tree==1, ]
cluster2_expr = diff_pro_expr[cut_tree==2, ]

cluster1_genelist <- unique(sapply(strsplit(rownames(cluster1_expr), "|", fixed = T), function(x){x[2]}))
cluster2_genelist <- unique(sapply(strsplit(rownames(cluster2_expr), "|", fixed = T), function(x){x[2]}))
write.table(cluster1_genelist, file="../data/2.prodata/diff_cluster1_genelist.txt", quote = F, row.names = F, col.names = F)
write.table(cluster2_genelist, file="../data/2.prodata/diff_cluster2_genelist.txt", quote = F, row.names = F, col.names = F)



cluster1_go <- enrichGO(cluster1_genelist, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                       ont = "ALL", 
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(cluster1_go, color = "p.adjust", showCategory=10)
cnetplot(overlap_go, showCategory = 10)

cluster2_go <- enrichGO(cluster2_genelist, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                        ont = "ALL", 
                        pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(cluster2_go, color = "p.adjust", showCategory=10)
cnetplot(overlap_go, showCategory = 10)



# GO
combined_genelist2 <- unique(unlist(strsplit(combined_genelist, split = "|", fixed = T))[seq(2, length(combined_genelist)*2, 2)])
write.table(combined_genelist2, file="../data/2.prodata/combined_genelist2.txt", quote = F, row.names = F, col.names = F)

overlap_go <- enrichGO(combined_genelist2, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                       ont = "ALL", 
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(overlap_go, color = "p.adjust", showCategory=10)
cnetplot(overlap_go, showCategory = 10)




















