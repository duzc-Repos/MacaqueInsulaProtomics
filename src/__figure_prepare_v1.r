setwd("C:/Users/86188/Desktop/work/brainnetome/project/help/01_Insular_Proteomics/1.analysis_init/src")


library(ggplot2)
library(clusterProfiler)
library(org.Mmu.eg.db)


posi_gene <- read.table("../res/use_v1/4.posi_cor_gene.txt")$V1
nega_gene <- read.table("../res/use_v1/4.nega_cor_gene.txt")$V1

bg_gene <- unique(sapply(strsplit(rownames(read.table("../data/2.prodata/protomics_expr_processed.tsv")), split = "|", fixed = TRUE),
                  function(x){x[2]}))

posi_go <- enrichGO(posi_gene, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                    ont = "ALL",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(posi_go, color = "p.adjust", showCategory=10)

nega_go <- enrichGO(nega_gene, OrgDb = org.Mmu.eg.db, keyType = "SYMBOL", 
                    ont = "ALL", 
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(nega_go, color = "p.adjust", showCategory=10)

write.table(as.data.frame(posi_go), "../res/use_v1/posi_cor_clusterprofile.tsv", row.names = TRUE, col.names = TRUE, quote = F, sep="\t")
write.table(as.data.frame(nega_go), "../res/use_v1/nega_cor_clusterprofile.tsv", row.names = TRUE, col.names = TRUE, quote = F, sep="\t")

