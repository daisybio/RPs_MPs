library(DESeq2)
library(EnhancedVolcano)
library(reshape2)
library(dplyr)
library(ggplot2)
library(stats)

# load input data for analysis
load("circRNA_analysis_input.RData")

### Script to re-create Fig. 5D
dds.circ <- DESeq2::DESeqDataSetFromMatrix(countData = round(circ_expr),
                                           colData = samplesheet,
                                           design = ~ condition)
dds.circ <- DESeq2::DESeq(dds.circ)
res.circ <- DESeq2::results(dds.circ)
# sort by p-value
res.circ <- res.circ[order(res.circ$padj),]
# create summary
DESeq2::summary(res.circ)
# define color map
ann_col = c("RPs"="#DE3654", "MPs"="#338ACC")

results <- res.circ
padj <- 0.05
log2FC <- 0
# filter for significance thresholds
signif.hits <- results[!is.na(results$padj) &
                         results$padj<as.double(padj) &
                         abs(results$log2FoldChange) > log2FC,]

host_symbols <- circ_RNAs[rownames(results),"gene_symbol"]
results$host_symbols <- host_symbols
# save DESeq2 results to file
write.csv(results, "DESeq2_results.csv", quote = F)

## Fig. 5D pca
deseq_vst <- DESeq2::vst(dds.circ, blind = F, nsub = 100)
PCA_data <- DESeq2::plotPCA(deseq_vst, intgroup = "condition", returnData = T)
percentVar <- round(100 * attr(PCA_data, "percentVar"))

p <- ggplot(PCA_data, aes(x = PC1, y = PC2, color = group)) + 
  geom_point(size = 3) +
  scale_color_manual(values = ann_col) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  ggtitle("PCA plot of circRNA expression")
ggsave(filename = "Fig_5D_PCA.pdf", dpi = 300, plot = p, width = 8, height = 8)

## Fig. 5D volcano
volcano <- EnhancedVolcano(results,
                           lab = host_symbols,
                           x = 'log2FoldChange',
                           y = 'padj',
                           FCcutoff = 1, pCutoff = 0.05,
                           title = "Volcano plot of circRNAs (labeled as host gene)",
                           subtitle = "Condition RPs vs MPs", 
                           pointSize = 3.0,
                           labSize = 3.0,
                           boxedLabels = T,
                           colAlpha = 4/5,
                           drawConnectors = T,
                           widthConnectors = 0.6,
                           colConnectors = 'black',
                           max.overlaps = 10,
                           # max.overlaps = Inf
)
# save plot to pdf
ggsave(filename = "Fig_5D_RP_vs_MP_volcano.pdf",
       dpi = 300, plot = volcano, width = 8, height = 8)


#### Supplementary plots
# heatmaps
library(pheatmap)
pc = 1
circ.counts = log2(counts(dds.circ, normalized=T)[rownames(signif.hits),]+pc)
ann_df = data.frame(samplesheet[,c(1,5)], row.names = 1)

# all de circRNAs
sh = signif.hits %>% as.data.frame %>% arrange(desc(abs(log2FoldChange)))
# top 25
n = 25
mt = circ.counts[rownames(sh %>% head(n)),]
p <- pheatmap(mt, show_rownames = F, treeheight_row = 0,
              annotation_col = ann_df, annotation_colors = list(condition = ann_col),
              cutree_cols = 2, show_colnames = F)
ggsave(filename = paste0("hm_", n, ".pdf"),
       dpi = 300, plot = p, width = 8, height = 8)

# top 50
n = 50
mt = circ.counts[rownames(sh %>% head(n)),]
p <- pheatmap(mt, show_rownames = F, treeheight_row = 0,
              annotation_col = ann_df, annotation_colors = list(condition = ann_col),
              cutree_cols = 2, show_colnames = F)
ggsave(filename = paste0("hm_", n, ".pdf"),
       dpi = 300, plot = p, width = 8, height = 8)
# top 100
n = 100
mt = circ.counts[rownames(sh %>% head(n)),]
p <- pheatmap(mt, show_rownames = F, treeheight_row = 0,
              annotation_col = ann_df, annotation_colors = list(condition = ann_col),
              cutree_cols = 2, show_colnames = F)
ggsave(filename = paste0("hm_", n, ".pdf"),
       dpi = 300, plot = p, width = 8, height = 8)
# all
n = "All"
mt = circ.counts[rownames(sh),]
p <- pheatmap(mt, show_rownames = F, treeheight_row = 0,
              annotation_col = ann_df, annotation_colors = list(condition = ann_col),
              cutree_cols = 2, show_colnames = F)
ggsave(filename = paste0("hm_", n, ".pdf"),
       dpi = 300, plot = p, width = 8, height = 8)

# extract detailed data
write.table(signif.hits, file = "signif_de.tsv", sep = "\t", quote = F, row.names = T)
RP_MP_down <- volcano$data %>% filter(Sig == "FC_P", log2FoldChange < 0) %>% arrange(log2FoldChange)
RP_MP_up <- volcano$data %>% filter(Sig == "FC_P", log2FoldChange > 0) %>% arrange(desc(log2FoldChange))

write.table(RP_MP_down, file = "RP_MP_down.tsv", sep = "\t", quote = F, row.names = T)
write.table(RP_MP_up, file = "RP_MP_up.tsv", sep = "\t", quote = F, row.names = T)

RP_MP_DE <- volcano$data %>% filter(Sig == "FC_P")
RP_MP_DE$ID <- rownames(RP_MP_DE)

## Supplementary Figure 6 (MA)
library(ggpubr)
results$ps <- ((log2(results$baseMean)/10) + abs(results$log2FoldChange) / 2) * 1.25

top = 5
top_labels <- rbind(RP_MP_up[1:top,],RP_MP_down[1:top,])
top_labels$lab <- paste0("bold('",top_labels$host_symbols,"')","~'",rownames(top_labels),"'")

results$lab <- paste0("bold('",results$host_symbols,"')~'",rownames(results),"'")

ma <- ggmaplot(results, fdr = 0.05, fc = 1, size = results$ps,
               main = "MA plot of circRNAs (labeled as host gene)",
               genenames = results$lab, legend = "top",
               label.rectangle = T, top = 0, ggtheme = theme_minimal()) +
  theme(text = element_text(size=14)) +
  geom_label_repel(data = as.data.frame(results[rp_spec,]),
                   aes(x = log2(baseMean+1), y = log2FoldChange, label = lab),
                   max.overlaps = Inf, box.padding = 1,
                   parse = T)
ggsave(filename = "Suppl_Fig_6.pdf",
       dpi = 300, plot = ma, width = 10, height = 10)

# save important input data for analysis
save(circ_RNAs, circ_expr, samples, samplesheet, rp_spec, file = "circRNA_analysis_input.RData")
