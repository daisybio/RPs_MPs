library(tidyr)
library(ggpubr)
library(ggtext)
library(EnhancedVolcano)
library(RColorBrewer)
library(scales)
library(ggnewscale)
library(gt)
library(gtExtras)
library(data.table)

get_specified_colors <- function(){
  return(as.list(c("col_upregulated" = "#B31B21",
              "col_downregulated" = "#0B589C",
              "col_mps" = "blue2",
              "col_rps" = "red3")))
}



myAgeDistrubutionPlot <- function(sample_annotation){
  sample_annotation_age <- sample_annotation[sample_annotation$RPs_MPs == "MPs", c("sample_name", "age")]
  sample_annotation_age$age <- as.numeric(sample_annotation_age$age)
  sample_annotation_age <- sample_annotation_age[order(sample_annotation_age$age, decreasing = FALSE),]
  sample_annotation_age$sample_name <- factor(sample_annotation_age$sample_name, levels = unique(sample_annotation_age$sample_name))
  
  plot <- ggplot(sample_annotation_age, aes(x=sample_name, y=age)) + 
    geom_point(stat="identity") +
    xlab("Patient") + 
    ylab("Age") +
    geom_hline(yintercept = mean(sample_annotation_age$age), color = "red") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  plot
}

myCountAnalysisPlot <- function(counts, data_type, plot_title){
  
  myPiePlot <- function(counts){
    expression_table <- data.frame(total_expression = rowSums(counts))
    expression_table$mp_expression <- rowSums(counts[, grepl('MPs', colnames(counts))])
    expression_table$rp_expression <- rowSums(counts[, grepl('RPs', colnames(counts))])
    
    
    not_expressed_all <- nrow(expression_table[expression_table$total_expression == "0",])
    not_expressed_mps <- nrow(expression_table[expression_table$mp_expression == "0" & expression_table$rp_expression != "0",])
    not_expressed_rps <- nrow(expression_table[expression_table$mp_expression != "0" & expression_table$rp_expression == "0",])
    expressed_all <- nrow(expression_table[expression_table$mp_expression != "0" & expression_table$rp_expression != "0",])
    not_expressed <- data.frame(
      class = factor(c("Expressed in no sample", 
                       "Expressed only in RPs", 
                       "Expressed only in MPs", 
                       "Expressed in MPs and RPs")),
      n = c(not_expressed_all,
            not_expressed_mps,
            not_expressed_rps,
            expressed_all),
      prop = c(round(not_expressed_all/nrow(expression_table)*100, digits = 1),
               round(not_expressed_mps/nrow(expression_table)*100, digits = 1),
               round(not_expressed_rps/nrow(expression_table)*100, digits = 1),
               round(expressed_all/nrow(expression_table)*100, digits = 1))
    )
    not_expressed <- not_expressed %>%
      arrange(desc(class)) %>%
      mutate(text_y = cumsum(prop) - prop/2)
    
    mycols <- c("#EFC000FF", "#868686FF", "blue2", "red3")
    plot <- ggplot(not_expressed, aes(x = "", y = prop, fill = class)) +
      geom_bar(width = 1, stat = "identity", color = "white") +
      coord_polar("y", start = 0)+
      scale_fill_manual(values = mycols) +
      geom_label_repel(aes(label = prop, y = text_y), size=5, show.legend = F, nudge_x = 1) +
      theme_void()
  }
  mySampleHist <- function(counts, data_type){
    expressed_samples <- data.frame(apply(counts, 1, function(x) sum(x != "0")))
    colnames(expressed_samples) <- "count_expressed_samples"
    expressed_samples$circRNA <- rownames(expressed_samples)
    ggplot(expressed_samples) + 
      geom_histogram(mapping = aes(x = count_expressed_samples), bins = ncol(counts)) + 
      theme_classic() + 
      xlab("Number of expressed samples") + 
      ylab(paste0("Number of ", data_type, "s"))
  }
  mySampleCountBoxplots <- function(counts){
    cols <- get_specified_colors()
    
    counts_temp <- data.table::copy(counts)
    counts_long <- melt(setDT(counts_temp, keep.rownames = TRUE), "rn")
    rm(counts_temp)
    
    colnames(counts_long) <- c("id", "sample_name", "counts")
    counts_long$RPs_MPs <- sapply(strsplit(as.character(counts_long$sample_name), "_"), "[", 2)
    counts_long$sample_number <- as.numeric(sapply(strsplit(as.character(counts_long$sample_name), "_"), "[", 1))
    counts_long <- counts_long[order(counts_long$sample_number, counts_long$RPs_MPs),]
    counts_long$sample_name <- factor(counts_long$sample_name, levels=unique(counts_long$sample_name))
    plot <- ggplot(counts_long) + 
      geom_boxplot(aes(x = sample_name, y = counts, fill = RPs_MPs), alpha = 0.65) +
      scale_fill_manual(values = c(cols$col_mps, cols$col_rps)) +
      scale_y_continuous(trans='log2') +
      labs(x = "", fill = "platelet type") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    return(plot)
  }
  
  pie_plot <- myPiePlot(counts)
  hist_plot <- mySampleHist(counts, data_type)
  box_plot <- mySampleCountBoxplots(counts)
  
  plot <- ggarrange(pie_plot, hist_plot, box_plot, ncol = 1, nrow = 3)
  plot <- annotate_figure(plot, top = text_grob(plot_title, 
                                                face = "bold", size = 16))
  return(plot)
}

myTopmiRNABarplot <- function(counts, num_tops = 15){
  counts <- as.data.frame(counts)
  count_sums <- data.table(mirna = rownames(counts), total_count = rowSums(counts))
  if(length(strsplit(count_sums$mirna, "_")[[1]]) == 6){
    count_sums <- count_sums  %>%
      separate(mirna, c("seq", "mir", "mism", "add", "t5", "t3"), "_") %>%
      select("mir", "total_count")
  } else if(length(strsplit(count_sums$mirna, "_")[[1]]) == 1){
    colnames(count_sums) <- c("mir", "total_count")
  }
  unique_count_sums <- data.table(mirna = c(), total_count = c())
  for(mirna in unique(count_sums$mir)){
    temp <- data.table(mirna = mirna, total_count = sum(count_sums[mir == mirna, total_count]))
    unique_count_sums <- rbind(unique_count_sums, temp)
  }
  unique_count_sums <- unique_count_sums[order(-total_count), ]
  other_count_sums <- unique_count_sums[(num_tops+1):nrow(unique_count_sums),]
  unique_count_sums <- unique_count_sums[0:num_tops,]
  unique_count_sums <- rbind(unique_count_sums, data.table(mirna = "others", total_count = sum(other_count_sums$total_count)))
  
  unique_count_sums <- unique_count_sums[order(total_count), ]
  unique_count_sums$mirna <- factor(unique_count_sums$mirna, levels = unique(unique_count_sums$mirna))
  total_expression <- sum(count_sums$total_count)
  
  my_colors <- rev(colorRampPalette(brewer.pal(9, "Set1"))(num_tops+1))
  options(scipen=10000)
  plot <- ggplot(unique_count_sums, aes(x = mirna, y = total_count, fill = mirna)) +
    geom_bar(stat = "identity", alpha = 1, show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = my_colors) +
    theme_classic() +
    labs(x = "", y = "Expression count") + 
    geom_text(aes(label=paste0(round(ceiling(total_count/total_expression * 10000)/100, digits = 2), "%")), size = 3) +
    theme(plot.title = element_text(hjust = 0.5))
}

myBiotypeBarplot <- function(counts, plot_title, num_tops = 10){
  
  get_biotype_table <- function(counts){
    attributes <- c("ensembl_gene_id", "gene_biotype")
    filters <- "ensembl_gene_id"
    mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    #attributes <- listAttributes(mart)
    #filters <- listFilters(mart)
    genes <- getBM(attributes = attributes,
                   filters = filters,
                   values = rownames(counts),
                   mart = mart)
    
    total_counts <- data.frame(ensembl_gene_id = rownames(counts), count = rowSums(counts))
    rownames(total_counts) <- c()
    
    #cat(paste0("Number of gene IDs in Ensembl annotation: ", length(genes$ensembl_gene_id), "\t(unique: ", uniqueN(genes$ensembl_gene_id), ")",
    #           "\nNumber of gene IDs in expression counts: ", length(total_counts$ensembl_gene_id), "\t(unique: ", uniqueN(total_counts$ensembl_gene_id), ")",
    #           "\nNumber of intersecting gene IDs: ", length(intersect(total_counts$ensembl_gene_id, genes$ensembl_gene_id)), "\t(unique: ", uniqueN(intersect(total_counts$ensembl_gene_id, genes$ensembl_gene_id)), ")"))
    
    not_annotated_genes <- data.frame(ensembl_gene_id = total_counts$ensembl_gene_id[!total_counts$ensembl_gene_id %in% genes$ensembl_gene_id])
    not_annotated_genes$gene_biotype <- rep("unknown", nrow(not_annotated_genes))
    genes <- rbind(genes, not_annotated_genes)
    #cat(paste0("Number of gene IDs in Ensembl annotation: ", length(genes$ensembl_gene_id), "\t(unique: ", uniqueN(genes$ensembl_gene_id), ")",
    #           "\nNumber of gene IDs in expression counts: ", length(total_counts$ensembl_gene_id), "\t(unique: ", uniqueN(total_counts$ensembl_gene_id), ")",
    #           "\nNumber of intersecting gene IDs: ", length(intersect(total_counts$ensembl_gene_id, genes$ensembl_gene_id)), "\t(unique: ", uniqueN(intersect(total_counts$ensembl_gene_id, genes$ensembl_gene_id)), ")"))
    
    genes <- merge(genes, total_counts, by="ensembl_gene_id")
    
    biotype_table <- data.frame(gene_biotype = c(), gene_biotype_frequency = c(), gene_expression_count = c())
    for(biotype in unique(genes$gene_biotype)){
      current_biotype_genes <- genes[genes$gene_biotype == biotype,]
      biotype_table <- rbind(biotype_table, data.frame(gene_biotype = biotype, gene_biotype_frequency = nrow(current_biotype_genes), gene_expression_count = sum(current_biotype_genes$count)))
    }
    
    return(biotype_table)
  }
  
  create_biotype_table_list <- function(biotype_table, num_tops){
    biotype_table_processed <- copy(biotype_table)
    biotype_table_processed <- biotype_table_processed[order(-biotype_table_processed$gene_expression_count),] %>%
      mutate(
        gene_biotype_frequency = round(gene_biotype_frequency/sum(biotype_table$gene_biotype_frequency)*100, digits = 5),
        gene_expression_count = round(gene_expression_count/sum(biotype_table$gene_expression_count)*100, digits = 5)
      )
    
    table <- gt(biotype_table_processed) %>%
      gt_theme_538() %>% 
      gt_plt_bar_pct(column = gene_biotype_frequency, scaled = TRUE, fill = "green", background = "lightgreen") %>%
      gt_plt_bar_pct(column = gene_expression_count, scaled = TRUE, fill = "blue", background = "lightblue") %>%
      cols_align("center", contains("scale")) %>%
      gt_highlight_rows(rows = (num_tops+1):nrow(biotype_table_processed), font_weight = "normal", fill = "lavender", alpha = 0.5)
    return(table)
  }
  
  create_biotype_barplot <- function(biotype_table, num_tops){
    picked_biotypes <- intersect(biotype_table$gene_biotype[order(biotype_table$gene_biotype_frequency, decreasing = TRUE)][0:10], 
                                 biotype_table$gene_biotype[order(biotype_table$gene_expression_count, decreasing = TRUE)][0:10])
    biotype_table_processed <- biotype_table[biotype_table$gene_biotype %in% picked_biotypes,]
    biotype_table_notpicked <- biotype_table[!biotype_table$gene_biotype %in% biotype_table_processed$gene_biotype,]
    biotype_table <- rbind(biotype_table_processed, biotype_table_notpicked)
    biotype_table <- biotype_table[order(-biotype_table$gene_expression_count),]
    biotype_table$gene_biotype <- factor(biotype_table$gene_biotype, level = unique(biotype_table$gene_biotype))
    
    biotype_table_processed <- rbind(biotype_table_processed, 
                                     data.table(gene_biotype = "other",
                                                gene_biotype_frequency = sum(biotype_table$gene_biotype_frequency[!biotype_table$gene_biotype %in% picked_biotypes]),
                                                gene_expression_count = sum(biotype_table$gene_expression_count[!biotype_table$gene_biotype %in% picked_biotypes])))
    biotype_table_processed <- biotype_table_processed[order(biotype_table_processed$gene_biotype_frequency),]
    biotype_table_processed$gene_biotype <- factor(biotype_table_processed$gene_biotype, level = unique(biotype_table_processed$gene_biotype))
    total_gene_expression <- sum(biotype_table_processed$gene_expression_count)
    
    plot <- ggplot(biotype_table_processed, aes(x = gene_biotype_frequency, y = gene_biotype, fill = gene_biotype)) +
      geom_bar(stat='identity', position = 'dodge', show.legend = FALSE) +
      scale_fill_brewer(palette="Spectral") +
      theme_classic() +
      theme(text = element_text(size = 10)) +
      labs(title = plot_title, x = "Number of genes", y = "") + 
      geom_text(aes(label=paste0(round(ceiling(gene_expression_count/total_gene_expression * 10000)/100, digits = 2), "%")), 
                size = 3)
  }
  
  create_biotype_piechart <- function(biotype_table, num_tops){
    picked_biotypes <- intersect(biotype_table$gene_biotype[order(biotype_table$gene_biotype_frequency, decreasing = TRUE)][0:10], 
                                 biotype_table$gene_biotype[order(biotype_table$gene_expression_count, decreasing = TRUE)][0:10])
    biotype_table_processed <- biotype_table[biotype_table$gene_biotype %in% picked_biotypes,]
    biotype_table_notpicked <- biotype_table[!biotype_table$gene_biotype %in% biotype_table_processed$gene_biotype,]
    biotype_table <- rbind(biotype_table_processed, biotype_table_notpicked)
    biotype_table <- biotype_table[order(-biotype_table$gene_expression_count),]
    biotype_table$gene_biotype <- factor(biotype_table$gene_biotype, level = unique(biotype_table$gene_biotype))
    
    biotype_table_processed <- rbind(biotype_table_processed, 
                                     data.table(gene_biotype = "other",
                                                gene_biotype_frequency = sum(biotype_table$gene_biotype_frequency[!biotype_table$gene_biotype %in% picked_biotypes]),
                                                gene_expression_count = sum(biotype_table$gene_expression_count[!biotype_table$gene_biotype %in% picked_biotypes])))
    total_gene_expression <- sum(biotype_table_processed$gene_expression_count)
    biotype_table_processed$prop <- round(biotype_table_processed$gene_expression_count/total_gene_expression, digits = 5)
    biotype_table_processed$text_y <- cumsum(biotype_table_processed$prop) -biotype_table_processed$prop/2
    biotype_table_processed$gene_biotype <- paste0(round(biotype_table_processed$prop*100, digits = 2), "% ", biotype_table_processed$gene_biotype)
    biotype_table_processed <- biotype_table_processed %>% arrange(desc(prop))
    biotype_table_processed$gene_biotype <- factor(biotype_table_processed$gene_biotype, level = unique(biotype_table_processed$gene_biotype))
    
    plot <- ggplot(biotype_table_processed, aes(x = "", y = prop, fill = gene_biotype)) +
      geom_bar(width = 1, stat = "identity", color = "white") +
      coord_polar("y", start = 0)+
      scale_fill_brewer(palette="Spectral") +
      theme_void()
    
    return(plot)
  }
  
  
  biotype_table <- get_biotype_table(counts)
  table <- create_biotype_table_list(biotype_table, num_tops)
  plot <- create_biotype_piechart(biotype_table, num_tops)
  
  return(List(biotype_table = table, biotype_plot = plot))
}


myMAplot <- function(deseq_results, fdr = 0.1, fc = 1.5) {
  cols <- get_specified_colors()
  sig <- rep(3, nrow(deseq_results))
  sig[which(deseq_results$padj <= fdr & deseq_results$log2FoldChange < 0 & abs(deseq_results$log2FoldChange) >= log2(fc))] = 2
  sig[which(deseq_results$padj <= fdr & deseq_results$log2FoldChange > 0 & abs(deseq_results$log2FoldChange) >= log2(fc))] = 1
  
  ggmaplot(deseq_results,
           fdr = fdr, fc = fc, size = 1,
           legend = "top", top = 0,
           font.legend = "bold",
           ggtheme = ggplot2::theme_minimal()) + 
    theme(axis.title = element_text(size = 13), legend.text = element_text(size = 11)) +
    geom_hline(yintercept = c(-fdr, fdr) , linetype = "dashed") + 
    scale_color_manual(labels = c(paste0("upregulated miRNAs: ", sum(sig == 1)), 
                                  paste0("downregulated miRNAs: ", sum(sig == 2)),
                                  ""),
                       values = c(cols$col_upregulated, cols$col_downregulated, "darkgrey"),
                       guide = guide_legend(override.aes = list(size = 4,
                                                                color = c(cols$col_upregulated, cols$col_downregulated, NA))))
}

myMAplot_lim <- function(deseq_results, fdr = 0.1, fc = 1.5, ylim) {
  cols <- get_specified_colors()
  sig <- rep(3, nrow(deseq_results))
  sig[which(deseq_results$padj <= fdr & deseq_results$log2FoldChange < 0 & abs(deseq_results$log2FoldChange) >= log2(fc))] = 2
  sig[which(deseq_results$padj <= fdr & deseq_results$log2FoldChange > 0 & abs(deseq_results$log2FoldChange) >= log2(fc))] = 1
  
  ggmaplot(deseq_results,
           fdr = fdr, fc = fc, size = 1,
           legend = "top", top = 0,
           font.legend = "bold",
           ylim = ylim,
           ggtheme = ggplot2::theme_minimal()) +
    theme(axis.title = element_text(size = 13), legend.text = element_text(size = 11)) +
    geom_hline(yintercept = c(-fdr, fdr) , linetype = "dashed") + 
    scale_color_manual(labels = c(paste0("upregulated miRNAs: ", sum(sig == 1)), 
                                  paste0("downregulated miRNAs: ", sum(sig == 2)),
                                  ""),
                       values = c(cols$col_upregulated, cols$col_downregulated, "darkgrey"),
                       guide = guide_legend(override.aes = list(size = 4,
                                                                color = c(cols$col_upregulated, cols$col_downregulated, NA))))
  
}

myPCAplot <- function (object, ...) {
  .local <- function (object, intgroup = "condition", ntop = 500, 
                      returnData = FALSE) {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object))
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df[intgroup, drop = FALSE], 1, paste, collapse = ":"))
    } else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:2]
      return(d)
    }
    
    color_set <- hue_pal()(length(unique(d$group)))
    cols <- get_specified_colors()
    if(length(unique(d$group)) == 4) {
      color_set <- c(cols$col_mps, "#00BFC4", cols$col_rps, "#F8766D")
      #color_set <- c(col_mps, col_rps, "#00BFC4", "#F8766D")
    } else if(length(unique(d$group)) == 2) {
      color_set <- c(cols$col_mps, cols$col_rps)
    } 
    
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point() + 
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
      coord_fixed() +
      scale_color_manual(values = color_set,
                         guide = guide_legend(override.aes = list(size = 3))) +
      theme_bw()
  }
  .local(object, ...)
}

myVolcanoPlot <- function(deseq_results, data_type, pCutoff = 0.01, FCcutoff = 0.5,
                          xlim = c(min(deseq_results[['log2FoldChange']], na.rm = TRUE) - 1.5, max(deseq_results[['log2FoldChange']], na.rm = TRUE) + 1.5), 
                          ylim = c(0, max(-log10(deseq_results[['padj']]), na.rm = TRUE) + 5),
                          lab_selection = NULL) {
  
  cols <- get_specified_colors()
  
  keyvals <- rep('grey25', nrow(deseq_results))
  names(keyvals) <- rep('NS', nrow(deseq_results))
  
  keyvals[which(abs(deseq_results$log2FoldChange) > FCcutoff & deseq_results$padj > pCutoff)] <- 'grey50'
  names(keyvals)[which(abs(deseq_results$log2FoldChange) > FCcutoff & deseq_results$padj > pCutoff)] <- 'Log2FC'
  
  keyvals[which(abs(deseq_results$log2FoldChange) < FCcutoff & deseq_results$padj < pCutoff)] <- 'grey75'
  names(keyvals)[which(abs(deseq_results$log2FoldChange) < FCcutoff & deseq_results$padj < pCutoff)] <- 'p-value'
  
  keyvals[which(deseq_results$log2FoldChange < -FCcutoff & deseq_results$padj < pCutoff)] <- cols$col_downregulated
  names(keyvals)[which(deseq_results$log2FoldChange < -FCcutoff & deseq_results$padj < pCutoff)] <- 'downregulated'
  
  keyvals[which(deseq_results$log2FoldChange > FCcutoff & deseq_results$padj < pCutoff)] <- cols$col_upregulated
  names(keyvals)[which(deseq_results$log2FoldChange > FCcutoff & deseq_results$padj < pCutoff)] <- 'upregulated'
  
  deseq_results_highlighted <- deseq_results[deseq_results$hgnc_symbol %in% lab_selection,]
  deseq_results_highlighted$padj <- -log10(deseq_results_highlighted$padj)
  
  #unique(keyvals)
  #unique(names(keyvals))
  
  if(data_type == "gene"){
    lab <- deseq_results$hgnc_symbol
  } else if(data_type == "miRNA"){
    lab <- rownames(deseq_results)
  }
  
  EnhancedVolcano(deseq_results,
                  lab = lab,
                  selectLab = lab_selection,
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = NULL,
                  subtitle = NULL,
                  xlim = xlim,
                  ylim = ylim,
                  pCutoff = pCutoff,
                  FCcutoff = FCcutoff,
                  pointSize = 2.0,
                  labSize = 3.0,
                  colCustom = keyvals,
                  colAlpha = 0.9,
                  legendPosition = 'top',
                  legendLabSize = 12,
                  legendIconSize = 4.5,
                  drawConnectors = TRUE,
                  widthConnectors = 0.2,
                  colConnectors = "green3"
                  #gridlines.major = TRUE,
                  #gridlines.minor = FALSE,
                  #border = 'partial',
                  #borderWidth = 1.5,
                  #borderColour = 'black'
  ) +
    geom_point(data = deseq_results_highlighted,
               mapping = aes(x=log2FoldChange, y = padj),
               color = "green3")
}



myGObarplot <- function(panther_result_up, panther_result_down, count_threshold, fdr_threshold, num_tops = 10, 
                        padj_threshold, background, gos_of_interest = c()){
  
  prepare_go_result <- function(panther_result_up, panther_result_down, 
                                count_threshold, fdr_threshold, gos_of_interest, num_tops){
    panther_result <- rbind(panther_result_up %>% mutate(enriched_platelet = rep("RPs", n())),
                            panther_result_down %>% mutate(enriched_platelet = rep("MPs", n()))) %>%
      separate(Term, c("GO_ID", "GO_Term"), sep = "~") %>%
      dplyr::select(Category, GO_ID, GO_Term, Count, PValue, `Fold Enrichment`, FDR, enriched_platelet) %>%
      rename(p_value = PValue, fold_enrichment = `Fold Enrichment`) %>%
      rename_all(tolower) %>%
      filter(count >= count_threshold,
             fdr <= fdr_threshold) %>%
      mutate(category = gsub("GOTERM_", "", category),
             category = gsub("_DIRECT", "", category),
             lab_color = ifelse(go_term %like% paste0(gos_of_interest, collapse = "|"), "red", "black"),
             log_pvalue = -log(p_value, base = 10),
             grouped_go_term = paste0(enriched_platelet, "_", go_term),
             grouped_go_term = factor(grouped_go_term, 
                                      levels = grouped_go_term[order(enriched_platelet, fold_enrichment)]),
             fold_enrichment = ifelse(enriched_platelet == "RPs", fold_enrichment, (fold_enrichment * -1)),
             rp_score = ifelse(enriched_platelet == "RPs", log_pvalue, rep(0, num_tops)),
             mp_score = ifelse(enriched_platelet == "RPs", rep(0, num_tops), log_pvalue),
             label_y = ifelse(enriched_platelet == "RPs", -0.2, 0.2),
             label_hjust = ifelse(enriched_platelet == "RPs", 1, 0))
    
    
    panther_result_tops <- setNames(data.table(matrix(nrow = 0, ncol = ncol(panther_result))), 
                                    colnames(panther_result))
    for(go_category in c("BP", "CC", "MF")){
      for(platelet in c("RPs", "MPs")){
        
        panther_result_subset <- panther_result %>%
          filter(category == go_category,
                 enriched_platelet == platelet) %>%
          slice_max(abs(fold_enrichment), n = num_tops)
        
        if(nrow(panther_result_subset) > num_tops){
          temp <- panther_result_subset[go_term %like% paste0(gos_of_interest, collapse = "|"), ]
          panther_result_subset <- panther_result_subset[!go_term %in% temp$go_term, ] %>%
            slice_max(abs(fold_enrichment), n = num_tops-nrow(temp)) %>%
            rbind(temp)
        }
        panther_result_tops <- rbind(panther_result_subset, panther_result_tops)
      }
    }
    
    return(panther_result_tops)
  }
  plot_single_go_barplot <- function(panther_result, rp_limit, mp_limit, go_category){
    terms <- data.frame(category_abb = c("BP", "CC", "MF"),
                        category = c("biological processes", "cellular components", "molecular functions"))
  
    plot <- ggplot(panther_result, aes(x = grouped_go_term, y = fold_enrichment)) + 
      geom_bar(stat = "identity") +

      geom_bar(stat = "identity", aes(fill = mp_score), 
               data = subset(panther_result, mp_score != 0)) +
      scale_fill_gradientn("GO enriched in MPs\nP value (-log10)", colors = brewer.pal(9, "Blues"),
                           limits=mp_limit,
                           guide = guide_colorbar(order = 2)) +
    
      new_scale("fill") +
      geom_bar(stat = "identity", aes(fill = rp_score), 
               data = subset(panther_result, rp_score != 0)) +
      scale_fill_gradientn("GO enriched in RPs\nP value (-log10)", colors = brewer.pal(9, "Reds"),
                           limits=rp_limit,
                           guide = guide_colorbar(order = 1)) +
      
      coord_flip(clip = "off") + 
      theme_classic(base_size = 12) +
      #theme(plot.margin = new_margin) +
      theme(legend.position = "none") +
    
      ggtitle(paste0("GO: ", terms$category[terms$category_abb == go_category])) + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "royalblue1", size = 16)) +
      
      labs(x = "", y = "Fold enrichment") +
      geom_text(aes(y = label_y, label = go_term, hjust = label_hjust, color = lab_color), show.legend = FALSE) +
      scale_colour_manual(values=c("black", "green4")) +
      
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.line = element_line(size = 1.2)) +
      geom_hline(yintercept = 0, size = 1.2)
      #geom_vline(xintercept = 0)
    
    return(plot)
  }
  get_shared_legend <- function(panther_result){
    plot <- ggplot(panther_result, aes(x = grouped_go_term, y = fold_enrichment)) +
      
      geom_bar(stat = "identity", aes(fill = mp_score), 
               data = subset(panther_result, mp_score != 0)) +
      scale_fill_gradientn("GO enriched in MPs\nP value (-log10)", colors = brewer.pal(9, "Blues"),
                           guide = guide_colorbar(order = 2)) +
      
      new_scale("fill") +
      geom_bar(stat = "identity", aes(fill = rp_score), 
               data = subset(panther_result, rp_score != 0)) +
      scale_fill_gradientn("GO enriched in RPs\nP value (-log10)", colors = brewer.pal(9, "Reds"),
                           guide = guide_colorbar(order = 1)) +
      
      theme_classic(base_size = 12) +
      theme(legend.position = "bottom")
    
    step1 <- ggplot_gtable(ggplot_build(plot))
    step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
    step3 <- step1$grobs[[step2]]
    return(step3)
  }
  

  
    panther_result <- prepare_go_result(panther_result_up = panther_result_up,
                                        panther_result_down = panther_result_down,
                                        count_threshold = count_threshold,
                                        fdr_threshold = fdr_threshold,
                                        gos_of_interest = gos_of_interest,
                                        num_tops = num_tops)

    plots <- List()
    mp_limit <- c(min(panther_result[enriched_platelet == "MPs", mp_score]),
                  max(panther_result[enriched_platelet == "MPs", mp_score]))
    rp_limit <- c(min(panther_result[enriched_platelet == "RPs", rp_score]),
                  max(panther_result[enriched_platelet == "RPs", rp_score]))
    for(term in c("BP", "CC", "MF")){
      plot <- plot_single_go_barplot(panther_result = panther_result[category == term,],
                                     mp_limit = mp_limit,
                                     rp_limit = rp_limit,
                                     go_category = term)
      plots[[term]] <- plot
    }
  
    plot_title <- paste0("Sign. padj = ", padj, ", count = ", count, ", fdr = ", fdr, ", Background = ", background)
    footnote <- paste0("Significant genes padj < ", padj_threshold,
                       "\nBackgroung for GO analysis: ", background,
                       "\nGene count per GO Term >= ", count_threshold, 
                       "\nFRD for GO Term <= ", fdr_threshold)
    shared_legend <- get_shared_legend(panther_result)
    plot <- ggarrange(plots$BP, plots$CC, plots$MF, nrow = 1, ncol = 3, legend.grob = shared_legend, legend = "bottom")
    plot <- annotate_figure(plot, 
                            top = text_grob(plot_title, color = "black", face = "bold", size = 18),
                            bottom = text_grob(footnote, color = "darkmagenta", face = "italic", size = 10, 
                                               hjust = 1, x = 1))
  
    return(plot)
}

myGObarplot_notgroupedGOterm <- function(panther_result_up, panther_result_down, count_threshold, fdr_threshold, num_tops = 10, 
                                         padj_threshold, background, gos_of_interest = c()){
  
  prepare_go_result <- function(panther_result, count_threshold, fdr_threshold, gos_of_interest){
    panther_result <- panther_result %>%
      separate(Term, c("GO_ID", "GO_Term"), sep = "~") %>%
      dplyr::select(Category, GO_ID, GO_Term, Count, PValue, `Fold Enrichment`, FDR) %>%
      rename(p_value = PValue, fold_enrichment = `Fold Enrichment`) %>%
      rename_all(tolower) %>%
      filter(count >= count_threshold,
             fdr <= fdr_threshold) %>%
      arrange(desc(fold_enrichment)) %>%
      mutate(category = gsub("GOTERM_", "", category),
             category = gsub("_DIRECT", "", category),
             lab_color = ifelse(go_term %like% paste0(gos_of_interest, collapse = "|"), "red", "black"),
             #go_term = factor(go_term, levels = unique(go_term)),
             log_pvalue = -log(p_value, base = 10))
    
    return(panther_result)
  }
  plot_single_go_barplot <- function(panther_result_up, panther_result_down, go_category, num_tops){
    
    subset_results <- function(panther_result, num_tops){
      panther_result_subset <- panther_result %>%
        filter(category == go_category) %>%
        top_n(num_tops, fold_enrichment) %>%
        arrange(fold_enrichment) %>%
        mutate(go_term = factor(go_term, levels = unique(go_term)),
               fold_enrichment = ifelse(enriched_platelet == "RPs", fold_enrichment, fold_enrichment * -1),
               rp_score = ifelse(enriched_platelet == "RPs", log_pvalue, rep(0, num_tops)),
               mp_score = ifelse(enriched_platelet == "RPs", rep(0, num_tops), log_pvalue),
               label_y = ifelse(enriched_platelet == "RPs", -0.2, 0.2),
               label_hjust = ifelse(enriched_platelet == "RPs", 1, 0))
      
      return(panther_result_subset)
    }
    
    
    terms <- data.frame(category_abb = c("BP", "CC", "MF"),
                        category = c("biological processes", "cellular components", "molecular functions"))
    panther_result <- rbind(subset_results(panther_result = panther_result_down, 
                                           num_tops = num_tops),
                            subset_results(panther_result = panther_result_up, 
                                           num_tops = num_tops))
    
    ggplot(panther_result, aes(x = go_term, y = fold_enrichment)) + 
      geom_bar(stat = "identity") +
      
      geom_bar(stat = "identity", aes(fill = mp_score), 
               data = subset(panther_result, mp_score != 0)) +
      scale_fill_gradientn("GO enriched in MPs\nP value (-log10)", colors = brewer.pal(9, "Blues"),
                           guide = guide_colorbar(order = 2)) +
      
      new_scale("fill") +
      geom_bar(stat = "identity", aes(fill = rp_score), 
               data = subset(panther_result, rp_score != 0)) +
      scale_fill_gradientn("GO enriched in RPs\nP value (-log10)", colors = brewer.pal(9, "Reds"),
                           guide = guide_colorbar(order = 1)) +
      
      
      coord_flip() + 
      theme_classic() +
      theme(legend.position = "bottom") +
      
      ggtitle(paste0("GO: ", terms$category[terms$category_abb == go_category])) + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "royalblue1", size = 16)) +
      
      labs(x = "", y = "Fold enrichment") +
      geom_text(aes(y = label_y, label = go_term, hjust = label_hjust, color = lab_color), size = 3.5, show.legend = FALSE) +
      scale_colour_manual(values=c("black", "green4")) +
      
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.line = element_line(size = 1.2)) +
      geom_hline(yintercept = 0, size = 1.2) + 
      geom_vline(xintercept = 0)
    
  }
  
  panther_result_up <- prepare_go_result(panther_result_up,
                                         count_threshold = count_threshold,
                                         fdr_threshold = fdr_threshold,
                                         gos_of_interest = gos_of_interest) %>%
    mutate(enriched_platelet = rep("RPs", n()))
  
  panther_result_down <- prepare_go_result(panther_result_down, 
                                           count_threshold = count_threshold,
                                           fdr_threshold = fdr_threshold,
                                           gos_of_interest = gos_of_interest) %>%
    mutate(enriched_platelet = rep("MPs", n()))
  
  
  plots <- List()
  for(term in c("BP", "CC", "MF")){
    plot <- plot_single_go_barplot(panther_result_up = panther_result_up, 
                                   panther_result_down = panther_result_down,
                                   go_category = term, 
                                   num_tops = num_tops)
    plots[[term]] <- plot
  }
  
  plot_title <- paste0("Sign. padj = ", padj, ", count = ", count, ", fdr = ", fdr, ", Background = ", background)
  footnote <- paste0("Significant genes padj < ", padj_threshold,
                     "\nBackgroung for GO analysis: ", background,
                     "\nGene count per GO Term >= ", count_threshold, 
                     "\nFRD for GO Term <= ", fdr_threshold)
  plot <- ggarrange(plots$BP, plots$CC, plots$MF, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")
  plot <- annotate_figure(plot, 
                          top = text_grob(plot_title, color = "black", face = "bold", size = 18),
                          bottom = text_grob(footnote, color = "darkmagenta", face = "italic", size = 10, 
                                             hjust = 1, x = 1))
  
  return(plot)
}

myGOsplitBarplot_bicolor <- function(panther_result_up, panther_result_down, count_threshold, fdr_threshold, num_tops = 10, 
                                     padj_threshold, background, gos_of_interest = c()){
  
  prepare_go_result <- function(panther_result, count_threshold, fdr_threshold, gos_of_interest){
    panther_result <- panther_result %>%
      separate(Term, c("GO_ID", "GO_Term"), sep = "~") %>%
      dplyr::select(Category, GO_ID, GO_Term, Count, PValue, `Fold Enrichment`, FDR) %>%
      rename(p_value = PValue, fold_enrichment = `Fold Enrichment`) %>%
      rename_all(tolower) %>%
      filter(count >= count_threshold,
             fdr <= fdr_threshold) %>%
      arrange(desc(fold_enrichment)) %>%
      mutate(category = gsub("GOTERM_", "", category),
             category = gsub("_DIRECT", "", category),
             lab_color = ifelse(go_term %like% paste0(gos_of_interest, collapse = "|"), "red", "black"),
             go_term = factor(go_term, levels = unique(go_term)),
             log_pvalue = -log(p_value, base = 10))
    
    return(panther_result)
  }
  plot_single_go_barplot <- function(panther_result, go_category, enriched_platelet, num_tops){
    panther_result_subset <- panther_result[category == go_category, ]
    panther_result_subset <- panther_result_subset[0:min(nrow(panther_result_subset), num_tops),]
    panther_result_subset <- panther_result_subset[order(panther_result_subset$fold_enrichment),]
    
    cols <- brewer.pal(9, "Reds")
    y_axis_position = "bottom"
    if(enriched_platelet == "MPs"){
      cols <- brewer.pal(9, "Blues")
      y_axis_position = "top"
      panther_result_subset$fold_enrichment <- panther_result_subset$fold_enrichment * (-1)
    }
    
    panther_result_subset$go_term <- factor(panther_result_subset$go_term,
                                            level = unique(panther_result_subset$go_term))
    
    plot <- ggplot(panther_result_subset, aes(x = go_term, y = fold_enrichment, fill = log_pvalue)) + 
      geom_bar(stat = "identity") +
      scale_fill_gradientn(colors = cols) +
      labs(x = "", y = "Fold enrichment")
    if(nrow(panther_result_subset) > 0){
      plot <- plot + 
        scale_x_discrete(position = y_axis_position)
    }
    plot <- plot +
      ggtitle(paste0("GO enriched in ", enriched_platelet)) + 
      theme(plot.title = element_text(hjust = 1)) +
      coord_flip() + 
      theme_classic() +
      theme(axis.text.y = ggtext::element_markdown(colour = panther_result_subset$lab_color), 
            legend.position = "bottom")
    
    return(plot)
  }
  
  panther_result_up <- prepare_go_result(panther_result_up, 
                                         count_threshold = count_threshold,
                                         fdr_threshold = fdr_threshold,
                                         gos_of_interest = gos_of_interest)
  panther_result_down <- prepare_go_result(panther_result_down, 
                                           count_threshold = count_threshold,
                                           fdr_threshold = fdr_threshold,
                                           gos_of_interest = gos_of_interest)
  
  
  plots <- List()
  terms <- data.frame(category_abb = c("BP", "CC", "MF"),
                      category = c("biological processes", "cellular components", "molecular functions"))
  
  for(term in terms$category_abb){
    plot_up <- plot_single_go_barplot(panther_result = panther_result_up, 
                                      go_category = term, 
                                      enriched_platelet = "RPs",
                                      num_tops = num_tops)
    plot_down <- plot_single_go_barplot(panther_result = panther_result_down, 
                                        go_category = term, 
                                        enriched_platelet = "MPs",
                                        num_tops = num_tops)
    plot <- ggarrange(plot_up, plot_down, nrow = 1, ncol = 2)
    plot <- annotate_figure(plot, left = text_grob(terms$category[terms$category_abb == term], 
                                                   color = "black", face = "bold", size = 16, rot = 90))
    plots[[term]] <- plot
  }
  
  plot_title <- paste0("Top ", num_tops, " enriched GOs")
  footnote <- paste0("Significant genes padj < ", padj_threshold,
                     "\nBackgroung for GO analysis: ", background,
                     "\nGene count per GO Term >= ", count_threshold, 
                     "\nFRD for GO Term <= ", fdr_threshold)
  plot <- ggarrange(plots$BP, plots$CC, plots$MF, nrow = 3, ncol = 1)
  plot <- annotate_figure(plot, 
                          top = text_grob(plot_title, color = "black", face = "bold", size = 18),
                          bottom = text_grob(footnote, color = "green4", face = "italic", size = 10, 
                                             hjust = 1, x = 1))
  
  return(plot)
}

myGOsplitBarplot_allblue <- function(panther_result_up, panther_result_down, count_threshold, fdr_threshold, num_tops = 10, 
                                     padj_threshold, background, gos_of_interest = c()){
  
  prepare_go_result <- function(panther_result, count_threshold, fdr_threshold, enriched_platelet, gos_of_interest){
    panther_result <- panther_result %>%
      separate(Term, c("GO_ID", "GO_Term"), sep = "~") %>%
      dplyr::select(Category, GO_ID, GO_Term, Count, PValue, `Fold Enrichment`, FDR) %>%
      rename(p_value = PValue, fold_enrichment = `Fold Enrichment`) %>%
      rename_all(tolower) %>%
      filter(count >= count_threshold,
             fdr <= fdr_threshold) %>%
      arrange(desc(fold_enrichment)) %>%
      mutate(category = gsub("GOTERM_", "", category),
             category = gsub("_DIRECT", "", category),
             lab_color = ifelse(go_term %like% paste0(gos_of_interest, collapse = "|"), "red", "black"),
             go_term = factor(go_term, levels = unique(go_term)),
             log_pvalue = -log(p_value, base = 10))
    
    return(panther_result)
  }
  plot_single_go_barplot <- function(panther_result, go_category, num_tops){
    panther_result_subset <- panther_result[category == go_category, ]
    panther_result_subset <- panther_result_subset[0:min(nrow(panther_result_subset), num_tops),]
    panther_result_subset <- panther_result_subset[order(panther_result_subset$fold_enrichment),]
    panther_result_subset$go_term <- factor(panther_result_subset$go_term,
                                            level = unique(panther_result_subset$go_term))
    
    plot <- ggplot(panther_result_subset, aes(x = go_term, y = fold_enrichment, fill = log_pvalue)) + 
      geom_bar(stat = "identity") +
      scale_fill_gradientn(colors = cols) +
      labs(x = "", y = "Fold enrichment")
      ggtitle(paste0("GO enriched in ", enriched_platelet)) + 
      theme(plot.title = element_text(hjust = 1)) +
      coord_flip() + 
      theme_classic() +
      theme(axis.text.y = ggtext::element_markdown(colour = panther_result_subset$lab_color), 
            legend.position = "bottom")
    
    return(plot)
  }
  
  panther_result_up <- prepare_go_result(panther_result_up, 
                                         count_threshold = count_threshold,
                                         fdr_threshold = fdr_threshold,
                                         gos_of_interest = gos_of_interest)
  panther_result_down <- prepare_go_result(panther_result_down, 
                                           count_threshold = count_threshold,
                                           fdr_threshold = fdr_threshold,
                                           gos_of_interest = gos_of_interest)
  
  
  plots <- List()
  terms <- data.frame(category_abb = c("BP", "CC", "MF"),
                      category = c("biological processes", "cellular components", "molecular functions"))
  
  for(term in terms$category_abb){
    plot_up <- plot_single_go_barplot(panther_result = panther_result_up, 
                                      go_category = term, 
                                      enriched_platelet = "RPs",
                                      num_tops = num_tops)
    plot_down <- plot_single_go_barplot(panther_result = panther_result_down, 
                                        go_category = term, 
                                        enriched_platelet = "MPs",
                                        num_tops = num_tops)
    plot <- ggarrange(plot_up, plot_down, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")
    plot <- annotate_figure(plot, left = text_grob(terms$category[terms$category_abb == term], 
                                                   color = "black", face = "bold", size = 16, rot = 90))
    plots[[term]] <- plot
  }
  
  plot_title <- paste0("Top ", num_tops, " enriched GOs")
  footnote <- paste0("Significant genes padj < ", padj_threshold,
                     "\nBackgroung for GO analysis: ", background,
                     "\nGene count per GO Term >= ", count_threshold, 
                     "\nFRD for GO Term <= ", fdr_threshold)
  plot <- ggarrange(plots$BP, plots$CC, plots$MF, nrow = 3, ncol = 1)
  plot <- annotate_figure(plot, 
                          top = text_grob(plot_title, color = "black", face = "bold", size = 18),
                          bottom = text_grob(footnote, color = "green4", face = "italic", size = 10, 
                                             hjust = 1, x = 1))
  
  return(plot)
}
