############# ------------- Overview Figure CyTOF ------------- #############

# Load packages
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(SingleCellExperiment)
  library(CATALYST)
  library(ggrepel)
  library(ggpubr)
  library(RColorBrewer)
  library(ggplotify)
  library(pheatmap)
  library(patchwork)
  source("cyanus_functions/prep_functions.R")
  source("cyanus_functions/de_functions.R")
  source("functions.R")
  library(BSDA)
})

path_to_data <- "/nfs/data/Bongiovanni-KrdIsar-platelets/Cyanus_RPsMPs/data/sce_objects"

#analysis_state <- "baseline" # or "stimulated"
analysis_state <- "stimulated" 


# For reproducibility
set.seed(1234)

# Fix colors
colors <- c("#FF1F5B", "#009ADE", "#C4C4C4")
names(colors) <- c("RP", "MP", "rest")


# Read data
mapping <- c(paste0("sce_", analysis_state, "_original_RPs_MPs_rest.rds"), 
             paste0("sce_", analysis_state, "_CD42b_RPs_MPs_rest.rds"), 
             paste0("sce_", analysis_state, "_DNA2_RPs_MPs_rest.rds"))
names(mapping) <- c('Original', 'CD42b', 'DNA2')

sce <- readRDS(file.path(path_to_data, mapping['Original']))
sce_CD42b <- readRDS(file.path(path_to_data, mapping['CD42b'])) 
sce_DNA2 <- readRDS(file.path(path_to_data, mapping['DNA2']))


######## Dimensionality Reduction ########

# tSNE colored by type (RPs, MPs, rest)
sce <- runDR(sce,
             dr = c("TSNE"),
             cells = 1000,
             features = "type",
             assay = "exprs")

tsne_RPs_MPs_rest_plot <- plotDR(sce,
                    dr = "TSNE",
                    color_by = "type") +
  scale_color_manual(name = "", values = colors) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(direction = "horizontal", title = "", override.aes = list(size = 5)))

ggsave(paste0("plots/TSNE_", analysis_state, "_RPs_MPs_rest.png"), width = 4, height = 4, dpi = 300)

# UMAP colored by type (RPs, MPs, rest)
sce <- runDR(sce,
             dr = c("UMAP"),
             cells = 1000,
             features = "type",
             assay = "exprs")

umap_RPs_MPs_rest_plot <- plotDR(sce,
                                 dr = "UMAP",
                                 color_by = "type") +
  scale_color_manual(name = "", values = colors) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(direction = "horizontal", title = "", override.aes = list(size = 5)))

ggsave(paste0("plots/UMAP_", analysis_state, "_RPs_MPs_rest.png"), width = 4, height = 4, dpi = 300)


# tSNE colored by expression, separated by RP and MP

sce_RPs_MPs <- sce[,sce$type != "rest"]

markers1 <- c("CD62P", "CD63", "GPVI", "PAR1", "CD40")
markers2 <- c("CD42a", "CD42b", "PEAR", "CD31", "PAC1")

tsne_expression_plot_1 <- plotDR(sce_RPs_MPs,
                               dr = c("TSNE"),
                               color_by = markers1,
                               facet_by = "type",
                               ncol = 4) +
  theme(legend.position = "bottom", legend.title = element_text(hjust = 1)) +
  guides(color = guide_colorbar(title = "Scaled Expression", direction = "horizontal", title.position = "left", title.vjust = 0.8))

tsne_expression_plot_2 <- plotDR(sce_RPs_MPs,
                                 dr = "TSNE",
                                 color_by = markers2,
                                 facet_by = "type",
                                 ncol = 4) +
  theme(legend.position = "bottom", legend.title = element_text(hjust = 1)) +
  guides(color = guide_colorbar(title = "Scaled Expression", direction = "horizontal", title.position = "left", title.vjust = 0.8))

tsne_expression_plots <- ggarrange(tsne_expression_plot_1, tsne_expression_plot_2, ncol = 1, common.legend = TRUE, labels = NULL, legend = "bottom")

ggsave(paste0("plots/TSNE_", analysis_state, "_marker_expression.png"), width = 12, height = 8, dpi = 300)

# UMAP colored by expression, separated by RP and MP
umap_expression_plot_1 <- plotDR(sce_RPs_MPs,
                                 dr = c("UMAP"),
                                 color_by = markers1,
                                 facet_by = "type",
                                 ncol = 4) +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(title = "Scaled Expression", direction = "horizontal", title.position = "left", title.vjust = 0.8))

umap_expression_plot_2 <- plotDR(sce_RPs_MPs,
                                 dr = "UMAP",
                                 color_by = markers2,
                                 facet_by = "type",
                                 ncol = 4) +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(title = "Scaled Expression", direction = "horizontal", title.position = "left", title.vjust = 0.8))

umap_expression_plots <- ggarrange(umap_expression_plot_1, umap_expression_plot_2, ncol = 1, common.legend = TRUE, labels = NULL, legend = "bottom")

ggsave(paste0("plots/UMAP_", analysis_state, "_marker_expression.png"), width = 12, height = 8, dpi = 300)


# Paired (patient-wise) analysis of marker expressions

# On original data
df_medians_original <- paired_boxes(sce, 'Original', paste0("tables/median_table_with_paired_results_", analysis_state, "_sign_test_original.csv"), method = "sign.test")

# Check differences
df_differences_original <- df_medians_original[, c("marker", "patient_id", "group", "Expression")]
# Remove gating markers
df_differences_original <- df_differences_original[!df_differences_original$marker %in% c("DNA1", "DNA2", "CD45"),]
df_differences_original <- dcast(df_differences_original, marker + patient_id ~ group, value.var = "Expression")
df_differences_original$differences <- df_differences_original$RP - df_differences_original$MP

medians_df <- df_differences_original %>%
  group_by(marker) %>%
  dplyr::summarise(median_diff = median(differences, na.rm = TRUE))

diff_density <- ggplot(df_differences_original, aes(x = differences)) + 
  geom_density(color = "red") + facet_wrap(~marker, scales="free_y") + geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_histogram(aes(y = ..density..), fill = "lightblue", alpha = 0.5, bins = 30) +
  geom_vline(data = medians_df, aes(xintercept = median_diff),
             color = "blue", linetype = "solid") +  # Median line +
  theme_bw() + labs(x = "Differences of RPs and MPs MSI", y = "Density")

diff_boxplot <- ggplot(df_differences_original, aes(x = marker, y = differences)) + 
  geom_boxplot(fill = "lightblue", alpha= 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") + theme_bw() + 
  labs(x = "Marker", y = "Differences of RPs and MPs MSI")

diff_density / diff_boxplot + plot_annotation(tag_levels = "A") + plot_layout(heights = c(1.5,1)) & theme(plot.tag = element_text(face = "bold"))
ggsave(paste0("plots/differences_", analysis_state, "_original.png"), width = 12, height = 10, dpi = 300)


symmetry_coeff_df <- df_differences_original %>%
  group_by(marker) %>%
  dplyr::summarise(median = median(differences, na.rm = TRUE), 
                   mean  = mean(differences, na.rm = TRUE), 
                   sd  = sd(differences, na.rm = TRUE),
                   q1  = quantile(differences, 0.25, na.rm = TRUE),
                   q3  = quantile(differences, 0.75, na.rm = TRUE))

symmetry_coeff_df$pearsons_skew <- (3 * (symmetry_coeff_df$mean - symmetry_coeff_df$median)) / symmetry_coeff_df$sd
symmetry_coeff_df$bowley <- (symmetry_coeff_df$q3 + symmetry_coeff_df$q1 - (2* symmetry_coeff_df$median)) / (symmetry_coeff_df$q3 - symmetry_coeff_df$q1)
symmetry_coeff_df <- symmetry_coeff_df[, c("marker", "pearsons_skew", "bowley")]
symmetry_coeff_df$pearsons_skew <- round(symmetry_coeff_df$pearsons_skew, 3)
symmetry_coeff_df$bowley <- round(symmetry_coeff_df$bowley, 3)
# Save the median table
write.csv(symmetry_coeff_df, paste0("tables/symmetry_coefficients_table_with_paired_results_", analysis_state, "_original.csv"), row.names = FALSE)


df <- df_medians_original[!df_medians_original$marker %in% c("CD45", "DNA1", "DNA2", "CD47"),]
violins_original <- ggplot(df[signif != ""], aes(x = group, y = Expression, color = group, fill = group))+
  geom_violin(alpha = 0.3)+
  geom_point()+
  geom_line(aes(group = patient_id), color = '#C4C4C4')+
  scale_color_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  scale_fill_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  facet_wrap(~marker_title, scales = 'free', ncol = 8)+
  theme_minimal()+
  theme(legend.position = 'none', axis.title.x=element_blank(), axis.text.x = element_blank(), strip.text.x = element_text(face = "bold"))
ggsave(paste0("plots/paired_boxes_", analysis_state, "_original.png"), width = 12, height = 4, dpi = 300)

df_medians_original_wilcoxon <- paired_boxes(sce, 'Original', paste0("tables/median_table_with_paired_results_", analysis_state, "_wilcoxon_original.csv"), method = "wilcoxon")


# Normalized by size
df_medians_CD42b <- paired_boxes(sce_CD42b, 'CD42b', paste0("tables/median_table_with_paired_results_", analysis_state, "_sign_test_CD42b.csv"), method = "sign.test")


df <- df_medians_CD42b[!df_medians_CD42b$marker %in% c("CD45", "DNA1", "DNA2"),]
violins_CD42b <- ggplot(df[signif != ""], aes(x = group, y = Expression, color = group, fill = group))+
  geom_violin(alpha = 0.3)+
  geom_point()+
  geom_line(aes(group = patient_id), color = '#C4C4C4')+
  scale_color_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  scale_fill_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  facet_wrap(~marker_title, scales = 'free', ncol = 4)+
  theme_minimal()+
  theme(legend.position = 'none', axis.title.x=element_blank(), axis.text.x = element_blank(), strip.text.x = element_text(face = "bold"))
ggsave(paste0("plots/paired_boxes_", analysis_state, "_CD42b.png"), width = 6, height = 5, dpi = 300)

df_medians_CD42b_wilcoxon <- paired_boxes(sce_CD42b, 'CD42b', paste0("tables/median_table_with_paired_results_", analysis_state, "_wilcoxon_CD42b.csv"), method = "wilcoxon")


# Normalized by RNA
df_medians_DNA2 <- paired_boxes(sce_DNA2, 'DNA2', paste0("tables/median_table_with_paired_results_", analysis_state, "_sign_test_DNA2.csv"), method = "sign.test")


df <- df_medians_DNA2[!df_medians_DNA2$marker %in% c("CD45"),]
violins_DNA2 <- ggplot(df[signif != ""], aes(x = group, y = Expression, color = group, fill = group))+
  geom_violin(alpha = 0.3)+
  geom_point()+
  geom_line(aes(group = patient_id), color = '#C4C4C4')+
  scale_color_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  scale_fill_manual(values = c("MP" = "#009ADE", "RP" = "#FF1F5B"))+
  facet_wrap(~marker_title, scales = 'free', ncol = 4)+
  theme_minimal()+
  theme(legend.position = 'none', axis.title.x=element_blank(), axis.text.x = element_blank(), strip.text.x = element_text(face = "bold"))
ggsave(paste0("plots/paired_boxes_", analysis_state, "_DNA2.png"), width = 6, height = 5, dpi = 300)

df_medians_DNA2_wilcoxon <- paired_boxes(sce_DNA2, 'DNA2', paste0("tables/median_table_with_paired_results_", analysis_state, "_wilcoxon_DNA2.csv"), method = "wilcoxon")


