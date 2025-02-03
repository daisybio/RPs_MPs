library(biomaRt)
library(data.table)
library(readr)
library(dplyr)
library(isomiRs)

#setwd("D:/Uni/Hiwi_ExBioLab/Rechts-der-Isar/Material/analyses_in_R/CCSonly_withoutSample6_filtered_restructured")


get_raw_gene_counts <- function(metadata){
  
  process_raw_counts <- function(metadata) {
    raw_counts <- read.table("../data/gene_counts/salmon.merged.gene_counts.tsv", header = TRUE, sep = "\t", stringsAsFactors = F)
    sample_names<- read.table("../data/sample_annotation/sample_names.tsv", header = TRUE, sep = "\t", stringsAsFactors = F)
    sample_names$sample <- gsub('-', '_', sample_names$sample)
    
    for( i in 2:ncol(raw_counts)){
      group_correspondence <- colnames(raw_counts)[i]
      sample_correspondence <- sample_names[sample_names$group == group_correspondence, 'sample']
      colnames(raw_counts)[i] <- sample_correspondence
    }
    
    row.names(raw_counts) <- raw_counts$gene_id
    raw_counts <- raw_counts[,-1]
    raw_counts <- subset(raw_counts, select = intersect(colnames(raw_counts), metadata$sample_name_huge))
    
    raw_counts <- raw_counts[, order(names(raw_counts))]
    metadata <- metadata[order(metadata$sample_name_huge), ]
    colnames(raw_counts) <- metadata$sample_name
    raw_counts <- raw_counts[, order(names(raw_counts))]
    metadata <- metadata[order(metadata$sample_number, metadata$RPs_MPs), ]
    
    
    # check if metadata and counts samples are in the same order
    all(rownames(metadata) %in% colnames(raw_counts))
    all(rownames(metadata) == colnames(raw_counts))
    return(raw_counts)
  }
  ensembl_gene_id_to_gene_symbol <- function(ensembl_ids){
    #fix to SSL error
    httr::set_config(httr::config(ssl_verifypeer = FALSE))
    
    mart <- useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl")
    
    genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"), 
                   filters = "ensembl_gene_id",
                   values = unique(ensembl_ids),
                   mart = mart)
    
    genes[which(genes$hgnc_symbol == ""),"hgnc_symbol"] <- genes[which(genes$hgnc_symbol == ""),"ensembl_gene_id"]
    
    result <- dplyr::left_join(data.frame(ensg = ensembl_ids), genes, by = c("ensg" = "ensembl_gene_id"))
    
    return(result)
  }
  
    raw_counts <- process_raw_counts(metadata)
    
    # replace gene ids with gene symbols
    gene_symbol <- ensembl_gene_id_to_gene_symbol(rownames(raw_counts))
  
    # concatenate gene symbols for duplicated gene ids
    duplicated_genes = gene_symbol$ensg[duplicated(gene_symbol$ensg)]
    id2symbol <- setDT(gene_symbol)[, list(hgnc_symbol = paste(hgnc_symbol, collapse = '/')), by = c('ensg', 'gene_biotype')]
  
    # set gene_name to symbol if available, otherwise to gene id
    id2symbol[id2symbol$hgnc_symbol == 'NA', 'hgnc_symbol'] <- id2symbol[id2symbol$hgnc_symbol == 'NA', 'ensg']
    counts <- round(raw_counts)
    
    return(List(counts = counts, id2symbol = id2symbol))
}

get_mirtop_mirna_counts <- function(metadata) {
  
  mirtop_data <- data.table(read_tsv("../data/mirtop_counts/mirtop_rawData.tsv") %>% dplyr::select(-starts_with("lib-smFRI")))
  mirtop_data <- mirtop_data %>% pivot_longer(cols = starts_with("sm"), names_sep = "_", names_to = c("sample", "r", "no"), values_to = "count")
  mirtop_data <- mirtop_data %>% group_by(seq, mir, mism, add, t5, t3, sample) %>% 
    summarize(count = sum(count)) %>%
    mutate(sample_name_huge = str_extract(sample, "DBO-[0-9]+"))
  mirtop_data$sample_name_huge <- gsub("-", "_", mirtop_data$sample_name_huge)
  
  mirtop_data <- merge(mirtop_data, metadata, by="sample_name_huge")
  mirtop_data <- mirtop_data[, c("seq", "mir", "mism", "add", "t5", "t3", "sample", "count", "sample_name")]
  mirtop_data <- mirtop_data %>% 
    pivot_wider(id_cols = c("seq", "mir", "mism", "add", "t5", "t3"),
                names_from = "sample_name",
                values_from = "count")
  
  isomir_data <- IsomirDataSeqFromMirtop(mirtop = mirtop_data, 
                                         coldata = metadata)
  dds <- isoDE(isomir_data, formula = ~patient + RPs_MPs)
  raw_counts <- data.frame(counts(dds))
  colnames(raw_counts) <- gsub("X", "", colnames(raw_counts))
  
  # check if metadata and counts samples are in the same order
  if(!all(rownames(metadata) %in% colnames(raw_counts))){
    cat("Samples in metadata and raw_counts do not match")
  }
  if(!all(rownames(metadata) == colnames(raw_counts))){
    raw_counts <- raw_counts[, order(names(raw_counts))]
  }
  return(raw_counts)
}








