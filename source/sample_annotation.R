library(tidyr)
library(dplyr)
library(stringr)
library(readxl)
library(data.table)

get_sample_annotation <- function(picked_disease = "all", remove_samples = c()){
  
  combine_sample_annotation <- function(){
    sample_annotation_old <- as.data.table(read.delim2("../data/sample_annotation/sample_annotation_old.tsv")) %>%
      mutate(sample_name_huge = str_extract(string = `Huge.ID`, pattern = "DBO-[0-9]*"),
             sample_name_new = New.ID,
             RPs_MPs = RPs.MPs,
             patient = str_extract(string = MUC.ID, pattern = "TO [0-9]+"))
    sample_annotation_old <- sample_annotation_old[,c("donor","RPs_MPs","disease","sample_name_huge","sample_name_new","patient", "gender", "treatment")]
    
    
    sample_annotation_new <- as.data.table(read_excel("../data/sample_annotation/200922_RNAsamples_inMilan.xlsx", 2)) %>% 
      mutate(sample_name_huge = str_extract(string = `HUGE ID`, pattern = "DBO-[0-9]*"),
             sample_name_new = `New ID`,
             patient = str_extract(string = `sample ID`, pattern = "TO [0-9]+|s[0-9]+"),
             donor = `Donor`,
             RPs.MPs = `type`,
             treatment = `medication`,
             RPs_MPs = `RPs.MPs`,
             disease = str_extract(string = `Disease`, pattern = "stable coronary artery disease")) %>%
      mutate(disease=replace(disease, disease=="stable coronary artery disease", "stable CAD"))
    sample_annotation_new <- sample_annotation_new[,c("donor","RPs_MPs","disease","sample_name_huge","sample_name_new","patient", "gender", "treatment", "age")]
    
    
    sample_annotation_additionals <- as.data.table(read.csv("../data/sample_annotation/Sampleinfo_RNAseqdata.csv", sep=";")) %>% 
      mutate(sample_name_huge = str_extract(string = HUGE.ID, pattern = "DBO-[0-9]*"),
             sample_name_new = new.name,
             patient = str_extract(string = sample.ID, pattern = "TO [0-9]+|s[0-9]+"),
             donor = Donor,
             RPs.MPs = type,
             treatment = medication,
             RPs_MPs = RPs.MPs,
             disease = str_extract(string = `Disease`, pattern = "stable coronary artery disease")) %>%
      mutate(disease=replace(disease, disease=="stable coronary artery disease", "stable CAD"))
    sample_annotation_additionals <- sample_annotation_additionals[,c("donor","RPs_MPs","disease","sample_name_huge","sample_name_new","patient", "gender", "treatment", "age")]
    
    
    # combine all sample_annotation fragments
    sample_annotation <- data.table()
    
    sample_annotation_additionals_new <- sample_annotation_additionals[sample_name_new %in% sample_annotation_new$sample_name_new,]
    if(all.equal(sample_annotation_new, sample_annotation_additionals_new)){
      sample_annotation <- copy(sample_annotation_additionals_new)
    }
    rm(sample_annotation_new, sample_annotation_additionals_new)
    
    sample_annotation_additionals_old <- sample_annotation_additionals[sample_name_new %in% sample_annotation_old$sample_name_new,]
    sample_annotation_old_additionals <- sample_annotation_old[sample_annotation_old$sample_name_new %in% sample_annotation_additionals$sample_name_new,]
    sample_annotation_old_additionals$age <- sample_annotation_additionals_old$age
    if(all.equal(sample_annotation_old_additionals, sample_annotation_additionals_old)){
      sample_annotation <- rbind(sample_annotation, sample_annotation_additionals_old)
    }
    rm(sample_annotation_additionals_old, sample_annotation_old_additionals)
    
    sample_annotation_old_missing <- sample_annotation_old[!sample_name_new %in% sample_annotation_additionals$sample_name_new,]
    sample_annotation_old_missing$age <- rep("NA", nrow(sample_annotation_old_missing))
    sample_annotation <- rbind(sample_annotation, sample_annotation_old_missing)
    rm(sample_annotation_old, sample_annotation_additionals, sample_annotation_old_missing)
    
    return(sample_annotation)
  }
  
  
    sample_annotation <- combine_sample_annotation()
    sample_annotation <- sample_annotation %>%
      mutate(disease=replace(disease, disease=="stable CAD", "CCS"),
             disease=replace(disease, disease=="MI", "ACS"))
    sample_annotation$sample_name_new <- gsub('-', '_', sample_annotation$sample_name_new)
    sample_annotation$sample_name_huge <- gsub('-', '_', sample_annotation$sample_name_huge)
    sample_annotation$patient <- gsub(' ', '_', sample_annotation$patient)
    sample_annotation$sample_name <- sample_annotation$sample_name_new
    sample_annotation <- sample_annotation[, sample_name_new := NULL]
    
    sample_annotation <- sample_annotation %>%
      separate(sample_name, c("sample_number", "type"), "_")
    sample_annotation <- as.data.table(sample_annotation)
    sample_annotation <- sample_annotation[order(sample_number, RPs_MPs),]
    sample_annotation$sample_name <- paste0(sample_annotation$sample_number, "_", sample_annotation$type)
    sample_annotation <- sample_annotation[, type := NULL]
    sample_annotation$sample_name <- factor(sample_annotation$sample_name, levels = unique(sample_annotation$sample_name))
    
  
    if(picked_disease == "all"){
      picked_disease = unique(sample_annotation$disease)
    }
  
    sample_annotation <- as.data.frame(sample_annotation[disease %in% picked_disease & !sample_name %in% remove_samples, ])
    rownames(sample_annotation) <- sample_annotation$sample_name
    return(sample_annotation)
}
