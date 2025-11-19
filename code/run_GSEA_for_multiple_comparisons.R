library(here)
library(limma)
library(DESeq2)
library(ggplot2)
library(stringr)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
source(here::here("code","p53_and_positional_genesets.R"))

#utility functions
get_sample_groupds <- function(samples, genotypes,comparison){
  #convert the "," separated list of genotypes to vectors and remove unneeded spaces
  group1_genotypes <- comparison$group1_genotypes
  group2_genotypes <- comparison$group2_genotypes
  group1_subtypes <- comparison$group1_subytypes
  group2_subtypes <- comparison$group2_subytypes
  group1_genotypes <- trimws(unlist(strsplit(group1_genotypes,",")))
  group2_genotypes <- trimws(unlist(strsplit(group2_genotypes,",")))
  group1 <- genotypes %in% group1_genotypes
  group2 <- genotypes %in% group2_genotypes
  #extract subtype groups, if they exist
  if(group1_subtypes != ''){
    group1_subtypes <- trimws(unlist(strsplit(group1_subtypes,",")))
    group2_subtypes <- trimws(unlist(strsplit(group2_subtypes,",")))
    group1 <- group1 & (subtypes %in% group1_subtypes)
    group2 <- group2 & (subtypes %in% group2_subtypes)
  }
  group1_samples <- samples[group1]
  group2_samples <- samples[group2]
  return(list(group1_samples = group1_samples,group2_samples = group2_samples))
}

run_DESeq2 <- function(group1_samples, group2_samples, count_matrix){
  set.seed(1234)
  count_matrix <- round(count_matrix)
  fillterd_count_matrix = count_matrix[,c(group1_samples, group2_samples)]
  
  #create DESeq2 comparison object and run
  cond_df <- data.frame(condition = c(rep('group1',length(group1_samples)),rep('group2',length(group2_samples))))
  rownames(cond_df) <- c(group1_samples, group2_samples)
  cond_df$condition <- factor(cond_df$condition, levels = c('group2','group1'))
  dds <- DESeqDataSetFromMatrix(countData = fillterd_count_matrix,
                                colData = cond_df,
                                design = ~ condition)
  dds <- DESeq(dds)
  return(data.frame(results(dds)))
}

run_limma <- function(group1_samples, group2_samples, log_intensity){
  set.seed(1234)
  
  fillterd_log_intensity = log_intensity[,c(group1_samples, group2_samples)]
  
  is_group_1 <- factor(colnames(fillterd_log_intensity) %in% group1_samples)
  limma.res <- lmFit(fillterd_log_intensity, model.matrix(~ is_group_1))
  limma.res <- eBayes(limma.res)
  limma.res <- topTable(limma.res, n=Inf)
  return(limma.res)
}

run_GSEA <- function(geneList, category, subcategory, geneset_name){
  if(subcategory == ''){
    subcategory = NULL
  }
  
  if(geneset_name == 'Positional'){
    m_t2g <- Positional
  }
  else if(geneset_name == 'p53'){
    m_t2g <- p53.genesets
  } else {
    m_t2g <- msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory) %>% 
      dplyr::select(gs_name, gene_symbol) 
  }
  set.seed(1234)
  return(GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1, seed = TRUE))
}

load_top10_hallmark_resutls <- function(dir){
  base.gsea.res <- readRDS(here::here(dir,'GSEA','Hallmarks_Del17p_vs_no_Del17p_all.rds'))
  top10 <- base.gsea.res@result$ID[1:10]
  top10_empty <- base.gsea.res@result[1:10,]
  top10_empty$setSize <- 0
  top10_empty$enrichmentScore <- 0
  top10_empty$NES <- 0
  top10_empty$pvalue <- 1
  top10_empty$p.adjust <- 1
  top10_empty$qvalue <- 1
}

#load data

#For METABRIC
# is_gene_counts <- FALSE
# 
# save_to <- "METABRIC"
# 
# expression_data <- read.csv(here::here("data","Gene_expression_METABRIC.csv"), header = F)
# 
# comparisons <- read.csv(here::here("data","Comparisons_METABRIC.csv"))

#For TCGA
is_gene_counts <- TRUE

save_to <- "TCGA"

expression_data <- read.csv(here::here("data","Gene_expression_TCGA.csv"), header = F)

comparisons <- read.csv(here::here("data","Comparisons_TCGA.csv"))

samples <- unname(unlist(as.vector(expression_data[1,-c(1,2)])))

genotypes <- unname(unlist(as.vector(expression_data[2,-c(1,2)])))

subtypes <- unname(unlist(as.vector(expression_data[3,-c(1,2)])))

genes <-  unname(unlist(as.vector(expression_data[-c(1,2,3),1])))

expression_data <- expression_data[-c(1,2,3),-c(1,2)]
colnames(expression_data) <- samples

expression_data <- expression_data[!duplicated(genes),]

expression_data <- apply(expression_data, 2, function(x) as.numeric(x))
expression_data <- as.data.frame(expression_data)

rownames(expression_data) <- genes[!duplicated(genes)]

genesets <- data.frame(name = c('Hallmarks','KEGG',"Reactome","Oncogenc_signatures","GOBP",'Positional','p53'),
                       category = c("H","C2","C2","C6","C5","",""),
                       subcategory = c("","CP:KEGG_LEGACY","CP:REACTOME","","GO:BP","",""))

#loop over all combinations
for(d in c('DEG','GSEA','plots')){
  if(!dir.exists(here::here(save_to,d))){
    dir.create(here::here(save_to,d), recursive = T)
  }
}

for(row in seq(1:nrow(comparisons))){
  tryCatch({
    print(paste0('Working in comparison number ',row))
    DEG_file_path <- here::here(save_to,'DEG',paste0(comparisons[row,'Comparison_name'],'.csv'))
    if(file.exists(DEG_file_path)){
      DEG <- read.csv(DEG_file_path, row.names = 1)
    } else {
      sample_groups <- get_sample_groupds(samples, genotypes, comparisons[row,])
      if(is_gene_counts){
        DEG <- run_DESeq2(sample_groups$group1_samples,sample_groups$group2_samples, expression_data)
      } else{
        DEG <- run_limma(sample_groups$group1_samples,sample_groups$group2_samples, expression_data)
      }
      write.csv(DEG, DEG_file_path)
    }
    if(is_gene_counts){
      DEG <- filter(DEG, !is.na(log2FoldChange)) %>% arrange(desc(log2FoldChange))
      geneList <- DEG$log2FoldChange
      names(geneList) <- rownames(DEG)
    } else {
      DEG <- DEG %>% arrange(desc(logFC))
      geneList <- DEG$logFC
      names(geneList) <- rownames(DEG)
    }
    
    for(j in seq(1:nrow(genesets))){
      category <- genesets[j,'category']
      subcategory <- genesets[j,'subcategory']
      geneset_name <- genesets[j,'name']
      
      gsea_file_path <- here::here(save_to,'GSEA',paste0(geneset_name,'_',comparisons[row,'Comparison_name'],'.rds'))
      if(file.exists(gsea_file_path)){
        gsea_res <- readRDS(gsea_file_path)
      } else {
        gsea_res <- run_GSEA(geneList, category, subcategory, geneset_name)
        saveRDS(gsea_res,gsea_file_path)
      }
      if((geneset_name == 'Hallmarks') & (comparisons$Comparison_name[j] != 'Del17p_vs_no_Del17p_all')){
        top10 <- load_top10_hallmark_resutls(save_to)
        gsea_res_fillterd = gsea_res %>% dplyr::filter(ID %in% top10)
        gsea_res_fillterd@result <- rbind(gsea_res_fillterd@result, dplyr::filter(top10_empty,! ID %in% gsea_res_fillterd@result$ID))
        gsea_res_fillterd@result$ID = factor(gsea_res_fillterd@result$ID, levels = top10)
      } else{
        gsea_res_fillterd <- dplyr::filter(gsea_res,  qvalue < 0.25)
      }
      if(nrow(gsea_res_fillterd) > 0){
        pdf_file_path <- here::here(save_to,'plots',paste0(geneset_name,'_',comparisons[row,'Comparison_name'],'.pdf'))
        if(! file.exists(pdf_file_path)){
          pdf(file=pdf_file_path,width=14,height=14)
          print(dotplot(gsea_res_fillterd, showCategory=10, x='NES', orderBy='ID',decreasing=F) + 
            scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + scale_fill_continuous(limit =c (0,1), low = "#e06663",high = "#327eba") +
            scale_x_continuous(limits = c(-3.5,3.5)))
          dev.off()
        }
      }
    }
  },error = function(cond){
    print(paste0('Error while working on row number ', row))
    print(cond)
  })
}

