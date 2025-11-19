# Load necessary libraries
library(survival)
library(ggplot2)
library(survminer)
library(ggpmisc)
library(dplyr)
library(tidyverse)

# Read the CSV file
METABRIC <- read.csv(here::here("data","Survival_METABRIC.csv"))
annotations = list("OS" = c("Overall survival"),"RFS" = c("Relapse free survival"))

# remove #N/A rows
METABRIC <- METABRIC [METABRIC$Genotype!='#N/A',]

for (annotation in names(annotations)) {
  annotation_df <- METABRIC[, c("PATIENT_ID", "Genotype", "Molecular_subtype", paste0(annotation, "_MONTHS"), paste0(annotation, "_STATUS"))]
  # Rename columns
  colnames(annotation_df) <- c("Sample_Name", "Genotype_Group", "Molecular_Subtype", "Survival_Months", "Status")
  print(head(annotation_df))
  annotation_df <- drop_na(annotation_df)
  
  annotation_df$Survival_Months <- as.numeric(annotation_df$Survival_Months)
  annotation_df <- drop_na(annotation_df)
  
  annotation_df <- annotation_df %>% mutate(Status = ifelse(Survival_Months <=120, Status, '0'))
  annotation_df <- annotation_df %>% mutate(Survival_Months = ifelse(Survival_Months <=120, Survival_Months, 120))
  
  annotation_df$Status <-sapply(annotation_df$Status, function(x) {
    val <- strsplit(x, ":")[[1]][1]
    val=="1"
    })

  
  
  subtypes <- c(unique(annotation_df$Molecular_Subtype), 'all')
  subtypes <- subtypes[subtypes!='NC']
  
  # Define the groups for comparison
  group1 <- c("Del17p/p53_WT", "Del17p/p53_mut")
  group2 <- c("WT_17p/p53_mut", "WT_17p/p53_WT")
  
  for (subtype in subtypes){
    print(subtype)
    if (subtype=='all'){
      subset_df <- annotation_df
    } else {
      subset_df <- annotation_df[annotation_df$Molecular_Subtype == subtype, ]
    }
    # Create a survival object for the subset
    surv_obj <- Surv(time = subset_df$Survival_Months, event = subset_df$Status)
    
    # Create a new column indicating the groups
    subset_df$Group <- ifelse(subset_df$Genotype_Group %in% group1, "Del17p", "Non-del17p")
    subset_df$Group <- factor(subset_df$Group, levels = c("Del17p", "Non-del17p"))
    
    # Fit survival curves
    fit <- survfit(surv_obj ~ Group, data = subset_df)
    
    # pairwise survdiff
    res <- pairwise_survdiff(Surv(Survival_Months, Status) ~ Group, data = subset_df)
    table <- res$p.value %>% as.data.frame() %>% mutate(across(where(is.numeric), ~ round(., digits = 5)))
    
    # Plot survival curves
    g <- ggsurvplot(fit, data = subset_df, risk.table = FALSE, pval = TRUE, 
                    pval.coord = c(0, .3), break.x.by=10,
                    palette = c("blue", "grey25"), legend.title = "")
    
    # Customize the plot
    plot_title <- paste(annotations[annotation], "METABRIC. \n", subtype, ifelse(subtype=='all',  "molecular subtypes", "molecular subtype"))
    g <- g + labs(title = plot_title,
                  x = "Time (months)",
                  y = "Survival probability",
                  caption = "p-value from Log-rank test") 
    g$plot <- g$plot + 
    theme(plot.title = element_text(hjust = 0.5, size = 20))
    
    
    # Print the plot
    # print(g)
    if (!file.exists(here::here("METABRIC", "survival"))) {
      # Create the directory
      dir.create(here::here("METABRIC", "survival"), recursive = T)
    }
    if (!file.exists(here::here("METABRIC", "survival", annotation))) {
      # Create the directory
      dir.create(here::here("METABRIC", "survival", annotation))
    }
    if (!file.exists(here::here("METABRIC", "survival", annotation, "Del_vs_non_del"))) {
      # Create the directory
      dir.create(here::here("METABRIC", "survival", annotation, "Del_vs_non_del"))
    }
    pdf(here::here("METABRIC", "survival", annotation, "Del_vs_non_del", paste0("Del_vs_Non_del_17p_", subtype, "_120_month.pdf")))
    print(g, newpage = FALSE)
    dev.off()
  }
  
  ###Genotypes
  
  # Define the groups for comparison
  group1 <- "Del17p/p53_WT"
  group2 <- "Del17p/p53_mut"
  group3 <- "WT_17p/p53_mut"
  group4 <- "WT_17p/p53_WT"
  
  for (subtype in subtypes){
    if (subtype=='all'){
      subset_df <- annotation_df
    } else {
      subset_df <- annotation_df[annotation_df$Molecular_Subtype == subtype, ]
    }
    # Create a survival object for the subset
    surv_obj <- Surv(time = subset_df$Survival_Months, event = subset_df$Status)
    
    # Create a new column indicating the groups
    subset_df$Group <- ifelse(subset_df$Genotype_Group == group1, "Del17p/p53_WT",
                              ifelse(subset_df$Genotype_Group == group2, "Del17p/p53_mut",
                                       ifelse(subset_df$Genotype_Group == group3, "WT_17p/p53_mut",
                                            ifelse(subset_df$Genotype_Group == group4, "WT_17p/p53_WT","other"))))
    
    # Convert Group to a factor
    subset_df$Group <- factor(subset_df$Group, levels = c("Del17p/p53_WT", "Del17p/p53_mut", "WT_17p/p53_mut", "WT_17p/p53_WT"))
    
    # Fit survival curves
    fit <- survfit(surv_obj ~ Group, data = subset_df)
    
    # pairwise survdiff
    res <- pairwise_survdiff(Surv(Survival_Months, Status) ~ Group, data = subset_df)
    table <- res$p.value %>% as.data.frame() %>% mutate(across(where(is.numeric), ~ round(., digits = 5)))
    
    custom_palette <- c("blue", "purple", "red", "gray")
    
    # Plot survival curves
    g <- ggsurvplot(fit, data = subset_df, risk.table = FALSE, pval = TRUE, 
                    palette =custom_palette,break.x.by=10,
                    pval.coord = c(0, .3), 
                    legend.title = "")
    
    # Customize the plot
    plot_title <- paste(annotations[annotation], "METABRIC. \n", subtype, ifelse(subtype=='all',  "molecular subtypes", "molecular subtype"))
    g <- g + labs(title = plot_title,
                  x = "Time (months)",
                  y = "Survival probability",
                  caption = "p-value from Log-rank test")
    g$plot <- g$plot + annotate(geom = "table", x = 0, y = 0, label = list(as.data.frame(table)), table.rownames = TRUE)+ 
      theme(plot.title = element_text(hjust = 0.5, size = 20))
    
    
    # Print the plot
    # print(g)
    if (!file.exists(here::here("METABRIC", "survival", annotation, "Genotypes"))) {
      # Create the directory
      dir.create(here::here("METABRIC", "survival", annotation, "Genotypes"))
    }
    pdf(here::here("METABRIC", "survival", annotation, "Genotypes", paste0("4_Genotypes_", subtype, "_120_month.pdf")))
    print(g, newpage = FALSE)
    dev.off()
  }
  
  
  # Define the groups for comparison
  group1 <- "Del17p/p53_WT"
  group2 <- "Del17p/p53_mut"
  group3 <- "WT_17p/p53_WT"
  subset_annotation_df <- annotation_df[annotation_df$Genotype_Group %in% c(group1, group2, group3), ]
  
  for (subtype in subtypes){
    if (subtype=='all'){
      subset_df <- subset_annotation_df
    } else {
      subset_df <- subset_annotation_df[subset_annotation_df$Molecular_Subtype == subtype, ]
    }
    # Create a survival object for the subset
    surv_obj <- Surv(time = subset_df$Survival_Months, event = subset_df$Status)
    
    # Create a new column indicating the groups
    subset_df$Group <- ifelse(subset_df$Genotype_Group == group1, "Del17p/p53_WT",
                              ifelse(subset_df$Genotype_Group == group2, "Del17p/p53_mut",
                                            ifelse(subset_df$Genotype_Group == group3, "WT_17p/p53_WT","other")))
    
    # Convert Group to a factor
    subset_df$Group <- factor(subset_df$Group, levels = c("Del17p/p53_WT", "Del17p/p53_mut", "WT_17p/p53_WT"))
    
    # Fit survival curves
    fit <- survfit(surv_obj ~ Group, data = subset_df)
    
    # pairwise survdiff
    res <- pairwise_survdiff(Surv(Survival_Months, Status) ~ Group, data = subset_df)
    table <- res$p.value %>% as.data.frame() %>% mutate(across(where(is.numeric), ~ round(., digits = 5)))
    
    custom_palette <- c("blue", "purple", "gray")
    
    # Plot survival curves
    g <- ggsurvplot(fit, data = subset_df, risk.table = FALSE, pval = TRUE, 
                    palette =custom_palette,break.x.by=10,
                    pval.coord = c(0, .3), 
                    legend.title = "")
    
    # Customize the plot
    plot_title <- paste(annotations[annotation], "METABRIC. \n", subtype, ifelse(subtype=='all',  "molecular subtypes", "molecular subtype"))
    g <- g + labs(title = plot_title,
                  x = "Time (months)",
                  y = "Survival probability",
                  caption = "p-value from Log-rank test")
    g$plot <- g$plot + annotate(geom = "table", x = 0, y = 0, label = list(as.data.frame(table)), table.rownames = TRUE)+ 
      theme(plot.title = element_text(hjust = 0.5, size = 20))
    
    # Print the plot
    # print(g)
    pdf(here::here("METABRIC", "survival", annotation, "Genotypes", paste0("3_Genotypes_", subtype, "_120_month.pdf")))
    print(g, newpage = FALSE)
    dev.off()
  }
  
}

