# Code for the "Chromosome-arm 17p loss renders breast cancer cells vulnerable to AURKB inhibition" paper. #
**To run the code, first downlad the TCGA and METABRIC expression files from https://doi.org/10.6084/m9.figshare.30655643.v1 and place them under the data directory** 

## Survival analysis ##

To recreate the TCGA/Metabric survival plots run the Survival_analysis_TCGA_120_month_limited.R or Survival_analysis_METABRIC_120_month_limited.R.
The results are saved to TCGA/survival or Metabric/survival directory.

## GSEA analysis ##
To recreate the GSEA resutlts run run_GSEA_for_multiple_comparisons.R. The results are saved to TCGA/plots or Metabric/plots directory.
the code has 4 settings the user needs to define (default values used for TCGA, values for METABRIC are commented out):
1. is_gene_counts - is the input count data or microarray data.
2. save_to - the name of the output directory.
3. expression_data - path to CSV with expression data. The first row contains sample IDs, the second sample genotypes, and the third subtype. The first column contains the gene Hugo_Symbol, and the second Entrez_Gene_Id (not used and can be ignored).
4. comparisons—path to CSV defining the comparisons. Each comparison has a unique name and a comma-separated list of genotypes and subtypes for each group.

### Comments ###
DESeq2 is used for RNA-Seq (TCGA) count data and limma for microarray data (METABRIC).

The geneset collections to use in the GSEA are hard-coded and can be changed in the code.
### Additional files ###
1. Gene_expression_TCGA.csv, Comparisons_TCGA.csv - input for TCGA BRCA comparisons, count data.
2. Gene_expression_METABRIC, Comparisons_METABRIC.csv - input for METABRIC comparisons, microarray data.
3. p53_and_positional_genesets.R - contain definitions for 2 aditinal gene set collections: Positional - arm level gene sets. p53 - 4 literature-based TP53-related signature, see references in the file.
