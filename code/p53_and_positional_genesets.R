library(dplyr)
library(msigdbr)
#p53 Pathways from two different papers


#https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.1057302/full"
TP53_induced <- c('DDB2', 'FAS', 'GADD45A', 'RPS27L', 'EDA2R', 'ACAD11', 'TRIM22', 'SPATA18', 'AEN', 'FDXR', 'MDM2', 'CDKN1A', 'PTCHD4', 'ZMAT3', 'PANK1', 'ALDH4A1', 'ESR1', 'RGCC', 'GADD45B', 'PHLDA3')
TP53_repressed <- c('CCNB1', 'PLK1', 'EED', 'CDK1', 'EZH2', 'CCNB2', 'E2F3', 'MYBL2', 'FOXM1', 'E2F2')


"https://www.biorxiv.org/content/10.1101/2022.07.28.501874v2.full#F7"
p53tru_DR <- c('CPE', 'BTG2', 'SCD', 'SESN1', 'SLC12A4', 'ZNF219', 'ATF3', 'ME3', 'RRAD', 'ME1', 'MAST4', 'TNFRSF10D', 'ABCB1', 'FAM13C', 'INKA2', 'SPATA18', 'EDA2R', 'TSPAN11', 'ACER2', 'VWCE', 'TLR3')
p53tru_UR <- c('PCNA', 'CDC25B', 'BAX', 'LAPTM5', 'TRAF4', 'POLD1', 'CKS2', 'SAC3D1', 'ENC1', 'RAP2B', 'GDF15', 'CCNB1', 'PRC1', 'ECT2', 'PLK1', 'PMAIP1', 'CDK1', 'ANLN', 'CCNB2', 'CDC20', 'RAD51', 'BIRC5', 'NEK2', 'GRHL3', 'CDC25C', 'ABCA12')

p53.genesets <- data.frame(
  gs_name = c(rep('TP53_induced',length(TP53_induced)), rep('TP53_repressed',length(TP53_repressed)),
              rep('p53tru_DR',length(p53tru_DR)), rep('p53tru_UR',length(p53tru_UR))),
  gene_symbol = c(TP53_induced, TP53_repressed, p53tru_DR, p53tru_UR)
)

#Positional geneset from Ron for arm level copy number
Positional <- msigdbr(species = "Homo sapiens", category = "C1") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::filter(gs_name != 'MT')

gsub('q.*','q',gsub('p.*','p',Positional$gs_name)) %>% unique() %>% sort()

Positional$gs_name <- gsub('q.*','q',gsub('p.*','p',Positional$gs_name))
