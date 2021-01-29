for( tot in 1:100)
{
  library("SAFEclustering")
  library('SC3')
  library(SingleCellExperiment)
  library('tidyverse')
  a<-read_rds('biase(changed).rds')
  print(a)
  data<-assay(a, "logcounts")
  counts(a) <-assay(a, "normcounts")
  library(SC3)
  deng <- sc3_prepare(a)
  gene_filter <- rowData(deng)$sc3_gene_filter
  matrix_filtered=exp(data[which(gene_filter),1:49])
  cluster.results <- individual_clustering(inputTags = matrix_filtered, mt_filter = FALSE, nGene_filter = FALSE, SC3 = TRUE, gene_filter = FALSE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nPC.seurat = NULL, resolution = 0.9, tSNE = TRUE, dimensions = 3, perplexity = 30, SEED = 123)
  cluster.ensemble <- SAFE(cluster_results = cluster.results, program.dir = "gpmetis_and_shmetis_for_Linux", k_min=3,k_max=3,
                           MCLA = TRUE, CSPA = TRUE, HGPA = TRUE, SEED = 123)
  library(cidr)
  cat(cluster.ensemble$optimal_clustering,file="SAFE_Biase.txt",append=TRUE)
  cat("\n",file="SAFE_Biase.txt",append=TRUE)
}