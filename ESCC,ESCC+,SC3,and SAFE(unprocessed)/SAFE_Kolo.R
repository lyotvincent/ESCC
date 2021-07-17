for( tot in 1:90)
{
  library("SAFEclustering")
  library('SC3')
  library(SingleCellExperiment)
  library('tidyverse')
  a<-read_rds('~/Desktop/bicluster-evolotionary/Gold (another copy)/Final/Result (5th copy)/kolodziejczyk.rds')
  print(a)
  data<-assay(a, "counts")
  library(SingleCellExperiment)
  library(SC3)
  sce <- SingleCellExperiment(assays = list(counts = counts(a),logcounts= log2(counts(a)+1)))
  rowData(sce)$feature_symbol=1:dim(a)[1]
  deng <- sc3_prepare(sce)
  gene_filter <- rowData(deng)$sc3_gene_filter
  matrix_filtered=data[which(gene_filter),1:704]
  cluster.results <- individual_clustering(inputTags = matrix_filtered, mt_filter = FALSE, nGene_filter = FALSE, SC3 = TRUE, gene_filter = FALSE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nPC.seurat = NULL, resolution = 0.9, tSNE = TRUE, dimensions = 3, perplexity = 30, SEED = 123)
  cluster.ensemble <- SAFE(cluster_results = cluster.results, program.dir = "~/Downloads/SAFEclustering-master/gpmetis_and_shmetis_for_Linux", k_min=3,k_max=3,
                           MCLA = TRUE, CSPA = TRUE, HGPA = TRUE, SEED = 123)
  library(cidr)
  cat(cluster.ensemble$optimal_clustering,file="~/Desktop/bicluster-evolotionary/Gold (another copy)/Final/Result (5th copy)/SAFE_Kolo.txt",append=TRUE)
  cat("\n",file="~/Desktop/bicluster-evolotionary/Gold (another copy)/Final/Result (5th copy)/SAFE_Kolo.txt",append=TRUE)
}