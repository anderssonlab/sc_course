# Sys.setenv(PATH = paste("/opt/software/gcc/11.2.0/bin", Sys.getenv("PATH"), sep = ":"))
# Sys.setenv(CXX = "/opt/software/gcc/11.2.0/bin/g++")
# 
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz")
# install.packages('Seurat')
# library(Seurat)
# 
# setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
# install.packages("Signac")
# library(Signac)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # BiocManager::install("AnnotationHub")
# install.packages("https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_2.3.4.tar.gz")
# BiocManager::install("ensembldb")
# BiocManager::install("biovizBase")
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# BiocManager::install("chromVAR")
# BiocManager::install("JASPAR2020")
# BiocManager::install("motifmatchr")
# install.packages('ggseqlogo')
# BiocManager::install("EnsDb.Hsapiens.v86")



# https://stuartlab.org/signac/articles/pbmc_vignette
# https://stuartlab.org/signac/articles/pbmc_multiomic

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)


counts_dir <- "input_data/filtered_feature_bc_matrix/"
frag_path <- "input_data/10k_PBMC_Multiome_nextgem_Chromium_X_atac_fragments.tsv.gz"
metadata_path <- "input_data/10k_PBMC_Multiome_nextgem_Chromium_X_per_barcode_metrics.csv"

counts <- Read10X(counts_dir)
metadata <- read.csv(
  file = metadata_path,
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = frag_path,
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

pbmc[['peaks']]

granges(pbmc)


peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]


# library(AnnotationHub)
# ah <- AnnotationHub()
# 
# # Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
# query(ah, "EnsDb.Hsapiens.v98")
# 
# ensdb_v98 <- ah[["AH75011"]]
# 
# # extract gene annotations from EnsDb
# annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
# 
# # change to UCSC style since the data was mapped to hg38
# seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# genome(annotations) <- "hg38"
# 
# 
# # add the gene information to the object
# Annotation(pbmc) <- annotations

library(EnsDb.Hsapiens.v86)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# Add the gene information to the object
Annotation(pbmc) <- annotation


# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc)

# add fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$atac_peak_region_fragments / pbmc$atac_fragments * 100

# add blacklist ratio
pbmc$blacklist_ratio <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 10000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 2 &
    TSS.enrichment > 4
)
pbmc


pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

DepthCor(pbmc)


pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()


gene.activities <- GeneActivity(pbmc)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)


DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("input_data/pbmc_10k_v3.rds")
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)


plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2


predicted_id_counts <- table(pbmc$predicted.id)

# Identify the predicted.id values that have more than 20 cells
major_predicted_ids <- names(predicted_id_counts[predicted_id_counts > 20])
pbmc <- pbmc[, pbmc$predicted.id %in% major_predicted_ids]

# change cell identities to the per-cell predicted labels
Idents(pbmc) <- pbmc$predicted.id


# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

# wilcox is the default option for test.use
da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'wilcox',
  min.pct = 0.1
)

head(da_peaks)


plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14+ Monocytes")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2


open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)

head(closest_genes_cd4naive)

head(closest_genes_cd14mono)

pbmc <- SortIdents(pbmc)

# find DA peaks overlapping gene of interest
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd4naive), LookupGeneCoords(pbmc, "CD4"))

CoveragePlot(
  object = pbmc,
  region = "CD4",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)

regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd14mono), LookupGeneCoords(pbmc, "LYZ"))

CoveragePlot(
  object = pbmc,
  region = "LYZ",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 5000
)


library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
pbmc <- AddMotifs(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)


top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])
DefaultAssay(pbmc) <- 'peaks'
enriched.motifs <- FindMotifs(
  object = pbmc,
  features = top.da.peak
)
MotifPlot(
  object = pbmc,
  motifs = head(rownames(enriched.motifs))
)


pbmc <- RunChromVAR(object = pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(pbmc) <- 'chromvar'


p1 <- DimPlot(pbmc, label = TRUE, pt.size = 0.1) + NoLegend()
p1
# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = pbmc,
  features = "MA0497.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2


differential.activity <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = pbmc,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)

saveRDS(pbmc,
        "output_data/pbmc.ATAC.rds")

