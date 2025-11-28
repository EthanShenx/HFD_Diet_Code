library(Seurat)
library(CellChat)
library(dplyr)
library(future)
library(biomaRt)

options(future.globals.maxSize = 200 * 1024^3)
plan(multicore, workers = 8)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/human_atlas/ST-patient-79_seurat.rds")

All <- NormalizeData(All)

colnames(All@meta.data)

DefaultAssay(All) <- "RNA"

All <- NormalizeData(All)

data.input <- Seurat::GetAssayData(All, slot = "data")
dim(data.input)

meta <- All@meta.data
meta$samples <- meta$sample_uuid

meta <- meta %>%
  dplyr::select(samples, cell_type)

colnames(meta) <- c("samples", "labels")

meta <- meta[colnames(data.input), , drop = FALSE]

coords <- All@meta.data[, c("array_row", "array_col")]
coords <- as.data.frame(coords)
colnames(coords) <- c("y", "x")  # names don't really matter, but this is clear
coords <- coords[colnames(data.input), , drop = FALSE]
rownames(coords) <- colnames(data.input)

head(coords)

spatial.factors <- data.frame(
  ratio = 1,
  tol   = 0.5
)
spatial.factors

d_spatial <- CellChat::computeCellDistance(
  coordinates = coords,
  ratio       = spatial.factors$ratio,
  tol         = spatial.factors$tol
)

summary(as.numeric(d_spatial[d_spatial > 0]))

ens_id <- rownames(data.input)
ens_id <- gsub("\\..*$", "", ens_id)
rownames(data.input) <- ens_id
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ens_ids <- rownames(data.input)
annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ens_ids,
  mart       = mart
)
annot <- annot %>%
  filter(hgnc_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)
symbol.vec <- annot$hgnc_symbol
names(symbol.vec) <- annot$ensembl_gene_id
data.mapped <- data.input[names(symbol.vec), ]
rownames(data.mapped) <- symbol.vec
expr.mat <- as.matrix(data.mapped)
expr.by.symbol <- rowsum(expr.mat, group = rownames(expr.mat))
data.input.symbol <- expr.by.symbol

cellchat <- CellChat::createCellChat(
  object          = data.input.symbol,
  meta            = meta,
  group.by        = "labels",
  datatype        = "spatial",
  coordinates     = coords,
  spatial.factors = spatial.factors
)

coords <- data.matrix(cellchat@images$coordinates)
storage.mode(coords) <- "double"
cellchat@images$coordinates <- coords
class(cellchat@images$coordinates)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- CellChat::subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(
  cellchat,
  variable.both = FALSE
)
cellchat <- CellChat::computeCommunProb(
  cellchat,
  type             = "truncatedMean",
  trim             = 0.1,
  distance.use     = TRUE,
  interaction.range = 80,
  scale.distance    = 1,
  contact.dependent = TRUE,
  contact.range     = 3
)