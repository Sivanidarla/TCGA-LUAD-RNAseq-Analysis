# TCGA-LUAD RNA-seq analysis pipeline
# Tumor-normal differential expression, TOP2A biomarker evaluation,
# and PP vs TRU subtype-level pathway interpretation.

# ------------------------------
# 0. Libraries
# ------------------------------

library(TCGAbiolinks)
library(DESeq2)
library(limma)
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pROC)
library(clusterProfiler)
library(ggplot2)

# ------------------------------
# 1. Download and prepare data
# ------------------------------

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, directory = "data/raw")

luad_data <- GDCprepare(query, directory = "data/raw")
saveRDS(luad_data, "data/luad_raw_se.rds")

counts <- assay(luad_data, "unstranded")
metadata <- as.data.frame(colData(luad_data))

metadata$condition <- factor(metadata$definition)
metadata$batch <- as.factor(substr(metadata$barcode, 6, 7)) # TSS proxy

# ------------------------------
# 2. Tumor-normal differential expression
# ------------------------------

keep <- rowSums(counts >= 10) >= 20
counts_filtered <- counts[keep, ]

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ condition
)

dds$condition <- relevel(dds$condition, ref = "Solid Tissue Normal")
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Primary solid Tumor", "Solid Tissue Normal")
)

# Batch-associated variation was adjusted only for visualization.
vsd <- vst(dds, blind = FALSE)
mat_corrected <- limma::removeBatchEffect(assay(vsd), batch = vsd$batch)
assay(vsd) <- mat_corrected

write.csv(as.data.frame(res), "results/Tumor_vs_Normal_DEGs.csv")

saveRDS(vsd, "data/vsd_corrected.rds")
saveRDS(dds, "data/dds_tumor_normal.rds")
saveRDS(res, "data/res_tumor_normal.rds")

# ------------------------------
# 3. Gene symbol mapping and volcano plot
# ------------------------------

ensembl_ids <- sub("\\..*$", "", rownames(res))

res$symbol <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

write.csv(res, "results/Tumor_vs_Normal_DEGs_with_Symbols.csv", row.names = TRUE)

png("plots/volcano_main_symbols.png", width = 1000, height = 1000, res = 120)

EnhancedVolcano(
  res,
  lab = res$symbol,
  x = "log2FoldChange",
  y = "padj",
  title = "Tumor vs Normal LUAD",
  subtitle = "Top genes labelled by gene symbol",
  pCutoff = 10e-50,
  FCcutoff = 2,
  pointSize = 2.0,
  labSize = 4.0,
  legendPosition = "right"
)

dev.off()

# ------------------------------
# 4. PP vs TRU subtype analysis
# ------------------------------

# Keep only PP and TRU samples.
sub_metadata <- metadata[
  make.names(metadata$paper_expression_subtype) %in% c("prox..prolif.", "TRU"),
]

sub_counts <- assay(vsd)[, rownames(sub_metadata)]

# Clean subtype names for model design.
sub_metadata$paper_expression_subtype <- make.names(
  sub_metadata$paper_expression_subtype
)

design_sub <- model.matrix(~ 0 + paper_expression_subtype, data = sub_metadata)

fit <- lmFit(sub_counts, design_sub)

# Explicit contrast: PP minus TRU.
cont_matrix <- makeContrasts(
  PP_vs_TRU = paper_expression_subtypeprox..prolif. - paper_expression_subtypeTRU,
  levels = design_sub
)

fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

res_sub <- topTable(fit2, coef = "PP_vs_TRU", number = Inf)

write.csv(res_sub, "results/Subtype_PP_vs_TRU_DEGs.csv", row.names = TRUE)

# ------------------------------
# 5. TOP2A biomarker evaluation
# ------------------------------

top2a_id <- rownames(res)[which(res$symbol == "TOP2A")][1]

gene_expr <- assay(vsd)[top2a_id, ]
labels <- ifelse(vsd$definition == "Solid Tissue Normal", 0, 1)

roc_obj <- roc(labels, gene_expr)

png("plots/ROC_TOP2A.png", width = 600, height = 600)

plot(
  roc_obj,
  main = paste("Biomarker Performance: TOP2A\nAUC =", round(auc(roc_obj), 3)),
  col = "#2c7fb8",
  lwd = 4
)

abline(a = 0, b = 1, lty = 2, col = "gray")

dev.off()

# ------------------------------
# 6. GO enrichment for PP-upregulated genes
# ------------------------------

rownames(res_sub) <- gsub("\\..*$", "", rownames(res_sub))
universe_ids <- rownames(res_sub)

# PP-up genes only.
sig_res <- res_sub[res_sub$adj.P.Val < 0.01 & res_sub$logFC > 1, ]
sig_ids <- rownames(sig_res)

ego <- enrichGO(
  gene = sig_ids,
  universe = universe_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

png("plots/pathway_enrichment_polished.png", width = 1000, height = 800, res = 130)

dotplot(ego, showCategory = 15) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  ggtitle("PP vs TRU LUAD Subtype: Top Enriched Biological Processes")

dev.off()

# ------------------------------
# 7. Curated target summary
# ------------------------------

targets <- data.frame(
  Gene = c("TOP2A", "CDK1", "AURKA", "CCNB1", "BUB1"),
  Mechanism = c(
    "DNA Topology",
    "Cell Cycle Checkpoint",
    "Mitotic Regulation",
    "Cyclin-dependent Signalling",
    "Spindle Assembly"
  ),
  Clinical_Status = c(
    "FDA Approved Target",
    "Clinical Trials",
    "Clinical Trials",
    "High-Value Biomarker",
    "Research Target"
  ),
  Pharma_Relevance = c(
    "Inhibition induces apoptosis",
    "Target for CDK inhibitors",
    "Inhibition halts mitosis",
    "Key proliferation marker",
    "Potential for immunotherapy synergy"
  )
)

write.csv(targets, "results/pharma_druggable_targets.csv", row.names = FALSE)

# End of script