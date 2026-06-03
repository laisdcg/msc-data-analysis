install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GEOquery")
library(GEOquery)
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("WGCNA")
library(WGCNA)
BiocManager::install("biomaRt")
library(biomaRt)
install.packages("umap")
library(umap)
BiocManager::install("edgeR")
library(edgeR)
install.packages("ggplot2")
library(ggplot2)
install.packages("dplyr")
library(dplyr)
BiocManager::install("gplots")
library(gplots)
install.packages("tidyr")
library(tidyr)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("enrichplot")
library(enrichplot)
install.packages("gprofiler2")
library(gprofiler2)
BiocManager::install("topGO")
library(topGO)
install.packages("factoextra")
library(factoextra)
install.packages("pheatmap")
library(pheatmap)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
BiocManager::install("ReactomePA")
library(ReactomePA)
BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)
BiocManager::install("GENIE3")
library(GENIE3)
install.packages("igraph")
library(igraph)
install.packages("DT")
library(DT)
BiocManager::install("KEGGREST")
library(KEGGREST)
BiocManager::install("pathview")
library(pathview)
BiocManager::install("GOstats")
library(GOstats)
BiocManager::install("TFBSTools")
library(TFBSTools)
BiocManager::install("GenomicRanges")
library(GenomicRanges)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene") # ou hg19 se você precisar
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
BiocManager::install("motifmatchr")
library(motifmatchr)
BiocManager::install("Biostrings")
library(Biostrings)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
BiocManager::install("JASPAR2020")
library(JASPAR2020)

BiocManager::install("GOSemSim")
library(GOSemSim)
BiocManager::install("rnextflow")
library(rnextflow)

packageVersion("GOSemSim")
packageVersion("enrichplot")


# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE221786", "file=GSE221786_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
gset<- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)

# Load the dataset
gse <- getGEO("GSE221786", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse <- gse[[1]]

# Extract phenotype data
pheno_data <- pData(gse)
pheno_data <- pheno_data[, c("title", "geo_accession", "disease state:ch1")]
print(pheno_data)

genes_of_interest <- c(
  "ABCB11", "ABCD2", "ACE", "ACKR2", "ACP5", "ACP7", "ADAMTS9", "ADCY7", "ADGRL2", "ADIPOQ", "AHR", "AIM2", "AKAP13", "ALB", "ANKRD55", "ANXA1",
  "APOL1", "APOL6", "AQP1", "AREG", "ARHGEF3", "ATG16L1", "ATG5", "ATOX1", "ATXN2L", "BEND2", "BGLAP", "BMP2", "BMP7", "BNIP3L", "BRAF", "BTG1",
  "C15ORF48", "C7ORF57", "CALD1", "CAPS2", "CARD9", "CASP1", "CASP10", "CAST", "CCL2", "CCL20", "CCL25", "CCL3", "CCL5", "CCND1", "CCND3", "CCR2",
  "CCR6", "CCR7", "CCRL2", "CD14", "CD160", "CD19", "CD274", "CD28", "CD36", "CD38", "CD3E", "CD4", "CD40", "CD40LG", "CD52", "CD63", "CD68",
  "CD69", "CD80", "CD83", "CD86", "CD8A", "CD8B", "CDC42BPB", "CEBPA", "CEBPG", "CENPK", "CLEC2B", "CLEC4D", "CLIC3", "CMAHP", "CMTM2",
  "COL1A1", "COX2", "CRB1", "CREM", "CRP", "CSF2", "CSMD1", "CSN3", "CSNK1A1", "CTLA4", "CTSB", "CTSK", "CUX1", "CX3CR1", "CXCL10", "CXCL13",
  "CXCL16", "CXCL2", "CXCL8", "CXCR2", "CYCS", "CYP1A1", "CYP4F22", "DDIT3", "DDX60", "DKK1", "DLAT", "DNAJA2", "DNAJB6", "DUSP4", "DYSF",
  "EDNRA", "EFCAB13", "EFCAB7", "EGF", "EGFR", "EGR3", "EIF5B", "ENO1", "EOMES", "ERAP1", "ERAP2", "ERN1", "ERP44", "EZH2", "FAXDC2", "FCGR1A",
  "FCGR2A", "FCGR3A", "FGB", "FNBP1", "FOS", "FOSL1", "FOSL2", "FOXP3", "FRZB", "FUT2", "GAPDH", "GATA3", "GBP1", "GBP3", "GBP5", "GEM", "GINS1",
  "GJB2", "GJB6", "GLUL", "GOLIM4", "GPR35", "GPT", "GZMA", "GZMB", "GZMK", "HCAR3", "HERC6", "HHAT", "HIF1A", "HK2", "HLA-A", "HLA-B", "HLA-C",
  "HLA-DQB1", "HLA-DRB1", "HSPA5", "HSPA6", "HYAL4", "ICAM1", "ICOS", "IFI16", "IFI6", "IFIH1", "IFIT3", "IFNA1", "IFNG", "IFNGR1", "IFNLR1",
  "IGF2", "IGF2-AS", "IL10", "IL12B", "IL13", "IL15", "IL17A", "IL17F", "IL17RA", "IL18", "IL1A", "IL1B", "IL1F10", "IL1R1", "IL1RN", "IL2",
  "IL21", "IL21-AS1", "IL22", "IL23A", "IL23R", "IL2RA", "IL33", "IL36RN", "IL37", "IL4", "IL5", "IL6", "IL6R", "IL7", "INS", "INS-IGF2", "IRAK1",
  "IRS1", "ITGA2B", "ITGAL", "ITGAM", "ITGAX", "JAK2", "JAK3", "JDP2", "JRKL", "JUN", "JUNB", "JUND", "KANK4", "KDM5B", "KIR2DS1", "KIR3DL1",
  "KIR3DL2", "KLRB1", "LAMP1", "LAMP2", "LEP", "LILRA5", "LILRB2", "LINC01185", "LINC01250", "LMO7", "LPAL2", "LRRK2", "LTA", "LURAP1L",
  "LURAP1L-AS1", "LYZ", "MAF", "MBL2", "MCAM", "MCL1", "MEFV", "MICA", "MIR146A", "MIR21", "MIX23", "MMP1", "MMP3", "MMP7", "MMP9", "MRPS23",
  "MSN", "MTHFR", "MUCL1", "MX1", "MYC", "MYNN", "NABP1", "NAMPT", "NCOA7", "NDUFS1", "NFKB1", "NKG7", "NLRP3", "NMI", "NOD2", "NOS2",
  "NOXRED1", "NPEPPS", "NRG1", "NT5C3A", "OLR1", "OSMR", "PDCD1", "PDLIM7", "PER1", "PFDN4", "PFDN5", "PFKL", "PGD", "PGK1", "PHEX", "PI3",
  "PIK3CD", "PINK1", "PLA2G4D", "PLCG1", "PLG", "PLIN5", "PLS1", "PPARD", "PPARG", "PPARGC1A", "PPARGC1B", "PRDM1", "PRF1", "PRTN3", "PSG2",
  "PSMC2", "PSMD7", "PSME2", "PTGS1", "PTGS2", "PTH", "PTPN22", "PTX3", "PYGL", "RAC1", "RBM45", "REL", "RETN", "RGPD6", "RIT1", "RORC", "RPL15",
  "RPL36AL", "RPL41", "RPL7", "RPS19", "RPS21", "RPS26", "RPS6KB1", "RPS7", "RSAD2", "RUNX2", "RUNX3", "S100A12", "S100A8", "S100A9", "S100P",
  "SAA1", "SAMD9", "SAR1A", "SCN1A", "SEC14L2", "SEC24B", "SELL", "SERPINA1", "SERPINB1", "SERPINE1", "SF3B1", "SF3B3", "SGK1", "SH3BGRL3",
  "SIAH1", "SLC1A2", "SLC2A3", "SLC51B", "SLC7A11", "SLC7A5", "SMAD3", "SMARCA4", "SMOX", "SOCS1", "SOD2", "SOST", "SOX4", "SP7", "SPCS3",
  "SPON2", "SPP1", "SSR1", "STAT1", "STAT3", "STIM1", "SUOX", "SYT1", "TALDO1", "TBX21", "TEK", "TFPI", "TGFA", "TGFB1", "TGFBR3", "TIMP1",
  "TLR2", "TLR3", "TLR4", "TLR9", "TMBIM6", "TMEM45A", "TMPRSS11B", "TNF", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF10A", "TNFRSF1A",
  "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF11", "TNFSF13B", "TNFSF15", "TNIP1", "TOMM5", "TOMM7", "TPST2", "TPT1", "TRAF2", "TRAF3IP2", "TRAF4",
  "TRAF5", "TRBV20OR9-2", "TRIM22", "TRIM69", "TTC39B", "TYK2", "TYMP", "UBE2L3", "UQCR10", "UTP11", "VCAM1", "VDR", "VEGFA", "VIM", "WDR1",
  "WNK1", "WWOX", "XAF1", "XBP1", "YOD1", "YWHAB", "ZC3H12A", "ZFP36", "ZMIZ1", "ZNF316", "ZNF415", "ZNF483"
)

# Carregue o pacote data.table se ainda não o fez
library(data.table)

# 2. Mantenha apenas as anotações para os genes de interesse
annot_filtered <- annot[annot$Symbol %in% genes_of_interest, ]

# 3. Filtre a matriz de contagens usando os GeneIDs como chave (MUITO MAIS SEGURO)
expr_data_filtered <- gset[rownames(gset) %in% annot_filtered$GeneID, ]

# Filter expression data for genes of interest
#expr_data_filtered <- gset[annot$Symbol %in% genes_of_interest, ]

print(expr_data_filtered)
express <- rownames(expr_data_filtered)

# Mapear Entrez IDs para símbolos de genes
gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(expr_data_filtered), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Substituir os rownames pelos símbolos de genes
rownames(expr_data_filtered) <- gene_symbols

# Visualizar o data frame atualizado
print("Data frame com rownames modificados para símbolos de genes:")
print(expr_data_filtered)

# sample selection
gsms <- "111111112222222222222222222"
sml <- strsplit(gsms, split="")[[1]]

# Filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
expr_data <- expr_data_filtered[, sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("CONTROL", "AS"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(expr_data))

# pre-filter low count genes

dsa <- DESeqDataSetFromMatrix(countData=expr_data, colData=sample_info, design= ~Group)
dsb <- DESeq(dsa, test="LRT", reduced = ~ 1)  # Use LRT for all-around gene ranking

# extract results for top genes table
r <- results(dsb, alpha=0.05, pAdjustMethod ="fdr")
#tT <- r[order(r$padj)[1:429], ]
#tT <- merge(as.data.frame(tT), expr_data, by = 0, sort = F)
##write.table(tT, file = stdout(), row.names = F, sep = "\t")

plotDispEsts(dsb, main="GSE221786 Dispersion Estimates")

# create histogram plot of p-values
hist(r$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white",
     xlab = "", ylab = "", main = "GSE221786 Frequencies of padj-values")

#Depois de filtrar quais os grupos apresentam algum genes expressos
cts <- list(c("Group",groups[2],groups[1]))

# Wald test to obtain contrast-specific results
dsc <- DESeq(dsa, test="Wald", sfType="poscount")
r <- results (dsc, contrast=cts[[1]], alpha=0.05, pAdjustMethod = "fdr")

# Venn diagram
library(gplots)
all_res <- list()

for (ct in cts) {
  i <- length(all_res)
  r <- results(dsc, contrast=ct, alpha=0.05, pAdjustMethod = "fdr")
  all_res[[i + 1]] <- rownames(r)[!is.na(r$padj) & r$padj < 0.05 & abs(r$log2FoldChange) >= 1]
  names(all_res)[i + 1] <- paste(ct, collapse="_")
}
venn(all_res)

#  Differential gene expression analysis
df <- as.data.frame(r)
df <- na.omit(df)

topGenes <- rownames(df)
heatmapData <- expr_data[topGenes, ]

# Carregue as bibliotecas necessárias
library(pheatmap)
library(RColorBrewer)

# 1. Calcular a média de expressão para cada gene
media_dos_genes <- rowMeans(heatmapData)

# 2. Ordenar as médias em ordem decrescente
media_ordenada <- sort(media_dos_genes, decreasing = TRUE)

# 3. Pegar os nomes dos 20 primeiros genes
top_20_genes_por_media <- names(head(media_ordenada, 20))

# 4. Filtrar a matriz original para manter apenas os 20 genes selecionados
heatmapData_top20 <- heatmapData[top_20_genes_por_media, ]

# 5. Criar anotação da coluna com grupos
annotation_data <- data.frame(
  group = sample_info$Group,
  row.names = rownames(sample_info)
)

# 6. Reordenar colunas por grupo
ordem_grupo <- order(annotation_data$group)  # Ordena por nome do grupo (alfabética)
heatmapData_top20 <- heatmapData_top20[, ordem_grupo]
annotation_data <- annotation_data[ordem_grupo, , drop = FALSE]

# Instalar se necessário
# install.packages("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

# 1. Matriz dos 20 genes
mat <- heatmapData_top20

# 2. Reescalonar por linha
mat_scaled <- t(scale(t(mat)))

# 3. Criar anotação de grupo
ha_col <- HeatmapAnnotation(
  Group = annotation_data$group,
  col = list(Group = c("AS" = "#d62728",
                       "CONTROL" = "green"))
)

# 4. Paleta de cores tipo RdBu invertido
cores <- colorRamp2(c(-2, 0, 2), rev(RColorBrewer::brewer.pal(3, "RdBu")))

# 5. Heatmap com clusterização dentro dos grupos
Heatmap(
  mat_scaled,
  name = "Z-score",
  top_annotation = ha_col,
  col = cores,
  cluster_columns = TRUE,              # << agora está ativado
  cluster_column_slices = TRUE,       # << cluster dentro de cada grupo
  column_split = annotation_data$group, # << separa as colunas por grupo
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_title = "Heatmap of the 20 most highly expressed genes (average) from GSE221786",
  heatmap_legend_param = list(title = "Z-score")
)


# UMAP plot (multi-dimensional scaling)
expr_data_umap <- expr_data[rowSums(expr_data) > 0, ] # Remover linhas com soma zero
u <- umap(t(expr_data_umap), n_neighbors=15, random_state=123)

plot(u$layout, main="GSE221786 UMAP", xlab="", ylab="", tcl=0.1, pch=19, col="blue")
text(u$layout, labels=colnames(expr_data_umap), cex=0.7, pos=3)

# Wald test to obtain contrast-specific results
dsd <- DESeq(dsa, test = "Wald", sfType = "poscount")
r <- results(dsd, contrast = c("Group", groups[2], groups[1]), alpha = 0.05, pAdjustMethod = "fdr")


# Filtrar apenas genes significativos
sig_genes <- subset(r, padj < 0.05 & abs(log2FoldChange) >= 1)


# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2,1), cex.main=1.5)
with(sig_genes, {
  plot(log2FoldChange, -log10(padj), main=paste(groups[2], "vs", groups[1]),
       xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
  text(log2FoldChange, -log10(padj), labels=rownames(r), cex=0.6, pos=4)
})
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 1),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)

# Plotar apenas genes significativos
#MD PLOT
par(mar=c(4,4,2,1), cex.main=1.5)
with(sig_genes, {
  plot(log10(baseMean), log2FoldChange,
       main=paste(groups[2], "vs", groups[1]),
       xlab="log10(mean of normalized counts)", ylab="log2FoldChange",
       pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1)
  text(log10(baseMean), log2FoldChange, labels=rownames(sig_genes), cex=0.6, pos=4)
})
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20, col=1:2)
abline(h=0)
palette(old.pal) # restaurar paleta

# Adicione os símbolos dos genes aos pontos do gráfico

# Instale os pacotes se ainda não tiver
# install.packages("ggplot2")
# install.packages("ggrepel")

library(ggplot2)
library(ggrepel)

plotVolcano <- function(res, title = "Volcano Plot") {
  res$group <- "NS"
  res$group[res$padj < 0.05 & res$log2FoldChange > 1] <- "Up"
  res$group[res$padj < 0.05 & res$log2FoldChange < -1] <- "Down"
  
  res$label <- ifelse(res$group != "NS", rownames(res), NA)
  
  ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
    scale_color_manual(values = c("Up" = "firebrick", "Down" = "dodgerblue", "NS" = "grey80")) +
    ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
    theme_minimal() +
    labs(title = title, x = "log2 Fold Change", y = "-log10 adjusted p-value", color = "Regulation")
}

# Exemplo de uso da função
plotVolcano(r, paste(groups[2], "vs", groups[1]))

# Filtrar genes diferencialmente expressos
de_genes <- rownames(subset(r, padj < 0.05 & abs(log2FoldChange) > 1))

# Exibir os símbolos dos genes antes da conversão
de_genes

# Filter expression data for genes of interest
expr_data <- gset[annot$Symbol %in% de_genes, ]
print(expr_data)
express <- rownames(expr_data)

# Mapear Entrez IDs para símbolos de genes
gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(expr_data), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Substituir os rownames pelos símbolos de genes
rownames(expr_data) <- gene_symbols

# Visualizar o data frame atualizado
print("Data frame com rownames modificados para símbolos de genes:")
print(expr_data)

# Convert expression data to long format for ggplot2
expr_datas <- as.data.frame(expr_data)
expr_datas$Gene <- rownames(expr_datas)
expr_datas <- pivot_longer(expr_datas, cols = -Gene, names_to = "Sample", values_to = "Expression")
expr_datas <- merge(expr_datas, pheno_data, by.x = "Sample", by.y = "geo_accession")
print(expr_datas)

# Check column names and a sample of the data
print(colnames(expr_datas))
head(expr_datas)

library(dplyr)

expr_datas_filtrado <- expr_datas %>%
  filter(`characteristics_ch1.6` != "phenotype: Rheumatoid arthritis")

expr_datas_filtrado <- expr_datas_filtrado %>%
  filter(`characteristics_ch1.6` != "phenotype: Rheumatoid or Psoriatic arthritis")

# Plot boxplot for gene expression
ggplot(expr_datas_filtrado, aes(x = `characteristics_ch1.6`, y = Expression, fill = `characteristics_ch1.6`)) +
  geom_boxplot() +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Gene Expression Boxplot", x = "Diagnosis", y = "Expression")

# Plot scatter plot for gene expression
ggplot(expr_datas_filtrado, aes(x = Sample, y = Expression, color = `characteristics_ch1.6`)) +
  geom_point() +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Gene Expression Boxplot", x = "Sample", y = "Expression")

# Aggregate data for bar plot
expr_datas <- expr_datas_filtrado %>%
  group_by(Gene, `characteristics_ch1.6`) %>%
  summarise(
    mean_expression = mean(Expression, na.rm = TRUE),
    sd_expression = sd(Expression, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot bar plot for gene expression
ggplot(expr_datas, aes(x = `characteristics_ch1.6`, y = mean_expression, fill = `characteristics_ch1.6`)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression),
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Gene Expression Bar Plot", x = "Disease", y = "Mean Expression")

sample_info <- sample_info %>%
  filter(`Group` != "lesion.Pso")

sample_info <- sample_info %>%
  filter(`Group` != "non.lesion.Pso")

# Differential gene expression analysis
df <- as.data.frame(expr_data)
df <- na.omit(df)

# Create a heatmap
topGenes <- rownames(df)
heatmapData <- expr_data[topGenes, ]
pheatmap(heatmapData, cluster_rows = TRUE, show_rownames = TRUE, show_colnames = TRUE,
         annotation_col = sample_info, main = "Top Differentially Expressed Genes Heatmap")

# Filtrar genes diferencialmente expressos
de_genes <- subset(r, padj < 0.05 & abs(log2FoldChange) > 1)

# Exibir os símbolos dos genes antes da conversão
de_genes

#Fazer um data.frame de de_genes
de_genes <- data.frame(de_genes)

# Suponha que o seu dataframe se chame df e você queira selecionar as colunas "coluna1", "coluna2" e "coluna3"
de_genes <- de_genes %>% select(log2FoldChange, padj)

# Supondo que rownames(df) sejam símbolos de genes
gene_symbols <- rownames(de_genes)

# Converta símbolos de genes para Entrez IDs
gene_entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Primeiro, converta os rownames de df em uma coluna
de_genes$Symbol <- rownames(de_genes)

# Realize a junção entre a tabela original e os IDs convertidos
de_genes <- merge(de_genes, gene_entrez_ids, by.x = "Symbol", by.y = "SYMBOL", all.x = TRUE)

# Ensure unique symbols
de_genes <- de_genes[!duplicated(de_genes$ENTREZID), ]

# Visualize os primeiros resultados para verificar a adição dos IDs
de_genes

BiocManager::install("STRINGdb")
library(STRINGdb)

#Workflow Alternativo com Integração do STRINGdb no R
# Inicializar o objeto STRINGdb para um organismo específico
string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 400, input_directory = "")

# Mapeamento de genes usando STRINGdb
mapped_genes <- string_db$map(de_genes, "Symbol", removeUnmappedRows = TRUE)

# Recuperar interações para os genes mapeados
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# Visualizar a rede usando igraph ou outras ferramentas de visualização
library(igraph)
g <- graph_from_data_frame(interactions)
plot(g)

# Realizar análise de enriquecimento GO
#Processo Biológico
enrich_result <- enrichGO(gene = mapped_genes$Symbol,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP", # Ontologia Biológica. Pode ser "BP", "MF", ou "CC"
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

barplot(enrich_result, showCategory = 10, title= "GO Enrichment Analysis - BP")

# Gráfico de rede dos termos enriquecidos
cnetplot(enrich_result, showCategory = 10)

# Mapa de calor
heatplot(enrich_result, showCategory = 10)

# Realizar análise de enriquecimento GO
#Função Molecular
enrich_result <- enrichGO(gene = mapped_genes$Symbol,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "MF", # Ontologia Biológica. Pode ser "BP", "MF", ou "CC"
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
barplot(enrich_result, showCategory = 10, title = "GO Enrichment Analysis - MF")

# Gráfico de rede dos termos enriquecidos
cnetplot(enrich_result, showCategory = 10)

# Mapa de calor
heatplot(enrich_result, showCategory = 10)

# Realizar análise de enriquecimento GO
#Componente Celular
enrich_result <- enrichGO(gene = mapped_genes$Symbol,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "CC", # Ontologia Biológica. Pode ser "BP", "MF", ou "CC"
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
barplot(enrich_result, showCategory = 10, title = "GO Enrichment Analysis - CC")

# Gráfico de rede dos termos enriquecidos
cnetplot(enrich_result, showCategory = 10)

# Mapa de calor
heatplot(enrich_result, showCategory = 10)

# Perform KEGG enrichment analysis
ekegg <- enrichKEGG(gene = mapped_genes$ENTREZID, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)
barplot(ekegg, showCategory = 10, title = "KEGG Enrichment Analysis")

# Gráfico de rede dos termos enriquecidos
cnetplot(ekegg, showCategory = 10)

# Mapa de calor
heatplot(ekegg, showCategory = 10)

BiocManager::install("ReactomePA")
library(ReactomePA)

# Perform Reactome pathway analysis
ereactome <- enrichPathway(gene = mapped_genes$ENTREZID, organism = "human", pAdjustMethod = "BH", qvalueCutoff = 0.05)
barplot(ereactome, showCategory = 10, title = "Reactome Pathway Analysis")

# Gráfico de rede dos termos enriquecidos
cnetplot(ereactome, showCategory = 10)

# Mapa de calor
heatplot(ereactome, showCategory = 10)

# Carregue os pacotes necessários
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(JASPAR2020)
library(TFBSTools)
library(SummarizedExperiment)

# Carregar a base de dados TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Extração de exons por gene
exons_by_gene <- exonsBy(txdb, by = "gene")

# Obter as coordenadas dos genes (tomando a primeira e última posição dos exons)
genes_info <- range(exons_by_gene)

# Visualizar as primeiras linhas
head(genes_info)

# Obter informações de transcrição
transcripts_info <- transcripts(txdb)

# Mapear genes de interesse para IDs Entrez (substitua 'de_genes' com a lista real de genes)
gene_entrez <- mapIds(org.Hs.eg.db, keys = mapped_genes$Symbol, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Filtrar os genes de interesse
promoters_info <- subset(genes_info, names(genes_info) %in% gene_entrez)

# Definir as regiões promotoras (2 kb upstream do TSS)
promoters <- promoters(promoters_info, upstream = 2000, downstream = 0)

# Filtrar os cromossomos principais
promoters <- keepStandardChromosomes(promoters, pruning.mode = "coarse")

# Remover quaisquer regiões que excedam os limites dos cromossomos
promoters <- trim(promoters)

# Converter promoters de CompressedGRangesList para GRanges
promoters_gr <- unlist(promoters)

# Carregar a base de dados de motivos JASPAR
motifs <- getMatrixSet(JASPAR2020, opts = list(species = "Homo sapiens"))

# Realizar o enriquecimento de motivos nas regiões promotoras
motifHits <- matchMotifs(motifs, promoters_gr, genome = BSgenome.Hsapiens.UCSC.hg38)

# Visualizar os resultados de motifHits
motifHits

# Contar o número de ocorrências dos motivos em cada região promotora
motif_counts <- countOverlaps(promoters_gr, motifHits)

# Adicionar as contagens ao GRanges com regiões promotoras
promoters_gr$motif_counts <- motif_counts

# Visualizar as primeiras linhas com contagens
head(promoters_gr)

# Resumo das contagens de motivos
summary(promoters_gr$motif_counts)

# Obter os nomes dos motivos e TFs
motif_names <- names(motifHits)

# Extrair a tabela de motivos do objeto RangedSummarizedExperiment
motif_data <- assays(motifHits)[[1]] # Supondo que a tabela de motivos esteja na primeira lista de assays

# Verificar as primeiras linhas da tabela de motivos
head(motif_data)

# Obter informações sobre os TFs
# O nome dos TFs pode estar nos colnames ou em um metadata associado
tf_names <- colnames(motif_data) # Ajuste conforme necessário

# Adicionar a contagem de motivos se disponível
motif_counts <- rowSums(motif_data) # Contar a ocorrência de motivos (ajuste conforme necessário)

# Verificar as dimensões do objeto motif_data
dim(motif_data) # Número de linhas e colunas

# Verificar o comprimento de tf_names e motif_counts
tf_names_length <- length(tf_names)
motif_counts_length <- length(motif_counts)

# Imprimir os comprimentos para diagnóstico
cat("Comprimento de tf_names:", tf_names_length, "\n")
cat("Comprimento de motif_counts:", motif_counts_length, "\n")

# Verificar se as dimensões do motif_data correspondem ao número de TFs e motivos
motif_data_dims <- dim(motif_data)
cat("Dimensões de motif_data (linhas, colunas):", motif_data_dims, "\n")

# Verificar se a quantidade de motivos é maior ou menor
if (tf_names_length > motif_counts_length) {
  # Verificar colunas adicionais
  tf_names <- tf_names[1:motif_counts_length] # Ajustar para o comprimento de motif_counts
} else {
  # Ajustar motif_counts para coincidir com tf_names
  motif_counts <- motif_counts[1:tf_names_length]
}

# Agora, crie o data frame com comprimentos ajustados
motif_summary <- data.frame(
  tf_name = tf_names,
  motif_count = motif_counts
)

# Visualizar os resultados
head(motif_summary)

# Verificar se há TFs ou motivos que estão faltando
head(tf_names) # Primeiros nomes de TFs
head(motif_counts) # Primeiras contagens de motivos

# Verificar a correspondência entre nomes de TFs e as colunas de motif_data
all(tf_names %in% colnames(motif_data)) # Deve retornar TRUE se todos os nomes de TFs estiverem presentes

# Visualizar os primeiros resultados
head(motif_summary)

# Ordenar e identificar principais TFs
tf_summary <- motif_summary[order(motif_summary$motif_count, decreasing = TRUE), ]
top_tf_summary <- head(tf_summary, 12) # Ajuste o número conforme necessário

# Visualizar os principais TFs
top_tf_summary

# Mapear Entrez IDs para símbolos de genes
gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(top_tf_summary), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Substituir os rownames pelos símbolos de genes
rownames(top_tf_summary) <- gene_symbols

# Visualizar o data frame atualizado
print(top_tf_summary)

# Encontrar o motivo pelo ID
motif <- motifs[["MA0066.1"]]
motif
