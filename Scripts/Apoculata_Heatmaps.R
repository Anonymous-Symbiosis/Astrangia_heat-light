#Astrangia poculata heat & light: Candidate circadian genes & heatmaps for DEGs associated with immuniry and cell cycle processes
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Supplementary Fig. S4-S7
#Note: In this script: apo = white corals, sym = brown corals

load("dds.RData")

library(ggplot2)
library(DESeq2)
library(dplyr)
library(stringr)
library(tibble)
library(pheatmap)


#----------Prepare LFC data for each treatment----------

##----------White corals----------
res_ADC <- results(dds_host.apo, contrast=c("treat", "NL_C", "CL_C"))
head(res_ADC)
summary(res_ADC)
#LFC > 0 (up)       : 237, 0.67%
#LFC < 0 (down)     : 82, 0.23%

res_AHC <- results(dds_host.apo, contrast=c("treat", "HL_C", "CL_C"))
head(res_AHC)
summary(res_AHC)
#LFC > 0 (up)       : 16, 0.045%
#LFC < 0 (down)     : 9, 0.025%

res_ADH <- results(dds_host.apo, contrast=c("treat", "NL_H", "CL_C"))
head(res_ADH)
summary(res_ADH)
#LFC > 0 (up)       : 824, 2.3%
#LFC < 0 (down)     : 1093, 3.1%

res_ACH <- results(dds_host.apo, contrast=c("treat", "CL_H", "CL_C"))
head(res_ACH)
summary(res_ACH)
#LFC > 0 (up)       : 136, 0.39%
#LFC < 0 (down)     : 155, 0.44%

res_AHH <- results(dds_host.apo, contrast=c("treat", "HL_H", "CL_C"))
head(res_AHH)
summary(res_AHH)
#LFC > 0 (up)       : 613, 1.7%
#LFC < 0 (down)     : 565, 1.6%

##----------Brown corals----------
res_SDC <- results(dds_host.sym, contrast=c("treat", "NL_C", "CL_C"))
head(res_SDC)
summary(res_SDC)
#LFC > 0 (up)       : 27, 0.077%
#LFC < 0 (down)     : 26, 0.074%

res_SHC <- results(dds_host.sym, contrast=c("treat", "HL_C", "CL_C"))
head(res_SHC)
summary(res_SHC)
#LFC > 0 (up)       : 15, 0.043%
#LFC < 0 (down)     : 9, 0.026%

res_SDH <- results(dds_host.sym, contrast=c("treat", "NL_H", "CL_C"))
head(res_SDH)
summary(res_SDH)
#LFC > 0 (up)       : 316, 0.9%
#LFC < 0 (down)     : 370, 1.1%

res_SCH <- results(dds_host.sym, contrast=c("treat", "CL_H", "CL_C"))
head(res_SCH)
summary(res_SCH)
#LFC > 0 (up)       : 173, 0.49%
#LFC < 0 (down)     : 260, 0.74%

res_SHH <- results(dds_host.sym, contrast=c("treat", "HL_H", "CL_C"))
head(res_SHH)
summary(res_SHH)
#LFC > 0 (up)       : 792, 2.3%
#LFC < 0 (down)     : 1010, 2.9%


#----------Circadian Genes----------

iso2gene = read.csv("../Data/apoculata_iso2gene.csv", header = TRUE) 
iso2gene$GN <- str_extract(iso2gene$GO_Description, "GN=\\w+")
iso2gene$GN <- gsub("GN=", "", iso2gene$GN)

#extract genes of interest
Circadian <- iso2gene[grepl("CLOCK|Cry|PARbZIP|Helt|Timeout|CK1|Hlf|Hebp2|Tef", iso2gene$GN, ignore.case = TRUE) &
    !iso2gene$GN %in% c("rbck1", "RACK1"),]

sig_ADH <- rownames(res_ADH)[which(res_ADH$padj < 0.1)]
sig_AHH <- rownames(res_AHH)[which(res_AHH$padj < 0.1)]
sig_SDH <- rownames(res_SDH)[which(res_SDH$padj < 0.1)]
sig_SHH <- rownames(res_SHH)[which(res_SHH$padj < 0.1)]
sig_ACH <- rownames(res_ACH)[which(res_ACH$padj < 0.1)]
sig_AHC <- rownames(res_AHC)[which(res_AHC$padj < 0.1)]
sig_ADC <- rownames(res_ADC)[which(res_ADC$padj < 0.1)]
sig_SCH <- rownames(res_SCH)[which(res_SCH$padj < 0.1)]
sig_SHC <- rownames(res_SHC)[which(res_SHC$padj < 0.1)]
sig_SDC <- rownames(res_SDC)[which(res_SDC$padj < 0.1)]

sig_union <- Reduce(union, list(sig_ADH, sig_AHH, sig_SDH, sig_SHH, sig_ADC, sig_ACH, sig_AHC, sig_SDC, sig_SCH, sig_SHC))
gene_use <- intersect(Circadian$gene, sig_union)

# LFC matrix
lfc_mat <- data.frame(
  gene = gene_use,
  ADC = res_ADC[gene_use, "log2FoldChange"],
  AHC = res_AHC[gene_use, "log2FoldChange"],
  ADH = res_ADH[gene_use, "log2FoldChange"],
  ACH = res_ACH[gene_use, "log2FoldChange"],
  AHH = res_AHH[gene_use, "log2FoldChange"],
  SDC = res_SDC[gene_use, "log2FoldChange"],
  SHC = res_SHC[gene_use, "log2FoldChange"],
  SDH = res_SDH[gene_use, "log2FoldChange"],
  SCH = res_SCH[gene_use, "log2FoldChange"],
  SHH = res_SHH[gene_use, "log2FoldChange"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# padj matrix
padj_mat <- data.frame(
  gene = gene_use,
  ADC = res_ADC[gene_use, "padj"],
  AHC = res_AHC[gene_use, "padj"],
  ADH = res_ADH[gene_use, "padj"],
  ACH = res_ACH[gene_use, "padj"],
  AHH = res_AHH[gene_use, "padj"],
  SDC = res_SDC[gene_use, "padj"],
  SHC = res_SHC[gene_use, "padj"],
  SDH = res_SDH[gene_use, "padj"],
  SCH = res_SCH[gene_use, "padj"],
  SHH = res_SHH[gene_use, "padj"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# make sure rownames/colnames match exactly
padj_mat <- padj_mat[rownames(lfc_mat), colnames(lfc_mat)]

lfc_long <- lfc_mat %>%
  as.data.frame() %>%
  rownames_to_column("GN") %>%
  pivot_longer(
    cols = -GN,
    names_to = "contrast",
    values_to = "log2FoldChange"
  )

padj_long <- padj_mat %>%
  as.data.frame() %>%
  rownames_to_column("GN") %>%
  pivot_longer(
    cols = -GN,
    names_to = "contrast",
    values_to = "padj"
  )

circadian_dot_df <- lfc_long %>%
  left_join(padj_long, by = c("GN", "contrast")) %>%
  mutate(
    sig = !is.na(padj) & padj < 0.1,
    neglog10padj = -log10(padj),
    contrast = factor(
      contrast,
      levels = c("ADC", "AHC", "ADH", "ACH", "AHH",
                 "SDC", "SHC", "SDH", "SCH", "SHH")
    ),
    GN = factor(GN, levels = rev(rownames(lfc_mat)))
  )

plot_df <- circadian_dot_df %>%
  filter(sig) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c("ADC", "AHC", "ADH", "ACH", "AHH",
                 "SDC", "SHC", "SDH", "SCH", "SHH")
    ),
    GN = factor(GN, levels = sort(unique(as.character(GN))))
  )

Circadian_dotplot <- ggplot(
  plot_df,
  aes(x = contrast, y = GN)
) +
  geom_point(
    aes(
      color = log2FoldChange,
      size = neglog10padj
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_color_gradient2(
    low = "turquoise4",
    mid = "grey85",
    high = "sienna3",
    midpoint = 0,
    limits = c(-8, 0),
    oob = squish,
    name = "log2FC"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = "-log10(padj)"
  ) +
  scale_y_discrete(
    limits = sort(unique(as.character(plot_df$GN))),
    position = "right"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

Circadian_dotplot
ggsave("../Figures/Circadian_dotplot.svg", Circadian_dotplot, width = 6, height = 4, units = "in", device = "svg")


#----------Immunity heatmaps----------

##----------Selected Genes----------

#extract genes of interest
TLR <- iso2gene[grepl("NFKB1|TRAF3|MAP3K7|TRAF6|MYD88_TLR|AP1G1|JUN|MAP2K7", iso2gene$GN, ignore.case = TRUE), ]
NFKB <- iso2gene[grepl("NFKB1|TRAF3|MAP3K7|TRAF6|MYD88_TLR|AP1G1|JUN|MAP2K7|DDX4|MALT1|CSNK1A1|CSNK1G2|CSNK1E", iso2gene$GN, ignore.case = TRUE), ]
TNF <- iso2gene[grepl("NFKB1|TRAF3|MAP3K7|AP1G1|JUN|MAP2K7|DNM1L|RPS6KA5|BCL3|CASP7|IFIH1|PIN1", iso2gene$GN, ignore.case = TRUE), ]
NOD <- iso2gene[grepl("NFKB1|TRAF3|MAP3K7|TRAF6|MYD88_TLR|AP1G1|PIN1|RNF31", iso2gene$GN, ignore.case = TRUE), ]
MAPK <- iso2gene[grepl("NFKB1|MAP3K7|TRAF6|MYD88_TLR|AP1G1|JUN|MAP2K7|MAP4K4|MYC|MYCBPAP|mycbp|RPS6KA5|GRB2|DAXX|FGFR4|FGFR2|FGF7|FGF20|FGF18|FGFR4|MAP3K13|MAP3K2|MAP3K20|MAP3K2|MAPK7|MKNK1|PLA2G4A|PRKCA|RASA1|TAOK1|Tie", iso2gene$GN, ignore.case = TRUE), ]
Ctype <- iso2gene[grepl("NFKB1|AP1G1|JUN|DDX4|MALT1|BCL3|CALM3", iso2gene$GN, ignore.case = TRUE), ]
TGFB <- iso2gene[grepl("MYC|MYCBPAP|mycbp|CREB1|ACVR1|CUL1", iso2gene$GN, ignore.case = TRUE), ]

TLR$Pathway  <- "TLR"
NFKB$Pathway <- "NFKB"
TNF$Pathway  <- "TNF"
NOD$Pathway  <- "NOD"
MAPK$Pathway <- "MAPK"
Ctype$Pathway <- "Ctype"
TGFB$Pathway <- "TGFB"

merged_pathways <- bind_rows(
  TLR, NFKB, TNF, NOD, MAPK, Ctype, TGFB
)

sig_ADH <- rownames(res_ADH)[which(res_ADH$padj < 0.1)]
sig_AHH <- rownames(res_AHH)[which(res_AHH$padj < 0.1)]
sig_SDH <- rownames(res_SDH)[which(res_SDH$padj < 0.1)]
sig_SHH <- rownames(res_SHH)[which(res_SHH$padj < 0.1)]
sig_ACH <- rownames(res_ACH)[which(res_ACH$padj < 0.1)]
sig_AHC <- rownames(res_AHC)[which(res_AHC$padj < 0.1)]
sig_ADC <- rownames(res_ADC)[which(res_ADC$padj < 0.1)]
sig_SCH <- rownames(res_SCH)[which(res_SCH$padj < 0.1)]
sig_SHC <- rownames(res_SHC)[which(res_SHC$padj < 0.1)]
sig_SDC <- rownames(res_SDC)[which(res_SDC$padj < 0.1)]

sig_union <- Reduce(union, list(sig_ADH, sig_AHH, sig_SDH, sig_SHH, sig_ADC, sig_ACH, sig_AHC, sig_SDC, sig_SCH, sig_SHC))
gene_use <- intersect(merged_pathways$gene, sig_union)

# LFC matrix
lfc_mat <- data.frame(
  gene = gene_use,
  ADC = res_ADC[gene_use, "log2FoldChange"],
  AHC = res_AHC[gene_use, "log2FoldChange"],
  ADH = res_ADH[gene_use, "log2FoldChange"],
  ACH = res_ACH[gene_use, "log2FoldChange"],
  AHH = res_AHH[gene_use, "log2FoldChange"],
  SDC = res_SDC[gene_use, "log2FoldChange"],
  SHC = res_SHC[gene_use, "log2FoldChange"],
  SDH = res_SDH[gene_use, "log2FoldChange"],
  SCH = res_SCH[gene_use, "log2FoldChange"],
  SHH = res_SHH[gene_use, "log2FoldChange"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# padj matrix
padj_mat <- data.frame(
  gene = gene_use,
  ADC = res_ADC[gene_use, "padj"],
  AHC = res_AHC[gene_use, "padj"],
  ADH = res_ADH[gene_use, "padj"],
  ACH = res_ACH[gene_use, "padj"],
  AHH = res_AHH[gene_use, "padj"],
  SDC = res_SDC[gene_use, "padj"],
  SHC = res_SHC[gene_use, "padj"],
  SDH = res_SDH[gene_use, "padj"],
  SCH = res_SCH[gene_use, "padj"],
  SHH = res_SHH[gene_use, "padj"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# make sure rownames/colnames match exactly
padj_mat <- padj_mat[rownames(lfc_mat), colnames(lfc_mat)]

# significance labels
sig_labels <- ifelse(
  !is.na(padj_mat) & padj_mat < 0.1,
  sprintf("%.2f", lfc_mat),
  ""
)
sig_labels[is.na(sig_labels)] <- ""
Immune_heatmap.selected<-pheatmap(
  lfc_mat,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "none",
  color = colorRampPalette(c("turquoise4", "white", "sienna3"))(100),
  breaks = seq(-3, 3, length.out = 101),
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  border_color = "grey70",
  fontsize_row = 10,
  display_numbers = sig_labels,
  number_color = "black",
  fontsize_number = 10,
)
Immune_heatmap.selected
ggsave("../Figures/Immune_heatmap.svg", Immune_heatmap.selected, width = 8, height = 7, units = "in", device = "svg")

##-----WGCNA Brown Module-----
Brown <- read.csv("../Data/GO_Output/BP_brown_fisher_0.2.csv", sep="",header = TRUE)

# choose GO terms related to immune functions
go_Immune <- c("GO:0002376", "GO:0002683", "GO:0006955", "GO:0045087","GO:0045088", "GO:0050776","")   

# extract genes from GO terms
genes_from_go <- Brown %>%
  filter(term %in% go_Immune, value > 0) %>%
  pull(seq) %>%
  unique()

# annotate with iso2gene
genes_from_go <- iso2gene %>%
  filter(gene %in% genes_from_go) %>%
  select(gene, GN, GO_Description)

gene_use.brown <- intersect(genes_from_go$gene, sig_union)

lfc_mat <- data.frame(
  gene = gene_use.brown,
  ADC = res_ADC[gene_use.brown, "log2FoldChange"],
  AHC = res_AHC[gene_use.brown, "log2FoldChange"],
  ADH = res_ADH[gene_use.brown, "log2FoldChange"],
  ACH = res_ACH[gene_use.brown, "log2FoldChange"],
  AHH = res_AHH[gene_use.brown, "log2FoldChange"],
  SDC = res_SDC[gene_use.brown, "log2FoldChange"],
  SHC = res_SHC[gene_use.brown, "log2FoldChange"],
  SDH = res_SDH[gene_use.brown, "log2FoldChange"],
  SCH = res_SCH[gene_use.brown, "log2FoldChange"],
  SHH = res_SHH[gene_use.brown, "log2FoldChange"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# padj matrix
padj_mat <- data.frame(
  gene = gene_use.brown,
  ADC = res_ADC[gene_use.brown, "padj"],
  AHC = res_AHC[gene_use.brown, "padj"],
  ADH = res_ADH[gene_use.brown, "padj"],
  ACH = res_ACH[gene_use.brown, "padj"],
  AHH = res_AHH[gene_use.brown, "padj"],
  SDC = res_SDC[gene_use.brown, "padj"],
  SHC = res_SHC[gene_use.brown, "padj"],
  SDH = res_SDH[gene_use.brown, "padj"],
  SCH = res_SCH[gene_use.brown, "padj"],
  SHH = res_SHH[gene_use.brown, "padj"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# make sure rownames/colnames match exactly
padj_mat <- padj_mat[rownames(lfc_mat), colnames(lfc_mat)]

# significance labels
sig_labels <- ifelse(
  !is.na(padj_mat) & padj_mat < 0.1,
  sprintf("%.2f", lfc_mat),
  ""
)
sig_labels[is.na(sig_labels)] <- ""
Immune_heatmap.brown<-pheatmap(
  lfc_mat,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "none",
  color = colorRampPalette(c("turquoise4", "white", "sienna3"))(100),
  border_color = "grey70",
  breaks = seq(-4, 4, length.out = 101),
  legend_breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
  fontsize_row = 10,
  display_numbers = sig_labels,
  number_color = "black",
  fontsize_number = 10,
)
Immune_heatmap.brown
ggsave("../Figures/Immune_heatmap_brown.svg", Immune_heatmap.brown, width = 8, height = 7, units = "in", device = "svg")


#----------Cell Cycle Heatmap----------

#extract genes of interest
genes <- iso2gene[grepl("c-MYC|E2F1|E2F2|E2F3|CDK4|CDK6|Dp-1|Dp-2|HDAC1|Mcm3|ORC3|Bub1|
                        CDC25B|CDC25C|CDK1|GADD45|Mdm2|14-3-3|cdc14A|Cdc20|Cyclin|Mob1|Mps1|Plk1", iso2gene$GN, ignore.case = TRUE), ]

gene_use <- intersect(genes$gene, sig_union)

# LFC matrix
lfc_mat <- data.frame(
  gene = gene_use,
  ADC = res_ADC[gene_use, "log2FoldChange"],
  AHC = res_AHC[gene_use, "log2FoldChange"],
  ADH = res_ADH[gene_use, "log2FoldChange"],
  ACH = res_ACH[gene_use, "log2FoldChange"],
  AHH = res_AHH[gene_use, "log2FoldChange"],
  SDC = res_SDC[gene_use, "log2FoldChange"],
  SHC = res_SHC[gene_use, "log2FoldChange"],
  SDH = res_SDH[gene_use, "log2FoldChange"],
  SCH = res_SCH[gene_use, "log2FoldChange"],
  SHH = res_SHH[gene_use, "log2FoldChange"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# padj matrix
padj_mat <- data.frame(
  gene = gene_use,
  ADC = res_ADC[gene_use, "padj"],
  AHC = res_AHC[gene_use, "padj"],
  ADH = res_ADH[gene_use, "padj"],
  ACH = res_ACH[gene_use, "padj"],
  AHH = res_AHH[gene_use, "padj"],
  SDC = res_SDC[gene_use, "padj"],
  SHC = res_SHC[gene_use, "padj"],
  SDH = res_SDH[gene_use, "padj"],
  SCH = res_SCH[gene_use, "padj"],
  SHH = res_SHH[gene_use, "padj"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# make sure rownames/colnames match exactly
padj_mat <- padj_mat[rownames(lfc_mat), colnames(lfc_mat)]

# significance labels
sig_labels <- ifelse(
  !is.na(padj_mat) & padj_mat < 0.1,
  sprintf("%.2f", lfc_mat),
  ""
)
sig_labels[is.na(sig_labels)] <- ""
CellCycle_heatmap.selected<-pheatmap(
  lfc_mat,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "none",
  color = colorRampPalette(c("turquoise4", "white", "sienna3"))(100),
  breaks = seq(-3, 3, length.out = 101),
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  border_color = "grey70",
  fontsize_row = 10,
  display_numbers = sig_labels,
  number_color = "black",
  fontsize_number = 10,
)

##-----WGCNA Yellow Module-----
Yellow <- read.csv("../Data/GO_Output/BP_yellow_fisher_0.2.csv", sep="",header = TRUE)

# choose GO terms related to immune functions
go_CellCycle <- c("GO:0007049", "GO:0007093", "GO:0000075", "GO:0007093;GO:0000075", "GO:0007095;GO:0044818", "GO:0007095", "GO:0044818",
                  "GO:0010564", "GO:0010972;GO:1902750", "GO:0010972", "GO:1902750", "GO:0022402", "GO:0033314;GO:0000076", "GO:0033314", "GO:0000076",
                  "GO:0045786;GO:0010948", "GO:0045786", "GO:0010948", "GO:0045930;GO:1901988", "GO:0045930", "GO:1901988", 
                  "GO:1901987;GO:1901990", "GO:1901987", "GO:1901990",
                  "GO:1901991","GO:1903047")

# extract genes from GO terms
genes_from_go <- Yellow %>%
  filter(term %in% go_CellCycle, value > 0) %>%
  pull(seq) %>%
  unique()

# annotate with iso2gene
genes_from_go <- iso2gene %>%
  filter(gene %in% genes_from_go) %>%
  select(gene, GN, GO_Description)

gene_use.yellow <- intersect(genes_from_go$gene, sig_union)

lfc_mat <- data.frame(
  gene = gene_use.yellow,
  ADC = res_ADC[gene_use.yellow, "log2FoldChange"],
  AHC = res_AHC[gene_use.yellow, "log2FoldChange"],
  ADH = res_ADH[gene_use.yellow, "log2FoldChange"],
  ACH = res_ACH[gene_use.yellow, "log2FoldChange"],
  AHH = res_AHH[gene_use.yellow, "log2FoldChange"],
  SDC = res_SDC[gene_use.yellow, "log2FoldChange"],
  SHC = res_SHC[gene_use.yellow, "log2FoldChange"],
  SDH = res_SDH[gene_use.yellow, "log2FoldChange"],
  SCH = res_SCH[gene_use.yellow, "log2FoldChange"],
  SHH = res_SHH[gene_use.yellow, "log2FoldChange"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# padj matrix
padj_mat <- data.frame(
  gene = gene_use.yellow,
  ADC = res_ADC[gene_use.yellow, "padj"],
  AHC = res_AHC[gene_use.yellow, "padj"],
  ADH = res_ADH[gene_use.yellow, "padj"],
  ACH = res_ACH[gene_use.yellow, "padj"],
  AHH = res_AHH[gene_use.yellow, "padj"],
  SDC = res_SDC[gene_use.yellow, "padj"],
  SHC = res_SHC[gene_use.yellow, "padj"],
  SDH = res_SDH[gene_use.yellow, "padj"],
  SCH = res_SCH[gene_use.yellow, "padj"],
  SHH = res_SHH[gene_use.yellow, "padj"]
) %>%
  left_join(
    iso2gene %>% select(gene, GN),
    by = "gene"
  ) %>%
  mutate(
    GN = ifelse(is.na(GN), gene, GN),
    GN = make.unique(GN)
  ) %>%
  column_to_rownames("GN") %>%
  select(-gene) %>%
  as.matrix()

# make sure rownames/colnames match exactly
padj_mat <- padj_mat[rownames(lfc_mat), colnames(lfc_mat)]

# significance labels
sig_labels <- ifelse(
  !is.na(padj_mat) & padj_mat < 0.1,
  sprintf("%.2f", lfc_mat),
  ""
)
sig_labels[is.na(sig_labels)] <- ""
CellCycle_heatmap.yellow<-pheatmap(
  lfc_mat,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "none",
  color = colorRampPalette(c("purple4", "white", "goldenrod2"))(100),
  border_color = "grey70",
  breaks = seq(-2,2, length.out = 101),
  legend_breaks = c(-2, -1, 0, 1, 2),
  fontsize_row = 10,
  display_numbers = sig_labels,
  number_color = "black",
  fontsize_number = 10,
)
CellCycle_heatmap.yellow
ggsave("../Figures/CellCycle_heatmap_yellow.svg", CellCycle_heatmap.yellow, width = 8, height = 7, units = "in", device = "svg")
