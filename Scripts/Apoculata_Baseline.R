#Astrangia poculata heat & light: baseline GE comparisons
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Fig. 2C in the manuscript

load("dds.Rdata")

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(plotly)


#----------Volcano Plot----------
BP.GO_Sym.control <- read.csv("../Data/GO_Output/MWU_BP_GO_Sym_control_signedP.csv", sep="")
BP.GO_Sym.control$significant <- "Not Sig"
BP.GO_Sym.control$significant[BP.GO_Sym.control$delta.rank > 0 & BP.GO_Sym.control$p.adj < 0.05] <- "Up"
BP.GO_Sym.control$significant[BP.GO_Sym.control$delta.rank < 0 & BP.GO_Sym.control$p.adj < 0.05] <- "Down"
BP.GO_Sym.control$negLog10Padj <- -log10(BP.GO_Sym.control$p.adj)

vol.GO.control <-ggplot(BP.GO_Sym.control, aes(x = delta.rank, y = negLog10Padj)) +
  geom_point(aes(fill = significant), shape = 21, color = "snow2", size = 3, stroke = 1) +
  scale_fill_manual(values = c("wheat2", "snow3", "tan4")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    x = "GO BP Delta Rank",
    y = "-log10 padj"
  ) +
  theme_minimal()+
  theme(
    plot.title=element_blank(),
    axis.ticks=element_line(linewidth=1),
    axis.title = element_text(size=14, face = "bold"),
    axis.text=element_text(size=12),
    panel.grid.major = element_line(color = scales::alpha("grey", 0.1)),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none")
vol.GO.control
ggsave("../Figures/volcano_symbiosis_control_GO.svg", vol.GO.control, width = 4.68, height =3.5, units = "in", device = "svg")

plot_ly(
  data = BP.GO_Sym.control,
  x = ~delta.rank,
  y = ~negLog10Padj,
  type = "scatter",
  mode = "markers",
  color = ~significant,
  colors = c("wheat2", "snow3", "tan4"),
  text = ~paste("Term:", term, "<br>GO:", name, "<br>padj:", negLog10Padj),
  hoverinfo = "text"
)

