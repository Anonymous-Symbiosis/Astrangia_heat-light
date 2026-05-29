#Astrangia poculata heat & light: Venn diagrams & LFC correlations of shared DEGs
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Fig. 3A, 4A, 4B, 5A & 5B in the manuscript

load("dds.RData")
load("valState.RData")

library(dplyr)
library(ggplot2)
library(eulerr)


#----------Extract DEGs from each treatment x symbiotic phenotype----------

##----------White corals----------
ADC_deg<-filter(as.data.frame(valState_ADC), padj<0.1)
ADC_deg.up<-filter(as.data.frame(valState_ADC), padj<0.1 & log2change>0)
ADC_deg.down<-filter(as.data.frame(valState_ADC), padj<0.1 & log2change<0)
ADC_deg$gene<-row.names(ADC_deg)
ADC_deg.up$gene<-row.names(ADC_deg.up)
ADC_deg.down$gene<-row.names(ADC_deg.down)

AHC_deg<-filter(as.data.frame(valState_AHC), padj<0.1)
AHC_deg.up<-filter(as.data.frame(valState_AHC), padj<0.1 & log2change>0)
AHC_deg.down<-filter(as.data.frame(valState_AHC), padj<0.1 & log2change<0)
AHC_deg$gene<-row.names(AHC_deg)
AHC_deg.up$gene<-row.names(AHC_deg.up)
AHC_deg.down$gene<-row.names(AHC_deg.down)

ADH_deg<-filter(as.data.frame(valState_ADH), padj<0.1)
ADH_deg.up<-filter(as.data.frame(valState_ADH), padj<0.1 & log2change>0)
ADH_deg.down<-filter(as.data.frame(valState_ADH), padj<0.1 & log2change<0)
ADH_deg$gene<-row.names(ADH_deg)
ADH_deg.up$gene<-row.names(ADH_deg.up)
ADH_deg.down$gene<-row.names(ADH_deg.down)

ACH_deg<-filter(as.data.frame(valState_ACH), padj<0.1)
ACH_deg.up<-filter(as.data.frame(valState_ACH), padj<0.1 & log2change>0)
ACH_deg.down<-filter(as.data.frame(valState_ACH), padj<0.1 & log2change<0)
ACH_deg$gene<-row.names(ACH_deg)
ACH_deg.up$gene<-row.names(ACH_deg.up)
ACH_deg.down$gene<-row.names(ACH_deg.down)

AHH_deg<-filter(as.data.frame(valState_AHH), padj<0.1)
AHH_deg.up<-filter(as.data.frame(valState_AHH), padj<0.1 & log2change>0)
AHH_deg.down<-filter(as.data.frame(valState_AHH), padj<0.1 & log2change<0)
AHH_deg$gene<-row.names(AHH_deg)
AHH_deg.up$gene<-row.names(AHH_deg.up)
AHH_deg.down$gene<-row.names(AHH_deg.down)

##----------Brown corals----------
SDC_deg<-filter(as.data.frame(valState_SDC), padj<0.1)
SDC_deg.up<-filter(as.data.frame(valState_SDC), padj<0.1 & log2change>0)
SDC_deg.down<-filter(as.data.frame(valState_SDC), padj<0.1 & log2change<0)
SDC_deg$gene<-row.names(SDC_deg)
SDC_deg.up$gene<-row.names(SDC_deg.up)
SDC_deg.down$gene<-row.names(SDC_deg.down)

SHC_deg<-filter(as.data.frame(valState_SHC), padj<0.1)
SHC_deg.up<-filter(as.data.frame(valState_SHC), padj<0.1 & log2change>0)
SHC_deg.down<-filter(as.data.frame(valState_SHC), padj<0.1 & log2change<0)
SHC_deg$gene<-row.names(SHC_deg)
SHC_deg.up$gene<-row.names(SHC_deg.up)
SHC_deg.down$gene<-row.names(SHC_deg.down)

SDH_deg<-filter(as.data.frame(valState_SDH), padj<0.1)
SDH_deg.up<-filter(as.data.frame(valState_SDH), padj<0.1 & log2change>0)
SDH_deg.down<-filter(as.data.frame(valState_SDH), padj<0.1 & log2change<0)
SDH_deg$gene<-row.names(SDH_deg)
SDH_deg.up$gene<-row.names(SDH_deg.up)
SDH_deg.down$gene<-row.names(SDH_deg.down)

SCH_deg<-filter(as.data.frame(valState_SCH), padj<0.1)
SCH_deg.up<-filter(as.data.frame(valState_SCH), padj<0.1 & log2change>0)
SCH_deg.down<-filter(as.data.frame(valState_SCH), padj<0.1 & log2change<0)
SCH_deg$gene<-row.names(SCH_deg)
SCH_deg.up$gene<-row.names(SCH_deg.up)
SCH_deg.down$gene<-row.names(SCH_deg.down)

SHH_deg<-filter(as.data.frame(valState_SHH), padj<0.1)
SHH_deg.up<-filter(as.data.frame(valState_SHH), padj<0.1 & log2change>0)
SHH_deg.down<-filter(as.data.frame(valState_SHH), padj<0.1 & log2change<0)
SHH_deg$gene<-row.names(SHH_deg)
SHH_deg.up$gene<-row.names(SHH_deg.up)
SHH_deg.down$gene<-row.names(SHH_deg.down)

#---------Compare DEG counts between phenotypes within each treatment----------
#Control Temperature, Dark
fit.DC <- euler(list(
  "White" = ADC_deg$gene,
  "Brown" = SDC_deg$gene
))
plot(fit.DC,
     fills = list(
       fill = c("wheat2", "tan4"),
       alpha = 1),
     edges = list(col = "black", lwd=2),
     labels = FALSE,
     quantities = FALSE)

#Control Temperature, High Light
fit.HC <- euler(list(
  "White" = AHC_deg$gene,
  "Brown" = SHC_deg$gene
))
plot(fit.HC,
     fills = list(
       fill = c("wheat2", "tan4"),
       alpha = 1),
     edges = list(col = "black", lwd=2),
     labels = FALSE,
     quantities = FALSE)

#Heated, Dark
fit.DH <- euler(list(
  "White" = ADH_deg$gene,
  "Brown" = SDH_deg$gene
))
plot(fit.DH,
     fills = list(
       fill = c("wheat2", "tan4"),
       alpha = 1),
     edges = list(col = "black", lwd=2),
     labels = FALSE,
     quantities = FALSE)

#Heated, Control Light
fit.CH <- euler(list(
  "White" = ACH_deg$gene,
  "Brown" = SCH_deg$gene
))

plot(fit.CH,
     fills = list(
       fill = c("wheat2", "tan4"),
       alpha = 1),
     edges = list(col = "black", lwd=2),
     labels = FALSE,
     quantities = FALSE)

#Heated, High Light
fit.HH <- euler(list(
  "White" = AHH_deg$gene,
  "Brown" = SHH_deg$gene
))
plot(fit.HH,
     fills = list(
       fill = c("wheat2", "tan4"),
       alpha = 1),
     edges = list(col = "black", lwd=2),
     labels = FALSE,
     quantities = FALSE)


#--------Extract DEGs from Venn diagrams for GO analyses-----------
##----------Control Temperature, Dark----------
#overlapping DEGs between ADC & SDC
DC.shared <- intersect(ADC_deg$gene, SDC_deg$gene)
ADC.shared <- ADC_deg[row.names(ADC_deg) %in% DC.shared, ]
SDC.shared <- SDC_deg[row.names(SDC_deg) %in% DC.shared,]
DC.shared <- merge(ADC.shared, SDC.shared, by="row.names", suffixes = c("_ADC", "_SDC"))
#plot to see whether there are antagonistic genes
ggplot(DC.shared, aes(x = log2change_SDC, y = log2change_ADC)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
  xlab("log2 Fold Change (Brown)") +
  ylab("log2 Fold Change (White)") +
  theme_minimal() 
#create input file for GO-MWU
#total DEGs shared between AHH & SHH
allGenes <- rownames(dds_host) 
rownames(DC.shared)=DC.shared$Row.names
DC.shared <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(DC.shared), 1, 0)
)
sum(DC.shared$fisher) #24
write.csv(DC.shared, "../Data/GO_Input/DC.shared_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs shared between ADC & SDC
DC.shared_up <- intersect(ADC_deg.up$gene, SDC_deg.up$gene)
DC.shared_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% DC.shared_up, 1, 0)
)
sum(DC.shared_up$fisher) #13
write.csv(DC.shared_up, "../Data/GO_Input/DC.shared_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs shared between ADC & SDC
DC.shared_down <- intersect(ADC_deg.down$gene, SDC_deg.down$gene)
DC.shared_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% DC.shared_down, 1, 0)
)
sum(DC.shared_down$fisher) #11
write.csv(DC.shared_down, "../Data/GO_Input/DC.shared_down_fisher.csv", quote = FALSE, row.names = FALSE)

#total DEGs unique to ADC
ADC.only <- setdiff(ADC_deg$gene, SDC_deg$gene)
ADC.only <- ADC_deg[row.names(ADC_deg) %in% ADC.only,]
allGenes <- rownames(dds_host) 
ADC.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(ADC.only), 1, 0)
)
sum(ADC.only$fisher) #295
write.csv(ADC.only, "../Data/GO_Input/ADC.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to ADC
ADC.only_up <- setdiff(ADC_deg.up$gene, SDC_deg$gene)
ADC.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% ADC.only_up, 1, 0)
)
sum(ADC.only_up$fisher) #224
write.csv(ADC.only_up, "../Data/GO_Input/ADC.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to ADC
ADC.only_down <- setdiff(ADC_deg.down$gene, SDC_deg$gene)
ADC.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% ADC.only_down, 1, 0)
)
sum(ADC.only_down$fisher) #71
write.csv(ADC.only_down, "../Data/GO_Input/ADC.only_down_fisher.csv", quote = FALSE, row.names = FALSE)

#DEGs unique to SDC
SDC.only <- setdiff(SDC_deg$gene, ADC_deg$gene)
SDC.only <- SDC_deg[row.names(SDC_deg) %in% SDC.only,]
SDC.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(SDC.only), 1, 0)
)
sum(SDC.only$fisher) #29
write.csv(SDC.only, "../Data/GO_Input/SDC.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to SDC
SDC.only_up <- setdiff(SDC_deg.up$gene, ADC_deg$gene)
SDC.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SDC.only_up, 1, 0)
)
sum(SDC.only_up$fisher) #14
write.csv(SDC.only_up, "../Data/GO_Input/SDC.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to SDC
SDC.only_down <- setdiff(SDC_deg.down$gene, ADC_deg$gene)
SDC.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SDC.only_down, 1, 0)
)
sum(SDC.only_down$fisher) #15
write.csv(SDC.only_down, "../Data/GO_Input/SDC.only_down_fisher.csv", quote = FALSE, row.names = FALSE)

##----------Control Temperature, High Light----------
#overlapping DEGs between AHC & SHC
HC.shared <- intersect(AHC_deg$gene, SHC_deg$gene)
AHC.shared <- AHC_deg[row.names(AHC_deg) %in% HC.shared, ]
SHC.shared <- SHC_deg[row.names(SHC_deg) %in% HC.shared,]
HC.shared <- merge(AHC.shared, SHC.shared, by="row.names", suffixes = c("_AHC", "_SHC"))
#plot to see whether there are antagonistic genes
ggplot(HC.shared, aes(x = log2change_SHC, y = log2change_AHC)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
  xlab("log2 Fold Change (Brown)") +
  ylab("log2 Fold Change (White)") +
  theme_minimal() 
#create input file for GO-MWU
#total DEGs shared between AHC & SHC
allGenes <- rownames(dds_host) 
rownames(HC.shared)=HC.shared$Row.names
HC.shared <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(HC.shared), 1, 0)
)
sum(HC.shared$fisher) #2
write.csv(HC.shared, "../Data/GO_Input/HC.shared_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs shared between AHC & SHC
HC.shared_up <- intersect(AHC_deg.up$gene, SHC_deg.up$gene)
HC.shared_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% HC.shared_up, 1, 0)
)
sum(HC.shared_up$fisher) #2
write.csv(HC.shared_up, "../Data/GO_Input/HC.shared_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs shared between AHC & SHC
HC.shared_down <- intersect(AHC_deg.down$gene, SHC_deg.down$gene)
HC.shared_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% HC.shared_down, 1, 0)
)
sum(HC.shared_down$fisher) #0

#total DEGs unique to AHC
AHC.only <- setdiff(AHC_deg$gene, SHC_deg$gene)
AHC.only <- AHC_deg[row.names(AHC_deg) %in% AHC.only,]
allGenes <- rownames(dds_host) 
AHC.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(AHC.only), 1, 0)
)
sum(AHC.only$fisher) #23
write.csv(AHC.only, "../Data/GO_Input/AHC.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to AHC
AHC.only_up <- setdiff(AHC_deg.up$gene, SHC_deg$gene)
AHC.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% AHC.only_up, 1, 0)
)
sum(AHC.only_up$fisher) #14
write.csv(AHC.only_up, "../Data/GO_Input/AHC.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to AHC
AHC.only_down <- setdiff(AHC_deg.down$gene, SHC_deg$gene)
AHC.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% AHC.only_down, 1, 0)
)
sum(AHC.only_down$fisher) #9
write.csv(AHC.only_down, "../Data/GO_Input/AHC.only_down_fisher.csv", quote = FALSE, row.names = FALSE)

#DEGs unique to SHC
SHC.only <- setdiff(SHC_deg$gene, AHC_deg$gene)
SHC.only <- SHC_deg[row.names(SHC_deg) %in% SHC.only,]
SHC.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(SHC.only), 1, 0)
)
sum(SHC.only$fisher) #22
write.csv(SHC.only, "../Data/GO_Input/SHC.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to SDC
SHC.only_up <- setdiff(SHC_deg.up$gene, AHC_deg$gene)
SHC.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SHC.only_up, 1, 0)
)
sum(SHC.only_up$fisher) #13
write.csv(SHC.only_up, "../Data/GO_Input/SHC.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to SDC
SHC.only_down <- setdiff(SHC_deg.down$gene, AHC_deg$gene)
SHC.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SHC.only_down, 1, 0)
)
sum(SHC.only_down$fisher) #9
write.csv(SHC.only_down, "../Data/GO_Input/SHC.only_down_fisher.csv", quote = FALSE, row.names = FALSE)

##----------Heated, Dark----------
#overlapping DEGs between ADH & SDH
DH.shared <- intersect(ADH_deg$gene, SDH_deg$gene)
ADH.shared <- ADH_deg[row.names(ADH_deg) %in% DH.shared, ]
SDH.shared <- SDH_deg[row.names(SDH_deg) %in% DH.shared,]
DH.shared <- merge(ADH.shared, SDH.shared, by="row.names", suffixes = c("_ADH", "_SDH"))
#plot to see whether there are antagonistic genes
fig_DH.shared <-ggplot(DH.shared, aes(x = log2change_SDH, y = log2change_ADH)) +
  geom_point(alpha = 0.7, shape=21, fill="deeppink4", color="deeppink4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted",linewidth = 0.8) +
  labs(
    x = expression(bold(log[2]~"FC (BDH)")),
    y = expression(bold(log[2]~"FC (WDH)"))
  ) +
  #xlab("log2 Fold Change (Brown)") +
  #ylab("log2 Fold Change (White)") +
  theme_minimal()+
  theme(
    plot.title=element_blank(),
    axis.ticks=element_line(linewidth=1),
    axis.title = element_text(size=14, face = "bold"),
    axis.text=element_text(size=12),
    panel.grid.major = element_line(color = scales::alpha("grey", 0.1)),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
fig_DH.shared
svg("../Figures/fig_DH.shared.svg", width =4.80, height = 4.68)
plot(fig_DH.shared)  
dev.off()
ggsave("../Figures/fig_DH.shared.svg", fig_DH.shared, width = 4.80, height = 4.68, units = "in", device = "svg")

#create input file for GO-MWU
#shared DEGs between ADH & SDH
allGenes <- rownames(dds_host) 
rownames(DH.shared)=DH.shared$Row.names
DH.shared <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(DH.shared), 1, 0)
)
sum(DH.shared$fisher) #351
write.csv(DH.shared, "../Data/GO_Input/DH.shared_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs shared between ADH & SDH
DH.shared_up <- intersect(ADH_deg.up$gene, SDH_deg.up$gene)
DH.shared_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% DH.shared_up, 1, 0)
)
sum(DH.shared_up$fisher) #153
write.csv(DH.shared_up, "../Data/GO_Input/DH.shared_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs shared between ADH & SDH
DH.shared_down <- intersect(ADH_deg.down$gene, SDH_deg.down$gene)
DH.shared_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% DH.shared_down, 1, 0)
)
sum(DH.shared_down$fisher) #198
write.csv(DH.shared_down, "../Data/GO_Input/DH.shared_down_fisher.csv", quote = FALSE, row.names = FALSE)

#DEGs unique to ADH
ADH.only <- setdiff(ADH_deg$gene, SDH_deg$gene)
ADH.only <- ADH_deg[row.names(ADH_deg) %in% ADH.only,]
ADH.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(ADH.only), 1, 0)
)
sum(ADH.only$fisher) #1566
write.csv(ADH.only, "../Data/GO_Input/ADH.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to ADH
ADH.only_up <- setdiff(ADH_deg.up$gene, SDH_deg$gene)
ADH.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% ADH.only_up, 1, 0)
)
sum(ADH.only_up$fisher) #671
write.csv(ADH.only_up, "../Data/GO_Input/ADH.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to ADH
ADH.only_down <- setdiff(ADH_deg.down$gene, SDH_deg$gene)
ADH.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% ADH.only_down, 1, 0)
)
sum(ADH.only_down$fisher) #895
write.csv(ADH.only_down, "../Data/GO_Input/ADH.only_down_fisher.csv", quote = FALSE, row.names = FALSE)

#Total DEGs unique to SDH
SDH.only <- setdiff(SDH_deg$gene, ADH_deg$gene)
SDH.only <- SDH_deg[row.names(SDH_deg) %in% SDH.only,]
SDH.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(SDH.only), 1, 0)
)
sum(SDH.only$fisher) #335
write.csv(SDH.only, "../Data/GO_Input/SDH.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to SDH
SDH.only_up <- setdiff(SDH_deg.up$gene, ADH_deg$gene)
SDH.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SDH.only_up, 1, 0)
)
sum(SDH.only_up$fisher) #163
write.csv(SDH.only_up, "../Data/GO_Input/SDH.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to SHH
SDH.only_down <- setdiff(SDH_deg.down$gene, ADH_deg$gene)
SDH.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SDH.only_down, 1, 0)
)
sum(SDH.only_down$fisher) #172
write.csv(SDH.only_down, "../Data/GO_Input/SDH.only_down_fisher.csv", quote = FALSE, row.names = FALSE)

##----------Heated, Control Light----------
#overlapping DEGs between ACH & SCH
CH.shared <- intersect(ACH_deg$gene, SCH_deg$gene)
ACH.shared <- ACH_deg[row.names(ACH_deg) %in% CH.shared, ]
SCH.shared <- SCH_deg[row.names(SCH_deg) %in% CH.shared,]
CH.shared <- merge(ACH.shared, SCH.shared, by="row.names", suffixes = c("_ACH", "_SCH"))
#plot to see whether there are antagonistic genes
ggplot(CH.shared, aes(x = log2change_SCH, y = log2change_ACH)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
  xlab("log2 Fold Change (Brown)") +
  ylab("log2 Fold Change (White)") +
  theme_minimal() 

#create input file for GO-MWU
allGenes <- rownames(dds_host) 
rownames(CH.shared)=CH.shared$Row.names
CH.shared <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(CH.shared), 1, 0)
)
sum(CH.shared$fisher) #113
write.csv(CH.shared, "../Data/GO_Input/CH.shared_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs shared between AHH & SHH
CH.shared_up <- intersect(ACH_deg.up$gene, SCH_deg.up$gene)
CH.shared_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% CH.shared_up, 1, 0)
)
sum(CH.shared_up$fisher) #65
write.csv(CH.shared_up, "../Data/GO_Input/CH.shared_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs shared between AHH & SHH
CH.shared_down <- intersect(ACH_deg.down$gene, SCH_deg.down$gene)
CH.shared_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% CH.shared_down, 1, 0)
)
sum(CH.shared_down$fisher) #48
write.csv(CH.shared_down, "../Data/GO_Input/CH.shared_down_fisher.csv", quote = FALSE, row.names = FALSE)

#total DEGs unique to ACH
ACH.only <- setdiff(ACH_deg$gene, SCH_deg$gene)
ACH.only <- ACH_deg[row.names(ACH_deg) %in% ACH.only,]
allGenes <- rownames(dds_host) 
ACH.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(ACH.only), 1, 0)
)
sum(ACH.only$fisher) #178
write.csv(ACH.only, "../Data/GO_Input/ACH.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to ACH
ACH.only_up <- setdiff(ACH_deg.up$gene, SCH_deg$gene)
ACH.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% ACH.only_up, 1, 0)
)
sum(ACH.only_up$fisher) #71
write.csv(ACH.only_up, "../Data/GO_Input/ACH.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to AHH
ACH.only_down <- setdiff(ACH_deg.down$gene, SCH_deg$gene)
ACH.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% ACH.only_down, 1, 0)
)
sum(ACH.only_down$fisher) #107
write.csv(ACH.only_down, "../Data/GO_Input/ACH.only_down_fisher.csv", quote = FALSE, row.names = FALSE)

#DEGs unique to SCH
SCH.only <- setdiff(SCH_deg$gene, ACH_deg$gene)
SCH.only <- SCH_deg[row.names(SCH_deg) %in% SCH.only,]
SCH.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(SCH.only), 1, 0)
)
sum(SCH.only$fisher) #320
write.csv(SCH.only, "../Data/GO_Input/SCH.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to SHH
SCH.only_up <- setdiff(SCH_deg.up$gene, ACH_deg$gene)
SCH.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SCH.only_up, 1, 0)
)
sum(SCH.only_up$fisher) #108
write.csv(SCH.only_up, "../Data/GO_Input/SCH.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to SCH
SCH.only_down <- setdiff(SCH_deg.down$gene, ACH_deg$gene)
SCH.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SCH.only_down, 1, 0)
)
sum(SCH.only_down$fisher) #212
write.csv(SCH.only_down, "../Data/GO_Input/SCH.only_down_fisher.csv", quote = FALSE, row.names = FALSE)

##----------Heated, High Light----------
#overlapping DEGs between AHH & SHH
HH.shared <- intersect(AHH_deg$gene, SHH_deg$gene)
AHH.shared <- AHH_deg[row.names(AHH_deg) %in% HH.shared, ]
SHH.shared <- SHH_deg[row.names(SHH_deg) %in% HH.shared,]
HH.shared <- merge(AHH.shared, SHH.shared, by="row.names", suffixes = c("_AHH", "_SHH"))
#plot to see whether there are antagonistic genes
fig_HH.shared <- ggplot(HH.shared, aes(x = log2change_SHH, y = log2change_AHH)) +
  geom_point(alpha = 0.7, shape=21, fill="rosybrown1", color="rosybrown4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted",linewidth = 0.8) +
  labs(
    x = expression(bold(log[2]~"FC (BHH)")),
    y = expression(bold(log[2]~"FC (WHH)"))
  ) +
  #xlab("log2 Fold Change (Brown)") +
  #ylab("log2 Fold Change (White)") +
  theme_minimal()+
  theme(
    plot.title=element_blank(),
    axis.ticks=element_line(linewidth=1),
    axis.title = element_text(size=14, face = "bold"),
    axis.text=element_text(size=12),
    panel.grid.major = element_line(color = scales::alpha("grey", 0.1)),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
fig_HH.shared
svg("../Figures/fig_HH.shared.svg", width =4.80, height = 4.68)
plot(fig_HH.shared)  
dev.off()
ggsave("../Figures/fig_HH.shared.svg", fig_HH.shared, width = 4.80, height = 4.68, units = "in", device = "svg")

HH.shared_diff <- HH.shared %>%
  dplyr::mutate(
    lfc_diff = log2change_AHH - log2change_SHH,
    abs_lfc_diff = abs(lfc_diff)
  ) %>%
  dplyr::filter(abs_lfc_diff > 1) %>%
  dplyr::arrange(desc(abs_lfc_diff))
more_white <- HH.shared_diff %>%
  dplyr::filter(lfc_diff > 0)
more_brown <- HH.shared_diff %>%
  dplyr::filter(lfc_diff < 0)

#create input file for GO-MWU
#total DEGs shared between AHH & SHH
allGenes <- rownames(dds_host) 
rownames(HH.shared)=HH.shared$Row.names
HH.shared <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(HH.shared), 1, 0)
)
sum(HH.shared$fisher) #602
write.csv(HH.shared, "../Data/GO_Input/HH.shared_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs shared between AHH & SHH
HH.shared_up <- intersect(AHH_deg.up$gene, SHH_deg.up$gene)
HH.shared_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% HH.shared_up, 1, 0)
)
sum(HH.shared_up$fisher) #268
write.csv(HH.shared_up, "../Data/GO_Input/HH.shared_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs shared between AHH & SHH
HH.shared_down <- intersect(AHH_deg.down$gene, SHH_deg.down$gene)
HH.shared_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% HH.shared_down, 1, 0)
)
sum(HH.shared_down$fisher) #333
write.csv(HH.shared_down, "../Data/GO_Input/HH.shared_down_fisher.csv", quote = FALSE, row.names = FALSE)

#total DEGs unique to AHH
AHH.only <- setdiff(AHH_deg$gene, SHH_deg$gene)
AHH.only <- AHH_deg[row.names(AHH_deg) %in% AHH.only,]
allGenes <- rownames(dds_host) 
AHH.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(AHH.only), 1, 0)
)
sum(AHH.only$fisher) #576
write.csv(AHH.only, "../Data/GO_Input/AHH.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to AHH
AHH.only_up <- setdiff(AHH_deg.up$gene, SHH_deg$gene)
AHH.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% AHH.only_up, 1, 0)
)
sum(AHH.only_up$fisher) #344
write.csv(AHH.only_up, "../Data/GO_Input/AHH.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to AHH
AHH.only_down <- setdiff(AHH_deg.down$gene, SHH_deg$gene)
AHH.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% AHH.only_down, 1, 0)
)
sum(AHH.only_down$fisher) #232
write.csv(AHH.only_down, "../Data/GO_Input/AHH.only_down_fisher.csv", quote = FALSE, row.names = FALSE)

#DEGs unique to SHH
SHH.only <- setdiff(SHH_deg$gene, AHH_deg$gene)
SHH.only <- SHH_deg[row.names(SHH_deg) %in% SHH.only,]
SHH.only <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% rownames(SHH.only), 1, 0)
)
sum(SHH.only$fisher) #1200
write.csv(SHH.only, "../Data/GO_Input/SHH.only_fisher.csv", quote = FALSE, row.names = FALSE)
#up-regulated DEGs unique to SHH
SHH.only_up <- setdiff(SHH_deg.up$gene, AHH_deg$gene)
SHH.only_up <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SHH.only_up, 1, 0)
)
sum(SHH.only_up$fisher) #524
write.csv(SHH.only_up, "../Data/GO_Input/SHH.only_up_fisher.csv", quote = FALSE, row.names = FALSE)
#down-regulated DEGs unique to SHH
SHH.only_down <- setdiff(SHH_deg.down$gene, AHH_deg$gene)
SHH.only_down <- data.frame(
  gene = allGenes,
  fisher = ifelse(allGenes %in% SHH.only_down, 1, 0)
)
sum(SHH.only_down$fisher) #676
write.csv(SHH.only_down, "../Data/GO_Input/SHH.only_down_fisher.csv", quote = FALSE, row.names = FALSE)
