#Astrangia poculata heat & light: differentially expressed genes (DEGs)
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to extract DEGs for downstream analyses
#Note: In this script: apo = white corals, sym = brown corals


library(DESeq2)
library(dplyr)
library(ggplot2)
library(grid)

load("dds.Rdata")


#----------Baseline comparison between symbiotic states----------
res_Sym <- results (dds_host.symstat_control, contrast=c("symstat", "Sym", "Apo"))
head(res_Sym)
summary(res_Sym)
#LFC > 0 (up)       : 70, 0.19%
#LFC < 0 (down)     : 96, 0.26%
valState_Sym=cbind(res_Sym$log2FoldChange, res_Sym$pvalue, res_Sym$padj)
head(valState_Sym)
colnames(valState_Sym)=c("log2change", "pval", "padj")
rownames(valState_Sym)<-rownames(res_Sym)
head(valState_Sym)

#prepare input files for GO-MWU analysis
GO_Sym<-mutate(as.data.frame(res_Sym),negP=-log(pvalue))
GO_Sym<-mutate(GO_Sym,signedlogP=case_when(log2FoldChange<0 ~ negP*-1, log2FoldChange>0 ~ negP))
GO_Sym$genes<-rownames(res_Sym)
GO_Sym_signedP<-dplyr::select(GO_Sym, genes, signedlogP)
write.csv(GO_Sym_signedP, file="../Data/GO_Input/GO_Sym_signedP.csv", row.names = F, quote = F)
GO_Sym_lfc<-dplyr::select(GO_Sym, genes,log2FoldChange)
write.csv(GO_Sym_lfc, file="../Data/GO_Input/GO_Sym_lfc.csv", row.names = F, quote = F)

res_Sym_control <- results (dds_host.symstat_control, contrast=c("symstat", "Sym", "Apo"))
head(res_Sym_control)
summary(res_Sym_control)
#LFC > 0 (up)       : 7, 0.022%
#LFC < 0 (down)     : 2, 0.0062%
valState_Sym_control=cbind(res_Sym_control$log2FoldChange, res_Sym_control$pvalue, res_Sym_control$padj)
head(valState_Sym_control)
colnames(valState_Sym_control)=c("log2change", "pval", "padj")
rownames(valState_Sym_control)<-rownames(res_Sym_control)
head(valState_Sym_control)

#prepare input files for GO-MWU analysis
GO_Sym_control<-mutate(as.data.frame(res_Sym_control),negP=-log(pvalue))
GO_Sym_control<-mutate(GO_Sym_control,signedlogP=case_when(log2FoldChange<0 ~ negP*-1, log2FoldChange>0 ~ negP))
GO_Sym_control$genes<-rownames(res_Sym_control)
GO_Sym_control_signedP<-dplyr::select(GO_Sym_control, genes, signedlogP)
write.csv(GO_Sym_control_signedP, file="../Data/GO_Input/GO_Sym_control_signedP.csv", row.names = F, quote = F)
GO_Sym_control_lfc<-dplyr::select(GO_Sym_control, genes,log2FoldChange)
write.csv(GO_Sym_control_lfc, file="../Data/GO_Input/GO_Sym_control_lfc.csv", row.names = F, quote = F)


#----------For subsequent comparative analysis of thermal responses----------
res_CH <- results(dds_host.treat, contrast=c("treat", "CL_H", "CL_C"))
head(res_CH)
summary(res_CH)
#LFC > 0 (up)       : 555, 1.5%
#LFC < 0 (down)     : 685, 1.8%
valState_CH=cbind(res_CH$log2FoldChange, res_CH$pvalue, res_CH$padj)
head(valState_CH)
colnames(valState_CH)=c("log2change", "pval", "padj")
rownames(valState_CH)<-rownames(res_CH)
head(valState_CH)

#prepare input files for GO-MWU analysis
GO_CH<-mutate(as.data.frame(res_CH),negP=-log(pvalue))
GO_CH<-mutate(GO_CH,signedlogP=case_when(log2FoldChange<0 ~ negP*-1, log2FoldChange>0 ~ negP))
GO_CH$genes<-rownames(GO_CH)
GO_CH_signedP<-dplyr::select(GO_CH, genes, signedlogP)
write.csv(GO_CH_signedP, file="../Data/GO_Input/GO_CH_signedP.csv", row.names = F, quote = F)
GO_CH_lfc<-dplyr::select(GO_CH, genes,log2FoldChange)
write.csv(GO_CH_lfc, file="../Data/GO_Input/GO_CH_lfc.csv", row.names = F, quote = F)

#----------For subsequent exploration of responses to constant darkness----------
res_DC <- results(dds_host.treat, contrast=c("treat", "NL_C", "CL_C"))
head(res_DC)
summary(res_DC)
#LFC > 0 (up)       : 249, 0.67%
#LFC < 0 (down)     : 132, 0.35%
valState_DC=cbind(res_DC$log2FoldChange, res_DC$pvalue, res_DC$padj)
head(valState_DC)
colnames(valState_DC)=c("log2change", "pval", "padj")
rownames(valState_DC)<-rownames(res_DC)
head(valState_DC)

#prepare input files for GO-MWU analysis
GO_DC<-mutate(as.data.frame(res_DC),negP=-log(pvalue))
GO_DC<-mutate(GO_DC,signedlogP=case_when(log2FoldDCange<0 ~ negP*-1, log2FoldDCange>0 ~ negP))
GO_DC$genes<-rownames(GO_DC)
GO_DC_signedP<-dplyr::select(GO_DC, genes, signedlogP)
write.csv(GO_DC_signedP, file="../Data/GO_Input/GO_DC_signedP.csv", row.names = F, quote = F)
GO_DC_lfc<-dplyr::select(GO_DC, genes,log2FoldDCange)
write.csv(GO_DC_lfc, file="../Data/GO_Input/GO_DC_lfc.csv", row.names = F, quote = F)


#----------Treatment comparisons within each symbiotic phenotype----------

##----------White corals----------
res_ADC <- results(dds_host.apo, contrast=c("treat", "NL_C", "CL_C"))
head(res_ADC)
summary(res_ADC)
#LFC > 0 (up)       : 237, 0.67%
#LFC < 0 (down)     : 82, 0.23%
valState_ADC=cbind(res_ADC$log2FoldChange, res_ADC$pvalue, res_ADC$padj)
head(valState_ADC)
colnames(valState_ADC)=c("log2change", "pval", "padj")
rownames(valState_ADC)<-rownames(res_ADC)
head(valState_ADC)

res_AHC <- results(dds_host.apo, contrast=c("treat", "HL_C", "CL_C"))
head(res_AHC)
summary(res_AHC)
#LFC > 0 (up)       : 16, 0.045%
#LFC < 0 (down)     : 9, 0.025%
valState_AHC=cbind(res_AHC$log2FoldChange, res_AHC$pvalue, res_AHC$padj)
head(valState_AHC)
colnames(valState_AHC)=c("log2change", "pval", "padj")
rownames(valState_AHC)<-rownames(res_AHC)
head(valState_AHC) 

res_ADH <- results(dds_host.apo, contrast=c("treat", "NL_H", "CL_C"))
head(res_ADH)
summary(res_ADH)
#LFC > 0 (up)       : 824, 2.3%
#LFC < 0 (down)     : 1093, 3.1%
valState_ADH=cbind(res_ADH$log2FoldChange, res_ADH$pvalue, res_ADH$padj)
head(valState_ADH)
colnames(valState_ADH)=c("log2change", "pval", "padj")
rownames(valState_ADH)<-rownames(res_ADH)
head(valState_ADH)

res_ACH <- results(dds_host.apo, contrast=c("treat", "CL_H", "CL_C"))
head(res_ACH)
summary(res_ACH)
#LFC > 0 (up)       : 136, 0.39%
#LFC < 0 (down)     : 155, 0.44%
valState_ACH=cbind(res_ACH$log2FoldChange, res_ACH$pvalue, res_ACH$padj)
head(valState_ACH)
colnames(valState_ACH)=c("log2change", "pval", "padj")
rownames(valState_ACH)<-rownames(res_ACH)
head(valState_ACH)
#prepare input file for GO-MWU
GO_ACH<-mutate(as.data.frame(res_ACH),negP=-log(pvalue))
GO_ACH<-mutate(GO_ACH,signedlogP=case_when(log2FoldChange<0 ~ negP*-1, log2FoldChange>0 ~ negP))
GO_ACH$genes<-rownames(GO_ACH)
GO_ACH_signedP<-dplyr::select(GO_ACH, genes, signedlogP)
write.csv(GO_ACH_signedP, file="../Data/GO_Input/GO_ACH_signedP.csv", row.names = F, quote = F)

res_AHH <- results(dds_host.apo, contrast=c("treat", "HL_H", "CL_C"))
head(res_AHH)
summary(res_AHH)
#LFC > 0 (up)       : 613, 1.7%
#LFC < 0 (down)     : 565, 1.6%
valState_AHH=cbind(res_AHH$log2FoldChange, res_AHH$pvalue, res_AHH$padj)
head(valState_AHH)
colnames(valState_AHH)=c("log2change", "pval", "padj")
rownames(valState_AHH)<-rownames(res_AHH)
head(valState_AHH)

#make data frame for white DEGs
host_degs.apo = data.frame(
  sample_id = c("ADC", "AHC", "ADH", "ACH", "AHH"),
  symstat = factor(c("Apo", "Apo", "Apo", "Apo", "Apo")),
  treat = factor(c("NL_C", "HL_C", "NL_H", "CL_H", "HL_H")),
  deg_up = c(237, 16, 824, 136, 613),
  deg_down = c(-82, -9, -1093, -155, -565)
)

##----------Brown corals----------
res_SDC <- results(dds_host.sym, contrast=c("treat", "NL_C", "CL_C"))
head(res_SDC)
summary(res_SDC)
#LFC > 0 (up)       : 27, 0.077%
#LFC < 0 (down)     : 26, 0.074%
valState_SDC=cbind(res_SDC$log2FoldChange, res_SDC$pvalue, res_SDC$padj)
head(valState_SDC)
colnames(valState_SDC)=c("log2change", "pval", "padj")
rownames(valState_SDC)<-rownames(res_SDC)
head(valState_SDC)

res_SHC <- results(dds_host.sym, contrast=c("treat", "HL_C", "CL_C"))
head(res_SHC)
summary(res_SHC)
#LFC > 0 (up)       : 15, 0.043%
#LFC < 0 (down)     : 9, 0.026%
valState_SHC=cbind(res_SHC$log2FoldChange, res_SHC$pvalue, res_SHC$padj)
head(valState_SHC)
colnames(valState_SHC)=c("log2change", "pval", "padj")
rownames(valState_SHC)<-rownames(res_SHC)
head(valState_SHC) 

res_SDH <- results(dds_host.sym, contrast=c("treat", "NL_H", "CL_C"))
head(res_SDH)
summary(res_SDH)
#LFC > 0 (up)       : 316, 0.9%
#LFC < 0 (down)     : 370, 1.1%
valState_SDH=cbind(res_SDH$log2FoldChange, res_SDH$pvalue, res_SDH$padj)
head(valState_SDH)
colnames(valState_SDH)=c("log2change", "pval", "padj")
rownames(valState_SDH)<-rownames(res_SDH)
head(valState_SDH)

res_SCH <- results(dds_host.sym, contrast=c("treat", "CL_H", "CL_C"))
head(res_SCH)
summary(res_SCH)
#LFC > 0 (up)       : 173, 0.49%
#LFC < 0 (down)     : 260, 0.74%
valState_SCH=cbind(res_SCH$log2FoldChange, res_SCH$pvalue, res_SCH$padj)
head(valState_SCH)
colnames(valState_SCH)=c("log2change", "pval", "padj")
rownames(valState_SCH)<-rownames(res_SCH)
head(valState_SCH)
#prepare input file for GO-MWU
GO_SCH<-mutate(as.data.frame(res_SCH),negP=-log(pvalue))
GO_SCH<-mutate(GO_SCH,signedlogP=case_when(log2FoldChange<0 ~ negP*-1, log2FoldChange>0 ~ negP))
GO_SCH$genes<-rownames(GO_SCH)
GO_SCH_signedP<-dplyr::select(GO_SCH, genes, signedlogP)
write.csv(GO_SCH_signedP, file="../Data/GO_Input/GO_SCH_signedP.csv", row.names = F, quote = F)

res_SHH <- results(dds_host.sym, contrast=c("treat", "HL_H", "CL_C"))
head(res_SHH)
summary(res_SHH)
#LFC > 0 (up)       : 792, 2.3%
#LFC < 0 (down)     : 1010, 2.9%
valState_SHH=cbind(res_SHH$log2FoldChange, res_SHH$pvalue, res_SHH$padj)
head(valState_SHH)
colnames(valState_SHH)=c("log2change", "pval", "padj")
rownames(valState_SHH)<-rownames(res_SHH)
head(valState_SHH)

#make data frame for brown DEGs
host_degs.sym = data.frame(
  sample_id = c("SDC", "SHC", "SDH", "SCH", "SHH"),
  symstat = factor(c("Sym", "Sym", "Sym", "Sym", "Sym")),
  treat = factor(c("NL_C", "HL_C", "NL_H", "CL_H", "HL_H")),
  deg_up = c(27, 15, 316, 173, 792),
  deg_down = c(-26, -9, -370, -260, -1010)
)

save(res_ADC, res_AHC, res_ADH, res_ACH, res_AHH,
     res_SDC, res_AHC, res_SDH, res_SCH, res_SHH,
     valState_ADC, valState_AHC, valState_ADH, valState_ACH, valState_AHH,
     valState_SDC, valState_SHC, valState_SDH, valState_SCH, valState_SHH,
     file="valState.RData")

#plot DEGs 
host_degs.all <- rbind(host_degs.apo, host_degs.sym)
host_degs.all$symstat <- recode(host_degs.all$symstat,
                                Apo = "White",
                                Sym = "Brown")
host_degs.all$treat <- factor(host_degs.all$treat,
                              levels = c("NL_C", "HL_C", "NL_H", "CL_H", "HL_H"),
                              labels = c("Control Temp, Dark",
                                         "Control Temp, High Light",
                                         "Heated, Dark",
                                         "Heated, Control Light",
                                         "Heated, High Light"))
host_degs.all$treat_code <- factor(c(rep(c("NL_C", "HL_C", "NL_H", "CL_H", "HL_H"), 2)),
                                   levels = c("NL_C", "HL_C", "NL_H", "CL_H", "HL_H"))
fig_degs <- ggplot(host_degs.all, aes(x = treat)) +
  geom_bar(aes(y = deg_up, fill = treat_code), stat = "identity", alpha = 1) +
  geom_bar(aes(y = deg_down, fill = treat_code), stat = "identity", alpha = 0.4) +
  facet_wrap(~ symstat) +
  scale_fill_manual(
    name = "Treatment",
    limits = c("NL_C", "HL_C", "NL_H", "CL_H", "HL_H"),
    values = c(
      "NL_C" = "navy",
      "HL_C" = "steelblue1",
      "NL_H" = "deeppink4",
      "CL_H" = "palevioletred1",
      "HL_H" = "rosybrown1"
    ),
    labels = c(
      "Control Temp, Dark",
      "Control Temp, High Light",
      "Heated, Dark",
      "Heated, Control Light",
      "Heated, High Light"
    )
  ) +
  labs(x = "Treatment", y = "Number of DEGs") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.ticks.length = unit(5, "pt"),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(color = "black")
  )
fig_degs
ggsave("../Supplementary/DEGs.svg", fig_degs, width = 7.58, height = 4.68, units = "in", device = "svg")


#----------Vst normalized count file with pvals----------
#get pvals
val_ADC=cbind(res_ADC$pvalue, res_ADC$padj)
head(val_ADC)
colnames(val_ADC)=c("pval.ADC", "padj.ADC")
length(val_ADC[,1])
table(complete.cases(val_ADC))

val_AHC=cbind(res_AHC$pvalue, res_AHC$padj)
head(val_AHC)
colnames(val_AHC)=c("pval.AHC", "padj.AHC")
length(val_AHC[,1])
table(complete.cases(val_AHC))

val_ADH=cbind(res_ADH$pvalue, res_ADH$padj)
head(val_ADH)
colnames(val_ADH)=c("pval.ADH", "padj.ADH")
length(val_ADH[,1])
table(complete.cases(val_ADH))

val_ACH=cbind(res_ACH$pvalue, res_ACH$padj)
head(val_ACH)
colnames(val_ACH)=c("pval.ACH", "padj.ACH")
length(val_ACH[,1])
table(complete.cases(val_ACH))

val_AHH=cbind(res_AHH$pvalue, res_AHH$padj)
head(val_AHH)
colnames(val_AHH)=c("pval.AHH", "padj.AHH")
length(val_AHH[,1])
table(complete.cases(val_AHH))

val_SDC=cbind(res_SDC$pvalue, res_SDC$padj)
head(val_SDC)
colnames(val_SDC)=c("pval.SDC", "padj.SDC")
length(val_SDC[,1])
table(complete.cases(val_SDC))

val_SHC=cbind(res_SHC$pvalue, res_SHC$padj)
head(val_SHC)
colnames(val_SHC)=c("pval.SHC", "padj.SHC")
length(val_SHC[,1])
table(complete.cases(val_SHC))

val_SDH=cbind(res_SDH$pvalue, res_SDH$padj)
head(val_SDH)
colnames(val_SDH)=c("pval.SDH", "padj.SDH")
length(val_SDH[,1])
table(complete.cases(val_SDH))

val_SCH=cbind(res_SCH$pvalue, res_SCH$padj)
head(val_SCH)
colnames(val_SCH)=c("pval.SCH", "padj.SCH")
length(val_SCH[,1])
table(complete.cases(val_SCH))

val_SHH=cbind(res_SHH$pvalue, res_SHH$padj)
head(val_SHH)
colnames(val_SHH)=c("pval.SHH", "padj.SHH")
length(val_SHH[,1])
table(complete.cases(val_SHH))

vst_host = vst(dds_host, blind = TRUE)
vst=assay(vst_host)
head(vst)
length(vst[,1]) #45670

vst.apo <- vst[,1:39]
head(vst.apo)
vst.sym <- vst[,40:77]
head(vst.sym)

vstpvals.apo = cbind(vst.apo, val_ADC, val_AHC, val_ADH, val_ACH, val_AHH)
head(vstpvals.apo)
dim(vstpvals.apo) #45670 49
table(complete.cases(vstpvals.apo))
#FALSE  TRUE 
#38678  6992 
write.csv(vstpvals.apo, "../Data/VSTandPVALS_apo.csv", quote=F)

vstpvals.sym = cbind(vst.sym, val_SDC, val_SHC, val_SDH, val_SCH, val_SHH)
head(vstpvals.sym)
dim(vstpvals.sym) #45670 48
table(complete.cases(vstpvals.sym))
#FALSE  TRUE 
#42716  2954 
write.csv(vstpvals.sym, "../Data/VSTandPVALS_sym.csv", quote=F)
