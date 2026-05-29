#Astrangia poculata heat & light: GE PCA & plasticity
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Fig. 2C & 2D in the manuscript
#Note: In this script: apo = white corals, sym = brown corals


library(DESeq2)
library(ggplot2)
library(dplyr)
library(vegan)
library(performance)
library(scales)


#----------GE PCA----------
load("dds.RData")

vst_host = vst(dds_host, blind = TRUE) #transform data
pcaData_host = plotPCA(vst_host, intgroup = c("symstat", "treat","heat","light"), returnData = TRUE)
percentVar_host <- 100 * round(attr(pcaData_host, "percentVar"), digits = 4)
pca_host <- stats::prcomp(t(assay(vst_host)), center = TRUE)

#statistics
host_pca.full <- vegan::adonis2(pca_host$x ~ vst_host$symstat
                                * vst_host$heat
                                * vst_host$light,
                                by = "terms", method = 'eu')
host_pca.full
#vst_host$symstat                               1    18460 0.01732 1.4024  0.008 ** 
#vst_host$heat                                  1    50697 0.04758 3.8513  0.001 ***
#vst_host$light                                 2    49304 0.04627 1.8728  0.001 ***
#vst_host$symstat:vst_host$heat                 1    11399 0.01070 0.8659  0.943    
#vst_host$symstat:vst_host$light                2    34262 0.03215 1.3014  0.001 ***
#vst_host$heat:vst_host$light                   2    23284 0.02185 0.8844  0.962    
#vst_host$symstat:vst_host$heat:vst_host$light  2    22584 0.02119 0.8579  0.996 

pca_pval.full <- host_pca.full[["Pr(>F)"]]
symbiosis_pca_pval.full <-substitute(italic(P[SP])==p, list(p = format(pca_pval.full[1], digits = 4)))
heat_pca_pval.full <-substitute(italic(P[H])==p, list(p = format(pca_pval.full[2], digits = 4)))
light_pca_pval.full <-substitute(italic(P[L])==p, list(p = format(pca_pval.full[3], digits = 4)))
sxh_pca_pval.full <-substitute(italic(P[SP~x~L])==p, list(p = format(pca_pval.full[4], digits = 4)))
sxl_pca_pval.full <-substitute(italic(P[SP~x~L])==p, list(p = format(pca_pval.full[5], digits = 4)))
hxl_pca_pval.full <-substitute(italic(P[H~x~L])==p, list(p = format(pca_pval.full[6], digits = 4)))
sxhxl_pca_pval.full <-substitute(italic(P[SP~x~H~x~L])==p, list(p = format(pca_pval.full[7], digits = 4)))

#plot GE PCA
shapes = c("Apo" = 1, "Sym" = 19)
GE.PCA_host <- ggplot(pcaData_host, aes(PC1, PC2)) + 
  geom_point(aes(colour = treat, shape = symstat), size = 3, stroke = 1.25) +
  stat_ellipse(geom = "polygon", level=0.95, alpha = 0.25, aes(fill = treat)) +
  scale_colour_manual(limits=c("NL_C","CL_C", "HL_C", "NL_H", "CL_H",  "HL_H"),
                      values=c("navy", "royalblue","steelblue1","deeppink4",  "palevioletred1",  "rosybrown1"),
                      labels=c("Control Temp, Dark","Control Temp, Control Light", "Control Temp, High Light","Heated, Dark", "Heated, Control Light",  "Heated, High Light"),
                      name="Combined Treatment") +
  scale_fill_manual(limits=c("NL_C","CL_C", "HL_C", "NL_H", "CL_H",  "HL_H"),
                    values=c("navy", "royalblue","steelblue1","deeppink4",  "palevioletred1",  "rosybrown1"),
                    labels=c("Control Temp, Dark","Control Temp, Control Light, ", "Control Temp, High Light","Heated, Dark", "Heated, Control Light",  "Heated, High Light"),
                    guide="none") +
  scale_shape_manual(values = shapes, breaks = c("Apo", "Sym"), labels = c("White", "Brown"), name = "Symbitoic Phenotype") +
  annotate("text", x = -25, y = 30, label = deparse(symbiosis_pca_pval.full), parse = TRUE, size = 4, hjust=0) +
  annotate("text", x = -25, y = 25, label = deparse(heat_pca_pval.full), parse = TRUE, size = 4, hjust=0) +
  annotate("text", x = -25, y = 20, label = deparse(light_pca_pval.full), parse = TRUE, size = 4, hjust=0) +
  annotate("text", x = -25, y = 15, label = deparse(sxl_pca_pval.full), parse = TRUE, size = 4, hjust=0) +
  xlab(paste0("PC1 (",percentVar_host[1],"%)")) + 
  ylab(paste0("PC2 (",percentVar_host[2],"%)")) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme_classic()+
  theme(
    plot.title=element_blank(),
    axis.ticks=element_line(linewidth=1),
    axis.title = element_text(size=14, face = "bold"),
    axis.text=element_text(size=12),
    legend.title = element_text(size=14, face = "bold"),
    legend.text = element_text(size=12),
    panel.grid.major = element_line(color = scales::alpha("grey", 0.1)),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
GE.PCA_host
ggsave("../Figures/GE.PCA_host.svg", GE.PCA_host, width = 7.58, height = 4.68, units = "in", device = "svg")


#----------GE Plasticity----------
source("PlasticityCustomFunctions_means.R")

#get all PCs using this function
plotPCA_allPCs <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE)
{
  rv <- rowVars(assay(object)) # Variance estimates for each row (column) in a matrix.
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))] # orders genes from largest -> smallest and pulls ntop (default 500) genes
  pca <- prcomp(t(assay(object)[select, ]), center = TRUE, scale. = FALSE) # performs PCA on ntop (default 500) from dataset
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(group = group, intgroup.df, name = colnames(object), pca$x)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *                                                                                                           100), "% variance")) + coord_fixed()
}

pcadata.plast = plotPCA_allPCs(vst_host, intgroup = c("stat_treat","symstat","treat","heat","light","genet", "tank"), returnData = TRUE)
pcadata.plast_sym = pcadata.plast %>%
  dplyr::filter(symstat == "Sym")
pcadata.plast_apo = pcadata.plast %>%
  dplyr::filter(symstat == "Apo")

pca.plast_sym <-  PCAplast(pca = pcadata.plast_sym[,c(10:ncol(pcadata.plast_sym))], # PCA dataframe containing the PCA eigenvalues
                           data = pcadata.plast_sym[,c(1:9)], # condition/treatment data 
                           sample_ID = "name", # name of column with unique ID per sample (if blank, will pull rownames for this)
                           num_pca = "2", # number of PCs to include in measure (default is 'all' if left blank)
                           control_col = "treat", # name of 'treatment' column
                           control_lvl = "CL_C") # control level of the treatment 
pca.plast_apo <-  PCAplast(pca = pcadata.plast_apo[,c(10:ncol(pcadata.plast_apo))], # PCA dataframe containing the PCA eigenvalues
                           data = pcadata.plast_apo[,c(1:9)], # condition/treatment data 
                           sample_ID = "name", # name of column with unique ID per sample (if blank, will pull rownames for this)
                           num_pca = "2", # number of PCs to include in measure (default is 'all' if left blank)
                           control_col = "treat", # name of 'treatment' column
                           control_lvl = "CL_C") # control level of the treatment 

# re-combine host and symbiont plasticity for plotting
all_plast = rbind(pca.plast_sym, pca.plast_apo) 
str(all_plast)
pca.plast_apo$treat <- factor(pca.plast_apo$treat,
                             levels = c("NL_C", "CL_C", "HL_C", "NL_H", "CL_H", "HL_H"))
pca.plast_sym$treat <- factor(pca.plast_sym$treat,
                              levels = c("NL_C", "CL_C", "HL_C", "NL_H", "CL_H", "HL_H"))

pd <- position_dodge(width = 0.8)
all_plast$light <- factor(
  all_plast$light,
  levels = c("Dark", "Control", "High")
)

fig_GE.Plast <- ggplot(
  all_plast,
  aes(x = light, y = dist,
      color = interaction(heat, light))
) +
  geom_boxplot(
    aes(fill = interaction(heat, light)),
    alpha = 0.6,
    outlier.colour = NA,
    linewidth = 1,
    position = pd
  ) +
  geom_point(
    aes(shape = symstat),
    size = 3,
    stroke = 1,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    )
  ) +
  facet_wrap(~ symstat, labeller = facet_labels) +
  ylab("Gene Expression Plasticity") +
  xlab("Light Treatment") +
  scale_color_manual(
    values = c(
      "Control.Control" = "royalblue",
      "Control.High" = "steelblue1",
      "Control.Dark" = "navy",
      "Heated.Control" = "palevioletred1",
      "Heated.High" = "rosybrown1",
      "Heated.Dark" = "deeppink4"
    )
  ) +scale_fill_manual(
    values = c(
      "Control.Control" = "royalblue",
      "Control.High" = "steelblue1",
      "Control.Dark" = "navy",
      "Heated.Control" = "palevioletred1",
      "Heated.High" = "rosybrown1",
      "Heated.Dark" = "deeppink4"
    )
  ) +
  scale_shape_manual(
    values = c("Apo" = 1, "Sym" = 19),
    labels = c("White", "Brown")
  )+
  scale_y_continuous(limits = c(0, 50)) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_line(linewidth = 1),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    legend.position = "none",
    strip.background = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.text = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(color = "gray90")
  )
fig_GE.Plast
ggsave("../Figures/fig_GE_Plast.svg", fig_GE.Plast, width = 7.58, height = 4.68, units = "in", device = "svg")

#statistics
lm_GE_Plast <- lmer(dist ~ symstat+(1|genet)+(1|tank), data = all_plast)
summary(lm_GE_Plast)
summary(glht(lm_GE_Plast, linfct = mcp(symstat = "Tukey")))  #N.S.
mod1 <- lmer(dist ~ heat*light+(1|genet)+(1|tank), data = pca.plast_apo, REML=FALSE)
mod2 <- lmer(dist ~ heat+light+(1|genet)+(1|tank), data = pca.plast_apo, REML=FALSE)
AICc(mod1, mod2) 
#df     AICc
#mod1  9 293.6474
#mod2  7 290.7550
lm_GE_Plast.apo <- lmer(dist ~ heat+light+(1|genet)+(1|tank), data = pca.plast_apo)
summary(lm_GE_Plast.apo)
summary(glht(lm_GE_Plast.apo, linfct = mcp(heat = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)    
#Heated - Control == 0   12.868      2.625   4.903 9.45e-07 ***
summary(glht(lm_GE_Plast.apo, linfct = mcp(light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)    
#Dark - Control == 0   14.266      3.760   3.794  0.00042 ***
#High - Control == 0    6.639      3.618   1.835  0.15805    
#High - Dark == 0      -7.627      3.879  -1.966  0.12046   
lm_GE_Plast.apo_treat <- lmer(dist ~ treat+(1|genet)+(1|tank), data = pca.plast_apo)
summary(lm_GE_Plast.apo_treat)
summary(glht(lm_GE_Plast.apo_treat, linfct = mcp(treat = "Tukey")))
posthoc.apo_treat <- glht(lm_GE_Plast.apo_treat,
                          linfct = mcp(treat = "Tukey"))
cld(posthoc.apo_treat)



mod1 <- lmer(dist ~ heat*light+(1|genet)+(1|tank), data = pca.plast_sym, REML=FALSE)
mod2 <- lmer(dist ~ heat+light+(1|genet)+(1|tank), data = pca.plast_sym, REML=FALSE)
AICc(mod1, mod2) 
#df     AICc
#mod1  9 281.5921
#mod2  7 276.7339
lm_GE_Plast.sym <- lmer(dist ~ heat+light+(1|genet)+(1|tank), data = pca.plast_sym)
summary(lm_GE_Plast.sym)
summary(glht(lm_GE_Plast.sym, linfct = mcp(heat = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)    
#Heated - Control == 0   14.338      2.961   4.842 1.29e-06 ***
summary(glht(lm_GE_Plast.sym, linfct = mcp(light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)
#Dark - Control == 0   5.2851     4.0588   1.302    0.394
#High - Control == 0   4.4139     4.0588   1.087    0.522
#High - Dark == 0     -0.8713     4.1500  -0.210    0.976

lm_GE_Plast.sym_treat <- lmer(dist ~ treat+(1|genet)+(1|tank), data = pca.plast_sym)
summary(lm_GE_Plast.sym_treat)
summary(glht(lm_GE_Plast.sym_treat, linfct = mcp(treat = "Tukey")))
posthoc.sym_treat <- glht(lm_GE_Plast.sym_treat,
                          linfct = mcp(treat = "Tukey"))
cld(posthoc.sym_treat)
