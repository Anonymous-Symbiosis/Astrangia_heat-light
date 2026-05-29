#Astrangia poculata heat & light: PCA and plasticity of physiological traits
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Fig. 2A & 2B in the manuscript
#Note: In this script: apo = white corals, sym = brown corals

load("Phys_df_PCA.Rdata")

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(vegan)
library(ggbiplot)
library(FactoMineR)
library(factoextra)
library(performance)
library(lme4)
library(multcomp)
library(MuMIn)


df_PCA<- df_PCA %>%
  dplyr::select("Sample_ID","Tank", "Genet","Treatment","State","Heat","Light","Stat_treat","Sym_cell_per_area", "Chla_per_area", "Chlc2_per_area", "Carb_host_per_area",
                "Carb_sym_per_area", "Protein_host_per_area", "Protein_sym_per_area", "PAM_d13", "Feeding_d13")
colnames(df_PCA)[9:17]<-c("sd", "chlA", "chlc2", "hc", "sc","hp", "sp", "Fv/Fm", "feeding", "color")
View(df_PCA)

PCA <- df_PCA %>%
  na.omit(df_PCA) %>%
  filter_all(all_vars(!is.infinite(.)))

#log transformed
PCA$sd = PCA$sd+1
PCA$chlA = PCA$chlA+1
PCA$chlc2 = PCA$chlc2+1
PCA$sc = PCA$sc+1
PCA$sp = PCA$sp +1
PCA$feeding = as.numeric(PCA$feeding)

PCA_log <- log(PCA[9:17])
str(PCA_log)

#calculate PCs
phys_pca_log<-prcomp(PCA_log, center = T, scale. = T)
summary(phys_pca_log)
#PC1=0.3588
#PC2=0.1902

var <- as.data.frame(phys_pca_log$rotation)
var$Variable <- rownames(var)
scaling_factor <- 5  # Adjust this factor as needed
var$PC1 <- var$PC1 * scaling_factor
var$PC2 <- var$PC2 * scaling_factor

phys_pcaData <- as.data.frame(phys_pca_log$x)
phys_pcaData$State <- PCA$State
phys_pcaData$Treatment <- PCA$Treatment
phys_pcaData$Heat <- PCA$Heat
phys_pcaData$Light <- PCA$Light

phys_pca<-prcomp(PCA[9:17], center = T, scale. = T)

#apo vs. sym
PCA_phys_symstat <- adonis2(phys_pca$x~PCA$State,
                            by="term",
                            method = 'eu')
PCA_phys_symstat
adonis2(formula = phys_pca$x ~ PCA$State, method = "eu", by = "term")
#Df SumOfSqs      R2     F Pr(>F)  
#PCA$State  1    22.52 0.04169 2.567  0.024 *
#Residual  59   517.48 0.95831               
#Total     60   540.00 1.00000 

pca_pval_symstat <- PCA_phys_symstat[["Pr(>F)"]]
symbiosis_pca_pval_symstat <-substitute(italic(P[Symbiosis])==p, list(p = format(pca_pval_symstat[1], digits = 4)))

phys_PCA_symstat <- ggplot(phys_pcaData, aes(PC1, PC2)) + 
  geom_point(aes(colour = State), size = 3, stroke = 1.25) +
  scale_colour_manual(limits=c("White", "Brown"),
                      values=c("wheat3","salmon4"),
                      labels=c("White", "Brown"),
                      name="Symbiotic Phenotype") +
  stat_ellipse(geom="polygon",type = "t", alpha = 0, 
               aes(color= State), show.legend = FALSE, linewidth=1)+
  annotate("text", x = 2.5, y = 3.5, label = deparse(symbiosis_pca_pval_symstat), parse = TRUE, size = 4.5, hjust=0) +
  labs(x = "PC1 (32.75%)", 
       y = "PC2 (17.05%)")+
  theme_classic()+
  geom_segment(data = var, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = var, aes(x = PC1, y = PC2, label = Variable), 
            size = 5,vjust = 1, hjust = 1)

phys_PCA_symstat+theme(
  plot.title=element_blank(),
  axis.ticks=element_line(linewidth=1),
  axis.title = element_text(size=14, face = "bold"),
  axis.text=element_text(size=12),
  legend.title = element_text(size=14, face = "bold"),
  legend.text = element_text(size=12))

#symstat x treatment
PCA_phys <- adonis2(phys_pca$x~PCA$State*
                      PCA$Heat*
                      PCA$Light,
                    by="term",
                    method = 'eu')
PCA_phys
#Df SumOfSqs      R2      F Pr(>F)    
#PCA$State                     1    22.52 0.04169 3.2284  0.006 ** 
#PCA$Heat                      1    12.58 0.02330 1.8038  0.106    
#PCA$Light                     2   114.43 0.21190 8.2036  0.001 ***
#PCA$State:PCA$Heat            1     6.16 0.01140 0.8828  0.478    
#PCA$State:PCA$Light           2    13.70 0.02537 0.9823  0.474    
#PCA$Heat:PCA$Light            2    22.08 0.04088 1.5828  0.080 .  
#PCA$State:PCA$Heat:PCA$Light  2     6.81 0.01260 0.4879  0.922    
#Residual                     49   341.74 0.63284    

pca_pval <- PCA_phys[["Pr(>F)"]]
symbiosis_pca_pval <-substitute(italic(P[SP])==p, list(p = format(pca_pval[1], digits = 4)))
heat_pca_pval <-substitute(italic(P[H])==p, list(p = format(pca_pval[2], digits = 4)))
light_pca_pval <-substitute(italic(P[L])==p, list(p = format(pca_pval[3], digits = 4)))
sxh_pca_pval <-substitute(italic(P[SP~x~L])==p, list(p = format(pca_pval[4], digits = 4)))
sxl_pca_pval <-substitute(italic(P[SP~x~L])==p, list(p = format(pca_pval[5], digits = 4)))
hxl_pca_pval <-substitute(italic(P[H~x~L])==p, list(p = format(pca_pval[6], digits = 4)))
sxhxl_pca_pval <-substitute(italic(P[SP~x~H~x~L])==p, list(p = format(pca_pval[7], digits = 4)))

shapes = c("White" = 1, "Brown" = 19)

phys_PCA <- ggplot(phys_pcaData, aes(PC1, PC2)) + 
  geom_point(aes(colour = Treatment, shape = State), size = 3, stroke = 1.25) +
  scale_colour_manual(limits=c("NL_C","CL_C", "HL_C", "NL_H", "CL_H",  "HL_H"),
                      values=c("navy", "royalblue","steelblue1","deeppink4",  "palevioletred1",  "rosybrown1"),
                      labels=c("Control Temp, Dark","Control Temp, Control Light", "Control Temp, High Light","Heated, Dark", "Heated, Control Light",  "Heated, High Light"),
                      name="Treatment") +
  stat_ellipse(geom="polygon",type = "t", alpha = 0, 
               aes(color= Treatment), show.legend = FALSE, linewidth=1)+
  scale_shape_manual(values = c(1,19), breaks = c("White", "Brown"), labels = c("White", "Brown"), name = "Symbitoic Phenotype") +
  annotate("text", x = -7, y = 4.5, label = deparse(symbiosis_pca_pval), parse = TRUE, size = 4.5, hjust=0) +
  annotate("text", x = -7, y = 3.5, label = deparse(light_pca_pval), parse = TRUE, size = 4.5, hjust=0) +
  labs(x = "PC1 (35.88%)", 
       y = "PC2 (19.02%)")+
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme_classic()+
  geom_segment(data = var, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = var, aes(x = PC1, y = PC2, label = Variable), 
            size = 5,vjust = 1, hjust = 1) + 
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
phys_PCA
svg("../Figures/phys_PCA.svg", width =7.58, height = 4.68)
plot(phys_PCA)  
dev.off()
ggsave("../Figures/phys_PCA.svg", phys_PCA, width = 7.58, height = 4.68, units = "in", device = "svg")


#---------- Phenome Plasticity -----------
#Compare White & Brown
df_PCs <- as.data.frame(phys_pca$x[, 1:2])
df_PCs$Sample_ID <- PCA$Sample_ID
df_PCs$symstat <- PCA$State
df_PCs$Treatment <- PCA$Treatment
df_PCs$stat_treat <- sub("_.*", "", df_PCs$Sample_ID)
df_PCs$Heat <- PCA$Heat
df_PCs$Light <- PCA$Light
df_PCs$Tank <- PCA$Tank
df_PCs$Genet <- PCA$Genet

str(df_PCs)
df_PCs$Treatment <- as.factor(df_PCs$Treatment)
df_PCs$stat_treat <- as.factor(df_PCs$stat_treat)
df_PCs$symstat <- as.factor(df_PCs$symstat)
df_PCs$Heat <- as.factor(df_PCs$Heat)
df_PCs$Light <- as.factor(df_PCs$Light)
df_PCs$Tank <- as.factor(PCA$Tank)
df_PCs$Genet <- as.factor(PCA$Genet)

#Calculate the PC center (mean) for ACC and SCC samples
PC.center_ACC <- df_PCs %>%
  filter(stat_treat == "ACC") %>%
  summarise(PC1.avg = mean(PC1), PC2.avg = mean(PC2))

PC.center_SCC <- df_PCs %>%
  filter(stat_treat == "SCC") %>%
  summarise(PC1.avg = mean(PC1), PC2.avg = mean(PC2))

#Assign baseline PC centers to samples starting with "A" and "S"
df_PCs <- df_PCs %>%
  mutate(
    baseline_PC1 = ifelse(grepl("^A", Sample_ID), PC.center_ACC$PC1.avg, 
                          ifelse(grepl("^S", Sample_ID), PC.center_SCC$PC1.avg, NA)),
    baseline_PC2 = ifelse(grepl("^A", Sample_ID), PC.center_ACC$PC2.avg, 
                          ifelse(grepl("^S", Sample_ID), PC.center_SCC$PC2.avg, NA))
  )

#Calculate Euclidean distance for each sample from the baseline PC center
df_PCs <- df_PCs %>%
  mutate(
    distance = sqrt((PC1 - baseline_PC1)^2 + (PC2 - baseline_PC2)^2)
  )

#df_PCs_nocontrols = df_PCs %>%
#  dplyr::filter(stat_treat != "ACC") %>%
#  dplyr::filter(stat_treat != "SCC")
#df_PCs_nocontrols$Treatment <- factor(df_PCs_nocontrols$Treatment,
#                                      levels = c("NL_C", "CL_C","HL_C",
#                                                 "NL_H", "CL_H", "HL_H"))

pd <- position_dodge(width = 0.8)
fig_Phys_Plast <- ggplot(
  df_PCs,
  aes(x = Light, y = distance,
      color = interaction(Heat, Light))
) +
  geom_boxplot(
    fill = NA,
    alpha = 0,
    outlier.colour = NA,
    linewidth = 1,
    position = position_dodge(width = 0.8)
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
  facet_wrap(~ symstat) +
  ylab("Phenome Plasticity") +
  xlab("Light Treatment") +
  scale_color_manual(
    values = c(
      "Control.Control" = "royalblue",
      "Control.High" = "steelblue1",
      "Control.Dark" = "navy",
      "Heat.Control" = "palevioletred1",
      "Heat.High" = "rosybrown1",
      "Heat.Dark" = "deeppink4"
    )
  ) +
  scale_shape_manual(values = c("White" = 1, "Brown" = 19)) +
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
fig_Phys_Plast
ggsave("../Figures/fig_Phys_Plast.svg", fig_Phys_Plast, width = 7.58, height = 4.68, units = "in", device = "svg")

df_PCs.apo <- df_PCs %>% filter(symstat=="White")
df_PCs.sym <- df_PCs %>% filter(symstat=="Brown")
lm_Phys_Plast <- lmer(distance ~ symstat+(1|Genet)+(1|Tank), data = df_PCs)
summary(lm_Phys_Plast)
summary(glht(lm_Phys_Plast, linfct = mcp(symstat = "Tukey")))  #N.S.

mod1 <- lmer(distance ~ Heat*Light+(1|Genet)+(1|Tank), data = df_PCs.apo, REML=FALSE)
mod2 <- lmer(distance ~ Heat+Light+(1|Genet)+(1|Tank), data = df_PCs.apo, REML=FALSE)
AICc(mod1, mod2) 
#     df      AICc
#mod1  9 100.51553
#mod2  7  93.26957

lm_Phys_Plast.apo <- lmer(distance ~ Heat+Light+(1|Genet)+(1|Tank), data = df_PCs.apo)
summary(lm_Phys_Plast.apo)
summary(glht(lm_Phys_Plast.apo, linfct = mcp(Heat = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)  
#Heat - Control == 0   0.8261     0.3510   2.354   0.0186 *
summary(glht(lm_Phys_Plast.apo, linfct = mcp(Light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)    
#Control - Dark == 0  -1.6837     0.4159  -4.049   <0.001 ***
#High - Dark == 0     -0.8620     0.3912  -2.203   0.0705 .  
#High - Control == 0   0.8217     0.4053   2.028   0.1055  
lm_Phys_Plast.apo_treat <- lmer(distance ~ Treatment+(1|Genet)+(1|Tank), data = df_PCs.apo)
summary(lm_Phys_Plast.apo_treat)
summary(glht(lm_Phys_Plast.apo_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.apo_treat <- glht(lm_Phys_Plast.apo_treat,
                    linfct = mcp(Treatment = "Tukey"))
cld(posthoc.apo_treat)

mod1 <- lmer(distance ~ Heat*Light+(1|Genet)+(1|Tank), data = df_PCs.sym, REML=FALSE)
mod2 <- lmer(distance ~ Heat+Light+(1|Genet)+(1|Tank), data = df_PCs.sym, REML=FALSE)
AICc(mod1, mod2) 
#df     AICc
#mod1  9 95.83576
#mod2  7 97.08718
lm_Phys_Plast.sym <- lmer(distance ~ Heat+Light+(1|Genet)+(1|Tank), data = df_PCs.sym)
summary(lm_Phys_Plast.sym)
summary(glht(lm_Phys_Plast.sym, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Phys_Plast.sym, linfct = mcp(Light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)   
#Control - Dark == 0 -1.80759    0.53989  -3.348  0.00244 **
#High - Dark == 0    -1.74094    0.56062  -3.105  0.00537 **
#High - Control == 0  0.06665    0.54513   0.122  0.99179
lm_Phys_Plast.sym_treat <- lmer(distance ~ Treatment+(1|Genet)+(1|Tank), data = df_PCs.sym)
summary(lm_Phys_Plast.sym_treat)
summary(glht(lm_Phys_Plast.sym_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.sym_treat <- glht(lm_Phys_Plast.sym_treat,
                          linfct = mcp(Treatment = "Tukey"))
cld(posthoc.sym_treat)

