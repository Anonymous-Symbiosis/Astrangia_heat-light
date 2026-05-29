#Astrangia poculata heat & light: Analyses of physiological traits
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Fig. 3A in the manuscript, Fig. S1 in the supplementary and prepare input data frame for PCA analysis
#Note: In this script: apo = white corals, sym = brown corals


data <- read.csv("../Data/Metadata_Astrangia_Clean_May1.csv")

library(dplyr)
library(ggplot2)
library(lme4)
library(multcomp)
library(performance)
library(MuMIn)
library(scales)
library(pavo)
library(gridExtra)
library(writexl)
library(tidyverse)
library(ordinal)
library(car)
library(emmeans)


#format the data frame
data$Stat_treat <- sub("_.*", "", data$Sample_ID)
data$Stat_treat = as.factor(data$Stat_treat)
data$Genet = as.factor(data$Genet)
data$Tank = as.factor(data$Tank)
data$State = factor(data$State, 
                    levels = c("apo", "sym"),
                    labels = c("White", "Brown"))
data$Light = factor(data$Light, 
                    levels = c("No Light", "Control", "High Light"),
                    labels = c("Dark", "Control", "High"))
data$Treatment = factor(data$Treatment,
                        levels = c("NL_C", "CL_C","HL_C",
                                   "NL_H", "CL_H", "HL_H"))
data$Heat = as.factor(data$Heat)

#----------Functions to remove outliers based on IQR method----------
detect_outlier <- function(x) {
  Quantile1 <- quantile(x, probs=.25)
  Quantile3 <- quantile(x, probs=.75)
  IQR = Quantile3-Quantile1
  x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}
remove_outlier <- function(dataframe,
                           columns=names(dataframe)) {
  for (col in columns) {
    dataframe <- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  print("Remove outliers")
  print(dataframe)
}


#----------Symbiont density----------
Sym.den <- data %>%
  mutate(Sym_cell_per_area = (Sym_density*Total_Slurry)/(Surface_area*0.01)) %>% #change surface area from mm^2 to cm^2
  na.omit(Sym.den) #remove NA/NaN
Sym.den.W <- Sym.den %>%
  filter(State == "White")
Sym.den.B <- Sym.den %>%
  filter(State == "Brown")

fig_Sym.den <- ggplot(Sym.den, aes(x = Light, y = Sym_cell_per_area, fill = interaction(Heat, Light), color = interaction(Heat, Light))) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, linewidth = 1, position = position_dodge(width = 0.8)) +
  geom_point(aes(shape = State), size = 3, stroke = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_wrap(~ State) +
  ylab(expression(Symbiont~Density~(cells/cm^2))) +
  xlab("Light Treatment") +
  scale_color_manual(values = c(
      "Control.Control" = "royalblue",
      "Control.High" = "steelblue1",
      "Control.Dark" = "navy",
      "Heat.Control" = "palevioletred1",
      "Heat.High" = "rosybrown1",
      "Heat.Dark" = "deeppink4")) +
  scale_fill_manual(values = c(
      "Control.Control" = "royalblue",
      "Control.High" = "steelblue1",
      "Control.Dark" = "navy",
      "Heat.Control" = "palevioletred1",
      "Heat.High" = "rosybrown1",
      "Heat.Dark" = "deeppink4")) +
  scale_shape_manual(values = c("White" = 1, "Brown" = 19)) +
  scale_y_continuous(labels = scientific_format()) +
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
fig_Sym.den
ggsave("../Figures_Factorial/Symbiont Density.svg", fig_Sym.den, width = 7, height = 4.68, units = "in", device = "svg")

#statistics 
lm_Sym.den <-  lmer(Sym_cell_per_area ~ State+(1|Genet)+(1|Tank), data = Sym.den)
summary(lm_Sym.den)
summary(glht(lm_Sym.den, linfct = mcp(State = "Tukey"))) 
#Brown - White 1.97e-05 ***
mod1 <- lmer(Sym_cell_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Sym.den.W, REML = FALSE)
mod2 <- lmer(Sym_cell_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Sym.den.W, REML = FALSE)
AICc(mod1, mod2)
#df     AICc
#mod1  9 1187.530
#mod2  7 1181.012
lm_Sym.den.W = lmer(Sym_cell_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Sym.den.W)
summary(lm_Sym.den.W)
summary(glht(lm_Sym.den.W, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Sym.den.W, linfct = mcp(Light = "Tukey"))) #N.S.
lm_Sym.den.W_treat <- lmer(Sym_cell_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Sym.den.W)
summary(lm_Sym.den.W_treat)
summary(glht(lm_Sym.den.W_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Sym.den.W_treat <- glht(lm_Sym.den.W_treat,
                          linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Sym.den.W_treat)

mod1 <- lmer(Sym_cell_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Sym.den.B, REML = FALSE)
mod2 <- lmer(Sym_cell_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Sym.den.B, REML = FALSE)
AICc(mod1, mod2)
#     df     AICc
#mod1  9 1234.569
#mod2  7 1229.573
lm_Sym.den.B = lmer(Sym_cell_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Sym.den.B)
summary(lm_Sym.den.B)
summary(glht(lm_Sym.den.B, linfct = mcp(Heat = "Tukey")))
#Heat - Control == 0  0.0197 *
summary(glht(lm_Sym.den.B, linfct = mcp(Light = "Tukey"))) #N.S.
lm_Sym.den.B_treat <- lmer(Sym_cell_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Sym.den.B)
summary(lm_Sym.den.B_treat)
summary(glht(lm_Sym.den.B_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Sym.den.B_treat <- glht(lm_Sym.den.B_treat,
                                linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Sym.den.B_treat)


#----------Chlorophyll a----------
Sym_chl <- Sym.den %>%
  mutate(Chla_per_area = ((Chl_a*5)/(Surface_area*0.01))) %>% #extrapolated to symbiont slurry volume (5 mL) #change surface area from mm^2 to cm^2
  mutate(Chlc2_per_area = ((Chl_c2*5)/(Surface_area*0.01))) %>%
  mutate(Chla_per_cell = Chla_per_area/Sym_cell_per_area) %>%
  mutate(Chlc2_per_cell = Chlc2_per_area/Sym_cell_per_area) %>%
  na.omit(Sym_chl) 

Chla_IQR <- remove_outlier(Sym_chl,'Chla_per_area')
Chla.cell_IQR <- remove_outlier(Sym_chl,'Chla_per_cell')
Chla_IQR.W <- Chla_IQR %>%
  filter(State == "White")
Chla_IQR.B <- Chla_IQR %>%
  filter(State == "Brown")
Chla.cell_IQR.W <- Chla.cell_IQR %>%
  filter(State == "White")
Chla.cell_IQR.B <- Chla.cell_IQR %>%
  filter(State == "Brown")

fig_Chla <- ggplot(Chla_IQR, aes(x = Light, y = Chla_per_area, fill = interaction(Heat, Light), color = interaction(Heat, Light))) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, linewidth = 1, position = position_dodge(width = 0.8)) +
  geom_point(aes(shape = State), size = 3, stroke = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_wrap(~ State) +
  ylab(expression(Chlorophyll~a~(μg/cm^2))) +
  xlab("Light Treatment") +
  scale_color_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
  scale_fill_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
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
fig_Chla
ggsave("../Figures_Factorial/Chla.svg", fig_Chla, width = 7, height = 4.68, units = "in", device = "svg")

#statistics 
lm_Chla = lmer(Chla_per_area ~ State+(1|Genet)+(1|Tank), data = Chla_IQR)
summary(lm_Chla)
summary(glht(lm_Chla, linfct = mcp(State = "Tukey"))) 
#Brown - White   2.703  0.0275 *
mod1 <- lmer(Chla_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Chla_IQR.W, REML = FALSE)
mod2 <- lmer(Chla_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Chla_IQR.W, REML = FALSE)
AICc(mod1, mod2)
#     df     AICc
#mod1  9 281.7604
#mod2  7 275.0212
lm_Chla.W = lmer(Chla_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Chla_IQR.W)
summary(lm_Chla.W)
summary(glht(lm_Chla.W, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Chla.W, linfct = mcp(Light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|) 
#Control - Dark == 0   16.897      3.822   4.421   <1e-04 ***
#High - Dark == 0      10.042      3.733   2.690   0.0194 *  
#High - Control == 0   -6.856      3.670  -1.868   0.1479 
lm_Chla.W_treat <- lmer(Chla_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Chla_IQR.W)
summary(lm_Chla.W_treat)
summary(glht(lm_Chla.W_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Chla.W_treat <- glht(lm_Chla.W_treat,
                                linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Chla.W_treat)

mod1 <- lmer(Chla_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Chla_IQR.B, REML = FALSE)
mod2 <- lmer(Chla_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Chla_IQR.B, REML = FALSE)
AICc(mod1, mod2)
#df     AICc
#mod1  9 307.7871
#mod2  7 301.7543
lm_Chla.B = lmer(Chla_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Chla_IQR.B)
summary(lm_Chla.B)
summary(glht(lm_Chla.B, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Chla.B, linfct = mcp(Light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)   
#Control - Dark == 0   19.992      5.808   3.442  0.00162 **
#High - Dark == 0      16.777      6.109   2.746  0.01668 * 
#High - Control == 0   -3.215      5.895  -0.545  0.84879  
lm_Chla.B_treat <- lmer(Chla_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Chla_IQR.B)
summary(lm_Chla.B_treat)
summary(glht(lm_Chla.B_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Chla.B_treat <- glht(lm_Chla.B_treat,
                             linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Chla.B_treat)


#----------Chlorophyll c2----------
Chlc2_IQR <- remove_outlier(Sym_chl,'Chlc2_per_area')
Chlc2.cell_IQR <- remove_outlier(Sym_chl,'Chlc2_per_cell')
Chlc2_IQR.W <- Chlc2_IQR %>%
  filter(State == "White")
Chlc2_IQR.B <- Chlc2_IQR %>%
  filter(State == "Brown")
Chlc2.cell_IQR.W <- Chlc2.cell_IQR %>%
  filter(State == "White")
Chlc2.cell_IQR.B <- Chlc2.cell_IQR %>%
  filter(State == "Brown")

fig_Chlc2 <- ggplot(Chlc2_IQR, aes(x = Light, y = Chlc2_per_area, fill = interaction(Heat, Light), color = interaction(Heat, Light))) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, linewidth = 1, position = position_dodge(width = 0.8)) +
  geom_point(aes(shape = State), size = 3, stroke = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_wrap(~ State) +
  ylab(expression(Chlorophyll~c2~(μg/cm^2))) +
  xlab("Light Treatment") +
  scale_color_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
  scale_fill_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
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
fig_Chlc2
ggsave("../Figures_Factorial/Chlc2.svg", fig_Chlc2, width = 7, height = 4.68, units = "in", device = "svg")

#statistics 
lm_Chlc2 = lmer(Chlc2_per_area ~ State+(1|Genet)+(1|Tank), data = Chlc2_IQR)
summary(lm_Chlc2)
summary(glht(lm_Chlc2, linfct = mcp(State = "Tukey"))) 
#Brown - White 0.0394 *
mod1 <- lmer(Chlc2_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Chlc2_IQR.W, REML = FALSE)
mod2 <- lmer(Chlc2_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Chlc2_IQR.W, REML = FALSE)
AICc(mod1, mod2)
#df     AICc
#mod1  9 236.2488
#mod2  7 229.7112
lm_Chlc2.W = lmer(Chlc2_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Chlc2_IQR.W)
summary(lm_Chlc2.W)
summary(glht(lm_Chlc2.W, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Chlc2.W, linfct = mcp(Light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)    
#Control - Dark == 0    8.061      2.037   3.957   <0.001 ***
#High - Dark == 0       5.012      1.990   2.519   0.0315 *  
#High - Control == 0   -3.049      1.956  -1.559   0.2635
lm_Chlc2.W_treat <- lmer(Chlc2_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Chlc2_IQR.W)
summary(lm_Chlc2.W_treat)
summary(glht(lm_Chlc2.W_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Chlc2.W_treat <- glht(lm_Chlc2.W_treat,
                             linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Chlc2.W_treat)

mod1 <- lmer(Chlc2_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Chlc2_IQR.B, REML = FALSE)
mod2 <- lmer(Chlc2_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Chlc2_IQR.B, REML = FALSE)
AICc(mod1, mod2)
#df     AICc
#mod1  9 231.9742
#mod2  7 226.9664
lm_Chlc2.B = lmer(Chlc2_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Chlc2_IQR.B)
summary(lm_Chlc2.B)
summary(glht(lm_Chlc2.B, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Chlc2.B, linfct = mcp(Light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)   
#Control - Dark == 0   10.128      2.984   3.394  0.00202 **
#High - Dark == 0       4.708      3.158   1.491  0.29511   
#High - Control == 0   -5.420      3.039  -1.784  0.17498 
lm_Chlc2.B_treat <- lmer(Chlc2_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Chlc2_IQR.B)
summary(lm_Chlc2.B_treat)
summary(glht(lm_Chlc2.B_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Chlc2.B_treat <- glht(lm_Chlc2.B_treat,
                              linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Chlc2.B_treat)


#----------Host carbohydrate----------
Carb <- data %>%
  mutate(Sym_cell_per_area = (Sym_density*Total_Slurry)/(Surface_area*0.01)) %>% #change surface area from mm^2 to cm^2
  mutate(Carb_host_per_area = ((Carb_host*Total_Slurry)/(Surface_area*0.01))) %>%  #change surface area from mm^2 to cm^2
  mutate(Carb_sym_per_area = ((Carb_sym*Total_Slurry)/(Surface_area*0.01))) %>%  #change surface area from mm^2 to cm^2
  na.omit() #remove NA/NaN
Carb_host_IQR <- remove_outlier(Carb,'Carb_host_per_area')
Carb_host_IQR.W <- Carb_host_IQR %>%
  filter(State == "White")
Carb_host_IQR.B <- Carb_host_IQR %>%
  filter(State == "Brown")

fig_Carb_host <- ggplot(Carb_host_IQR, aes(x = Light, y = Carb_host_per_area, fill = interaction(Heat, Light), color = interaction(Heat, Light))) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, linewidth = 1, position = position_dodge(width = 0.8)) +
  geom_point(aes(shape = State), size = 3, stroke = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_wrap(~ State) +
  ylab(expression(Host~Carbohydrate~(μg/cm^2))) +
  xlab("Light Treatment") +
  scale_color_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
  scale_fill_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
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
fig_Carb_host
ggsave("../Figures_Factorial/Carb_host.svg", fig_Carb_host, width = 7, height = 4.68, units = "in", device = "svg")

#statistics 
lm_Carb_host = lmer(Carb_host_per_area ~ State+(1|Genet)+(1|Tank), data = Carb_host_IQR)
summary(lm_Carb_host)
summary(glht(lm_Carb_host, linfct = mcp(State = "Tukey")))  #N.S.
mod1 <- lmer(Carb_host_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Carb_host_IQR.W, REML=FALSE)
mod2 <- lmer(Carb_host_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Carb_host_IQR.W, REML=FALSE)
AICc(mod1, mod2)
#     df       AICc
#mod1  9  -7.179354
#mod2  7 -13.263495
lm_Carb_host.W = lmer(Carb_host_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Carb_host_IQR.W)
summary(lm_Carb_host.W)
summary(glht(lm_Carb_host.W, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Carb_host.W, linfct = mcp(Light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|) 
#High - Dark == 0     0.16850    0.06529   2.581   0.0266 *
lm_Carb_host.W_treat <- lmer(Carb_host_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Carb_host_IQR.W)
summary(lm_Carb_host.W_treat)
summary(glht(lm_Carb_host.W_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Carb_host.W_treat <- glht(lm_Carb_host.W_treat,
                              linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Carb_host.W_treat)

mod1 <- lmer(Carb_host_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Carb_host_IQR.B, REML=FALSE)
mod2 <- lmer(Carb_host_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Carb_host_IQR.B, REML=FALSE)
AICc(mod1, mod2)
#     df      AICc
#mod1  9 -26.26534
#mod2  7 -32.33164
lm_Carb_host.B = lmer(Carb_host_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Carb_host_IQR.B)
summary(lm_Carb_host.B)
summary(glht(lm_Carb_host.B, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Carb_host.B, linfct = mcp(Light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)  
#High - Dark == 0     0.14960    0.05931   2.522   0.0313 *
lm_Carb_host.B_treat <- lmer(Carb_host_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Carb_host_IQR.B)
summary(lm_Carb_host.B_treat)
summary(glht(lm_Carb_host.B_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Carb_host.B_treat <- glht(lm_Carb_host.B_treat,
                                  linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Carb_host.B_treat)


#----------Symbiont carbohydrate----------
Carb_sym_IQR <- remove_outlier(Carb,'Carb_sym_per_area')
Carb_sym_IQR.W <- Carb_sym_IQR %>%
  filter(State == "White")
Carb_sym_IQR.B <- Carb_sym_IQR %>%
  filter(State == "Brown")

fig_Carb_sym <- ggplot(Carb_sym_IQR, aes(x = Light, y = Carb_sym_per_area, fill = interaction(Heat, Light), color = interaction(Heat, Light))) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, linewidth = 1, position = position_dodge(width = 0.8)) +
  geom_point(aes(shape = State), size = 3, stroke = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_wrap(~ State) +
  ylab(expression(Symbiont~Carbohydrate~(μg/cm^2))) +
  xlab("Light Treatment") +
  scale_color_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
  scale_fill_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
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
fig_Carb_sym
ggsave("../Figures_Factorial/Carb_symbiont.svg", fig_Carb_sym, width = 7, height = 4.68, units = "in", device = "svg")

#statistics 
lm_Carb_sym = lmer(Carb_sym_per_area ~ State+(1|Genet)+(1|Tank), data = Carb_sym_IQR)
summary(lm_Carb_sym)
summary(glht(lm_Carb_sym, linfct = mcp(State = "Tukey")))  #N.S.
mod1 <- lmer(Carb_sym_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Carb_sym_IQR.W, REML = FALSE)
mod2 <- lmer(Carb_sym_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Carb_sym_IQR.W, REML = FALSE)
AICc(mod1, mod2)
#     df      AICc
#mod1  9 -85.54187
#mod2  7 -82.88756
lm_Carb_sym.W = lmer(Carb_sym_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Carb_sym_IQR.W)
summary(lm_Carb_sym.W)
summary(glht(lm_Carb_sym.W, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Carb_sym.W, linfct = mcp(Light = "Tukey"))) #N.S.
lm_Carb_sym.W_treat <- lmer(Carb_sym_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Carb_sym_IQR.W)
summary(lm_Carb_sym.W_treat)
summary(glht(lm_Carb_sym.W_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Carb_sym.W_treat <- glht(lm_Carb_sym.W_treat,
                                  linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Carb_sym.W_treat)
mod1 <- lmer(Carb_sym_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Carb_sym_IQR.B, REML = FALSE)
mod2 <- lmer(Carb_sym_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Carb_sym_IQR.B, REML = FALSE)
AICc(mod1, mod2)
#     df      AICc
#mod1  9 -97.96913
#mod2  7 -97.23336
lm_Carb_sym.B = lmer(Carb_sym_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Carb_sym_IQR.B)
summary(lm_Carb_sym.B)
summary(glht(lm_Carb_sym.B, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Carb_sym.B, linfct = mcp(Light = "Tukey"))) #N.S.
lm_Carb_sym.B_treat <- lmer(Carb_sym_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Carb_sym_IQR.B)
summary(lm_Carb_sym.B_treat)
summary(glht(lm_Carb_sym.B_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Carb_sym.B_treat <- glht(lm_Carb_sym.B_treat,
                                 linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Carb_sym.B_treat)


#----------Host protein----------
Protein <- data %>%
  mutate(Sym_cell_per_area = (Sym_density*Total_Slurry)/(Surface_area*0.01)) %>% #change surface area from mm^2 to cm^2
  mutate(Protein_host_per_area = (((Protein_host/0.08)*Total_Slurry)/(Surface_area*0.01))) %>% #dilution=0.08;change surface area from mm^2 to cm^2
  mutate(Protein_sym_per_area = (((Protein_sym/0.08)*Total_Slurry)/(Surface_area*0.01))) %>% #dilution=0.08;change surface area from mm^2 to cm^2
  na.omit(Protein) #remove NA/NaN

Protein_host_IQR <- remove_outlier(Protein,'Protein_host_per_area')
Protein_host_IQR.W <- Protein_host_IQR %>%
  filter(State == "White")
Protein_host_IQR.B <- Protein_host_IQR %>%
  filter(State == "Brown")

fig_Protein_host <- ggplot(Protein_host_IQR, aes(x = Light, y = Protein_host_per_area, fill = interaction(Heat, Light), color = interaction(Heat, Light))) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, linewidth = 1, position = position_dodge(width = 0.8)) +
  geom_point(aes(shape = State), size = 3, stroke = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_wrap(~ State) +
  ylab(expression(Host~Protein~(μg/cm^2))) +
  xlab("Light Treatment") +
  scale_color_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
  scale_fill_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
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
fig_Protein_host
ggsave("../Figures_Factorial/Protein_host.svg", fig_Protein_host, width = 7, height = 4.68, units = "in", device = "svg")

#statistics 
lm_Protein_host = lmer(Protein_host_per_area ~ State+(1|Genet)+(1|Tank), data = Protein_host_IQR)
summary(lm_Protein_host)
summary(glht(lm_Protein_host, linfct = mcp(State = "Tukey")))  #N.S.
mod1 <- lmer(Protein_host_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Protein_host_IQR.W, REML = FALSE)
mod2 <- lmer(Protein_host_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Protein_host_IQR.W, REML = FALSE)
AICc(mod1, mod2)
#     df     AICc
#mod1  9 39.89133
#mod2  7 34.56212
lm_Protein_host.W = lmer(Protein_host_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Protein_host_IQR.W)
summary(lm_Protein_host.W)
summary(glht(lm_Protein_host.W, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Protein_host.W, linfct = mcp(Light = "Tukey"))) #N.S.
lm_Protein_host.W_treat <- lmer(Protein_host_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Protein_host_IQR.W)
summary(lm_Protein_host.W_treat)
summary(glht(lm_Protein_host.W_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Protein_host.W_treat <- glht(lm_Protein_host.W_treat,
                                 linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Protein_host.W_treat)
mod1 <- lmer(Protein_host_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Protein_host_IQR.B, REML = FALSE)
mod2 <- lmer(Protein_host_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Protein_host_IQR.B, REML = FALSE)
AICc(mod1, mod2)
#     df      AICc
#mod1  9 2.7978882
#mod2  7 0.3068078
lm_Protein_host.B = lmer(Protein_host_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Protein_host_IQR.B)
summary(lm_Protein_host.B)
summary(glht(lm_Protein_host.B, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Protein_host.B, linfct = mcp(Light = "Tukey"))) #N.S.
lm_Protein_host.B_treat <- lmer(Protein_host_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Protein_host_IQR.B)
summary(lm_Protein_host.B_treat)
summary(glht(lm_Protein_host.B_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Protein_host.B_treat <- glht(lm_Protein_host.B_treat,
                                     linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Protein_host.B_treat)


#----------Symbiont protein----------
Protein_sym_IQR <- remove_outlier(Protein,'Protein_sym_per_area')
Protein_sym_IQR.W <- Protein_sym_IQR %>%
  filter(State == "White")
Protein_sym_IQR.B <- Protein_sym_IQR %>%
  filter(State == "Brown")

fig_Protein_symbiont <- ggplot(Protein_sym_IQR, aes(x = Light, y = Protein_sym_per_area, fill = interaction(Heat, Light), color = interaction(Heat, Light))) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, linewidth = 1, position = position_dodge(width = 0.8)) +
  geom_point(aes(shape = State), size = 3, stroke = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_wrap(~ State) +
  ylab(expression(Symbiont~Protein~(μg/cm^2))) +
  xlab("Light Treatment") +
  scale_color_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
  scale_fill_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
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
fig_Protein_symbiont
ggsave("../Figures_Factorial/Protein_symbiont.svg", fig_Protein_symbiont, width = 7, height = 4.68, units = "in", device = "svg")

#statistics 
lm_Protein_sym = lmer(Protein_sym_per_area ~ State+(1|Genet)+(1|Tank), data = Protein_sym_IQR)
summary(lm_Protein_sym)
summary(glht(lm_Protein_sym, linfct = mcp(State = "Tukey")))  #N.S.
mod1 <- lmer(Protein_sym_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Protein_sym_IQR.W, REML = FALSE)
mod2 <- lmer(Protein_sym_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Protein_sym_IQR.W, REML = FALSE)
AICc(mod1, mod2)
#     df      AICc
#mod1  9 -81.46383
#mod2  7 -88.25595
lm_Protein_sym.W = lmer(Protein_sym_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Protein_sym_IQR.W)
summary(lm_Protein_sym.W)
summary(glht(lm_Protein_sym.W, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Protein_sym.W, linfct = mcp(Light = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)  
#High - Control == 0 -0.07172    0.02479  -2.893   0.0107 *
lm_Protein_sym.W_treat <- lmer(Protein_sym_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Protein_sym_IQR.W)
summary(lm_Protein_sym.W_treat)
summary(glht(lm_Protein_sym.W_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Protein_sym.W_treat <- glht(lm_Protein_sym.W_treat,
                                     linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Protein_sym.W_treat)
mod1 <- lmer(Protein_sym_per_area ~ Heat*Light+(1|Genet)+(1|Tank), data = Protein_sym_IQR.B, REML = FALSE)
mod2 <- lmer(Protein_sym_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Protein_sym_IQR.B, REML = FALSE)
AICc(mod1, mod2)
#     df      AICc
#mod1  9 -82.75503
#mod2  7 -86.99614
lm_Protein_sym.B = lmer(Protein_sym_per_area ~ Heat+Light+(1|Genet)+(1|Tank), data = Protein_sym_IQR.B)
summary(lm_Protein_sym.B)
summary(glht(lm_Protein_sym.B, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_Protein_sym.B, linfct = mcp(Light = "Tukey"))) #N.S.
lm_Protein_sym.B_treat <- lmer(Protein_sym_per_area ~ Treatment+(1|Genet)+(1|Tank), data = Protein_sym_IQR.B)
summary(lm_Protein_sym.B_treat)
summary(glht(lm_Protein_sym.B_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.Protein_sym.B_treat <- glht(lm_Protein_sym.B_treat,
                                    linfct = mcp(Treatment = "Tukey"))
cld(posthoc.Protein_sym.B_treat)


#----------Photochemical efficiency----------
PAM <- data%>%
  mutate(PAM_change=(PAM_d13-PAM_d6)/PAM_d6)%>%
  na.omit(PAM)%>%
  filter_all(all_vars(!is.infinite(.)))
PAM.W <- PAM %>%
  filter(State == "White")
PAM.B <- PAM %>%
  filter(State == "Brown")

fig_PAM_d6 <- ggplot(PAM, aes(x = Light, y = PAM_d6, fill = interaction(Heat, Light), color = interaction(Heat, Light))) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, linewidth = 1, position = position_dodge(width = 0.8)) +
  geom_point(aes(shape = State), size = 3, stroke = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_wrap(~ State) +
  ylab("Photochemical Efficiency on day 6 (Fv/Fm)") +
  xlab("Light Treatment") +
  scale_color_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
  scale_fill_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
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
fig_PAM_d6
ggsave("../Figures_Factorial/PAM_d6.svg", fig_PAM_d6, width = 7, height = 4.68, units = "in", device = "svg")

#statistics 
lm_PAM_d6 = lmer(PAM_d6 ~ State+(1|Genet)+(1|Tank), data = PAM)
summary(lm_PAM_d6)
summary(glht(lm_PAM_d6, linfct = mcp(State = "Tukey")))  #N.S.
mod1 <- lmer(PAM_d6 ~ Heat*Light+(1|Genet)+(1|Tank), data = PAM.W, REML = FALSE)
mod2 <- lmer(PAM_d6 ~ Heat+Light+(1|Genet)+(1|Tank), data = PAM.W, REML = FALSE)
AICc(mod1, mod2)
#df      AICc
#mod1  9 -75.88629
#mod2  7 -82.33813
lm_PAM_d6.W = lmer(PAM_d6 ~ Heat+Light+(1|Genet)+(1|Tank), data = PAM.W)
summary(lm_PAM_d6.W)
summary(glht(lm_PAM_d6.W, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_PAM_d6.W, linfct = mcp(Light = "Tukey"))) 
#Estimate Std. Error z value Pr(>|z|)    
#Control - Dark == 0  0.19295    0.04007   4.816   <1e-05 ***
#High - Control == 0 -0.24457    0.03988  -6.133   <1e-05 ***
lm_PAM_d6.W_treat <- lmer(PAM_d6 ~ Treatment+(1|Genet)+(1|Tank), data = PAM.W)
summary(lm_PAM_d6.W_treat)
summary(glht(lm_PAM_d6.W_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.PAM_d6.W_treat <- glht(lm_PAM_d6.W_treat,
                                    linfct = mcp(Treatment = "Tukey"))
cld(posthoc.PAM_d6.W_treat)
mod1 <- lmer(PAM_d6 ~ Heat*Light+(1|Genet)+(1|Tank), data = PAM.B, REML = FALSE)
mod2 <- lmer(PAM_d6 ~ Heat+Light+(1|Genet)+(1|Tank), data = PAM.B, REML = FALSE)
AICc(mod1, mod2)
#     df      AICc
#mod1  9 -21.37323
#mod2  7 -27.68407
lm_PAM_d6.B = lmer(PAM_d6 ~ Heat+Light+(1|Genet)+(1|Tank), data = PAM.B)
summary(lm_PAM_d6.B)
summary(glht(lm_PAM_d6.B, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_PAM_d6.B, linfct = mcp(Light = "Tukey"))) 
#Estimate Std. Error z value Pr(>|z|)    
#Control - Dark == 0  0.216903   0.057671   3.761 0.000545 ***
#High - Control == 0 -0.212208   0.057671  -3.680 0.000666 ***
lm_PAM_d6.B_treat <- lmer(PAM_d6 ~ Treatment+(1|Genet)+(1|Tank), data = PAM.B)
summary(lm_PAM_d6.B_treat)
summary(glht(lm_PAM_d6.B_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.PAM_d6.B_treat <- glht(lm_PAM_d6.B_treat,
                               linfct = mcp(Treatment = "Tukey"))
cld(posthoc.PAM_d6.B_treat)
  
fig_PAM_d13 <- ggplot(PAM, aes(x = Light, y = PAM_d13, fill = interaction(Heat, Light), color = interaction(Heat, Light))) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, linewidth = 1, position = position_dodge(width = 0.8)) +
  geom_point(aes(shape = State), size = 3, stroke = 1, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  facet_wrap(~ State) +
  ylab("Photochemical Efficiency on Day 13 (Fv/Fm)") +
  xlab("Light Treatment") +
  scale_color_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
  scale_fill_manual(values = c(
    "Control.Control" = "royalblue",
    "Control.High" = "steelblue1",
    "Control.Dark" = "navy",
    "Heat.Control" = "palevioletred1",
    "Heat.High" = "rosybrown1",
    "Heat.Dark" = "deeppink4")) +
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
fig_PAM_d13
ggsave("../Figures_Factorial/PAM_d13.svg", fig_PAM_d13, width = 7, height = 4.68, units = "in", device = "svg")

lm_PAM_d13 = lmer(PAM_d13 ~ State+(1|Genet)+(1|Tank), data = PAM)
summary(lm_PAM_d13)
summary(glht(lm_PAM_d13, linfct = mcp(State = "Tukey")))
#Estimate Std. Error z value Pr(>|z|)  
#Brown - White == 0  0.07876    0.03343   2.356   0.0185 *
lm_PAM_d13.W_treat <- lmer(PAM_d13 ~ Treatment+(1|Genet)+(1|Tank), data = PAM.W)
summary(lm_PAM_d13.W_treat)
summary(glht(lm_PAM_d13.W_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.PAM_d13.W_treat <- glht(lm_PAM_d13.W_treat,
                               linfct = mcp(Treatment = "Tukey"))
cld(posthoc.PAM_d13.W_treat)
mod1 <- lmer(PAM_d13 ~ Heat*Light+(1|Genet)+(1|Tank), data = PAM.W, REML = FALSE)
mod2 <- lmer(PAM_d13 ~ Heat+Light+(1|Genet)+(1|Tank), data = PAM.W, REML = FALSE)
AICc(mod1, mod2)
#     df      AICc
#mod1  9 -53.40083
#mod2  7 -57.67087
lm_PAM_d13.W = lmer(PAM_d13 ~ Heat+Light+(1|Genet)+(1|Tank), data = PAM.W)
summary(lm_PAM_d13.W)
summary(glht(lm_PAM_d13.W, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_PAM_d13.W, linfct = mcp(Light = "Tukey"))) #N.S.
mod1 <- lmer(PAM_d13 ~ Heat*Light+(1|Genet)+(1|Tank), data = PAM.B, REML = FALSE)
mod2 <- lmer(PAM_d13 ~ Heat+Light+(1|Genet)+(1|Tank), data = PAM.B, REML = FALSE)
AICc(mod1, mod2)
#     df      AICc
#mod1  9 -35.51814
#mod2  7 -40.66875
lm_PAM_d13.B = lmer(PAM_d13 ~ Heat+Light+(1|Genet)+(1|Tank), data = PAM.B)
summary(lm_PAM_d13.B)
summary(glht(lm_PAM_d13.B, linfct = mcp(Heat = "Tukey"))) #N.S.
summary(glht(lm_PAM_d13.B, linfct = mcp(Light = "Tukey"))) 
#Estimate Std. Error z value Pr(>|z|)    
#Control - Dark == 0  0.26795    0.05361   4.998   <0.001 ***
#High - Dark == 0     0.14558    0.05615   2.593   0.0259 *  
#High - Control == 0 -0.12237    0.05361  -2.283   0.0583 . 
lm_PAM_d13.B_treat <- lmer(PAM_d13 ~ Treatment+(1|Genet)+(1|Tank), data = PAM.B)
summary(lm_PAM_d13.B_treat)
summary(glht(lm_PAM_d13.B_treat, linfct = mcp(Treatment = "Tukey")))
posthoc.PAM_d13.B_treat <- glht(lm_PAM_d13.B_treat,
                                linfct = mcp(Treatment = "Tukey"))
cld(posthoc.PAM_d13.B_treat)


#----------Feeding behavior----------
Feeding <- data %>%
  dplyr::select(Sample_ID, Genet, Tank, State, Treatment, Light, Heat,
                Feeding_d2, Feeding_d6, Feeding_d9, Feeding_d13) %>%
  tidyr::pivot_longer(cols = starts_with("Feeding_d"), names_to = "Day", values_to = "Feeding") %>%
  dplyr::mutate(Day = dplyr::case_when(
    Day == "Feeding_d2" ~ "2",
    Day == "Feeding_d6" ~ "6",
    Day == "Feeding_d9" ~ "9",
    Day == "Feeding_d13" ~ "13"),
    Day = factor(Day, levels = c("2", "6", "9", "13")),
    Feeding = as.numeric(as.character(Feeding)),
    Treatment = factor(Treatment,
                       levels = c("NL_C", "CL_C", "HL_C", "NL_H", "CL_H", "HL_H")),
    Heat = factor(Heat, levels = c("Control", "Heat"))
  )

fig_Feeding <- ggplot(Feeding, aes(x = Day, y = Feeding, color = Treatment, group = Treatment,shape = State)) +
  geom_jitter(width = 0.12, height = 0.08, size = 2.5, alpha = 0.8, stroke=1) +
  stat_summary(fun = median, geom = "line", linewidth = 0.8) +
  facet_grid(State ~ Heat) +
  scale_color_manual(name = "Treatment",
                     limits = c("NL_C", "CL_C", "HL_C", "NL_H", "CL_H", "HL_H"),
                     values = c("navy","royalblue","steelblue1","deeppink4","palevioletred1","rosybrown1"),
                     labels = c("Control Temp, Dark", "Control Temp, Control Light", "Control Temp, High Light",
                                "Heated, Dark", "Heated, Control Light", "Heated, High Light")) +
  scale_y_continuous(breaks = 1:5) +
  scale_shape_manual(values = c("White" = 1, "Brown" = 19))+
  labs(x = "Day", y = "Feeding behavior score", color = "Treatment") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_line(linewidth = 1),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right",
    strip.background = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.text = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(color = "gray90")
  )
fig_Feeding
ggsave("../Figures/Feeding.svg", fig_Feeding, width = 7, height = 4.68, units = "in", device = "svg")

Feeding$Feeding <- as.factor(Feeding$Feeding)
Feeding$Treatment <- factor(Feeding$Treatment,
                            levels = c("NL_C", "CL_C", "HL_C", "NL_H", "CL_H", "HL_H"))
Feeding.W <- Feeding %>%
  filter(State == "White")
mod_W_null <- clmm(Feeding ~ Day + (1 | Genet) + (1|Tank),data = Feeding.W)
mod_W_treat <- clmm(Feeding ~ Day + Heat * Light + (1 | Genet) + (1|Tank),
                    data = Feeding.W)
mod_W_treat.2 <- clmm(Feeding ~ Day + Heat + Light + (1 | Genet) + (1|Tank),
                    data = Feeding.W)
AICc(mod_W_null, mod_W_treat, mod_W_treat.2)
#df     AICc
#mod_W_null     9 275.1449
#mod_W_treat   14 253.3140
#mod_W_treat.2 12 257.2327
emmeans(mod_W_treat.2, pairwise ~ Heat) #N.S.
emmeans(mod_W_treat.2, pairwise ~ Light)
#contrast       estimate    SE  df z.ratio p.value
#Dark - Control   -3.178 0.523 Inf  -6.081 <0.0001
#Dark - High      -2.991 0.542 Inf  -5.514 <0.0001
#Control - High    0.187 0.580 Inf   0.323  0.9442
mod_W_treat.group <- clmm(Feeding ~ Day + Treatment + (1 | Genet) + (1|Tank),
                             data = Feeding.W)
emmeans(mod_W_treat.group, pairwise ~ Treatment)
emm_W <- emmeans(mod_W_treat.group, ~ Treatment)
cld(emm_W,adjust = "tukey",sort = FALSE,Letters = letters)

Feeding.B <- Feeding %>%
  filter(State == "Brown")
mod_B_null <- clmm(Feeding ~ Day + (1 | Genet) + (1|Tank),data = Feeding.B)
mod_B_treat <- clmm(Feeding ~ Day + Heat + Light + (1 | Genet) + (1|Tank),
                    data = Feeding.B)
anova(mod_B_null, mod_B_treat)
#no.par    AIC  logLik LR.stat df Pr(>Chisq)    
#mod_B_null       9 255.41 -118.71                          
#mod_B_treat     12 236.77 -106.38  24.646  3  1.831e-05 ***
emmeans(mod_B_treat, pairwise ~ Heat)
#Control - Heat    -1.09 0.51 Inf  -2.139  0.0324 *
emmeans(mod_B_treat, pairwise ~ Light)
#Control - High Light     0.716 0.861 Inf   0.832  0.6833
#Control - No Light       4.557 0.956 Inf   4.769 <0.0001
#High Light - No Light    3.841 0.900 Inf   4.267 <0.0001
mod_B_treat.group <- clmm(Feeding ~ Day + Treatment + (1 | Genet) + (1|Tank),
                          data = Feeding.B)
emmeans(mod_B_treat.group, pairwise ~ Treatment)
emm_B <- emmeans(mod_B_treat.group, ~ Treatment)
cld(emm_B,adjust = "tukey",sort = FALSE,Letters = letters)


#----------Prepare data for PCA---------
df_PCA <- data %>%
  dplyr::select("Sample_ID", "Tank", "Genet",  "State", "Treatment", "Light", "Heat", "Stat_treat")
df_PCA <- df_PCA %>%
  left_join(Sym.den %>% dplyr::select(Sample_ID, Sym_cell_per_area), by = "Sample_ID")%>%
  left_join(Chla_IQR %>% dplyr::select(Sample_ID, Chla_per_area), by = "Sample_ID")%>%
  left_join(Chlc2_IQR %>% dplyr::select(Sample_ID, Chlc2_per_area), by = "Sample_ID")%>%
  left_join(Carb_host_IQR %>% dplyr::select(Sample_ID, Carb_host_per_area), by = "Sample_ID")%>%
  left_join(Carb_sym_IQR %>% dplyr::select(Sample_ID, Carb_sym_per_area), by = "Sample_ID")%>%
  left_join(Protein_host_IQR %>% dplyr::select(Sample_ID, Protein_host_per_area), by = "Sample_ID")%>%
  left_join(Protein_sym_IQR %>% dplyr::select(Sample_ID, Protein_sym_per_area), by = "Sample_ID")%>%
  left_join(PAM %>% dplyr::select(Sample_ID, PAM_d13), by = "Sample_ID")%>%
  left_join(data %>% dplyr::select(Sample_ID, Feeding_d13), by="Sample_ID")
save(df_PCA,file="Phys_df_PCA.Rdata")


#----------Baseline spectroscopic traits----------
Astrangia_symbiotic <- read.csv("../Data/Astrangia_sym_spectro.csv")
Astrangia_Aposymbiotic <- read.csv("../Data/Astrangia_apo_spectro.csv")
dim(Astrangia_symbiotic); summary(Astrangia_symbiotic); head(Astrangia_symbiotic)#Get data information 
dim(Astrangia_Aposymbiotic); summary(Astrangia_Aposymbiotic); head(Astrangia_Aposymbiotic)#Get data information 

#'Data frames must be converted to rspec objects to use in further pavo functions. Wavelenght range needs to be selected for data analysis. The wavelenght range 400-750nm is selected because it is the spectral range of solar radiation that photosynthetic organisms, such as corals, are able to use in the process of photosynthesis.NOTE: The  exact values of the selected wavelenght range must be introduced in the lim argument.
#Convert data frames to rspec objects with specific wl range
Astrangia_symbiotic <-as.rspec(Astrangia_symbiotic,whichwl = 1,interp = FALSE,lim = c(450,750)); 
Astrangia_Aposymbiotic <-as.rspec(Astrangia_Aposymbiotic,whichwl = 1,interp = FALSE,lim = c(450,750))
is.rspec(Astrangia_symbiotic); is.rspec(Astrangia_Aposymbiotic) #Corroborate if data frames are rspec objects


##----------Reflectance----------
#'Reflectance data obtained from spectrophotometer requires further processing before analysis. Processing involves quality control and averaging data.
#'Values between 725 and 750 nm can be used as “quality control” for reflectance spectra. Deviation from 100% of reflectance in the range between 725 and 750 nm highlights the amount of uncontrolled residual scattering of this determination. For a reflectance spectrum to be accepted as correct, the amount of residual scattering must be lower than 20% (Reflectance range between 725-750 nm must be higher than 80%). To visualize and corroborate this condition it is suggested to plot the spectra. Note: reflectance spectrum with more than 20% of residual scattering should not be considered. 
par(mfrow=c(1,2))#Format plot
par(mar=c(4,4,2,2))#Format plot margins
plot(Astrangia_symbiotic,ylab="Reflectance (R)",main="Astrangia symbiotic (R)", ylim=c(0,100)); 
plot(Astrangia_Aposymbiotic,ylab="Reflectance (R)",main="Astrangia aposymbiotic (R)", ylim=c(0,100)) #Plot reflectance spectra

#Averaging spectra
par(mfrow = c(1, 1), cex.axis = 1.4, cex.lab = 1.6, family = "sans", font.lab = 2, mar = c(5.5, 5.5, 2, 2))
aggplot(Astrangia_symbiotic, FUN.center = mean, FUN.error = function(x) sd(x), lcol = "tan4", shadecol = "tan4", xlab = "Wavelength (nm)", ylab = "Reflectance (R)", ylim = c(0, 100), lwd = 2.5)
par(new = TRUE)
aggplot(Astrangia_Aposymbiotic, FUN.center = mean, FUN.error = function(x) sd(x), lcol = "wheat2", shadecol = "wheat3", xlab = "Wavelength (nm)", ylab = "Reflectance (R)", ylim = c(0, 100), lwd = 2.5)


##----------Absorptance----------
#'Estimated absorbance spectra according to Shibata (1969) and Enriquez et al (2005) can be calculated from reflectance data as D=log (1/R). Absorbance along the wavelenght range of 720 − 750 nm reflects primarily the amount of light losses through scattering. Subtracting absorbance at 750 nm to the whole spectrum allows the correction for residual scattering and non-pigmented absorption. In other words, at 750 nm there is no absorption from photosynthetic tissue. Data frame conversion to rspec object allows to efficiently perform perform mathematical operations to the whole spectra for data corrections.
D_Astrangia_symbiotic<-Astrangia_symbiotic; D_Astrangia_symbiotic[,2:6]<-(log10(1/(Astrangia_symbiotic/100)[,2:6]))#Calculate absorbance D=log(1/R)
D_Astrangia_Aposymbiotic<-Astrangia_Aposymbiotic; D_Astrangia_Aposymbiotic[,2:6]<-(log10(1/(Astrangia_Aposymbiotic/100)[,2:6]))#Calculate absorbance D=log(1/R)

par(mfrow = c(1, 1), cex.axis = 1.0, cex.lab = 1.2, family = "sans", font.lab = 2, mar = c(5.5, 5.5, 2, 2))
aggplot(D_Astrangia_symbiotic, FUN.center = mean, FUN.error = function(x) sd(x), lcol = "tan4", shadecol = "tan4", xlab = "Wavelength (nm)", ylab = "log(1/R)", ylim = c(0, 1.3), lwd = 2.5)
par(new = TRUE)
aggplot(D_Astrangia_Aposymbiotic, FUN.center = mean, FUN.error = function(x) sd(x), lcol = "wheat2", shadecol = "wheat3", xlab = "Wavelength (nm)", ylab = "log(1/R)", ylim = c(0, 1.3), lwd = 2.5)

#'Physicists developed the term absorptance (A) to denote the fraction of incident light absorbed by a particular surface. For physiological purposes, absorptance describes the relative amount of solar energy that has the potential to be used in photosynthesis for organic carbon fixation. Absorptance (A) can be derived from transmitted (T) and reflected (R) light measurements, as: A=1-T-R. For organisms with an skeleton, which transmits an insignificant fraction of light (for corals with skeleton thickness >3 mm; see Enríquez et al 2005), absorptance can be calculated as the average value of: A=1-R.
A_Astrangia_symbiotic<-Astrangia_symbiotic; A_Astrangia_symbiotic[,2:6]<-(100-(Astrangia_symbiotic)[,2:6])#Calculate absorbance A=1-R
A_Astrangia_Aposymbiotic<-Astrangia_Aposymbiotic; A_Astrangia_Aposymbiotic[,2:6]<-(100-(Astrangia_Aposymbiotic)[,2:6])#Calculate absorbance A=1-R

par(mfrow = c(1, 1), cex.axis = 1.4, cex.lab = 1.6, family = "sans", font.lab = 2, mar = c(5.5, 5.5, 2, 2))
aggplot(A_Astrangia_symbiotic, FUN.center = mean, FUN.error = function(x) sd(x), lcol = "tan4", shadecol = "tan4", xlab = "Wavelength (nm)", ylab = "Absorptance (x100)", ylim = c(0, 100), lwd = 2.5)
par(new = TRUE)
aggplot(A_Astrangia_Aposymbiotic, FUN.center = mean,FUN.error = function(x) sd(x), lcol = "wheat2", shadecol = "wheat3",xlab = "Wavelength (nm)", ylab = "Absorptance (x100)", ylim = c(0, 100), lwd = 2.5)

#Average absorptance
write_xlsx(A_Astrangia_symbiotic, "../Data/A_Astrangia_symbiotic.xlsx")
write_xlsx(A_Astrangia_Aposymbiotic, "../Data/A_Astrangia_Aposymbiotic.xlsx")
