#Astrangia poculata heat & light: Comparative GE analysis of thermal responses
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Fig. 3C in the manuscript


library(tidyverse)
library(ggplot2)
library(dplyr)


#----------Base plot for camparative analysis----------
#load GO-MWU result files from each dataset for comparative analysis
#A. poculata samples under heated and control light condition in this study
CH_GoBP <- read.csv("../Data/GO_Output/MWU_BP_GO_CH_signedP.csv", sep="", header=TRUE)
#Acropora hyacinthus from Barshis et al. 2013
Barshis_GoBP <- read.csv("../Data/GO_Output/MWU_BP_L_Barshis_bleachResillience_PRJNA177515_For_MWU.csv", sep="", header=TRUE)
#Acropora hyacinthus from Palumbi et al. 2014
Palumbi_GoBP <- read.csv("../Data/GO_Output/MWU_BP_k1_Palumbi_lab_heat_resilience_PRJNA274410_For_MWU.csv", sep="", header=TRUE)
#Acropora millepora from the heat experiment done in Dixon et al. 2020
Dixon2020_GoBP <- read.csv("../Data/GO_Output/MWU_BP_j1_thisStudy_PRJNA559404_For_MWU.csv", sep="", header=TRUE)
#A. poculata from Wuitchik et al. 2024
Wuitchik_GoBP <- read.csv("../Data/GO_Output/MWU_BP_heat_results_modified_pvalues_Astrangia.csv", sep="", header=TRUE)
#Oculina arbuscula from Aichelman et al. 2024
Aichelman_GoBP <- read.csv("../Data/GO_Output/MWU_BP_heat_results_modified_pvalues_Oculina.csv", sep="", header=TRUE)

CH <- CH_GoBP %>% select(name, delta.rank, p.adj) %>% rename(delta.rank_CH = "delta.rank", p.adj_CH = "p.adj") #%>% filter(p.adj_SCH <0.05)
Barshis <- Barshis_GoBP %>% select(name, delta.rank #, p.adj
) %>% rename(delta.rank_Barshis = "delta.rank") #%>% filter(p.adj <0.1)
Palumbi <- Palumbi_GoBP %>% select(name, delta.rank #, p.adj
) %>% rename(delta.rank_Palumbi = "delta.rank") #%>% filter(p.adj<0.1)
Dixon2020 <- Dixon2020_GoBP %>% select(name, delta.rank #, p.adj
) %>% rename(delta.rank_Dixon2020 = "delta.rank") #%>% filter(p.adj<0.1)
Wuitchik <- Wuitchik_GoBP %>% select(name, delta.rank #, p.adj
) %>% rename(delta.rank_Wuitchik = "delta.rank") #%>% filter(p.adj<0.1)
Aichelman <- Aichelman_GoBP %>% select(name, delta.rank #, p.adj
) %>% rename(delta.rank_Aichelman = "delta.rank") #%>% filter(p.adj<0.1)

GO <- list(CH, Barshis, Palumbi, Dixon2020, Wuitchik, Aichelman) %>%
  reduce(full_join, by = "name") %>%
  select(!p.adj_CH) %>%
  rename(CH="delta.rank_CH", Barshis="delta.rank_Barshis",
         Palumbi="delta.rank_Palumbi", Dixon2020="delta.rank_Dixon2020", Wuitchik="delta.rank_Wuitchik", Aichelman="delta.rank_Aichelman")

GO_CH <- GO %>%
  pivot_longer(
    cols = c(Barshis, Palumbi, Dixon2020, Wuitchik, Aichelman),
    names_to = "project",
    values_to = "delta.rank"
  )

project.colors <- c(
  "Barshis" = "#E69F00",      # orange
  "Palumbi" = "#708090",      # light blue
  "Dixon2020" = "#009E73",    # teal green
  "Wuitchik" = "#8FB8DE",     # muted sky blue
  "Aichelman" = "#D4A6C8"     # soft dusty rose
)


CH_comparison <- ggplot(GO_CH, aes(x = CH, y = delta.rank, color = project, text = name)) +
  geom_point(size=2, alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = project.colors) +
  labs(
    x = "Heated, Control Light Delta Rank",
    y = "Coral Project Delta Rank",
    color = "BioProject"
  ) +
  theme_minimal()
CH_comparison


#----------Highlight selected GO terms with mean ± SE----------
CH.shared_up_GoBP <- read.csv("../Data/GO_Output/MWU_BP_CH.shared_up_fisher.csv", sep="", header=TRUE)
SCH.only_up_GoBP <- read.csv("../Data/GO_Output/MWU_BP_SCH.only_up_fisher.csv", sep="", header=TRUE)
CH.shared_up <- CH.shared_up_GoBP %>%
  select(name, p.adj) %>% filter(p.adj < 0.05)
SCH.only_up <- SCH.only_up_GoBP %>%
  select(name, p.adj) %>% filter(p.adj < 0.05)

highlight_terms <- union(CH.shared_up$name, SCH.only_up$name)

highlight_terms_df <- GO_CH %>%
  filter(name %in% highlight_terms) %>%
  mutate(
    group = case_when(
      name %in% CH.shared_up$name ~ "Shared",
      name %in% SCH.only_up$name ~ "BCH Only"
    )
  )
GO.highlight <- highlight_terms_df %>%
  group_by(name, group) %>%
  summarize(
    x = mean(CH, na.rm = TRUE),
    y = mean(delta.rank, na.rm = TRUE),
    se_y = sd(delta.rank, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
highlight.fill.colors <- c(
  "Shared" = "#b98b52",
  "BCH Only"  = "tan4"
)
GO_CH$project <- factor(GO_CH$project, levels = c("Wuitchik", "Aichelman", "Barshis", "Palumbi", "Dixon2020"))

CH_comparison <- CH_comparison +
  geom_errorbar(data = GO.highlight, aes(x = x, ymin = y - se_y, ymax = y + se_y),
                inherit.aes = FALSE, width = 0.2, color = "black") +
  geom_point(
    data = GO.highlight,
    aes(x = x, y = y, fill = group),
    inherit.aes = FALSE,
    shape = 21,
    color = "black",
    size = 4,
    stroke = 0.7,
    alpha = 1
  ) +
  scale_fill_manual(values = highlight.fill.colors, name = "Highlighted Terms")
CH_comparison <- CH_comparison+theme(
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

CH_comparison
svg("../Figures/fig_CH_comparison.svg", width =7.58, height = 4.68)
ggsave("../Figures/fig_CH_comparison.svg", CH_comparison, width = 7.58, height = 4.98, units = "in", device = "svg")

