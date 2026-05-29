#Astrangia poculata heat & light: Bubble plots for GO term comparisons
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Fig. 4C & 5C in the manuscript


library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)


#----------Heated, High Light----------
BP_HH.shared_up<-read.csv("../Data/GO_Output/MWU_BP_HH.shared_up_fisher.csv", sep = "")
BP_HH.shared_up <- BP_HH.shared_up %>%
  filter(p.adj <= 0.05)
BP_HH.shared_up$category <- "BP"
deg2go_HH.shared_up<-read.csv("../Data/GO_Output/BP_HH.shared_up_fisher.csv", sep="")
deg2go_HH.shared_up <- subset(deg2go_HH.shared_up, value==1)
deg2go_HH.shared_up <- deg2go_HH.shared_up %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_HH.shared_up)
BP_HH.shared_up <- left_join(BP_HH.shared_up, deg2go_HH.shared_up, by = "term")
identical(BP_HH.shared_up$name.x, BP_HH.shared_up$name.y)  #TRUE
BP_HH.shared_up <- BP_HH.shared_up %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

BP_HH.shared_down<-read.csv("../Data/GO_Output/MWU_BP_HH.shared_down_fisher.csv", sep = "")
BP_HH.shared_down<- BP_HH.shared_down %>%
  filter(p.adj <= 0.05) #0
BP_HH.shared_down$DEG_overlap <- integer(0)

BP_AHH.only_up<-read.csv("../Data/GO_Output/MWU_BP_AHH.only_up_fisher.csv", sep = "")
BP_AHH.only_up <- BP_AHH.only_up %>%
  filter(p.adj <= 0.05)
BP_AHH.only_up$category <- "BP"
deg2go_AHH.only_up<-read.csv("../Data/GO_Output/BP_AHH.only_up_fisher.csv", sep="")
deg2go_AHH.only_up <- subset(deg2go_AHH.only_up, value==1)
deg2go_AHH.only_up <- deg2go_AHH.only_up %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_AHH.only_up)
BP_AHH.only_up <- left_join(BP_AHH.only_up, deg2go_AHH.only_up, by = "term")
identical(BP_AHH.only_up$name.x, BP_AHH.only_up$name.y)  #TRUE
BP_AHH.only_up <- BP_AHH.only_up %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

BP_AHH.only_down<-read.csv("../Data/GO_Output/MWU_BP_AHH.only_down_fisher.csv", sep = "")
BP_AHH.only_down <- BP_AHH.only_down %>%
  filter(p.adj <= 0.05)
BP_AHH.only_down$category <- "BP"
deg2go_AHH.only_down<-read.csv("../Data/GO_Output/BP_AHH.only_down_fisher.csv", sep="")
deg2go_AHH.only_down <- subset(deg2go_AHH.only_down, value==1)
deg2go_AHH.only_down <- deg2go_AHH.only_down %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_AHH.only_down)
BP_AHH.only_down <- left_join(BP_AHH.only_down, deg2go_AHH.only_down, by = "term")
identical(BP_AHH.only_down$name.x, BP_AHH.only_down$name.y)  #TRUE
BP_AHH.only_down <- BP_AHH.only_down %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

BP_SHH.only_up<-read.csv("../Data/GO_Output/MWU_BP_SHH.only_up_fisher.csv", sep = "")
BP_SHH.only_up <- BP_SHH.only_up %>%
  filter(p.adj <= 0.05)
BP_SHH.only_up$category <- "BP"
deg2go_SHH.only_up<-read.csv("../Data/GO_Output/BP_SHH.only_up_fisher.csv", sep="")
deg2go_SHH.only_up <- subset(deg2go_SHH.only_up, value==1)
deg2go_SHH.only_up <- deg2go_SHH.only_up %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_SHH.only_up)
BP_SHH.only_up <- left_join(BP_SHH.only_up, deg2go_SHH.only_up, by = "term")
identical(BP_SHH.only_up$name.x, BP_SHH.only_up$name.y)  #TRUE
BP_SHH.only_up <- BP_SHH.only_up %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

BP_SHH.only_down<-read.csv("../Data/GO_Output/MWU_BP_SHH.only_down_fisher.csv", sep = "")
BP_SHH.only_down <- BP_SHH.only_down %>%
  filter(p.adj <= 0.05)
BP_SHH.only_down$category <- "BP"
deg2go_SHH.only_down<-read.csv("../Data/GO_Output/BP_SHH.only_down_fisher.csv", sep="")
deg2go_SHH.only_down <- subset(deg2go_SHH.only_down, value==1)
deg2go_SHH.only_down <- deg2go_SHH.only_down %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_SHH.only_down)
BP_SHH.only_down <- left_join(BP_SHH.only_down, deg2go_SHH.only_down, by = "term")
identical(BP_SHH.only_down$name.x, BP_SHH.only_down$name.y)  #TRUE
BP_SHH.only_down <- BP_SHH.only_down %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

all_go_terms.HH <- unique(c(
  BP_HH.shared_up$name,
  BP_HH.shared_down$name,
  BP_AHH.only_up$name,
  BP_AHH.only_down$name,
  BP_SHH.only_up$name,
  BP_SHH.only_down$name
))

all_sources.HH <- c("Shared_up", "Shared_down",
                    "White_up", "White_down",
                    "Brown_up", "Brown_down")
full_grid.HH <- expand.grid(name = all_go_terms.HH, source = all_sources.HH, stringsAsFactors = FALSE)

add_source <- function(df, label) {
  df %>%
    filter(p.adj <= 0.05) %>%
    mutate(source = label) %>%
    select(name, DEG_overlap, p.adj, source)
}

BP_combined.HH <- bind_rows(
  add_source(BP_HH.shared_up, "Shared_up"),
  add_source(BP_HH.shared_down, "Shared_down"),
  add_source(BP_AHH.only_up, "White_up"),
  add_source(BP_AHH.only_down, "White_down"),
  add_source(BP_SHH.only_up, "Brown_up"),
  add_source(BP_SHH.only_down, "Brown_down")
)

dummy_row.HH <- data.frame(
  name = "No significant GO terms",
  DEG_overlap = 0,
  p.adj = 1,
  source = "HH.shared.down",
  direction = "Downregulated"
)
BP_combined.HH <- bind_rows(BP_combined.HH, dummy_row.HH)

BP_complete.HH <- full_grid.HH %>%
  left_join(BP_combined.HH, by = c("name", "source"))

BP_complete.HH <- BP_complete.HH %>%
  mutate(direction = ifelse(grepl("up", source, ignore.case = TRUE),
                            "Upregulated", "Downregulated"),
         direction = factor(direction, levels = c("Upregulated", "Downregulated")))

#Manually pick GO terms related to cell cycle, immunity, metabolism, and transport

HH_terms <- c("immune system process","cell death","defense response","immune response","cell cycle","I-kappaB kinase/NF-kappaB signaling","regulation of mitotic cell cycle","regulation of cell cycle process","cell cycle process",
              "regulation of protein metabolic process","negative regulation of protein metabolic process",
              "positive regulation of protein metabolic process","regulation of cell death","negative regulation of cell cycle","negative regulation of cell cycle phase transition",
              "regulation of defense response to virus","regulation of cell cycle","response to cell cycle checkpoint signaling",
              "regulation of cell cycle phase transition","regulation of mitotic cell cycle phase transition","negative regulation of mitotic cell cycle phase transition",
              "mitotic cell cycle process","regulation of cellular protein catabolic process","cellular protein metabolic process","proteolysis involved in cellular protein catabolic process","cell cycle checkpoint signaling",
              "protein metabolic process","regulation of protein catabolic process","positive regulation of cell cycle",
              "positive regulation of cell cycle process","polysaccharide catabolic process","polysaccharide metabolic process","aminoglycan metabolic process","aminoglycan catabolic process","amino sugar metabolic process","peptide biosynthetic process","peptide metabolic process","lipid transport","carbohydrate catabolic process",
              "amide biosynthetic process","amino sugar catabolic process","carbohydrate derivative catabolic process")
HH_pattern <- paste(HH_terms, collapse = "|")
BP_plot.HH <- BP_complete.HH %>%
  filter(grepl(HH_pattern, name, ignore.case = TRUE))%>%
  filter(!name %in%  c("meiotic cell cycle process", "meiotic cell cycle", "meiosis I cell cycle process"))

BP_plot.HH <- BP_plot.HH %>%
  mutate(group = case_when(
    grepl("cell cycle", name, ignore.case = TRUE) ~ "Cell Cycle",
    grepl("immune system process|cell death|defense response|immune response|I-kappaB kinase/NF-kappaB signaling|regulation of cell death|regulation of defense response to virus", name, ignore.case = TRUE) ~ "Immunity",
    grepl("polysaccharide|carbohydrate", name, ignore.case = TRUE) ~ "Carbohydrate",
    grepl("aminoglycan|amino sugar|peptide|amide", name, ignore.case =TRUE) ~ "Nitrogen",
    grepl("protein", name, ignore.case = TRUE) ~ "Protein",
    grepl("transport", name, ignore.case = TRUE) ~ "Transport",
    TRUE ~ "Other"
  ))

BP_plot.HH$source <- factor(BP_plot.HH$source, levels = c(
  "Shared_up", "Shared_down",
  "White_up", "White_down",
  "Brown_up", "Brown_down"
))
BP_plot.HH$group <- factor(BP_plot.HH$group, levels = c("Transport","Protein","Nitrogen", "Carbohydrate", "Immunity", "Cell Cycle"))
BP_plot.HH <- BP_plot.HH %>%
  arrange(group, name) %>%
  mutate(name_ordered = factor(name, levels = unique(name)))
BP_plot.HH <- BP_plot.HH %>%
  mutate(
    end = word(name, -2, -1)
  ) %>%
  arrange(group, end, name) %>%
  mutate(name_ordered = factor(name, levels = rev(unique(name))))
fig_HH.bubble <-ggplot(BP_plot.HH, aes(x = source, y = forcats::fct_rev(name_ordered), size = DEG_overlap, color = -log10(p.adj))) +
  geom_point(alpha = 0.8, na.rm = TRUE) +
  scale_color_gradient(low = "lightblue", high = "darkblue", na.value = "transparent") +
  scale_size_continuous(
    range = c(0, 15),
    limits = c(0, 80),
    breaks = c(20, 40, 60, 80),
    name = "DEG Count"
  ) +
  facet_wrap(~direction, scales = "free_x") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(size=10,angle = 45, hjust = 1,face="bold"),
    axis.text.y = element_text(size=10),    
    strip.background = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.text = element_text(size=14,face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(color = "gray90"),
    axis.title = element_blank(),
    legend.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=8),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
fig_HH.bubble
ggsave("../Figures/HH_bubble.svg", fig_HH.bubble, width = 8.4, height = 6.18, units = "in", device = "svg")


#----------Heated, Dark----------
BP_DH.shared_up<-read.csv("../Data/GO_Output/MWU_BP_DH.shared_up_fisher.csv", sep = "")
BP_DH.shared_up <- BP_DH.shared_up %>%
  filter(p.adj <= 0.05)
BP_DH.shared_up$category <- "BP"
deg2go_DH.shared_up<-read.csv("../Data/GO_Output/BP_DH.shared_up_fisher.csv", sep="")
deg2go_DH.shared_up <- subset(deg2go_DH.shared_up, value==1)
deg2go_DH.shared_up <- deg2go_DH.shared_up %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_DH.shared_up)
BP_DH.shared_up <- left_join(BP_DH.shared_up, deg2go_DH.shared_up, by = "term")
identical(BP_DH.shared_up$name.x, BP_DH.shared_up$name.y)  #TRUE
BP_DH.shared_up <- BP_DH.shared_up %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

BP_DH.shared_down<-read.csv("../Data/GO_Output/MWU_BP_DH.shared_down_fisher.csv", sep = "")
BP_DH.shared_down<- BP_DH.shared_down %>%
  filter(p.adj <= 0.05)
BP_DH.shared_down$category <- "BP"
deg2go_DH.shared_down<-read.csv("../Data/GO_Output/BP_DH.shared_down_fisher.csv", sep="")
deg2go_DH.shared_down <- subset(deg2go_DH.shared_down, value==1)
deg2go_DH.shared_down <- deg2go_DH.shared_down %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_DH.shared_down)
BP_DH.shared_down <- left_join(BP_DH.shared_down, deg2go_DH.shared_down, by = "term")
identical(BP_DH.shared_down$name.x, BP_DH.shared_down$name.y)  #TRUE
BP_DH.shared_down <- BP_DH.shared_down %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

BP_ADH.only_up<-read.csv("../Data/GO_Output/MWU_BP_ADH.only_up_fisher.csv", sep = "")
BP_ADH.only_up <- BP_ADH.only_up %>%
  filter(p.adj <= 0.05)
BP_ADH.only_up$category <- "BP"
deg2go_ADH.only_up<-read.csv("../Data/GO_Output/BP_ADH.only_up_fisher.csv", sep="")
deg2go_ADH.only_up <- subset(deg2go_ADH.only_up, value==1)
deg2go_ADH.only_up <- deg2go_ADH.only_up %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_ADH.only_up)
BP_ADH.only_up <- left_join(BP_ADH.only_up, deg2go_ADH.only_up, by = "term")
identical(BP_ADH.only_up$name.x, BP_ADH.only_up$name.y)  #TRUE
BP_ADH.only_up <- BP_ADH.only_up %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

BP_ADH.only_down<-read.csv("../Data/GO_Output/MWU_BP_ADH.only_down_fisher.csv", sep = "")
BP_ADH.only_down <- BP_ADH.only_down %>%
  filter(p.adj <= 0.05)
BP_ADH.only_down$category <- "BP"
deg2go_ADH.only_down<-read.csv("../Data/GO_Output/BP_ADH.only_down_fisher.csv", sep="")
deg2go_ADH.only_down <- subset(deg2go_ADH.only_down, value==1)
deg2go_ADH.only_down <- deg2go_ADH.only_down %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_ADH.only_down)
BP_ADH.only_down <- left_join(BP_ADH.only_down, deg2go_ADH.only_down, by = "term")
identical(BP_ADH.only_down$name.x, BP_ADH.only_down$name.y)  #TRUE
BP_ADH.only_down <- BP_ADH.only_down %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

BP_SDH.only_up<-read.csv("../Data/GO_Output/MWU_BP_SDH.only_up_fisher.csv", sep = "")
BP_SDH.only_up <- BP_SDH.only_up %>%
  filter(p.adj <= 0.05) #0
BP_SDH.only_up$DEG_overlap <- integer(0)

BP_SDH.only_down<-read.csv("../Data/GO_Output/MWU_BP_SDH.only_down_fisher.csv", sep = "")
BP_SDH.only_down <- BP_SDH.only_down %>%
  filter(p.adj <= 0.05)
BP_SDH.only_down$category <- "BP"
deg2go_SDH.only_down<-read.csv("../Data/GO_Output/BP_SDH.only_down_fisher.csv", sep="")
deg2go_SDH.only_down <- subset(deg2go_SDH.only_down, value==1)
deg2go_SDH.only_down <- deg2go_SDH.only_down %>%
  group_by(term, name) %>%
  summarise(DEG_overlap = n(), .groups = "drop")
head(deg2go_SDH.only_down)
BP_SDH.only_down <- left_join(BP_SDH.only_down, deg2go_SDH.only_down, by = "term")
identical(BP_SDH.only_down$name.x, BP_SDH.only_down$name.y)  #TRUE
BP_SDH.only_down <- BP_SDH.only_down %>%
  select(-name.y) %>%
  dplyr::rename(name = name.x)

all_go_terms.DH <- unique(c(
  BP_DH.shared_up$name,
  BP_DH.shared_down$name,
  BP_ADH.only_up$name,
  BP_ADH.only_down$name,
  BP_SDH.only_up$name,
  BP_SDH.only_down$name
))

all_sources.DH <- c("Shared_up", "Shared_down",
                    "White_up", "White_down",
                    "Brown_up", "Brown_down")
full_grid.DH <- expand.grid(name = all_go_terms.DH, source = all_sources.DH, stringsAsFactors = FALSE)

add_source <- function(df, label) {
  df %>%
    filter(p.adj <= 0.05) %>%
    mutate(source = label) %>%
    select(name, DEG_overlap, p.adj, source)
}

# Merge all treatments
BP_combined.DH <- bind_rows(
  add_source(BP_DH.shared_up, "Shared_up"),
  add_source(BP_DH.shared_down, "Shared_down"),
  add_source(BP_ADH.only_up, "White_up"),
  add_source(BP_ADH.only_down, "White_down"),
  add_source(BP_SDH.only_up, "Brown_up"),
  add_source(BP_SDH.only_down, "Brown_down")
)

dummy_row.DH <- data.frame(
  name = "No significant GO terms",
  DEG_overlap = 0,
  p.adj = 1,
  source = "HH.shared.down",
  direction = "Downregulated"
)
BP_combined.DH <- bind_rows(BP_combined.DH, dummy_row.DH)

BP_complete.DH <- full_grid.DH %>%
  left_join(BP_combined.DH, by = c("name", "source"))

BP_complete.DH <- BP_complete.DH %>%
  mutate(direction = ifelse(grepl("up", source, ignore.case = TRUE),
                            "Upregulated", "Downregulated"),
         direction = factor(direction, levels = c("Upregulated", "Downregulated")))
#Manually pick GO terms related to cell cycle, immunity, metabolism, and transport

DH_terms <- c("immune system process","regulation of immune system process","negative regulation of immune system process",
              "immune response-regulating signaling pathway","cellular protein metabolic process","cell death","defense response","inflammatory response","immune response","protein metabolic process","regulation of protein metabolic process",
              "negative regulation of protein metabolic process","positive regulation of protein metabolic process","defense response to protozoan","regulation of cell death","innate immune response",
              "defense response to virus","lipid metabolic process","steroid metabolic process","steroid metabolic process","sterol transport","pattern recognition receptor signaling pathway","toll-like receptor signaling pathway",
              "immune system process","regulation of immune system process","immune response-regulating signaling pathway","peptide biosynthetic process","cellular protein metabolic process","peptide metabolic process","cellular amino acid metabolic process","glutamine metabolic process","cell death",
              "defense response","immune response","cell cycle","regulation of mitotic cell cycle","glutamine family amino acid metabolic process","protein metabolic process","regulation of protein metabolic process","negative regulation of protein metabolic process",
              "regulation of cell death","positive regulation of cell death","negative regulation of cell death","regulation of I-kappaB kinase/NF-kappaB signaling","positive regulation of I-kappaB kinase/NF-kappaB signaling","cellular amide metabolic process","amide biosynthetic process","innate immune response","regulation of defense response to virus",
              "regulation of immune response","regulation of cell cycle","sulfur amino acid metabolic process","peptide metabolic process","cellular modified amino acid metabolic process","amino acid transport","regulation of oligopeptide transport","cellular modified amino acid catabolic process","peptide catabolic process","cellular amide metabolic process","prostaglandin biosynthetic process",
              "lipid metabolic process","fatty acid biosynthetic process","unsaturated fatty acid biosynthetic process","prostaglandin metabolic process","lipid biosynthetic process"
)
DH_pattern <- paste(DH_terms, collapse = "|")
BP_plot.DH <- BP_complete.DH %>%
  filter(grepl(DH_pattern, name, ignore.case = TRUE))%>%
  filter(!name %in%  c("meiotic cell cycle process", "meiotic cell cycle", "meiosis I cell cycle process","cytoplasmic pattern recognition receptor signaling pathway"))

BP_plot.DH <- BP_plot.DH %>%
  mutate(group = case_when(
    grepl("cell cycle", name, ignore.case = TRUE) ~ "Cell Cycle",
    grepl("immune system process|regulation of immune system process|negative regulation of immune system process|immune response-regulating signaling pathway|cell death|defense response|inflammatory response|immune response|defense response to protozoan|regulation of cell death|innate immune response|defense response to virus|pattern recognition receptor signaling pathway|toll-like receptor signaling pathway|immune system process|regulation of immune system process|immune response-regulating signaling pathway|regulation of I-kappaB kinase/NF-kappaB signaling|positive regulation of I-kappaB kinase/NF-kappaB signaling|innate immune response|regulation of defense response to virus|regulation of immune response", name, ignore.case = TRUE) ~ "Immunity",
    grepl("lipid|steroid metabolic process|steroid metabolic process|prostaglandin|fatty acid", name, ignore.case = TRUE) ~ "Carbohydrate",
    grepl("peptide metabolic|peptide catabolic|peptide biosynthetic|amino acid metabolic|amino acid catabolic|glutamine|amide", name, ignore.case =TRUE) ~ "Nitrogen",
    grepl("protein", name, ignore.case = TRUE) ~ "Protein",
    grepl("transport", name, ignore.case = TRUE) ~ "Transport",
    TRUE ~ "Other"
  ))

BP_plot.DH$source <- factor(BP_plot.DH$source, levels = c(
  "Shared_up", "Shared_down",
  "White_up", "White_down",
  "Brown_up", "Brown_down"
))
BP_plot.DH$group <- factor(BP_plot.DH$group, levels = c("Transport","Protein","Nitrogen", "Carbohydrate", "Immunity", "Cell Cycle"))
BP_plot.DH <- BP_plot.DH %>%
  arrange(group, name) %>%
  mutate(name_ordered = factor(name, levels = unique(name)))
BP_plot.DH <- BP_plot.DH %>%
  mutate(
    end = word(name, -2, -1)
  ) %>%
  arrange(group, end, name) %>%
  mutate(name_ordered = factor(name, levels = rev(unique(name))))
fig_DH.bubble <-ggplot(BP_plot.DH, aes(x = source, y = forcats::fct_rev(name_ordered), size = DEG_overlap, color = -log10(p.adj))) +
  geom_point(alpha = 0.8, na.rm = TRUE) +
  scale_color_gradient(low = "lightblue", high = "darkblue", na.value = "transparent") +
  scale_size_continuous(
    range = c(0, 15),
    limits = c(0, 80),
    breaks = c(20, 40, 60, 80),
    name = "DEG Count"
  ) +
  facet_wrap(~direction, scales = "free_x") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(size=10,angle = 45, hjust = 1,face="bold"),
    axis.text.y = element_text(size=10),    
    strip.background = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.text = element_text(size=14,face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(color = "gray90"),
    axis.title = element_blank(),
    legend.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=8),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
fig_DH.bubble
ggsave("../Figures/DH_bubble.svg", fig_DH.bubble, width = 9.4, height = 7.2, units = "in", device = "svg")
