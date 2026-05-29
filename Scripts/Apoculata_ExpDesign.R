#Astrangia poculata heat & light: Experimental design
#Light environment shapes divergent coral thermal responses across symbiotic states
#This script was used to generate Fig. 1B in the manuscript

library(ggplot2)
library(dplyr)


Tank.temp <- read.csv("../Data/Water_Temp_tanks.csv", header=TRUE)
Tank.temp$Tank <- as.factor(Tank.temp$Tank)
Tank.temp$Treatment <- as.factor(Tank.temp$Treatment)
str(Tank.temp)

#calculate mean and standard deviation
temp_summary <- Tank.temp %>%
  group_by(Day, Treatment) %>%
  summarise(mean_temp = mean(Temp, na.rm = TRUE),
            sd_temp = sd(Temp, na.rm = TRUE),
            target_temp = mean(Target_Temp, na.rm = TRUE),
            .groups = "drop")
#create summary table with target temps for each treatment
target.temp <- temp_summary %>%
  dplyr::select(Day, Treatment, target_temp) %>%
  distinct()

#create dummy data for sampling events
sampling <- data.frame(x = c(2, 5.95, 9, 12.95, 6.05, 10, 13.05, 14),
                       type = c("Feeding behavior", "Feeding behavior", "Feeding behavior", "Feeding behavior",
                                "PAM", "PAM", "PAM",                                                             
                                "RNA and other phys"))                                                            

ExpDesign <- ggplot(Tank.temp, aes(x = Day, y = Temp, color = Treatment)) +
  geom_vline(data = sampling,
             aes(xintercept = x, linetype = type),
             color = "black",
             alpha = 0.3,
             linewidth = 0.8,
             show.legend = TRUE) +
  scale_linetype_manual(name = "Sampling events",
                        values = c("Feeding behavior" = "dotted",
                                   "PAM" = "longdash",
                                   "RNA and other phys" = "solid")) +
  geom_step(data = target.temp,
            aes(x = Day, y = target_temp, group = Treatment, color = Treatment),
            alpha = 0.2, linetype = "solid", linewidth = 3, inherit.aes = FALSE) +
  geom_jitter(alpha = 0.5, size = 2, width = 0.2, height = 0) +
  geom_errorbar(data = temp_summary,
                aes(x = Day, y = mean_temp, ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp),
                width = 0, linewidth = 0.8, color = "black", inherit.aes = FALSE) +
  geom_point(data = temp_summary,
             aes(x = Day, y = mean_temp, fill = Treatment),
             size = 3, shape = 21, color = "black", stroke = 1, inherit.aes = FALSE) +
  scale_x_continuous(breaks = 1:14, limits = c(1, 14)) +
  scale_y_continuous(breaks = seq(10, 35, 5), minor_breaks = seq(10, 35, 1)) +
  scale_color_manual(values = c("Control_Temp" = "royalblue", "Heated" = "palevioletred1"),
                     labels = c("Control_Temp" = "control", "Heated" = "heated")) +
  scale_fill_manual(values = c("Control_Temp" = "royalblue", "Heated" = "palevioletred1"),
                    labels = c("Control_Temp" = "control", "Heated" = "heated")) +
  labs(x = "Experimental Day", y = "Temperature (°C)", color = "Treatment") +
  theme_bw(base_size = 14)

ExpDesign <-ExpDesign+ theme(
  plot.title=element_blank(),
  axis.ticks=element_line(linewidth=1),
  axis.title = element_text(size=14, face = "bold"),
  axis.text=element_text(size=12),
  legend.title = element_text(size=14, face = "bold"),
  legend.text = element_text(size=12),
  panel.grid.major = element_line(color = scales::alpha("grey", 0.1)),
  panel.grid.minor = element_line(color = scales::alpha("grey", 0.1)),
  strip.text = element_text(size =12, face = "bold"),
  strip.background = element_rect(fill = "grey90", color = "black")
)
ExpDesign
ggsave("../Figures/ExpDesign.svg", ExpDesign, width = 10, height = 6, units = "in", device = "svg")