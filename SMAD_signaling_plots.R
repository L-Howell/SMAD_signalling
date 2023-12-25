#Libraries ####
#libraries ####
library(svDialogs)
library(ggplot2)
library(tidyr)
library(plyr)
library(gganimate)
library(transformr)
library(dplyr)
library(tidyverse)
library(broom)
library(magick)
library(patchwork)
library(lubridate)
library(zoo)
library(extrafont)
library(showtext)
library(rstatix)
library(changepoint)
library(ggsignif)
library(purrr)
library(RColorBrewer)
library(FSA)

#Load font
showtext_auto()
font_add("Arial", "C:/Windows/Fonts/arial.ttf")

#Figure 1 ####
#Panel A ####
fig1a<- read.csv("Fig1B.csv")

fig1a_long <- fig1a %>% 
  pivot_longer(cols = everything(),names_to = "condition", values_to = "Intensity")


# Calculate means and standard deviations
fig1a_long_summary <- fig1a_long %>%
  group_by(condition) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity))

# Convert 'condition' to a factor and specify the order of the levels
fig1a_long_summary$condition <- factor(fig1a_long_summary$condition, levels = c("Un.", "TGF.B", "VACV", "ECTV", "HSV"))

# Create a bar plot with error bars
p<-ggplot(fig1a_long_summary, aes(x = condition, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(), fill="grey50", color="black", size=0.2) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(.9), size=0.2) +
  labs(x = "",   
       y = expression(atop("CAGA"[12] ~ "luciferase activity", "(Normalised RLU)")), 
       title = "") +
  scale_x_discrete(labels = c("Un." = "Untreated",  "TGF.B" = expression("TGF-" * beta), "VACV" = "VACV", "HSV" = "HSV-1", "ECTV" = "ECTV")) +
  scale_y_continuous(
    limits = c(NA, 125),  # Extend beyond 100 for annotations
    breaks = c(0, 25, 50, 75, 100),  # Primary axis setup
    sec.axis = sec_axis(~., breaks = NULL, labels = NULL)  # Invisible secondary axis
  ) +
  theme_minimal()+
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"))+  # Hide y-axis ticks 
  geom_signif(
    comparisons = list(
      c("Un.", "TGF.B"),
      c("Un.", "VACV"),
      c("Un.", "ECTV"),
      c("Un.", "HSV")
    ),
    annotations = c("****", "****", "ns", "****"),
    map_signif_level=TRUE,
    y_position = c(85, 95, 105, 115),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )


ggsave("Fig1A.pdf",width=60, height=70, units="mm", p)


#Panel B ####
fig1b<- read.csv("Fig1B_2.csv")

fig1b_long <- fig1b %>% 
  pivot_longer(cols = everything(),names_to = "condition", values_to = "density")


# Calculate means and standard deviations
fig1b_long_summary <- fig1b_long %>%
  group_by(condition) %>%
  summarise(Mean = mean(density), SD = sd(density))

# Convert 'condition' to a factor and specify the order of the levels
fig1b_long_summary$condition <- factor(fig1b_long_summary$condition, levels = c("Un.", "TGF.B", "VACV", "HSV"))

# Create a bar plot with error bars
p<-ggplot(fig1b_long_summary, aes(x = condition, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(), fill="grey50", color="black", size=0.2) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(.9), size=0.2) +
  labs(x = "",   
       y = "Relative abundance (AU)", 
       title = "") +
  scale_x_discrete(labels = c("Un." = "Untreated",  "TGF.B" = expression("TGF-" * beta), "VACV" = "VACV", "HSV" = "HSV-1", "ECTV" = "ECTV")) +
  scale_y_continuous(
    limits = c(NA, 125),  # Extend beyond 100 for annotations
    breaks = c(0, 25, 50, 75, 100),  # Primary axis setup
    sec.axis = sec_axis(~., breaks = NULL, labels = NULL)  # Invisible secondary axis
  ) +
  theme_minimal()+
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"))+ 
  geom_signif(
    comparisons = list(
      c("Un.", "TGF.B"),
      c("Un.", "VACV"),
      c("Un.", "HSV")
    ),
    annotations = c("***", "***", "***"),
    map_signif_level=TRUE,
    y_position = c(90, 100, 110),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )

ggsave("Fig1B.pdf",width=45, height=70, units="mm", p)

#Panel C ####
fig1c<- read.csv("Fig1C.csv")

# Pivot longer then create new columns for Condition and SB treatment status
fig1c_long <- fig1c %>% 
  pivot_longer(cols = everything(),names_to = "condition", values_to = "Intensity")%>%
  mutate(Treatment = ifelse(grepl("_SB", condition), "+", "-"),
         Condition = gsub("_SB", "", condition),
         Condition = factor(Condition, levels = c("Un.", "TGF.B", "VACV", "HSV")),
         Treatment = factor(Treatment, levels = c("-", "+")))


# Calculate means and standard deviations
fig1c_long_summary <- fig1c_long %>%
  group_by(Condition, Treatment) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity))


# Convert 'condition' to a factor and specify the order of the levels
fig1c_long_summary$Condition <- factor(fig1c_long_summary$Condition, levels = c("Un.", "TGF.B", "VACV", "HSV"))

# Define the levels of Condition and their corresponding y-positions for the annotations
condition_levels <- unique(fig1c_long_summary$Condition)

# Create a data frame for annotations with one row per condition
annotations <- data.frame(
  Condition = condition_levels,
  x = 1,  # Start of the bar (for SB-)
  xend = 2,  # End of the bar (for SB+)
  y = c(100,100,100,100),  # Unique y-position for each condition
  annotation = c("ns", "****", "ns", "****")  # Example annotations, replace with your actual data
)

# Create a bar plot with error bars
p<-ggplot(fig1c_long_summary, aes(x = Treatment, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(), color="black", size=0.2, aes(fill=Treatment)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(width = 0.9), size=0.2) +
  labs(x = "", y = expression(atop("CAGA"[12] ~ "luciferase activity", "(Normalised RLU)")), title = "") +
  scale_fill_manual(values = c("-" = "grey70", "+" = "grey40")) +
  geom_text(aes(label = Treatment, y = Inf), vjust = -4, size = 3) +# Add custom text labels for SB treatment
  facet_wrap(~ Condition, nrow = 1, scales = "free_x", strip.position = "bottom",labeller = labeller(Condition = as_labeller(c("Un." = "Untreated", 
                                                                                                                               "TGF.B" = "TGF-β", 
                                                                                                                               "VACV" = "VACV", 
                                                                                                                               "HSV" = "HSV-1")))) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 8,face="bold"),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        strip.text.x = element_text(size = 5),
        strip.placement = "outside",  # Place strips outside the plot
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,105))

# Add significance annotations for each condition
for (i in 1:nrow(annotations)) {
  annotation_row <- annotations[i, ]
  p <- p + geom_signif(
    data = subset(fig1c_long_summary, Condition == annotation_row$Condition),
    xmin = 1,  # Start of the bar (for SB-)
    xmax = 2,  # End of the bar (for SB+)
    y_position = annotation_row$y,
    annotation = annotation_row$annotation,
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25,
    vjust = -0.5)
}

print(p)

ggsave("Fig1C.pdf",width=80, height=60, units="mm", p)

#Panel D ####
fig1d<- read.csv("Fig1F.csv")
# Pivot longer then create new columns for Condition and SB treatment status
fig1d_long <- fig1d %>% 
  pivot_longer(cols = everything(),names_to = "condition", values_to = "Intensity")%>%
  mutate(Treatment = ifelse(grepl("_SB", condition), "+", "-"),
         Condition = gsub("_SB", "", condition),
         Condition = factor(Condition, levels = c("Un.", "TGF.B", "VACV", "HSV")),
         Treatment = factor(Treatment, levels = c("-", "+")))


# Calculate means and standard deviations
fig1d_long_summary <- fig1d_long %>%
  group_by(Condition, Treatment) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity))


# Convert 'condition' to a factor and specify the order of the levels
fig1d_long_summary$Condition <- factor(fig1d_long_summary$Condition, levels = c("Un.", "TGF.B", "VACV", "HSV"))

# Define the levels of Condition and their corresponding y-positions for the annotations
condition_levels <- unique(fig1d_long_summary$Condition)

# Create a data frame for annotations with one row per condition
annotations <- data.frame(
  Condition = condition_levels,
  x = 1,  # Start of the bar (for SB-)
  xend = 2,  # End of the bar (for SB+)
  y = 135,  # Unique y-position for each condition
  annotation = c("ns", "****", "ns", "****")  # Example annotations, replace with your actual data
)
annotations$Condition <- factor(annotations$Condition, levels = c("Un.", "TGF.B", "VACV", "HSV"))

# Create a bar plot with error bars
p<-ggplot(fig1d_long_summary, aes(x = Treatment, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(), aes(fill=Treatment), color="black", size=0.2) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(width = 0.9), size=0.2) +
  labs(x = "", y = "Mean PAI-1/cell (NFI)", title = "") +
  scale_fill_manual(values = c("-" = "grey70", "+" = "grey40")) +
  geom_text(aes(label = Treatment, y = Inf), vjust = -4, size = 3) +# Add custom text labels for SB treatment
  facet_wrap(~ Condition, nrow = 1, scales = "free_x", strip.position = "bottom",labeller = labeller(Condition = as_labeller(c("Un." = "Untreated", 
                                                                                                                               "TGF.B" = "TGF-β", 
                                                                                                                               "VACV" = "VACV", 
                                                                                                                               "HSV" = "HSV-1")))) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 8,face="bold"),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        strip.text.x = element_text(size = 5),
        strip.placement = "outside",  # Place strips outside the plot
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,140))

# Add significance annotations for each condition
for (i in 1:nrow(annotations)) {
  annotation_row <- annotations[i, ]
  p <- p + geom_signif(
    data = subset(fig1d_long_summary, Condition == annotation_row$Condition),
    xmin = 1,  # Start of the bar (for SB-)
    xmax = 2,  # End of the bar (for SB+)
    y_position = annotation_row$y,
    annotation = annotation_row$annotation,
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25,
    vjust = -0.5)
}

ggsave("Fig1D.pdf",width=80, height=60, units="mm", p)

#Figure 2 ####
#Panel A ####
fig2a<- read.csv("Fig2A.csv")

fig2a_long <- fig2a %>% 
  pivot_longer(cols = everything(),names_to = "condition", values_to = "density")


# Calculate means and standard deviations
fig2a_long_summary <- fig2a_long %>%
  group_by(condition) %>%
  summarise(Mean = mean(density), SD = sd(density))

# Convert 'condition' to a factor and specify the order of the levels
fig2a_long_summary$condition <- factor(fig2a_long_summary$condition, levels = c("HaCaTWT", "Clone1.3", "Clone2.3"))

KO_palette <- c('HaCaTWT' = '#F9A63A', 'Clone1.3' = '#E89F95', 'Clone2.3' = '#EC7B67')
# Create a bar plot with error bars
p<-ggplot(fig2a_long_summary, aes(x = condition, y = Mean, fill=condition)) +
  geom_bar(stat = "identity", position = position_dodge(), color="black", size=0.2) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(.9), size=0.2) +
  scale_fill_manual(values = KO_palette) +
  labs(x = "",   
       y = "Relative SMAD2 abundance (AU)", 
       title = "") +
  scale_x_discrete(labels = c("HaCaTWT" = "WT",  "Clone1.3" = "1.3", "Clone2.3" = "2.3")) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 7),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"))+ 
  geom_signif(
    comparisons = list(
      c("HaCaTWT", "Clone1.3"),
      c("HaCaTWT", "Clone2.3")
    ),
    annotations = c("***", "****"),
    map_signif_level=TRUE,
    y_position = c(11, 12),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )

ggsave("Fig2A.pdf",width=30, height=50, units="mm", p)

#Panel B  ####
fig2b<- read.csv("Fig2B.csv")

fig2b_long <- fig2b %>% 
  pivot_longer(cols = everything(),names_to = "condition", values_to = "density")


# Calculate means and standard deviations
fig2b_long_summary <- fig2b_long %>%
  group_by(condition) %>%
  summarise(Mean = mean(density), SD = sd(density)) %>% 
  filter(condition!="Clone14.5")

# Convert 'condition' to a factor and specify the order of the levels
fig2b_long_summary$condition <- factor(fig2b_long_summary$condition, levels = c("HaCaTWT", "Clone2.11","Clone14.7"))

KO_palette <- c('HaCaTWT' = '#F9A63A', 'Clone2.11' = '#9FD3D3', 'Clone14.7' = '#449596')
# Create a bar plot with error bars
p<-ggplot(fig2b_long_summary, aes(x = condition, y = Mean/10, fill=condition)) +
  geom_bar(stat = "identity", position = position_dodge(), color="black", size=0.2) +
  geom_errorbar(aes(ymin = (Mean/10) - (SD/10), ymax = (Mean/10) + (SD/10)), width = 0.2, position = position_dodge(.9), size=0.2) +
  scale_fill_manual(values = KO_palette) +
  labs(x = "",   
       y = "Relative SMAD3 abundance (AU)", 
       title = "") +
  scale_x_discrete(labels = c("HaCaTWT" = "WT",  "Clone2.11" = "2.11", "Clone14.7"="14.7")) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 7),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text=element_text(family ="Arial"))+ 
  geom_signif(
    comparisons = list(
      c("HaCaTWT", "Clone2.11"),
      c("HaCaTWT", "Clone14.7")
    ),
    annotations = c("****", "****"),
    map_signif_level=TRUE,
    y_position = c(10.5, 11.5, 12.5),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )+
  coord_cartesian(ylim = c(0,12.5))

ggsave("Fig2B.pdf",width=30, height=50, units="mm", p)

#Panel C  ####
fig2c<- read.csv("Fig2C.csv")

fig2c_long <- fig2c %>% # Pivoting longer, renaming, and cleaning up the cell_type column
  pivot_longer(
    cols = -X,
    names_to = "cell_type",
    values_to = "Intensity"
  ) %>%
  mutate(
    cell_type = sub("\\.\\d+$", "", cell_type), # Remove trailing ".number"
    Condition = X # Rename X to Condition
  ) %>%
  select(-X) # Remove the old X column


# Calculate means and standard deviations
fig2c_long_summary <- fig2c_long %>%
  group_by(Condition, cell_type) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity))


# Convert 'condition' to a factor and specify the order of the levels
fig2c_long_summary$Condition <- factor(fig2c_long_summary$Condition, levels = c("Un.", "TGF-β", "VACV-WR", "ECTV","HSV-1"))
fig2c_long_summary$cell_type <- factor(fig2c_long_summary$cell_type, levels = c("HaCaTWT", "HaCaTSMAD2KO", "HaCaTSMAD3KO"))
# Define the levels of Condition and their corresponding y-positions for the annotations
condition_levels <- unique(fig2c_long_summary$Condition)

# Create a data frame for annotations with three rows per condition (for three comparisons)
annotations <- expand.grid(
  Condition = condition_levels,
  Comparison = 1:3
)

# Assign x, xend, y, and annotation based on the comparison
# Update these values based on your actual data and desired annotations
annotations$x <- ifelse(annotations$Comparison == 1, 1, ifelse(annotations$Comparison == 2, 2, 1))
annotations$xend <- ifelse(annotations$Comparison == 1, 2, ifelse(annotations$Comparison == 2, 3, 3))
annotations$y <- c(110, 125, 140)  # Adjust these values as needed

#Signifiance annotations for each level of Condition
condition_annotations <- list(
  "Un." = c("ns", "***", "**"),
  "TGF-β" = c("****", "****", "ns"),
  "VACV-WR" = c("****", "****", "ns"),
  "ECTV" = c("**", "ns", "ns"),
  "HSV-1" = c("****", "****", "ns")
)

# Apply the specific annotations for each condition
annotations$annotation <- mapply(function(cond, comp) condition_annotations[[cond]][comp], 
                                 annotations$Condition, annotations$Comparison)

KO_palette <- c('HaCaTWT' = '#F9A63A', 'HaCaTSMAD2KO' = '#EC7B67', 'HaCaTSMAD3KO' = '#449596')
# Create a bar plot with error bars
p<-ggplot(fig2c_long_summary, aes(x = cell_type, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), aes(fill=cell_type), color="black", size=0.1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(width = 0.9), size=0.2) +
  labs(x = "", y = "ARE luciferase activity \n(Normalised RLU)", title = "",
       fill = "") +
  scale_fill_manual(values = KO_palette,
                    labels = c('HaCaTWT' = 'WT', 'HaCaTSMAD2KO' = 'SMAD2 KO', 'HaCaTSMAD3KO' = 'SMAD3 KO')) +
  geom_text(aes(label = cell_type, y = Inf), vjust = -4, size = 3) +# Add custom text labels for SB treatment
  facet_wrap(~ Condition, nrow = 1, scales = "free_x", strip.position = "bottom",labeller = labeller(Condition = as_labeller(c("Un." = "Untreated", 
                                                                                                                               "TGF-β" = "TGF-β", 
                                                                                                                               "VACV-WR" = "VACV", 
                                                                                                                               "ECTV"="ECTV",
                                                                                                                               "HSV-1" = "HSV-1")))) +
  scale_y_continuous(
    limits = c(NA, 145),  # Extend beyond 100 for annotations
    breaks = c(0, 25, 50, 75, 100),  # Primary axis setup
    sec.axis = sec_axis(~., breaks = NULL, labels = NULL)  # Invisible secondary axis
  ) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        strip.text.x = element_text(size = 5),
        legend.box.margin = margin(t = -20, unit = "pt"),
        #strip.placement = "outside",  # Place strips outside the plot
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,145))

# Add significance annotations for each condition
for (i in 1:nrow(annotations)) {
  annotation_row <- annotations[i, ]
  p <- p + geom_signif(
    data = subset(fig2c_long_summary, Condition == annotation_row$Condition),
    xmin = annotation_row$x,
    xmax = annotation_row$xend,
    y_position = annotation_row$y,
    annotation = annotation_row$annotation,
    tip_length = 0.0125,
    textsize = 1.5,
    size = 0.25,
    vjust = -0.5)
}


ggsave("Fig2C.pdf",width=80, height=60, units="mm", p)

#Panel D  ####
fig2d<- read.csv("Fig2D.csv")

fig2d_long <- fig2d %>% # Pivoting longer, renaming, and cleaning up the cell_type column
  pivot_longer(
    cols = -X,
    names_to = "cell_type",
    values_to = "Intensity"
  ) %>%
  mutate(
    cell_type = sub("\\.\\d+$", "", cell_type), # Remove trailing ".number"
    Condition = X # Rename X to Condition
  ) %>%
  select(-X) # Remove the old X column


# Calculate means and standard deviations
fig2d_long_summary <- fig2d_long %>%
  group_by(Condition, cell_type) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity))


# Convert 'condition' to a factor and specify the order of the levels
fig2d_long_summary$Condition <- factor(fig2d_long_summary$Condition, levels = c("Un.", "TGF-B", "VACV", "ECTV","HSV"))
fig2d_long_summary$cell_type <- factor(fig2d_long_summary$cell_type, levels = c("HaCaTWT", "HaCaTSMAD2KO", "HaCaTSMAD3KO"))
# Define the levels of Condition and their corresponding y-positions for the annotations
condition_levels <- unique(fig2d_long_summary$Condition)

# Create a data frame for annotations with three rows per condition (for three comparisons)
annotations <- expand.grid(
  Condition = condition_levels,
  Comparison = 1:3
)

# Assign x, xend, y, and annotation based on the comparison
# Update these values based on your actual data and desired annotations
annotations$x <- ifelse(annotations$Comparison == 1, 1, ifelse(annotations$Comparison == 2, 2, 1))
annotations$xend <- ifelse(annotations$Comparison == 1, 2, ifelse(annotations$Comparison == 2, 3, 3))
annotations$y <- c(110, 125, 140)  # Adjust these values as needed

#Signifiance annotations for each level of Condition
condition_annotations <- list(
  "Un." = c("ns", "ns", "ns"),
  "TGF-B" = c("ns", "****", "***"),
  "VACV" = c("ns", "****", "****"),
  "ECTV" = c("ns", "ns", "*"),
  "HSV" = c("ns", "****", "****")
)

# Apply the specific annotations for each condition
annotations$annotation <- mapply(function(cond, comp) condition_annotations[[cond]][comp], 
                                 annotations$Condition, annotations$Comparison)

KO_palette <- c('HaCaTWT' = '#F9A63A', 'HaCaTSMAD2KO' = '#EC7B67', 'HaCaTSMAD3KO' = '#449596')
# Create a bar plot with error bars
p<-ggplot(fig2d_long_summary, aes(x = cell_type, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), aes(fill=cell_type), color="black", size=0.1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(width = 0.9), size=0.2) +
  labs(x = "", y = expression(atop("SBR"[6] ~ "luciferase activity", "(Normalised RLU)")), title = "",
       fill = "") +
  scale_fill_manual(values = KO_palette,
                    labels = c('HaCaTWT' = 'WT', 'HaCaTSMAD2KO' = 'SMAD2 KO', 'HaCaTSMAD3KO' = 'SMAD3 KO')) +
  geom_text(aes(label = cell_type, y = Inf), vjust = -4, size = 3) +# Add custom text labels for SB treatment
  facet_wrap(~ Condition, nrow = 1, scales = "free_x", strip.position = "bottom",labeller = labeller(Condition = as_labeller(c("Un." = "Untreated", 
                                                                                                                               "TGF-B" = "TGF-β", 
                                                                                                                               "VACV" = "VACV", 
                                                                                                                               "ECTV"="ECTV",
                                                                                                                               "HSV" = "HSV-1")))) +
  scale_y_continuous(
    limits = c(NA, 145),  # Extend beyond 100 for annotations
    breaks = c(0, 25, 50, 75, 100),  # Primary axis setup
    sec.axis = sec_axis(~., breaks = NULL, labels = NULL)  # Invisible secondary axis
  ) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        strip.text.x = element_text(size = 5),
        legend.box.margin = margin(t = -20, unit = "pt"),
        #strip.placement = "outside",  # Place strips outside the plot
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,145))

# Add significance annotations for each condition
for (i in 1:nrow(annotations)) {
  annotation_row <- annotations[i, ]
  subset_data <- subset(fig2d_long_summary, Condition == annotation_row$Condition)
  
  if (nrow(subset_data) > 0) {
    p <- p + geom_signif(
      data = subset_data,
      xmin = annotation_row$x,
      xmax = annotation_row$xend,
      y_position = annotation_row$y,
      annotation = annotation_row$annotation,
      tip_length = 0.0125,
      textsize = 1.5,
      size = 0.25,
      vjust = -0.5)
  }
}

print(p)
ggsave("Fig2D.pdf",width=80, height=60, units="mm", p)

#Panel E ####
fig2e<- read.csv("Fig2E.csv")

fig2e_long <- fig2e %>% # Pivoting longer, renaming, and cleaning up the cell_type column
  pivot_longer(
    cols = -X,
    names_to = "cell_type",
    values_to = "Intensity"
  ) %>%
  mutate(
    cell_type = sub("\\.\\d+$", "", cell_type), # Remove trailing ".number"
    Condition = X # Rename X to Condition
  ) %>%
  select(-X) # Remove the old X column


# Calculate means and standard deviations
fig2e_long_summary <- fig2e_long %>%
  group_by(Condition, cell_type) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity))


# Convert 'condition' to a factor and specify the order of the levels
fig2e_long_summary$Condition <- factor(fig2e_long_summary$Condition, levels = c("Un.", "TGF-β", "VACV-WR", "ECTV","HSV-1"))
fig2e_long_summary$cell_type <- factor(fig2e_long_summary$cell_type, levels = c("HaCaTWT", "HaCaTSMAD2KO", "HaCaTSMAD3KO"))
# Define the levels of Condition and their corresponding y-positions for the annotations
condition_levels <- unique(fig2e_long_summary$Condition)

# Create a data frame for annotations with three rows per condition (for three comparisons)
annotations <- expand.grid(
  Condition = condition_levels,
  Comparison = 1:3
)

# Assign x, xend, y, and annotation based on the comparison
# Update these values based on your actual data and desired annotations
annotations$x <- ifelse(annotations$Comparison == 1, 1, ifelse(annotations$Comparison == 2, 2, 1))
annotations$xend <- ifelse(annotations$Comparison == 1, 2, ifelse(annotations$Comparison == 2, 3, 3))
annotations$y <- c(110, 125, 140)  # Adjust these values as needed

#Signifiance annotations for each level of Condition
condition_annotations <- list(
  "Un." = c("ns", "ns", "ns"),
  "TGF-β" = c("ns", "****", "****"),
  "VACV-WR" = c("ns", "****", "****"),
  "ECTV" = c("ns", "ns", "ns"),
  "HSV-1" = c("ns", "****", "****")
)

# Apply the specific annotations for each condition
annotations$annotation <- mapply(function(cond, comp) condition_annotations[[cond]][comp], 
                                 annotations$Condition, annotations$Comparison)

KO_palette <- c('HaCaTWT' = '#F9A63A', 'HaCaTSMAD2KO' = '#EC7B67', 'HaCaTSMAD3KO' = '#449596')
# Create a bar plot with error bars
p<-ggplot(fig2e_long_summary, aes(x = cell_type, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), aes(fill=cell_type), color="black", size=0.1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(width = 0.9), size=0.2) +
  labs(x = "", y = expression(atop("CAGA"[12] ~ "luciferase activity", "(Normalised RLU)")), title = "",
       fill = "") +
  scale_fill_manual(values = KO_palette,
                    labels = c('HaCaTWT' = 'WT', 'HaCaTSMAD2KO' = 'SMAD2 KO', 'HaCaTSMAD3KO' = 'SMAD3 KO')) +
  geom_text(aes(label = cell_type, y = Inf), vjust = -4, size = 3) +# Add custom text labels for SB treatment
  facet_wrap(~ Condition, nrow = 1, scales = "free_x", strip.position = "bottom",labeller = labeller(Condition = as_labeller(c("Un." = "Untreated", 
                                                                                                                               "TGF-β" = "TGF-β", 
                                                                                                                               "VACV-WR"="VACV",
                                                                                                                               "ECTV"="ECTV",
                                                                                                                               "HSV-1" = "HSV-1")))) +
  scale_y_continuous(
    limits = c(NA, 145),  # Extend beyond 100 for annotations
    breaks = c(0, 25, 50, 75, 100),  # Primary axis setup
    sec.axis = sec_axis(~., breaks = NULL, labels = NULL)  # Invisible secondary axis
  ) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        strip.text.x = element_text(size = 5),
        legend.box.margin = margin(t = -20, unit = "pt"),
        #strip.placement = "outside",  # Place strips outside the plot
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,145))

# Add significance annotations for each condition
for (i in 1:nrow(annotations)) {
  annotation_row <- annotations[i, ]
  subset_data <- subset(fig2e_long_summary, Condition == annotation_row$Condition)
  
  if (nrow(subset_data) > 0) {
    p <- p + geom_signif(
      data = subset_data,
      xmin = annotation_row$x,
      xmax = annotation_row$xend,
      y_position = annotation_row$y,
      annotation = annotation_row$annotation,
      tip_length = 0.0125,
      textsize = 1.5,
      size = 0.25,
      vjust = -0.5)
  }
}

ggsave("Fig2E.pdf",width=80, height=60, units="mm", p)
#Panel F ####
fig2f<- read.csv("Fig2F.csv")

fig2f_long <- fig2f %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("cell_type", "Condition"),
    names_sep = "_",
    values_to = "Intensity"
  ) %>% 
  # Separate the condition from the cell type prefix
  mutate(cell_type = sub("(.*)_", "", cell_type),
         Condition = sub("(.*)_", "", Condition))

# Calculate means and standard deviations
fig2f_long_summary <- fig2f_long %>%
  group_by(Condition, cell_type) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity))


# Convert 'condition' to a factor and specify the order of the levels
fig2f_long_summary$Condition <- factor(fig2f_long_summary$Condition, levels = c("Un.", "TGF.B", "VACV", "HSV"))
fig2f_long_summary$cell_type <- factor(fig2f_long_summary$cell_type, levels = c("WT", "S2KO", "S3KO"))
# Define the levels of Condition and their corresponding y-positions for the annotations
condition_levels <- unique(fig2f_long_summary$Condition)

# Create a data frame for annotations with three rows per condition (for three comparisons)
annotations <- expand.grid(
  Condition = condition_levels,
  Comparison = 1:3
)

# Assign x, xend, y, and annotation based on the comparison
# Update these values based on your actual data and desired annotations
annotations$x <- ifelse(annotations$Comparison == 1, 1, ifelse(annotations$Comparison == 2, 2, 1))
annotations$xend <- ifelse(annotations$Comparison == 1, 2, ifelse(annotations$Comparison == 2, 3, 3))
annotations$y <- c(110, 125, 140)  # Adjust these values as needed

#Signifiance annotations for each level of Condition
condition_annotations <- list(
  "Un." = c("*", "***", "ns"),
  "TGF.B" = c("ns", "****", "****"),
  "VACV" = c("ns", "****", "****"),
  "HSV" = c("***", "****", "****")
)

# Apply the specific annotations for each condition
annotations$annotation <- mapply(function(cond, comp) condition_annotations[[cond]][comp], 
                                 annotations$Condition, annotations$Comparison)

KO_palette <- c('WT' = '#F9A63A', 'S2KO' = '#EC7B67', 'S3KO' = '#449596')
# Create a bar plot with error bars
p<-ggplot(fig2f_long_summary, aes(x = cell_type, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), aes(fill=cell_type), color="black", size=0.1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(width = 0.9), size=0.2) +
  labs(x = "", y = "Mean PAI-1/cell (NFI)", title = "",
       fill = "") +
  scale_fill_manual(values = KO_palette,
                    labels = c('WT' = 'WT', 'S2KO' = 'SMAD2 KO', 'S3KO' = 'SMAD3 KO')) +
  geom_text(aes(label = cell_type, y = Inf), vjust = -4, size = 3) +# Add custom text labels for SB treatment
  facet_wrap(~ Condition, nrow = 1, scales = "free_x", strip.position = "bottom",labeller = labeller(Condition = as_labeller(c("Un." = "Untreated", 
                                                                                                                               "TGF.B" = "TGF-β", 
                                                                                                                               "VACV" = "VACV", 
                                                                                                                               "ECTV"="ECTV",
                                                                                                                               "HSV" = "HSV-1")))) +
  scale_y_continuous(
    limits = c(NA, 145),  # Extend beyond 100 for annotations
    breaks = c(0, 25, 50, 75, 100),  # Primary axis setup
    sec.axis = sec_axis(~., breaks = NULL, labels = NULL)  # Invisible secondary axis
  ) +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        strip.text.x = element_text(size = 5),
        legend.box.margin = margin(t = -20, unit = "pt"),
        #strip.placement = "outside",  # Place strips outside the plot
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,145))

# Add significance annotations for each condition
for (i in 1:nrow(annotations)) {
  annotation_row <- annotations[i, ]
  subset_data <- subset(fig2f_long_summary, Condition == annotation_row$Condition)
  
  if (nrow(subset_data) > 0) {
    p <- p + geom_signif(
      data = subset_data,
      xmin = annotation_row$x,
      xmax = annotation_row$xend,
      y_position = annotation_row$y,
      annotation = annotation_row$annotation,
      tip_length = 0.0125,
      textsize = 1.5,
      size = 0.25,
      vjust = -0.5)
  }
}

ggsave("Fig2F.pdf",width=80, height=60, units="mm", p)

#Figure 3 ####
#Panel A ####
fig3a<- read.csv("Fig3A.csv")

fig3a_long <- fig3a %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Condition", "cell_type"), # Specify the new column names
    names_sep = "_", # Define the separator used in the column names
    values_to = "Intensity"
  ) %>%
  mutate(
    # If further cleaning of cell_type is needed, apply it here
    cell_type = sub("\\.\\d+$", "", cell_type)
  )


# Calculate means and standard deviations
fig3a_long_summary <- fig3a_long %>%
  group_by(Condition, cell_type) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity))


# Convert 'condition' to a factor and specify the order of the levels
fig3a_long_summary$Condition <- factor(fig3a_long_summary$Condition, levels = c("Un.", "TGF.B", "VACV"))
fig3a_long_summary$cell_type <- factor(fig3a_long_summary$cell_type, levels = c("noRes", "SMAD2SSXS", "SMAD2Res"))
# Define the levels of Condition and their corresponding y-positions for the annotations
condition_levels <- unique(fig3a_long_summary$Condition)

# Create a data frame for annotations with three rows per condition (for three comparisons)
annotations <- expand.grid(
  Condition = condition_levels,
  Comparison = 1:3
)

# Assign x and xend based on the comparison
annotations$x <- ifelse(annotations$Comparison == 1, 1, ifelse(annotations$Comparison == 2, 2, 1))
annotations$xend <- ifelse(annotations$Comparison == 1, 2, ifelse(annotations$Comparison == 2, 3, 3))

# Assign y based on Comparison
y_values <- c(105, 115, 125)  # Specify the y values for each comparison
annotations$y <- y_values[annotations$Comparison] # Adjust these values as needed

#Signifiance annotations for each level of Condition
condition_annotations <- list(
  "Un." = c("ns", "****", "****"),
  "TGF.B" = c("****", "****", "****"),
  "VACV" = c("****", "**", "****")
)

# Apply the specific annotations for each condition
annotations$annotation <- mapply(function(cond, comp) condition_annotations[[cond]][comp], 
                                 annotations$Condition, annotations$Comparison)

SMAD2rescue_palette <- c('noRes' = '#EC7B67', 'SMAD2SSXS' = '#E0A39B', 'SMAD2Res' = '#D6C9C8')

# Create a bar plot with error bars
p<-ggplot(fig3a_long_summary, aes(x = cell_type, y = Mean, fill=cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), color="black", size=0.1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(width = 0.9), size=0.2) +
  labs(x = "", y = "ARE luciferase activity \n(Normalised RLU)", title = "",
       fill = "") +
  scale_fill_manual(values = SMAD2rescue_palette,
                    labels = c(
                      "noRes" = "Control",
                      "SMAD2SSXS" = expression("+SMAD2"^paste(Delta, "SSXS")),
                      "SMAD2Res" = expression("+SMAD2")
                    )) +
  geom_text(aes(label = cell_type, y = Inf), vjust = -4, size = 3) +# Add custom text labels for SB treatment
  facet_wrap(~ Condition, nrow = 1, scales = "free_x", strip.position = "bottom",labeller = labeller(Condition = as_labeller(c("Un." = "Untreated", 
                                                                                                                               "TGF.B" = "TGF-β", 
                                                                                                                               "VACV" = "VACV")))) +
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        strip.text.x = element_text(size = 5),
        legend.box.margin = margin(t = -20, unit = "pt"),
        #strip.placement = "outside",  # Place strips outside the plot
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,130))

# Add significance annotations for each condition
for (i in 1:nrow(annotations)) {
  annotation_row <- annotations[i, ]
  subset_data <- subset(fig3a_long_summary, Condition == annotation_row$Condition)
  
  if (nrow(subset_data) > 0) {
    p <- p + geom_signif(
      data = subset_data,
      xmin = annotation_row$x,
      xmax = annotation_row$xend,
      y_position = annotation_row$y,
      annotation = annotation_row$annotation,
      tip_length = 0.0125,
      textsize = 1.5,
      size = 0.25,
      vjust = -0.5)
  }
}

ggsave("Fig3A.pdf",width=85, height=60, units="mm", p)

#Panel B ####
fig3b<- read.csv("Fig3B.csv")

fig3b_long <- fig3b %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Condition", "cell_type"), # Specify the new column names
    names_sep = "_", # Define the separator used in the column names
    values_to = "Intensity"
  ) %>%
  mutate(
    # If further cleaning of cell_type is needed, apply it here
    cell_type = sub("\\.\\d+$", "", cell_type)
  )


# Calculate means and standard deviations
fig3b_long_summary <- fig3b_long %>%
  group_by(Condition, cell_type) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity))


# Convert 'condition' to a factor and specify the order of the levels
fig3b_long_summary$Condition <- factor(fig3b_long_summary$Condition, levels = c("Un.", "TGF.B", "VACV"))
fig3b_long_summary$cell_type <- factor(fig3b_long_summary$cell_type, levels = c("noRes", "SMAD3SSXS", "SMAD3Res"))
# Define the levels of Condition and their corresponding y-positions for the annotations
condition_levels <- unique(fig3b_long_summary$Condition)

# Create a data frame for annotations with three rows per condition (for three comparisons)
annotations <- expand.grid(
  Condition = condition_levels,
  Comparison = 1:3
)

# Assign x and xend based on the comparison
annotations$x <- ifelse(annotations$Comparison == 1, 1, ifelse(annotations$Comparison == 2, 2, 1))
annotations$xend <- ifelse(annotations$Comparison == 1, 2, ifelse(annotations$Comparison == 2, 3, 3))

# Assign y based on Comparison
y_values <- c(75, 85, 95)  # Specify the y values for each comparison
annotations$y <- y_values[annotations$Comparison] # Adjust these values as needed

#Signifiance annotations for each level of Condition
condition_annotations <- list(
  "Un." = c("ns", "ns", "ns"),
  "TGF.B" = c("*", "**", "***"),
  "VACV" = c("**", "***", "****")
)

# Apply the specific annotations for each condition
annotations$annotation <- mapply(function(cond, comp) condition_annotations[[cond]][comp], 
                                 annotations$Condition, annotations$Comparison)

SMAD3rescue_palette <- c('noRes' = '#449596', 'SMAD3SSXS' = '#78BABA', 'SMAD3Res' = '#AFD8D8')
# Create a bar plot with error bars
p<-ggplot(fig3b_long_summary, aes(x = cell_type, y = Mean, fill=cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), color="black", size=0.1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(width = 0.9), size=0.2) +
  labs(x = "", y = expression(atop("CAGA"[12] ~ "luciferase activity", "(Normalised RLU)")), title = "",
       fill = "") +
  scale_fill_manual(values = SMAD3rescue_palette,
                    labels = c(
                      "noRes" = "Control",
                      "SMAD3SSXS" = expression("+SMAD3"^paste(Delta, "SSXS")),
                      "SMAD3Res" = expression("+SMAD3")
                    ))  +
  geom_text(aes(label = cell_type, y = Inf), vjust = -4, size = 3) +# Add custom text labels for SB treatment
  facet_wrap(~ Condition, nrow = 1, scales = "free_x", strip.position = "bottom",labeller = labeller(Condition = as_labeller(c("Un." = "Untreated", 
                                                                                                                               "TGF.B" = "TGF-β", 
                                                                                                                               "VACV" = "VACV")))) +
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        strip.text.x = element_text(size = 5),
        legend.box.margin = margin(t = -20, unit = "pt"),
        #strip.placement = "outside",  # Place strips outside the plot
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,95))

# Add significance annotations for each condition
for (i in 1:nrow(annotations)) {
  annotation_row <- annotations[i, ]
  subset_data <- subset(fig3b_long_summary, Condition == annotation_row$Condition)
  
  if (nrow(subset_data) > 0) {
    p <- p + geom_signif(
      data = subset_data,
      xmin = annotation_row$x,
      xmax = annotation_row$xend,
      y_position = annotation_row$y,
      annotation = annotation_row$annotation,
      tip_length = 0.0125,
      textsize = 1.5,
      size = 0.25,
      vjust = -0.5)
  }
}

ggsave("Fig3B.pdf",width=85, height=60, units="mm", p)
#Figure 4 ####
#Panel B ####
fig4b<- read.csv("Fig4B.csv")

fig4b_long <- fig4b %>% 
  pivot_longer(cols = -X,names_to = "cell_type", values_to = "count") %>% 
  mutate(
    cell_type = sub("\\.\\d+$", "", cell_type), # Remove trailing ".number"
    condition = X) %>%  # Rename X to Condition
  dplyr::select(-X) %>% 
  filter(cell_type!="Smad2Hyp")


# Calculate means and standard deviations
fig4b_long_summary <- fig4b_long %>%
  group_by(condition, cell_type) %>%
  summarise(Mean = mean(count), SD = sd(count))

# Convert 'condition' to a factor and specify the order of the levels
fig4b_long_summary$condition <- factor(fig4b_long_summary$condition, levels = c("VACV", "ECTV", "HSV"))
fig4b_long_summary$cell_type <- factor(fig4b_long_summary$cell_type, levels = c("WT", "Smad2KO", "Smad3KO"))
# Define the levels of Condition and their corresponding y-positions for the annotations
condition_levels <- unique(fig4b_long_summary$condition)

# Create a data frame for annotations with three rows per condition (for three comparisons)
annotations <- expand.grid(
  condition = condition_levels,
  Comparison = 1:3
)

# Assign x and xend based on the comparison
annotations$x <- ifelse(annotations$Comparison == 1, 1, ifelse(annotations$Comparison == 2, 2, 1))
annotations$xend <- ifelse(annotations$Comparison == 1, 2, ifelse(annotations$Comparison == 2, 3, 3))

# Assign y based on Comparison
y_values <- c(145, 160, 175)  # Specify the y values for each comparison
annotations$y <- y_values[annotations$Comparison] # Adjust these values as needed

#Signifiance annotations for each level of Condition
condition_annotations <- list(
  "VACV" = c("***", "*", "ns"),
  "ECTV" = c("*", "***", "*"),
  "HSV" = c("****", "ns", "*")
)

# Apply the specific annotations for each condition
annotations$annotation <- mapply(function(cond, comp) condition_annotations[[cond]][comp], 
                                 annotations$condition, annotations$Comparison)

KO_palette <- c('WT' = '#F9A63A', 'Smad2KO' = '#EC7B67', 'Smad3KO' = '#449596')

# Create a bar plot with error bars
p<-ggplot(fig4b_long_summary, aes(x = cell_type, y = Mean, fill=cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), color="black", size=0.15) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(width = 0.9), size=0.2) +
  labs(x = "", y = "Mean plaque count", title = "",
       fill = "") +
  scale_fill_manual(values = KO_palette,
                    labels = c('WT' = 'WT', 'Smad2KO' = 'SMAD2KO', 'Smad3KO' = 'SMAD3KO')) +
  geom_text(aes(label = cell_type, y = Inf), vjust = -4, size = 3) +# Add custom text labels for SB treatment
  facet_wrap(~ condition, nrow = 1, scales = "free_x", strip.position = "bottom",labeller = labeller(condition = as_labeller(c("VACV" = "VACV", 
                                                                                                                               "ECTV" = "ECTV", 
                                                                                                                               "HSV" = "HSV-1")))) +
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        strip.text.x = element_text(size = 5),
        legend.box.margin = margin(t = -20, unit = "pt"),
        #strip.placement = "outside",  # Place strips outside the plot
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,178))

# Add significance annotations for each condition
for (i in 1:nrow(annotations)) {
  annotation_row <- annotations[i, ]
  subset_data <- subset(fig4b_long_summary, condition == annotation_row$condition)
  
  if (nrow(subset_data) > 0) {
    p <- p + geom_signif(
      data = subset_data,
      xmin = annotation_row$x,
      xmax = annotation_row$xend,
      y_position = annotation_row$y,
      annotation = annotation_row$annotation,
      tip_length = 0.0125,
      textsize = 1.5,
      size = 0.25,
      vjust = -0.5)
  }
}

print(p)
ggsave("Fig4B.pdf",width=75, height=60, units="mm", p)

#Panel C ####
fig4c<- read.csv("Fig4C.csv")

fig4c_long <- fig4c %>% 
  pivot_longer(cols = everything(),names_to = "cell_type", values_to = "area") %>% 
  dplyr::filter(cell_type!="SMAD2Hyp")


# Calculate means and standard deviations
fig4c_long_summary <- fig4c_long %>%
  group_by(cell_type) %>%
  summarise(Mean = mean(area), SD = sd(area))

# Convert 'condition' to a factor and specify the order of the levels
fig4c_long_summary$cell_type <- factor(fig4c_long_summary$cell_type, levels = c("HaCaTWT", "HaCaTSMAD2KO", "HaCaTSMAD3KO"))
fig4c_long$cell_type <- factor(fig4c_long$cell_type, levels = c("HaCaTWT", "HaCaTSMAD2KO", "HaCaTSMAD3KO"))

KO_palette <- c('HaCaTWT' = '#F9A63A', 'HaCaTSMAD2KO' = '#EC7B67', 'HaCaTSMAD3KO' = '#449596')

mean_sd <- function(x) {
  return(c(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x)))
}
# Create a bar plot with error bars
p <- ggplot(fig4c_long, aes(x = cell_type, y = area)) +
  geom_jitter(stat = "identity", width = 0.3, aes(fill = cell_type), size=1, color="black",stroke=0.15, shape = 21) +
  scale_fill_manual(values = KO_palette) +
  labs(x = "",   
       y = expression("Plaque area (mm"^2*")"), 
       title = "") +
  scale_x_discrete(labels = c("HaCaTWT" = "WT",  "HaCaTSMAD2KO" = "SMAD2KO", "HaCaTSMAD3KO" = "SMAD3KO")) +
  scale_y_continuous(
    limits = c(NA, 20),  # Extend beyond 100 for annotations
    breaks = c(0, 5, 10, 15),  # Primary axis setup
    sec.axis = sec_axis(~., breaks = NULL, labels = NULL)  # Invisible secondary axis
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial")) +
  stat_summary(
    fun.data = mean_sd, 
    geom = "errorbar", 
    position = position_dodge(width = 0.3), 
    width = 0.2,
    color = "black",
    alpha=0.35,
    size=0.25
  ) +
  stat_summary(
    fun = mean, 
    geom = "crossbar", 
    width = 0.8,
    size = 0.125, # Adjust size to match the line width in your plot
    position = position_dodge(width = 0.3)
  )+
  geom_signif(
    comparisons = list(
      c("HaCaTWT", "HaCaTSMAD2KO"),
      c("HaCaTWT", "HaCaTSMAD3KO"),
      c("HaCaTSMAD2KO", "HaCaTSMAD3KO")
    ),
    annotations = c("****", "****", "*****"),
    map_signif_level=TRUE,
    y_position = c(14.5, 16.5, 18.5),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )

ggsave("Fig4C.pdf",width=57, height=55, units="mm", p)

#Panel D ####
fig4d<- read.csv("Fig4D.csv")

fig4d_long <- fig4d %>% 
  pivot_longer(cols = everything(),names_to = "cell_type", values_to = "area") %>% 
  dplyr::filter(cell_type!="Hyp")


# Calculate means and standard deviations
fig4d_long_summary <- fig4d_long %>%
  group_by(cell_type) %>%
  summarise(Mean = mean(area), SD = sd(area))

# Convert 'condition' to a factor and specify the order of the levels
fig4d_long_summary$cell_type <- factor(fig4d_long_summary$cell_type, levels = c("WT", "SMAD2KO", "SMAD3KO"))
fig4d_long$cell_type <- factor(fig4d_long$cell_type, levels = c("WT", "SMAD2KO", "SMAD3KO"))

KO_palette <- c('WT' = '#F9A63A', 'SMAD2KO' = '#EC7B67', 'SMAD3KO' = '#449596')

mean_sd <- function(x) {
  return(c(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x)))
}
# Create a bar plot with error bars
p <- ggplot(fig4d_long, aes(x = cell_type, y = area)) +
  geom_jitter(stat = "identity", width = 0.3, aes(fill = cell_type), size=1, color="black",stroke=0.15, shape = 21) +
  scale_fill_manual(values = KO_palette) +
  labs(x = "",   
       y = expression("Plaque area (mm"^2*")"), 
       title = "") +
  scale_x_discrete(labels = c("WT" = "WT",  "SMAD2KO" = "SMAD2KO", "SMAD3KO" = "SMAD3KO")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial")) +
  stat_summary(
    fun.data = mean_sd, 
    geom = "errorbar", 
    position = position_dodge(width = 0.3), 
    width = 0.2,
    color = "black",
    alpha=0.35,
    size=0.25
  ) +
  stat_summary(
    fun = mean, 
    geom = "crossbar", 
    width = 0.8,
    size = 0.125, # Adjust size to match the line width in your plot
    position = position_dodge(width = 0.3)
  )+
  geom_signif(
    comparisons = list(
      c("WT", "SMAD2KO"),
      c("WT", "SMAD3KO"),
      c("SMAD2KO", "SMAD3KO")
    ),
    annotations = c("ns", "****", "*"),
    map_signif_level=TRUE,
    y_position = c(6.5, 7.5, 8.5),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )

ggsave("Fig4D.pdf",width=57, height=55, units="mm", p)

#Panel E ####
fig4e<- read.csv("Fig4E.csv")

fig4e_long <- fig4e %>% 
  pivot_longer(cols = everything(),names_to = "cell_type", values_to = "area") %>% 
  dplyr::filter(cell_type!="Hyp")


# Calculate means and standard deviations
fig4e_long_summary <- fig4e_long %>%
  group_by(cell_type) %>%
  summarise(Mean = mean(area), SD = sd(area))

# Convert 'condition' to a factor and specify the order of the levels
fig4e_long_summary$cell_type <- factor(fig4e_long_summary$cell_type, levels = c("WT", "SMAD2KO", "SMAD3KO"))
fig4e_long$cell_type <- factor(fig4e_long$cell_type, levels = c("WT", "SMAD2KO", "SMAD3KO"))

KO_palette <- c('WT' = '#F9A63A', 'SMAD2KO' = '#EC7B67', 'SMAD3KO' = '#449596')

mean_sd <- function(x) {
  return(c(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x)))
}
# Create a bar plot with error bars
p <- ggplot(fig4e_long, aes(x = cell_type, y = area)) +
  geom_jitter(stat = "identity", width = 0.3, aes(fill = cell_type), size=1, color="black",stroke=0.15, shape = 21) +
  scale_fill_manual(values = KO_palette) +
  labs(x = "",   
       y = expression("Plaque area (mm"^2*")"), 
       title = "") +
  scale_x_discrete(labels = c("WT" = "WT",  "SMAD2KO" = "SMAD2KO", "SMAD3KO" = "SMAD3KO")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial")) +
  stat_summary(
    fun.data = mean_sd, 
    geom = "errorbar", 
    position = position_dodge(width = 0.3), 
    width = 0.2,
    color = "black",
    alpha=0.35,
    size=0.25
  ) +
  stat_summary(
    fun = mean, 
    geom = "crossbar", 
    width = 0.8,
    size = 0.125, # Adjust size to match the line width in your plot
    position = position_dodge(width = 0.3)
  )+
  geom_signif(
    comparisons = list(
      c("WT", "SMAD2KO"),
      c("WT", "SMAD3KO"),
      c("SMAD2KO", "SMAD3KO")
    ),
    annotations = c("****", "ns", "****"),
    map_signif_level=TRUE,
    y_position = c(17, 19, 21),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )
print(p)
ggsave("Fig4E.pdf",width=57, height=55, units="mm", p)

#Figure 5 ####
#Panel A ####
fig5a<- read.csv("Fig5A.csv")

fig5a_long <- fig5a %>% # Pivoting longer, renaming, and cleaning up the cell_type column
  pivot_longer(
    cols = -Time,
    names_to = "cell_type",
    values_to = "titre"
  ) %>%
  mutate(
    cell_type = sub("\\.\\d+$", "", cell_type), # Remove trailing ".number"
  )

# Calculate means and standard deviations
fig5a_long_summary <- fig5a_long %>%
  group_by(cell_type, Time) %>%
  summarise(Mean = mean(titre), 
            SD = sd(titre),
            SEM = SD / sqrt(n()))  # Calculate SEM

fig5a_long_summary$cell_type <- factor(fig5a_long_summary$cell_type, levels = c("HaCaTWT", "HaCaTSMAD2KO", "HaCaTSMAD3KO"))
fig5a_long$cell_type <- factor(fig5a_long$cell_type, levels = c("HaCaTWT", "HaCaTSMAD2KO", "HaCaTSMAD3KO"))

KO_palette <- c('HaCaTWT' = '#F9A63A', 'HaCaTSMAD2KO' = '#EC7B67', 'HaCaTSMAD3KO' = '#449596')

p <- ggplot(fig5a_long_summary, aes(x = Time, y = Mean, group = cell_type, color = cell_type)) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM, color = cell_type), width = 2, position = position_dodge(0.2))+
  scale_color_manual(values=KO_palette,
                     labels = c("HaCaTWT" = "WT",  "HaCaTSMAD2KO" = "SMAD2KO", "HaCaTSMAD3KO" = "SMAD3KO"))+
  geom_path(size=0.5) +
  geom_point(aes(shape=cell_type), show.legend=FALSE)+
  theme_minimal() +
  labs(x = "Time (HPI)", y = "Titre (PFU/mL)", color = "")+
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial"))+
  coord_cartesian(ylim=c(0,6900000))

ggsave("Fig5A.pdf",width=85, height=55, units="mm", p)

#Panel C ####

#SMAD2KO
center_x <- 1215
center_y <- 1147
x_range <- 1600
y_range <- 1600
p<- data_uniqueID %>%
  filter(grepl("20221004_13", unique_ID)) %>%   # Controls which plaque is displayed
  ggplot(aes(x = AreaShape_Center_Y, y = AreaShape_Center_X, group = unique_ID, colour = HPI)) +
  geom_path(size = 0.2, lineend = "round") +
  scale_color_gradientn(colours = c("#D43031", "#F0BA1A", "#FEFEFE", "#252973"), na.value = "transparent") +
  labs(x = "", y = "",  # Remove axis labels
       title = "", colour = "Time \n(HPI)") +
  theme(legend.position = "none",
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        panel.background = element_rect(fill = 'black', colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(center_y - x_range/2, center_y + x_range/2), ylim = c(center_x - y_range/2, center_x + y_range/2)) +
  scale_x_continuous(breaks = NULL, labels = NULL) +
  scale_y_continuous(breaks = NULL, labels = NULL)

ggsave("Fig5B_SMAD2KO.pdf",width=54.5, height=60, units="mm", p)
#SMAD2KO
p <-  data_uniqueID_means_per_time %>% 
  filter(Metadata_Condition=="SMAD2KO") %>%   
  ggplot(aes(x=HPI, y=mean_integrateddistance, group = plaque_ID)) +
  geom_path(size=0.25, color="#EC7B67") +
  scale_x_continuous()+
  labs(x = "Time (HPI)", 
       y = expression(paste("Mean total distance traveled (",mu,"m)")),
       title = "", colour = "cell type")+
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial"))+
  coord_cartesian(ylim=c(0,70))

ggsave("Fig5B_SMAD2KO_dist.pdf",width=50, height=40, units="mm", p)
  
#SMAD3KO
center_x <- 1166.5
center_y <- 1145.5
x_range <- 1600
y_range <- 1600
p<-data_uniqueID %>%
  filter(grepl("20221004_22", unique_ID)) %>%  # Controls which plaque is displayed
  ggplot(aes(x = AreaShape_Center_Y, y = AreaShape_Center_X, group = unique_ID, colour = HPI)) +
  geom_path(size = 0.2, lineend = "round") +
  scale_color_gradientn(colours = c("#D43031", "#F0BA1A", "#FEFEFE", "#252973"), na.value = "transparent") +
  labs(x = "", y = "",  # Remove axis labels
       title = "", colour = "Time \n(HPI)") +
  theme(legend.position = "none",
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        panel.background = element_rect(fill = 'black', colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(center_y - x_range/2, center_y + x_range/2), ylim = c(center_x - y_range/2, center_x + y_range/2)) +
  scale_x_continuous(breaks = NULL, labels = NULL) +
  scale_y_continuous(breaks = NULL, labels = NULL)

ggsave("Fig5B_SMAD3KO.pdf",width=54.5, height=60, units="mm", p)
#SMAD3KO
p <-  data_uniqueID_means_per_time %>% 
  filter(Metadata_Condition=="SMAD3KO") %>%   
  ggplot(aes(x=HPI, y=mean_integrateddistance, group = plaque_ID)) +
  geom_path(size=0.25, color="#449596") +
  scale_x_continuous()+
  labs(x = "Time (HPI)", 
       y = expression(paste("Mean total distance traveled (",mu,"m)")),
       title = "", colour = "cell type")+
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial"))+
  coord_cartesian(ylim=c(0,70))

ggsave("Fig5B_SMAD3KO_dist.pdf",width=50, height=40, units="mm", p)
#WT
center_x <- 1197
center_y <- 870
x_range <- 1600
y_range <- 1600
p<-data_uniqueID %>%
  filter(grepl("20221004_20", unique_ID)) %>%  # Controls which plaque is displayed
  ggplot(aes(x = AreaShape_Center_Y, y = AreaShape_Center_X, group = unique_ID, colour = HPI)) +
  geom_path(size = 0.2, lineend = "round") +
  scale_color_gradientn(colours = c("#D43031", "#F0BA1A", "#FEFEFE", "#252973"), na.value = "transparent") +
  labs(x = "", y = "",  # Remove axis labels
       title = "", colour = "Time \n(HPI)") +
  theme(legend.position = "none",
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        panel.background = element_rect(fill = 'black', colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(center_y - x_range/2, center_y + x_range/2), ylim = c(center_x - y_range/2, center_x + y_range/2)) +
  scale_x_continuous(breaks = NULL, labels = NULL) +
  scale_y_continuous(breaks = NULL, labels = NULL)

ggsave("Fig5B_WT.pdf",width=54.5, height=60, units="mm", p)
#WT
p <-  data_uniqueID_means_per_time %>% 
  filter(Metadata_Condition=="WT") %>%   
  ggplot(aes(x=HPI, y=mean_integrateddistance, group = plaque_ID)) +
  geom_path(size=0.25, color="#F9A63A") +
  scale_x_continuous()+
  labs(x = "Time (HPI)", 
       y = expression(paste("Mean total distance traveled (",mu,"m)")),
       title = "", colour = "cell type")+
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial"))+
  coord_cartesian(ylim=c(0,70))

ggsave("Fig5B_WT_dist.pdf",width=50, height=40, units="mm", p)

#Panel D ####
data_uniqueID_CP_metrics_summary<- data_uniqueID_CP_metrics %>% 
  group_by(Metadata_Condition) %>% 
  summarise(mean=mean(TrackObjects_IntegratedDistance, na.rm=TRUE),
            sd=sd(TrackObjects_IntegratedDistance, na.rm=TRUE))

data_uniqueID_CP_metrics$Metadata_Condition <- factor(data_uniqueID_CP_metrics$Metadata_Condition, levels = c("WT", "SMAD2KO", "SMAD3KO"))

# Define your color palette
KO_palette <- c('WT' = '#F9A63A', 'SMAD2KO' = '#EC7B67', 'SMAD3KO' = '#449596')

p <- ggplot(data_uniqueID_CP_metrics, aes(y = TrackObjects_IntegratedDistance, x = Metadata_Condition, fill = Metadata_Condition)) +
  geom_violin(size = 0.15) +
  scale_fill_manual(values = KO_palette) +
  labs(title = "", x = "", y = "Total distance traveled (µm)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial"))+
  geom_signif(
    comparisons = list(
      c("WT", "SMAD2KO"),
      c("WT", "SMAD3KO"),
      c("SMAD2KO", "SMAD3KO")
    ),
    annotations = c("****", "ns", "****"),
    map_signif_level=TRUE,
    y_position = c(300, 330, 360),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )

# Print the plot
ggsave("Fig5C.pdf",width=50, height=65, units="mm", p)

#Panel E ####

directionality_df_summary<- directionality_df %>% 
  group_by(Metadata_Condition.y) %>% 
  summarise(mean=mean(directionality, na.rm=TRUE),
            sd=sd(directionality, na.rm=TRUE))

directionality_df$Metadata_Condition.y <- factor(directionality_df$Metadata_Condition.y, levels = c("WT", "SMAD2KO", "SMAD3KO"))

p <- ggplot(directionality_df, aes(y = directionality, x = Metadata_Condition.y, fill = Metadata_Condition.y)) +
  geom_violin(size = 0.15) +
  scale_fill_manual(values = KO_palette) +
  labs(title = "", x = "", y = "Radial migration bias") +
  scale_y_continuous(
    limits = c(NA, 1.6),  # Extend beyond 100 for annotations
    breaks = c(-1, -0.5, 0, 0.5, 1),  # Primary axis setup
    sec.axis = sec_axis(~., breaks = NULL, labels = NULL)  # Invisible secondary axis
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial"))+
  geom_signif(
    comparisons = list(
      c("WT", "SMAD2KO"),
      c("WT", "SMAD3KO"),
      c("SMAD2KO", "SMAD3KO")
    ),
    annotations = c("****", "ns", "****"),
    map_signif_level=TRUE,
    y_position = c(1, 1.2, 1.4),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )

# Print the plot
ggsave("Fig5D.pdf",width=50, height=65, units="mm", p)

#Panel F ####
data_uniqueID_per_cell_summary<- data_uniqueID_per_cell %>% 
  group_by(Metadata_Condition) %>% 
  summarise(mean=mean(radial_velocity, na.rm=TRUE),
            sd=sd(radial_velocity, na.rm=TRUE))

#Plot average total distance traveled per cell
data_uniqueID_per_cell$Metadata_Condition <- factor(data_uniqueID_per_cell$Metadata_Condition, levels = c("WT", "SMAD2KO", "SMAD3KO"))

# Plot
p <- ggplot(data_uniqueID_per_cell, aes(y = radial_velocity, x = Metadata_Condition, fill = Metadata_Condition)) +
  geom_violin(size = 0.15) +
  scale_fill_manual(values = KO_palette) +
  labs(title = "", x = "", y = "Radial velocity (µm/hr)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 5),
        legend.key.height = unit(0.25, "cm"),
        text = element_text(family = "Arial"))+
  geom_signif(
    comparisons = list(
      c("WT", "SMAD2KO"),
      c("WT", "SMAD3KO"),
      c("SMAD2KO", "SMAD3KO")
    ),
    annotations = c("****", "ns", "****"),
    map_signif_level=TRUE,
    y_position = c(18, 21, 24),
    tip_length = 0.0125,
    textsize = 2,
    size = 0.25 
  )

# Print the plot
ggsave("Fig5E.pdf",width=50, height=65, units="mm", p)



#Supplementary figures ####
