################################################################################
#
# Evaluate ensemble model outputs
#
###############################################################################

rm(list = ls())

source("scripts/util/__Util__MASTER.R")
file_name <- "Non50-50Mixes_SupplFigure"

####################
# Function to process and plot
####################
# Function to process data
process_mix_data <- function(task_distribution_matrix) {
  # Process raw data
  task_dist <- task_dist %>% 
    group_by(n) %>% 
    mutate(set = paste0(Mix, "-", replicate)) %>% 
    mutate(set = factor(set, 
                        levels = c(sort(apply(expand.grid(unique(Mix), unique(replicate)), 1, paste, collapse = "-")))))
  
  # Calculate mix means
  task_dist_summary <- task_dist %>% 
    group_by(Mix, set) %>% 
    summarise(Task1 = mean(Task1),
              Task2 = mean(Task2)) %>% 
    group_by(Mix) %>% 
    summarise(Task1_mean = mean(Task1),
              Task1_SD = sd(Task1),
              Task1_SE = sd(Task1) / sqrt(length(Task1)),
              Task2_mean = mean(Task2),
              Task2_SD = sd(Task2),
              Task1_SE = sd(Task2) / sqrt(length(Task2))) %>% 
    filter(!Mix %in% c(0, 1)) %>% 
    mutate(Line = "Mixed",
           Group_mean = TRUE)
  
  # Calculate Mix X Line means
  task_dist_lines <- task_dist %>% 
    group_by(set, Mix, Line) %>% 
    summarise(Task1 = mean(Task1),
              Task2 = mean(Task2)) %>% 
    group_by(Mix, Line) %>% 
    summarise(Task1_mean = mean(Task1),
              Task1_SD = sd(Task1),
              Task1_SE = sd(Task1) / sqrt(length(Task1)),
              Task2_mean = mean(Task2),
              Task2_SD = sd(Task2),
              Task1_SE = sd(Task2) / sqrt(length(Task2))) %>% 
    mutate(Group_mean = Mix %in% c(0, 1))
  
  #Bind
  task_dist_summary <- task_dist_summary %>% 
    bind_rows(task_dist_lines) %>% 
    mutate(Group_mean = as.factor(Group_mean))
  
  # Null hypothesis lines
  null_hypothesis <- task_dist_summary %>% 
    filter(Mix %in% c(0, 1))
  
  # Make list and return
  return_item <- list()
  return_item[[1]] <- task_dist_summary
  return_item[[2]] <- null_hypothesis
  return(return_item)
}

# Function to plot data
plot_mix_data <- function(task_distribution_data, null_hypothesis_data) {
  gg_dist_sum <- ggplot(data = task_distribution_data, aes(y = Task1_mean, x = Mix, color = Line)) +
    geom_line(data = null_hypothesis_data, 
              aes(y = Task1_mean, x = Mix, group = NA),
              color = "grey60",
              linetype = "dashed",
              size = 0.3) +
    geom_errorbar(aes(ymin = Task1_mean - Task1_SE, ymax = Task1_mean + Task1_SE), 
                  position = position_dodge(width = 0.05),
                  width = 0.035,
                  size = 0.2) +
    geom_point(aes(size = Group_mean),
               position = position_dodge(width = 0.05)) +
    theme_classic() +
    labs(x = "Fraction of A individuals in colony",
         y = "Frequency Task 1, mean \u00B1 s.e.") +
    scale_color_manual(values = c("#E52521", "#2B4B9B", "#7C217F")) +
    scale_size_manual(values = c(0.2, 1), 
                          guide = FALSE) +
    scale_y_continuous(limits = c(0, 0.50), breaks = seq(0, 1, 0.1)) +
    theme_ctokita() +
    theme(legend.position = "none",
          legend.background = element_blank())
  return(gg_dist_sum)
}

####################
# Load, process, and plot data
####################
# Panel A (linear mix)
load("output/Rdata/Mix_AThreshM_14.00_14.00_BThreshM_10.00_10.00_deltas_0.60_0.60_threshSlope_7_Aalpha_6.00_6.00_Balpha_2.00_2.00_quitP_0.20.Rdata")
panelA <- process_mix_data(task_dist)
panelA_data <- panelA[[1]]
panelA_null <- panelA[[2]]
gg_panelA <- plot_mix_data(panelA_data, panelA_null)
ggsave(gg_panelA, filename = "output/Task_dist/svg_files/non-5050-panelA.svg", width = 90, height = 45, units = "mm")

# Panel B (convex, high demand)
load("output/Rdata/Mix_AThreshM_10.00_10.00_BThreshM_10.00_10.00_deltas_1.50_1.50_threshSlope_7_Aalpha_6.00_6.00_Balpha_2.00_2.00_quitP_0.20.Rdata")
panelB <- process_mix_data(task_dist)
panelB_data <- panelB[[1]]
panelB_null <- panelB[[2]]
gg_panelB <- plot_mix_data(panelB_data, panelB_null)
ggsave(gg_panelB, filename = "output/Task_dist/svg_files/non-5050-panelB.svg", width = 90, height = 45, units = "mm")

# Panel C (convex, low demand)
load("output/Rdata/Mix_AThreshM_10.00_10.00_BThreshM_10.00_10.00_deltas_0.60_0.60_threshSlope_7_Aalpha_6.00_6.00_Balpha_2.00_2.00_quitP_0.20.Rdata")
panelC <- process_mix_data(task_dist)
panelC_data <- panelC[[1]]
panelC_null <- panelC[[2]]
gg_panelC <- plot_mix_data(panelC_data, panelC_null)
ggsave(gg_panelC, filename = "output/Task_dist/svg_files/non-5050-panelC.svg", width = 90, height = 45, units = "mm")

# Panel D (concave)
load("output/Rdata/Mix_AThreshM_20.00_20.00_BThreshM_10.00_10.00_deltas_0.60_0.60_threshSlope_7_Aalpha_6.00_6.00_Balpha_2.00_2.00_quitP_0.20.Rdata")
panelD <- process_mix_data(task_dist)
panelD_data <- panelD[[1]]
panelD_null <- panelD[[2]]
gg_panelD <- plot_mix_data(panelD_data, panelD_null)
ggsave(gg_panelD, filename = "output/Task_dist/svg_files/non-5050-panelD.svg", width = 90, height = 45, units = "mm")
gg_panelD_legend <- gg_panelD +
  theme(legend.position = "right")
ggsave(gg_panelD_legend, filename = "output/Task_dist/svg_files/non-5050-panelDlegend.svg", width = 90, height = 45, units = "mm")


