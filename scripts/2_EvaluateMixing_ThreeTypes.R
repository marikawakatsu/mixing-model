################################################################################
#
# Evaluate ensemble model outputs
#
################################################################################
rm(list = ls())

source("scripts/util/__Util__MASTER.R")

load("output/Rdata/Mix_ThreeTypes_AThreshM_10.00_10.00_BThreshM_15.00_15.00_CThreshM_7.50_7.50_deltas_0.60_0.60_threshSlope_7_Aalpha_2.00_2.00_Balpha_3.00_3.00_Calpha_1.50_1.50_quitP_0.20.Rdata")

file_name <- "Mix_threetypes-Scenario2"

####################
# Final task distributions
####################
###### Prepare ######
# Raw data
task_dist <- task_dist %>% 
  group_by(n) %>% 
  mutate(set = paste0(Mix, "-", replicate))

# Order sets for appropriate plotting
plot_order <- c()
mixes <- c("A", "B", "C", "AB", "AC", "BC", "ABC")
for (mix in mixes) {
  mix_reps <- paste(mix, seq(1, 100), sep = "-")
  plot_order <- c(plot_order, mix_reps)
}
task_dist <- task_dist %>% 
  mutate(set = factor(set, levels = plot_order))

# Mix means
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
  mutate(Line = as.character(Mix),
         Group_mean = TRUE)


# Mix x Line means
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
  mutate(Group_mean = FALSE)

#Bind
task_dist_summary <- task_dist_summary %>% 
  bind_rows(task_dist_lines) %>% 
  mutate(Group_mean = as.numeric(Group_mean),
         Line = factor(Line, levels = mixes))


###### Plot ######
line_cols <- c("#2B4B9B", "#EE751C", "#009640")
mix_col <- c(line_cols, "#8D615B", "#2B706E", "#77862E", "black")

# Plot raw data
gg_dist <- ggplot(data = task_dist, aes(y = Task1, x = set, color = Line)) +
  geom_point(size = 0.3) +
  theme_classic() +
  labs(x = "Replicate",
       y = "Frequency Task 1") +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = line_cols) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_ctokita() +
  theme(axis.text.x = element_blank(), 
        strip.background = element_rect(color = NA, fill = "grey85")) +
  facet_grid(n~.)
gg_dist

ggsave(filename = paste0("output/Task_dist/", file_name, ".png"), width = 5, height = 2, dpi = 400)

# Plot summary points
pure_A <- task_dist_summary$Task1_mean[task_dist_summary$Mix == "A" & task_dist_summary$Group_mean == 1]
pure_B <- task_dist_summary$Task1_mean[task_dist_summary$Mix == "B" & task_dist_summary$Group_mean == 1]
pure_C <- task_dist_summary$Task1_mean[task_dist_summary$Mix == "C" & task_dist_summary$Group_mean == 1]
gg_dist_sum <- ggplot(data = task_dist_summary, aes(y = Task1_mean, x = Mix, color = Line)) +
  geom_hline(yintercept = pure_A, 
             size = 0.1,
             linetype = "dotted",
             color = line_cols[1]) +
  geom_hline(yintercept = pure_B, 
             size = 0.1,
             linetype = "dotted",
             color = line_cols[2]) +
  geom_hline(yintercept = pure_C, 
             size = 0.1,
             linetype = "dotted",
             color = line_cols[3]) +
  geom_errorbar(aes(ymin = Task1_mean - Task1_SE, ymax = Task1_mean + Task1_SE), 
                position = position_dodge(width = 0.05),
                width = 0.035,
                size = 0.2) +
  geom_point(aes(size = Group_mean),
             position = position_dodge(width = 0.05)) +
  theme_classic() +
  labs(x = "Mix",
       y = "Frequency Task 1, mean \u00B1 s.e.") +
  scale_color_manual(values = mix_col) +
  scale_size_continuous(range = c(0.1, 1.5), 
                        guide = FALSE) +
  # scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 1, 0.1)) +
  theme_ctokita() +
  theme(legend.background = element_blank())
gg_dist_sum

ggsave(filename = paste0("output/Task_dist/", file_name, "_Means.png"), width = 4, height = 2, dpi = 400)
