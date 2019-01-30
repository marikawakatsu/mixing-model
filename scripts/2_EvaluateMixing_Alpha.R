################################################################################
#
# Evaluate ensemble model outputs
#
################################################################################
rm(list = ls())

source("scripts/util/__Util__MASTER.R")

load("output/Rdata/Mix_AThreshM_10.00_10.00_BThreshM_10.00_10.00_deltas_0.60_0.10_threshSlope_7_Aalpha_2.00_2.00_Balpha_0.50_0.50_quitP_0.20.Rdata")

file_name <- "Mix_Alphas"

####################
# Final task distributions
####################
# Prepare
task_dist <- task_dist %>% 
  group_by(n) %>% 
  mutate(set = paste0(Mix, "-", replicate)) %>% 
  mutate(set = factor(set, 
                      levels = c(sort(apply(expand.grid(unique(Mix), unique(replicate)), 1, paste, collapse = "-")))))

# Plot
gg_dist <- ggplot(data = task_dist, aes(y = Task1, x = set, color = Line)) +
  geom_point(size = 0.3) +
  theme_classic() +
  labs(x = "Replicate",
       y = "Frequency Task 1") +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = c("#ca0020", "#0571b0")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_ctokita() +
  theme(axis.text.x = element_blank(), 
        strip.background = element_rect(color = NA, fill = "grey85")) +
  facet_grid(n~.)
gg_dist

ggsave(filename = paste0("output/Task_dist/", file_name, ".png"), width = 3, height = 3, dpi = 400)

