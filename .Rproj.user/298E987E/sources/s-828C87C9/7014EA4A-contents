################################################################################
#
# Evaluate ensemble model outputs
#
################################################################################
rm(list = ls())

source("scripts/util/__Util__MASTER.R")

load("output/Rdata/AlphaDiff_OneHighOneLow.Rdata")

file_name <- "AlphaDiff_OneHighOneLow"

####################
# Final task distributions
####################
# Prepare
task_dist <- task_dist %>% 
  group_by(n) %>% 
  mutate(set = paste0(Mix, "-", replicate)) %>% 
  mutate(set = factor(set, 
                      levels = c(paste("A", unique(replicate), sep = "-"),  
                                 paste("B", unique(replicate), sep = "-"),
                                 paste("AB", unique(replicate), sep = "-")) ))

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

