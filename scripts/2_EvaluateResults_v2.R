################################################################################
#
# Evaluate ensemble model outputs
#
################################################################################
rm(list = ls())

source("scripts/util/__Util__MASTER.R")

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
Ns             <- c(4, 16) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 10000 #number of generations to run simulation 
corrStep       <- 200 #number of time steps for calculation of correlation 
reps           <- 10 #number of replications per simulation (for ensemble) !!Change!!

# Threshold Parameters
mixes          <- c("A", "B", "AB")
A_ThreshM      <- c(10, 10) #population threshold means for clone line A !!Change!!
A_ThreshSD     <- A_ThreshM * 0.1 #population threshold standard deviations for clone line A !!Change!!
B_ThreshM      <- c(10, 10) #population threshold means for clone line B !!Change!!
B_ThreshSD     <- B_ThreshM * 0.1 #population threshold standard deviations for clone line B !!Change!!
InitialStim    <- c(0, 0) #intital vector of stimuli
deltas         <- c(0.6, 0.6) #vector of stimuli increase rates  
threshSlope    <- 7 #exponent parameter for threshold curve shape
alpha          <- m
A_alpha        <- c(m, m*3) #efficiency of task performance
B_alpha        <- c(m, m*3)
quitP          <- c(0.2, 0.2) #probability of quitting task once active !!Change!!

file_name <- sprintf("AThreshM_%1.2f_%1.2f_AThreshSD_%1.2f_%1.2f_BThreshM_%1.2f_%1.2f_BThreshSD_%1.2f_%1.2f_deltas_%1.2f_%1.2f_threshSlope_%d_%d_Aalpha_%1.2f_%1.2f_Balpha_%1.2f_%1.2f_quitP_%1.2f_%1.2f",
                     A_ThreshM[1], A_ThreshM[2], A_ThreshSD[1]/A_ThreshM[1], A_ThreshSD[2]/A_ThreshM[2], 
                     B_ThreshM[1], B_ThreshM[2], B_ThreshSD[1]/B_ThreshM[1], B_ThreshSD[2]/B_ThreshM[2],
                     deltas[1], deltas[2], threshSlope, threshSlope, A_alpha[1], A_alpha[2], 
                     B_alpha[1], B_alpha[2], quitP[1], quitP[2])

load(paste0("output/Rdata/", file_name, ".Rdata"))

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

