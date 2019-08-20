################################################################################
#
# Evaluate ensemble model outputs
# Updated 06/05/19: 
#       Plots the behavioral variation + specialization in colonies
#
################################################################################

rm(list = ls())

# Fig. 3a-b
# params <- matrix(c(2, 2, 1, 1, 0.6,	0.6, 10, 10, 10, 10), nrow = 1, ncol = 10, byrow = TRUE) # 3a
params <- matrix(c(2, 2, 1, 1, 0.4,	0.4, 10, 10, 10, 10), nrow = 1, ncol = 10, byrow = TRUE) # 3b
                 
# Fig. S1 (1b and 1c require manual changes below)
# params <- matrix(c(2, 2, 2, 2, 0.6,	0.6, 10, 10, 20, 20), nrow = 1, ncol = 10, byrow = TRUE) # S1a
# params <- matrix(c(2, 2, 2, 2, 0.6,	0.6, 10, 10, 10, 10), nrow = 1, ncol = 10, byrow = TRUE) # S1b and S1c


# Plotting
ymax <- 0.5 # max y for plotting
yinc <- 0.1 # y-axis increments
figH <- 1.5 # figure height for printing; default width is 3


for (INDEX in 1:nrow(params)){
  source("scripts/util/__Util__MASTER.R")
  
  ####################
  # Set global variables
  ####################
  # Initial paramters: Free to change
  # Base parameters
  Ns             <- c(16) #vector of number of individuals to simulate
  m              <- 2 #number of tasks
  gens           <- 10000 #number of generations to run simulation 
  corrStep       <- 200 #number of time steps for calculation of correlation 
  reps           <- 10 #number of replications per simulation (for ensemble) !!Change!!
  
  # Threshold Parameters
  mixes          <- c("A", "B", "AB")
  A_ThreshM      <- c(params[INDEX,7], params[INDEX,8]) #population threshold means for clone line A !!Change!!
  A_ThreshSD     <- A_ThreshM * 0.1 #population threshold standard deviations for clone line A !!Change!!
  B_ThreshM      <- c(params[INDEX,9], params[INDEX,10]) #population threshold means for clone line B !!Change!!
  B_ThreshSD     <- B_ThreshM * 0.1 #population threshold standard deviations for clone line B !!Change!!
  InitialStim    <- c(0, 0) #intital vector of stimuli
  deltas         <- c(params[INDEX,5], params[INDEX,6]) #vector of stimuli increase rates  
  threshSlope    <- 7 #exponent parameter for threshold curve shape
  alpha          <- m
  # A_alpha        <- c(m, m*3) #efficiency of task performance
  # B_alpha        <- c(m*3, m)
  A_alpha        <- c(params[INDEX,1], params[INDEX,2])
  B_alpha        <- c(params[INDEX,3], params[INDEX,4])
  quitP          <- c(0.2, 0.2) #probability of quitting task once active !!Change!!
  
  file_name1 <- sprintf("N16only_AThreshM_%1.2f_%1.2f_BThreshM_%1.2f_%1.2f_deltas_%1.2f_%1.2f_threshSlope_%d_Aalpha_%1.2f_%1.2f_Balpha_%1.2f_%1.2f_quitP_%1.2f",
                        A_ThreshM[1], A_ThreshM[2], B_ThreshM[1], B_ThreshM[2], deltas[1], deltas[2], threshSlope, 
                        A_alpha[1], A_alpha[2], B_alpha[1], B_alpha[2], quitP[1])  # for quitp[1] = quitP[2]
  
  file_name2 <- sprintf("N16only_AThreshM_%1.2f_%1.2f_AThreshSD_%1.2f_%1.2f_BThreshM_%1.2f_%1.2f_BThreshSD_%1.2f_%1.2f_deltas_%1.2f_%1.2f_threshSlope_%d_%d_Aalpha_%1.2f_%1.2f_Balpha_%1.2f_%1.2f_quitP_%1.2f_%1.2f",
                        A_ThreshM[1], A_ThreshM[2], A_ThreshSD[1]/A_ThreshM[1], A_ThreshSD[2]/A_ThreshM[2], 
                        B_ThreshM[1], B_ThreshM[2], B_ThreshSD[1]/B_ThreshM[1], B_ThreshSD[2]/B_ThreshM[2],
                        deltas[1], deltas[2], threshSlope, threshSlope, A_alpha[1], A_alpha[2], 
                        B_alpha[1], B_alpha[2], quitP[1], quitP[2])
  
  file_name <- file_name2
  rm(file_name1, file_name2)
  
  # load(paste0("output/Rdata/", file_name, "reps_100.Rdata"))
  load(paste0("output/Rdata/", file_name, ".Rdata"))
  
  
  ####################
  # Final task distributions
  ####################
  # Prepare data
  task_dist <- task_dist %>% 
    group_by(n) %>% 
    mutate(set = paste0(Mix, "-", replicate)) %>% 
    mutate(set = factor(set, 
                        levels = c(paste("A", unique(replicate), sep = "-"),  
                                   paste("B", unique(replicate), sep = "-"),
                                   paste("AB", unique(replicate), sep = "-")) ))
  
  # Prepare means and SDs
  task_VarMean_byrepbyLine <- task_dist %>% # stats by replicate and Line
    group_by(n, Mix, replicate, Line, set) %>%
    summarise(SD1 = sd(Task1),
              SD2 = sd(Task2),
              Mean1 = mean(Task1),
              Mean2 = mean(Task2))
  
  task_VarMean_byrepABonly <- task_dist %>% # stats by replicate, mixed colonies only
    filter(Mix == "AB") %>%
    group_by(n, Mix, replicate, set) %>%
    summarise(SD1 = sd(Task1),
              SD2 = sd(Task2),
              Mean1 = mean(Task1),
              Mean2 = mean(Task2)) %>%
    mutate(Line = "Mixed")
  
  task_VarMean_byrep <- rbind(task_VarMean_byrepbyLine, task_VarMean_byrepABonly) %>% group_by(n, Mix)
  
  task_VarMean_byMixbyLine <- task_dist %>% # stats by Mix
    group_by(n, Mix, replicate, Line, set) %>%
    summarise(Mean1rep = mean(Task1),
              Mean2rep = mean(Task2)) %>%
    group_by(n, Mix, Line) %>%
    summarise(SD1 = sd(Mean1rep),
              SD2 = sd(Mean2rep),
              SE1 = sd(Mean1rep) / sqrt(length(Mean1rep)),
              SE2 = sd(Mean2rep) / sqrt(length(Mean2rep)),
              Mean1 = mean(Mean1rep),
              Mean2 = mean(Mean2rep))
  
  task_VarMean_byMixABonly <- task_dist %>% # stats by Mix, mixed colonies only
    filter(Mix == "AB") %>%
    group_by(n, Mix, replicate, set) %>%
    summarise(Mean1rep = mean(Task1),
              Mean2rep = mean(Task2)) %>%
    group_by(n, Mix) %>%
    summarise(SD1 = sd(Mean1rep),
              SD2 = sd(Mean2rep),
              SE1 = sd(Mean1rep) / sqrt(length(Mean1rep)),
              SE2 = sd(Mean2rep) / sqrt(length(Mean2rep)),
              Mean1 = mean(Mean1rep),
              Mean2 = mean(Mean2rep)) %>%
    mutate(Line = "Mixed")
  
  task_VarMean_byMix <- rbind(task_VarMean_byMixbyLine, task_VarMean_byMixABonly) %>% group_by(n, Mix)
  
  # !!! NEW !!! change 
  # 1) "Line" -> "Group"
  colnames(task_VarMean_byMix)[colnames(task_VarMean_byMix) == "Line"] <- "Group"
  colnames(task_VarMean_byrep)[colnames(task_VarMean_byrep) == "Line"] <- "Group"
  colnames(task_dist)[colnames(task_dist) == "Line"] <- "Group"
  
  # 2) Change category names
  task_VarMean_byMix$Group[task_VarMean_byMix$Group == "A"] <- "Group A"
  task_VarMean_byMix$Group[task_VarMean_byMix$Group == "B"] <- "Group B"
  task_VarMean_byMix$Group[task_VarMean_byMix$Group == "AB"] <- "Mixed"
  
  task_VarMean_byMix$Mix <- as.character(task_VarMean_byMix$Mix)
  task_VarMean_byMix$Mix[task_VarMean_byMix$Mix == "A"] <- "Group A"
  task_VarMean_byMix$Mix[task_VarMean_byMix$Mix == "B"] <- "Group B"
  task_VarMean_byMix$Mix[task_VarMean_byMix$Mix == "AB"] <- "Mixed"
  
  task_VarMean_byrep$Group[task_VarMean_byrep$Group == "A"] <- "Group A"
  task_VarMean_byrep$Group[task_VarMean_byrep$Group == "B"] <- "Group B"
  task_VarMean_byrep$Group[task_VarMean_byrep$Group == "AB"] <- "Mixed"
  
  task_VarMean_byrep$Mix <- as.character(task_VarMean_byrep$Mix)
  task_VarMean_byrep$Mix[task_VarMean_byrep$Mix == "A"] <- "Group A"
  task_VarMean_byrep$Mix[task_VarMean_byrep$Mix == "B"] <- "Group B"
  task_VarMean_byrep$Mix[task_VarMean_byrep$Mix == "AB"] <- "Mixed"
  
  task_dist$Group[task_dist$Group == "A"] <- "Group A"
  task_dist$Group[task_dist$Group == "B"] <- "Group B"
  task_dist$Group[task_dist$Group == "AB"] <- "Mixed"
  
  task_dist$Mix <- as.character(task_dist$Mix)
  task_dist$Mix[task_dist$Mix == "A"] <- "Group A"
  task_dist$Mix[task_dist$Mix == "B"] <- "Group B"
  task_dist$Mix[task_dist$Mix == "AB"] <- "Mixed"
  
  # Means of means
  gg_dist3 <- ggplot(data = task_VarMean_byrep, aes(y = Mean1, x = Mix, colour = Group)) +
    geom_point(size = 0.3, alpha = 0.2, stroke = 0, 
               position = position_dodge(width = 1)) +
    labs(x = "",
         y = "Frequency task 1, mean \u00B1 s.e.") +
    scale_color_manual(values = c("#E52521","#2B4B9B","#7C217F")) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, yinc)) +
    theme_mk() +
    theme(axis.text.x = element_text(colour = c("#E52521","#2B4B9B","#7C217F"))) +
    geom_hline( yintercept = mean( task_VarMean_byMix[task_VarMean_byMix$Mix != "Mixed",]$Mean1 ),
                lty = 1, size = 0.1, color = "gray30" ) +
    geom_point(data = task_VarMean_byMix, aes(x = Mix, y = Mean1),
               size = 0.8, alpha = 1, stroke = 0.2, 
               # shape = 21, fill = NA, 
               position = position_dodge(width = 1)) +
    geom_errorbar(data = task_VarMean_byMix, 
                  aes(x = Mix, ymin = Mean1 - SE1, ymax = Mean1 + SE1),
                  size = 0.2, width = 0.6, 
                  position = position_dodge(width = 1)) +
    theme(legend.position="none")
  
  gg_dist3
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_Task1Summary_SE_nolegend.png"), width = figH, height = figH*1.15, dpi = 800)
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_Task1Summary_SE.png"), width = 2.25, height = figH, dpi = 800)
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_reps_100_Task1Summary_SE.png"), width = 2.25, height = figH, dpi = 800)
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_reps_100_Task1Summary_SE_nolegend.png"), width = figH, height = figH*1.15, dpi = 800)
  
  gg_dist4 <- ggplot(data = task_VarMean_byrep, aes(y = Mean2, x = Mix, colour = Group)) +
    geom_point(size = 0.3, alpha = 0.2, stroke = 0,
               position = position_dodge(width = 1)) +
    labs(x = "",
         y = "Frequency task 2, mean \u00B1 s.e.") +
    scale_color_manual(values = c("#E52521", "#2B4B9B", "#7C217F")) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, yinc)) +
    theme_mk() +
    theme(axis.text.x = element_text(colour = c("#E52521","#2B4B9B","#7C217F"))) +
    # theme(axis.text.x = Mix) +
    geom_point(data = task_VarMean_byMix, aes(x = Mix, y = Mean2),
               size = 0.8, alpha = 1, stroke = 0.2,
               # shape = 21, fill = NA,
               position = position_dodge(width = 1)) +
    geom_errorbar(data = task_VarMean_byMix,
                  aes(x = Mix, ymin = Mean2 - SE2, ymax = Mean2 + SE2),
                  size = 0.2, width = 0.4,
                  position = position_dodge(width = 1))

  gg_dist4
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_reps_100_Task2BehavVar.png"), width = 3, height = figH, dpi = 400)
  
  ####################
  # Task Rank Correlation
  ####################
  # Compute mean
  task_corr <- task_corr %>% 
    mutate(TaskMean = (Task1 + Task2) / 2)
  
  # Calculate means and SE
  task_corr_VarMean <- task_corr %>% 
    group_by(n, Mix) %>% 
    summarise(SpecMean = mean(TaskMean),
              SpecSE = sd(TaskMean) / sqrt(length(TaskMean)),
              SpecCI = 1.96 * SpecSE)
  
  # Change label names
  task_corr$Mix <- as.character(task_corr$Mix)
  task_corr$Mix[task_corr$Mix == "A"] <- "Group A"
  task_corr$Mix[task_corr$Mix == "B"] <- "Group B"
  task_corr$Mix[task_corr$Mix == "AB"] <- "Mixed"
  
  task_corr_VarMean$Mix <- as.character(task_corr_VarMean$Mix)
  task_corr_VarMean$Mix[task_corr_VarMean$Mix == "A"] <- "Group A"
  task_corr_VarMean$Mix[task_corr_VarMean$Mix == "B"] <- "Group B"
  task_corr_VarMean$Mix[task_corr_VarMean$Mix == "AB"] <- "Mixed"
  
  # Plot
  gg_corr <- ggplot(data = task_corr_VarMean, aes(x = Mix, y = SpecMean, colour = Mix)) +
    geom_point(data = task_corr, aes(x = Mix, y = TaskMean, colour = Mix),
               size = 0.3, alpha = 0.2, stroke = 0, 
               position = position_dodge(width = 1)) +
    theme_mk() +
    theme(legend.position = "none",
          axis.text.x = element_text(colour = c("#E52521","#2B4B9B","#7C217F"))) +
    labs(x = "",
         y = "Specialization, mean \u00B1 s.e.") +
    scale_color_manual(values = c("#E52521", "#2B4B9B", "#7C217F")) +
    scale_y_continuous(limits = c(-0.1, 1), breaks = seq(-1, 1, 0.2)) +
    # Mean and SE portion of plot
    geom_errorbar(aes(x = Mix, ymin = SpecMean - SpecSE, ymax = SpecMean + SpecSE),
                  size = 0.2, width = 0.6) +
    geom_point(size = 0.8, alpha = 1, stroke = 0.2)
  
  gg_corr
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_Spec_nolegend.png"), width = figH, height = figH*1.15, dpi = 800)
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_reps_100_Spec.png"), width = 1.5, height = figH, dpi = 800)
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_reps_100_Spec_nolegend.png"), width = figH, height = figH*1.15, dpi = 800)
  
  ####################
  # Task variance by group size
  ####################
  # Calculate average SD by mix
  
  # WITHOUT per-group lines
  # task_VarMean_byrep <- task_VarMean_byrep[task_VarMean_byrep$Mix == task_VarMean_byrep$Group,] %>%
  #   mutate(SD = (SD1 + SD2)/2)
  # 
  # task_VarMean_SD <- task_VarMean_byrep %>%
  #   group_by(n, Mix) %>%
  #   summarise(SDMean = mean(SD),
  #             SDSE = sd(SD) / sqrt(length(SD)))
  
  # WITH per-group lines
  task_VarMean_byrep <- task_VarMean_byrep %>%
    mutate(SD = (SD1 + SD2)/2)
  
  task_VarMean_SD <- task_VarMean_byrep %>%
    group_by(n, Mix, Group) %>%
    summarise(SDMean = mean(SD),
              SDSE = sd(SD) / sqrt(length(SD)))
  
  # # Plot behavioral variation by mix -- WITHOUT per-group line in the Mixed case
  # gg_var <- ggplot(data = task_VarMean_SD, aes(x = Mix, y = SDMean, colour = Mix)) +
  #   geom_point(data = task_VarMean_byrep, aes(x = Mix, y = SD),
  #              size = 0.3, alpha = 0.2, stroke = 0, 
  #              position = position_dodge(width = 1)) +
  #   theme_mk() +
  #   theme(legend.position = "none",
  #         axis.text.x = element_text(colour = c("#E52521","#2B4B9B","#7C217F"))) +
  #   xlab("Mix") +
  #   ylab("Behavioral variation") +
  #   scale_color_manual(values = c("#E52521", "#2B4B9B", "#7C217F")) +
  #   scale_y_continuous(limits = c(0, 0.25), breaks = seq(-1, 1, 0.05)) +
  #   # Mean and SE portion of plot
  #   geom_errorbar(aes(x = Mix, ymin = SDMean - SDSE, ymax = SDMean + SDSE, colour = Mix),
  #                 size = 0.2, width = 0.4) +
  #   geom_point(size = 0.8, alpha = 1, stroke = 0.2)
  
  # Plot behavioral variation by mix -- WITH per-group line in the Mixed case
  gg_var <- ggplot(data = task_VarMean_SD, aes(x = Mix, y = SDMean, colour = Group)) +
    geom_point(data = task_VarMean_byrep, aes(x = Mix, y = SD),
               size = 0.3, alpha = 0.2, stroke = 0, 
               position = position_dodge(width = 1)) +
    theme_mk() +
    theme(legend.position = "none",
          axis.text.x = element_text(colour = c("#E52521","#2B4B9B","#7C217F"))) +
    xlab("Mix") +
    labs(x = "",
         y = "Behavioral variation, mean \u00B1 s.e.") +
    scale_color_manual(values = c("#E52521", "#2B4B9B", "#7C217F")) +
    scale_y_continuous(limits = c(0, 0.2), breaks = seq(-1, 1, 0.05)) +
    # Mean and SE portion of plot
    geom_errorbar(aes(x = Mix, ymin = SDMean - SDSE, ymax = SDMean + SDSE),
                  size = 0.2, width = 0.6, 
                  position = position_dodge(width = 1)) +
    geom_point(size = 0.8, alpha = 1, stroke = 0.2, 
               position = position_dodge(width = 1))
  
  gg_var
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_reps_100_Var.png"), width = 1.5, height = figH, dpi = 800)
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_Var_Sep_nolegend.png"), width = figH, height = figH*1.15, dpi = 800)
  # ggsave(filename = paste0("output/Task_dist/", file_name, "_reps_100_Var_Sep_nolegend.png"), width = figH, height = figH*1.15, dpi = 800)

  
}

