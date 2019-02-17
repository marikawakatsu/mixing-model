################################################################################
#
# Evaluate ensemble model outputs
# Updated 02/16/19:
#     Plot mixes with means and predictions by line
#       
################################################################################

rm(list = ls())
# base
# params <- matrix(c(2, 2, 2, 2, 0.6,	0.6, 10, 10, 10, 10),
#                  nrow = 1, ncol = 10, byrow = TRUE)

# vary alphas only
# params <- matrix(c(2, 2, 6, 6, 0.6,	0.6, 10, 10, 10, 10,
#                    2, 6, 6, 2, 0.6,	0.6, 10, 10, 10, 10,
#                    2, 2, 1, 1, 0.6,	0.6, 10, 10, 10, 10,
#                    2, 1, 1, 2, 0.6,	0.6, 10, 10, 10, 10,
#                    3, 3, 1, 1, 0.6,	0.6, 10, 10, 10, 10,
#                    3, 1, 1, 3, 0.6,	0.6, 10, 10, 10, 10),
#                  nrow = 6, ncol = 10, byrow = TRUE)

# vary deltas only
# params <- matrix(c(2, 2, 2, 2, 0.6,	1.0, 10, 10, 10, 10,
#                    2, 2, 2, 2, 0.6,	0.2, 10, 10, 10, 10),
#                  nrow = 2, ncol = 10, byrow = TRUE)

# vary both alphas and deltas
# params <- matrix(c(2, 2, 6, 6, 0.6,	1.0, 10, 10, 10, 10,
#                    2, 6, 6, 2, 0.6,	1.0, 10, 10, 10, 10,
#                    2, 2, 6, 6, 1.0,	0.6, 10, 10, 10, 10,
#                    2, 6, 6, 2, 1.0,	0.6, 10, 10, 10, 10,
#                    2, 2, 1, 1, 0.6,	1.0, 10, 10, 10, 10,
#                    2, 1, 1, 2, 0.6,	1.0, 10, 10, 10, 10,
#                    2, 2, 1, 1, 1.0,	0.6, 10, 10, 10, 10,
#                    2, 1, 1, 2, 1.0,	0.6, 10, 10, 10, 10),
#                  nrow = 8, ncol = 10, byrow = TRUE)

# vary both alphas and deltas, v2
# params <- matrix(c(2, 2, 6, 6, 1.0,	1.0, 10, 10, 10, 10,
#                    2, 6, 6, 2, 1.0,	1.0, 10, 10, 10, 10,
#                    2, 2, 1, 1, 1.0,	1.0, 10, 10, 10, 10,
#                    2, 1, 1, 2, 1.0,	1.0, 10, 10, 10, 10),
#                  nrow = 4, ncol = 10, byrow = TRUE)
# 
# params <- matrix(c(2, 2, 6, 6, 1,	1, 10, 10, 10, 10),
#                  nrow = 1, ncol = 10, byrow = TRUE)

params <- matrix(c(2, 6, 6, 2, 0.6,	0.6, 10, 10, 10, 10,
                   2, 6, 6, 2, 1.0,	1.0, 10, 10, 10, 10,
                   2, 1, 1, 2, 0.6,	0.6, 10, 10, 10, 10,
                   2, 1, 1, 2, 1.0,	1.0, 10, 10, 10, 10),
                 nrow = 4, ncol = 10, byrow = TRUE)
# 
# params <- matrix(c(2, 2, 2, 2, 0.6,	0.6, 10, 15, 15, 10),
#                  nrow = 1, ncol = 10, byrow = TRUE)

for (INDEX in 1:nrow(params)){
  # rm(list = ls())
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
  A_ThreshSD     <- A_ThreshM * 0. #1 #population threshold standard deviations for clone line A !!Change!!
  B_ThreshM      <- c(params[INDEX,9], params[INDEX,10]) #population threshold means for clone line B !!Change!!
  B_ThreshSD     <- B_ThreshM * 0. #1 #population threshold standard deviations for clone line B !!Change!!
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
  
  # # NEW 02/16/19
  task_dist <- as.data.frame(task_dist)

  # Prepare means and SDs
  task_VarMean_byrepbyLine <- task_dist %>% # stats by replicate and Line
    # mutate(Line = paste0(Line," simulation")) %>% 
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
    # mutate(Line = "Mixed simulation")
  
  task_VarMean_byrep <- rbind(task_VarMean_byrepbyLine, task_VarMean_byrepABonly) %>% group_by(n, Mix) 
  
  task_VarMean_byMixbyLine <- task_dist %>% # stats by Mix
    # mutate(Line = paste0(Line," simulation")) %>% 
    group_by(n, Mix, replicate, Line, set) %>%
    summarise(Mean1rep = mean(Task1),
              Mean2rep = mean(Task2)) %>%
    group_by(n, Mix, Line) %>%
    summarise(SD1 = sd(Mean1rep),
              SD2 = sd(Mean2rep),
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
              Mean1 = mean(Mean1rep),
              Mean2 = mean(Mean2rep)) %>%
    mutate(Line = "Mixed")
    # mutate(Line = "Mixed simulation")
  
  task_VarMean_byMix <- rbind(task_VarMean_byMixbyLine, task_VarMean_byMixABonly) %>% group_by(n, Mix)
  
  # NEW: Analytical predictions
  if( A_ThreshM[[1]]==B_ThreshM[[2]] & A_ThreshM[[1]]==B_ThreshM[[2]] & A_alpha[[1]]==A_alpha[[2]] & A_alpha[[1]]==B_alpha[[2]]){
    print("Varying mu")
    a <- A_ThreshM[[1]]
    b <- A_ThreshM[[2]]
    s <- ((1/2)*( -(a^threshSlope + b^threshSlope) + sqrt( (a^threshSlope + b^threshSlope)^2 + (a^threshSlope * b^threshSlope) * 8*deltas[[1]]*quitP[[1]] / (A_alpha[[1]] - 2*deltas[[1]]*(1+quitP[[1]])) ) ))^(1/threshSlope)
    n1A_pred <- (1 / quitP[[1]])*(s^threshSlope / (s^threshSlope + a^threshSlope))*(2 - (s^threshSlope / (s^threshSlope + b^threshSlope)) )*(1/2 - deltas[[1]]/A_alpha[[1]])
    n2A_pred <- (1 / quitP[[1]])*(s^threshSlope / (s^threshSlope + b^threshSlope))*(2 - (s^threshSlope / (s^threshSlope + a^threshSlope)) )*(1/2 - deltas[[1]]/A_alpha[[1]])
    n1B_pred <- n2A_pred
    n2B_pred <- n1A_pred
    
    }else if( A_ThreshM[[1]]==A_ThreshM[[2]] & A_ThreshM[[1]]==B_ThreshM[[2]] ){
      print("Varying delta and alpha")
      n1A_pred <- 2*deltas[[1]]/(A_alpha[[1]]+B_alpha[[1]])
      n2A_pred <- 2*deltas[[1]]/(A_alpha[[2]]+B_alpha[[2]])
      n1B_pred <- n1A_pred
      n2B_pred <- n2A_pred
  
    }else{
      print("Check your conditions!")
      break
    }
  task_Mean_byMix_pred <- data.frame(
    Mean1 = c(deltas[[1]]/A_alpha[[1]], deltas[[1]]/B_alpha[[1]], n1A_pred, n1B_pred, 2*deltas[[1]]/(A_alpha[[1]]+B_alpha[[1]])),
    Mean2 = c(deltas[[2]]/A_alpha[[2]], deltas[[2]]/B_alpha[[2]], n2A_pred, n2B_pred, 2*deltas[[2]]/(A_alpha[[2]]+B_alpha[[2]])),
    Line = c("A","B","A","B","Mixed"),
    Mix = c("A","B","AB","AB","AB"),
    n = 16,
    SD1 = NA,
    SD2 = NA
  )
  
  # NEW: only for the varying delta and alpha type
  if( A_ThreshM[[1]]==A_ThreshM[[2]] & A_ThreshM[[1]]==B_ThreshM[[2]] ){
    task_VarMean_byrep <- task_VarMean_byrep[task_VarMean_byrep$Line == "Mixed" | (task_VarMean_byrep$Mix %in% c("A","B")), ]
    task_VarMean_byMix <- task_VarMean_byMix[task_VarMean_byMix$Line == "Mixed" | (task_VarMean_byMix$Mix %in% c("A","B")), ]
    task_Mean_byMix_pred <- task_Mean_byMix_pred[task_Mean_byMix_pred$Line == "Mixed" | (task_Mean_byMix_pred$Mix %in% c("A","B")), ]
  }
  
  task_Mean_byMix_pred <- task_Mean_byMix_pred %>% mutate("Type" = "Prediction")
  task_VarMean_byMix <- task_VarMean_byMix %>% mutate("Type" = "Simulation")
  task_VarMean_comb <- merge(task_VarMean_byMix, task_Mean_byMix_pred, all=TRUE)
  
  # Means of means
  gg_dist3 <- 
    ggplot(data = task_VarMean_comb, aes(y = Mean1, x = Mix, colour = Line, shape = Type)) +
    geom_point(data = task_VarMean_comb, aes(x = Mix, y = Mean1),
               size = 1, alpha = 1, stroke = 0.3, position = position_dodge(width = 0.7)) +
    theme_classic() +
    labs(x = "Mix",
         y = "Frequency Task 1") +
    scale_color_manual(values = c("#ca0020", "#0571b0","#80007F")) +
    scale_shape_manual(values = rep(c(1,16),3)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_ctokita() +
    # theme(axis.text.x = Mix) +
    geom_errorbar(data = task_VarMean_comb, aes(x = Mix, ymin = Mean1 - SD1, ymax = Mean1 + SD1),
                  size = 0.05, width = 0.4, position = position_dodge(width = 0.7)) 
  
  gg_dist3
  ggsave(filename = paste0("output/Task_dist/vs_analytical/", file_name, "_Task1_comp.png"), width = 3, height = 1.5, dpi = 400)
  
  gg_dist4 <- 
    ggplot(data = task_VarMean_comb, aes(y = Mean2, x = Mix, colour = Line, shape = Type)) +
    geom_point(data = task_VarMean_comb, aes(x = Mix, y = Mean2),
               size = 1, alpha = 1, stroke = 0.3, position = position_dodge(width = 0.7)) +
    theme_classic() +
    labs(x = "Mix",
         y = "Frequency Task 2") +
    scale_color_manual(values = c("#ca0020", "#0571b0","#80007F")) +
    scale_shape_manual(values = rep(c(1,16),3)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_ctokita() +
    # theme(axis.text.x = Mix) +
    geom_errorbar(data = task_VarMean_comb, aes(x = Mix, ymin = Mean2 - SD2, ymax = Mean2 + SD2),
                  size = 0.05, width = 0.4, position = position_dodge(width = 0.7)) 
  
  gg_dist4
  ggsave(filename = paste0("output/Task_dist/vs_analytical/", file_name, "_Task2_comp.png"), width = 3, height = 1.5, dpi = 400)
  
}
