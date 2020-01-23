################################################################################
#
# Evaluate ensemble model outputs -- quick check for robustness fig
#
################################################################################

rm(list = ls())

#######################
# Set parameter range #
#######################

# Robustness check
# mu_sweep    <- seq(10, 20, by = 1) # range of AThreshM
# alpha_sweep <- seq(2, 7, by = 0.5)  # range of Aalpha
# mu_sweep    <- seq(10, 20, by = 1) # range of AThreshM
# alpha_sweep <- seq(2, 7, by = 0.5)  # range of Aalpha
###### NEW 012120
mu_sweep    <- seq(6, 20, by = 1) # range of AThreshM
alpha_sweep <- seq(1, 7, by = 0.5)  # range of Aalpha

params        <- matrix(data = NA, nrow = length(mu_sweep)*length(alpha_sweep), ncol = 10)
params[,3:4]  <- 2    # efficiency of A
params[,5:6]  <- 0.6  # demand rate
params[,9:10] <- 10   # mean threshold of A

# Fill in params matrix
for(i in 1:nrow(params)){
  params[i,7:8] <- mu_sweep[(i-1)%/%length(alpha_sweep)+1]
  params[i,1:2] <- alpha_sweep[((i-1)%%length(alpha_sweep)+1)]
}

# ...or select individual param set
params <- matrix(c(5, 5, 2, 2, 0.6,	0.6, 10, 10, 10, 10), nrow = 1, ncol = 10, byrow = TRUE)  #4a
params <- matrix(c(1.5, 1.5, 2, 2, 0.6,	0.6, 7, 7, 10, 10), nrow = 1, ncol = 10, byrow = TRUE)  #4c
params <- matrix(c(3, 3, 2, 2, 0.6,	0.6, 15, 15, 10, 10), nrow = 1, ncol = 10, byrow = TRUE)  #4d

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
  load(paste0("output/Rdata/", file_name, "_robust_50.Rdata"))
  
  
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
  # 1) "Line" -> "Type"
  colnames(task_VarMean_byMix)[colnames(task_VarMean_byMix) == "Line"] <- "Type"
  colnames(task_VarMean_byrep)[colnames(task_VarMean_byrep) == "Line"] <- "Type"
  colnames(task_dist)[colnames(task_dist) == "Line"] <- "Type"
  
  # 2) Change category names
  task_VarMean_byMix$Type[task_VarMean_byMix$Type == "A"] <- "Type X"
  task_VarMean_byMix$Type[task_VarMean_byMix$Type == "B"] <- "Type Y"
  task_VarMean_byMix$Type[task_VarMean_byMix$Type == "AB"] <- "Mixed"
  
  task_VarMean_byMix$Mix <- as.character(task_VarMean_byMix$Mix)
  task_VarMean_byMix$Mix[task_VarMean_byMix$Mix == "A"] <- "Type X"
  task_VarMean_byMix$Mix[task_VarMean_byMix$Mix == "B"] <- "Type Y"
  task_VarMean_byMix$Mix[task_VarMean_byMix$Mix == "AB"] <- "Mixed"
  
  task_VarMean_byrep$Type[task_VarMean_byrep$Type == "A"] <- "Type X"
  task_VarMean_byrep$Type[task_VarMean_byrep$Type == "B"] <- "Type Y"
  task_VarMean_byrep$Type[task_VarMean_byrep$Type == "AB"] <- "Mixed"
  
  task_VarMean_byrep$Mix <- as.character(task_VarMean_byrep$Mix)
  task_VarMean_byrep$Mix[task_VarMean_byrep$Mix == "A"] <- "Type X"
  task_VarMean_byrep$Mix[task_VarMean_byrep$Mix == "B"] <- "Type Y"
  task_VarMean_byrep$Mix[task_VarMean_byrep$Mix == "AB"] <- "Mixed"
  
  task_dist$Type[task_dist$Type == "A"] <- "Type X"
  task_dist$Type[task_dist$Type == "B"] <- "Type Y"
  task_dist$Type[task_dist$Type == "AB"] <- "Mixed"
  
  task_dist$Mix <- as.character(task_dist$Mix)
  task_dist$Mix[task_dist$Mix == "A"] <- "Type X"
  task_dist$Mix[task_dist$Mix == "B"] <- "Type Y"
  task_dist$Mix[task_dist$Mix == "AB"] <- "Mixed"
  
  # Adjust x label order - 9/19/19
  task_VarMean_byrep$Mix <- factor(task_VarMean_byrep$Mix, levels = c("Type X","Type Y","Mixed"))
  task_VarMean_byMix$Mix <- factor(task_VarMean_byMix$Mix, levels = c("Type X","Type Y","Mixed"))
  task_VarMean_byrep$Type <- factor(task_VarMean_byrep$Type, levels = c("Type X","Type Y","Mixed"))
  task_VarMean_byMix$Type <- factor(task_VarMean_byMix$Type, levels = c("Type X","Type Y","Mixed"))
  
  # Fig title based on pattern - 1/22/20
  # # extract relevant data
  # X_pure <- task_VarMean_byMix$Mean1[task_VarMean_byMix$Mix=="Type X" & task_VarMean_byMix$Group=="Type X"]
  # Y_pure <- task_VarMean_byMix$Mean1[task_VarMean_byMix$Mix=="Type Y" & task_VarMean_byMix$Group=="Type Y"]
  # X_mix <- task_VarMean_byMix$Mean1[task_VarMean_byMix$Mix=="Mixed" & task_VarMean_byMix$Group=="Type X"]
  # Y_mix <- task_VarMean_byMix$Mean1[task_VarMean_byMix$Mix=="Mixed" & task_VarMean_byMix$Group=="Type Y"]
  # X_pure_SE <- task_VarMean_byMix$SE1[task_VarMean_byMix$Mix=="Type X" & task_VarMean_byMix$Group=="Type X"]
  # Y_pure_SE <- task_VarMean_byMix$SE1[task_VarMean_byMix$Mix=="Type Y" & task_VarMean_byMix$Group=="Type Y"]
  # X_mix_SE <- task_VarMean_byMix$SE1[task_VarMean_byMix$Mix=="Mixed" & task_VarMean_byMix$Group=="Type X"]
  # Y_mix_SE <- task_VarMean_byMix$SE1[task_VarMean_byMix$Mix=="Mixed" & task_VarMean_byMix$Group=="Type Y"]
  # 
  # # compute and store deg of behavioral amplification / contagion ("amp")
  # # amp <- abs(Y_mix - X_mix) - abs(Y_pure - X_pure)
  # amp <- (Y_mix - X_mix) - (Y_pure - X_pure)  # new definition
  # robustcheck[INDEX,'amps'] <- amp
  # 
  # # if there is a 'flip' in order, record NaN; will appear as gray
  # CIfactor <- 0
  # 
  # if( X_pure < Y_pure ){
  #   if( (X_mix - X_mix_SE*CIfactor) > (Y_mix + Y_mix_SE*CIfactor) ){
  #     robustcheck[INDEX,'amps'] <- NaN
  #     print("CASE1")
  #     print(paste0(INDEX))
  #     print(paste0(params[INDEX,]))
  #   }
  # }else if( Y_pure < X_pure ){
  #   if( (Y_mix - Y_mix_SE*CIfactor) > (X_mix + X_mix_SE*CIfactor) ){
  #     robustcheck[INDEX,'amps'] <- NaN
  #     print("CASE2")
  #     print(paste0(INDEX))
  #     print(paste0(params[INDEX,]))
  #   }
  # }
  
  # Means of means
  gg_dist3 <- ggplot(data = task_VarMean_byrep, aes(y = Mean1, x = Mix, colour = Type)) +
    geom_point(size = 0.3, alpha = 0.2, stroke = 0, 
               position = position_dodge(width = 1)) +
    labs(x = "",
         y = "Task 1 performance, mean \u00B1 s.e.") +
    scale_color_manual(values = c("#E52521","#2B4B9B","#7C217F")) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, yinc)) +
    theme_mk() +
    theme(legend.position = "none",
          axis.text.x     = element_text(colour = c("#E52521","#2B4B9B","#7C217F")),
          axis.title.x    = element_text(size=0)) +
    geom_hline( yintercept = mean( task_VarMean_byMix[task_VarMean_byMix$Mix != "Mixed",]$Mean1 ),
                lty = 1, size = 0.1, color = "gray30" ) +
    geom_point(data = task_VarMean_byMix, aes(x = Mix, y = Mean1),
               size = 0.8, alpha = 1, stroke = 0.2, 
               # shape = 21, fill = NA, 
               position = position_dodge(width = 1)) +
    geom_errorbar(data = task_VarMean_byMix, 
                  aes(x = Mix, ymin = Mean1 - SE1, ymax = Mean1 + SE1),
                  size = 0.2, width = 0.6, 
                  position = position_dodge(width = 1))
  
  gg_dist3
  
  ggsave(filename = paste0("output/Parameter_exp/", file_name, "_robust_50.png"), width = figH, height = figH*1.15, dpi = 800)

  
}

