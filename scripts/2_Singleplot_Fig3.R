################################################################################
#
# Script for Fig. 3 in MS
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")

#################
### Functions ###
#################

prepdata_dist <- function(taskdist){
  
  # Prepare data
  taskdist <- taskdist %>% 
    group_by(n) %>% 
    mutate(set = paste0(Mix, "-", replicate)) %>% 
    mutate(set = factor(set, 
                        levels = c(paste("A", unique(replicate), sep = "-"),  
                                   paste("B", unique(replicate), sep = "-"),
                                   paste("AB", unique(replicate), sep = "-")) ))
  
  # Prepare means and SDs
  taskVarMeanbyrepbyLine <- taskdist %>% # stats by replicate and Line
    group_by(n, Mix, replicate, Line, set) %>%
    summarise(SD1 = sd(Task1),
              SD2 = sd(Task2),
              Mean1 = mean(Task1),
              Mean2 = mean(Task2))
  
  taskVarMeanbyrepABonly <- taskdist %>% # stats by replicate, mixed colonies only
    filter(Mix == "AB") %>%
    group_by(n, Mix, replicate, set) %>%
    summarise(SD1 = sd(Task1),
              SD2 = sd(Task2),
              Mean1 = mean(Task1),
              Mean2 = mean(Task2)) %>%
    mutate(Line = "Mixed")
  
  taskVarMeanbyrep <- rbind(taskVarMeanbyrepbyLine, taskVarMeanbyrepABonly) %>% group_by(n, Mix)
  
  taskVarMeanbyMixbyLine <- taskdist %>% # stats by Mix
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
  
  taskVarMeanbyMixABonly <- taskdist %>% # stats by Mix, mixed colonies only
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
  
  taskVarMeanbyMix <- rbind(taskVarMeanbyMixbyLine, taskVarMeanbyMixABonly) %>% group_by(n, Mix)
  
  # !!! NEW !!! change 
  # 1) "Line" -> "Group"
  colnames(taskVarMeanbyMix)[colnames(taskVarMeanbyMix) == "Line"] <- "Group"
  colnames(taskVarMeanbyrep)[colnames(taskVarMeanbyrep) == "Line"] <- "Group"
  colnames(taskdist)[colnames(taskdist) == "Line"] <- "Group"
  
  # 2) Change category names
  taskVarMeanbyMix$Group[taskVarMeanbyMix$Group == "A"] <- "Group A"
  taskVarMeanbyMix$Group[taskVarMeanbyMix$Group == "B"] <- "Group B"
  taskVarMeanbyMix$Group[taskVarMeanbyMix$Group == "AB"] <- "Mixed"
  
  taskVarMeanbyMix$Mix <- as.character(taskVarMeanbyMix$Mix)
  taskVarMeanbyMix$Mix[taskVarMeanbyMix$Mix == "A"] <- "Group A"
  taskVarMeanbyMix$Mix[taskVarMeanbyMix$Mix == "B"] <- "Group B"
  taskVarMeanbyMix$Mix[taskVarMeanbyMix$Mix == "AB"] <- "Mixed"
  
  taskVarMeanbyrep$Group[taskVarMeanbyrep$Group == "A"] <- "Group A"
  taskVarMeanbyrep$Group[taskVarMeanbyrep$Group == "B"] <- "Group B"
  taskVarMeanbyrep$Group[taskVarMeanbyrep$Group == "AB"] <- "Mixed"
  
  taskVarMeanbyrep$Mix <- as.character(taskVarMeanbyrep$Mix)
  taskVarMeanbyrep$Mix[taskVarMeanbyrep$Mix == "A"] <- "Group A"
  taskVarMeanbyrep$Mix[taskVarMeanbyrep$Mix == "B"] <- "Group B"
  taskVarMeanbyrep$Mix[taskVarMeanbyrep$Mix == "AB"] <- "Mixed"
  
  taskdist$Group[taskdist$Group == "A"] <- "Group A"
  taskdist$Group[taskdist$Group == "B"] <- "Group B"
  taskdist$Group[taskdist$Group == "AB"] <- "Mixed"
  
  taskdist$Mix <- as.character(taskdist$Mix)
  taskdist$Mix[taskdist$Mix == "A"] <- "Group A"
  taskdist$Mix[taskdist$Mix == "B"] <- "Group B"
  taskdist$Mix[taskdist$Mix == "AB"] <- "Mixed"
  
  return(list(taskVarMeanbyrep, taskVarMeanbyMix))
  
}

##############
# Make plots #
##############

params <- matrix(c(2, 2, 2, 2, 0.6,	0.6, 10, 10, 10, 10),
                 nrow = 1, ncol = 10, byrow = TRUE)

# Plotting
ymax <- 0.5 # max y for plotting
yinc <- 0.1 # y-axis increments
figH <- 1.5 # figure height for printing; default width is 3


make_dd_dist3 <- function(Aalpha1, Aalpha2, Balpha1, Balpha2, delta1, delta2, AThreshM1, AThreshM2, AThreshM2, BThreshM1, BThreshM2)
for (INDEX in 1:nrow(params)){
  
  ########################
  # Set global variables #
  ########################
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
  B_ThreshSD     <- B_ThreshM * 0.5 #population threshold standard deviations for clone line B !!Change!!
  InitialStim    <- c(0, 0) #intital vector of stimuli
  deltas         <- c(params[INDEX,5], params[INDEX,6]) #vector of stimuli increase rates  
  threshSlope    <- 7 #exponent parameter for threshold curve shape
  alpha          <- m
  # A_alpha        <- c(m, m*3) #efficiency of task performance
  # B_alpha        <- c(m*3, m)
  A_alpha        <- c(params[INDEX,1], params[INDEX,2])
  B_alpha        <- c(params[INDEX,3], params[INDEX,4])
  quitP          <- c(0.2, 0.2) #probability of quitting task once active !!Change!!
  
  #############
  # Load data #
  #############
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
  
  load(paste0("output/Rdata/", file_name, "reps_100.Rdata"))
  
  ############################
  # Prepare data for gg_dist #
  ############################
  task_dist_byrep <- prepdata_dist(task_dist)[[1]]
  task_dist_byMix <- prepdata_dist(task_dist)[[2]]

  }

###############
# Make panels #
###############

gg_dist3_a <- make_gg_dist3(params[1,])
gg_dist3_b <- make_gg_dist3(params[1,])
gg_dist3_c <- make_gg_dist3(params[1,])
gg_dist3_d <- make_gg_dist3(params[1,])


