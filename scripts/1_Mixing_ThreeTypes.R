################################################################################
#
# Modified version of 1_MixingTest_v2.R, allowing for mixes of three types of individuals
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters 
Ns             <- c(18) #vector of number of individuals to simulate NOTE: We chose 18 to allow 6 of each type.
m              <- 2 #number of tasks
gens           <- 10000 #number of generations to run simulation 
corrStep       <- 200 #number of time steps for calculation of correlation 
reps           <- 100 #number of replications per simulation (for ensemble) !!Change!!

# Threshold Parameters
mix_ratios     <- data.frame("A" = c(1, 0, 0, 0.5, 0.5, 0,   1/3),
                             "B" = c(0, 1, 0, 0.5, 0,   0.5, 1/3),
                             "C" = c(0, 0, 1, 0,   0.5, 0.5, 1/3)) # %of each type
A_ThreshM      <- c(10, 10)        # population threshold means for clone line A
A_ThreshSD     <- A_ThreshM * 0.1  # population threshold standard deviations for clone line A (DON'T change)
B_ThreshM      <- c(15, 15)        # population threshold means for clone line B 
B_ThreshSD     <- B_ThreshM * 0.1  # population threshold standard deviations for clone line B (DON'T change)
C_ThreshM      <- c(7.5, 7.5)        # population threshold means for clone line A
C_ThreshSD     <- C_ThreshM * 0.1  # population threshold standard deviations for clone line C (DON'T change)
InitialStim    <- c(0, 0)          # intital vector of stimuli
deltas         <- c(0.6, 0.6)      # vector of stimuli increase rates  
threshSlope    <- 7                # exponent parameter for threshold curve shape (DON'T change)
A_alpha        <- c(2, 2)          # efficiency of task performance for A type
B_alpha        <- c(3, 3)          # efficiency of task performance for B type
C_alpha        <- c(1.5, 1.5)          # efficiency of task performance for C type
quitP          <- c(0.2, 0.2, 0.2) # probability of quitting task once active for each line A, B, C (DON'T change)

file_name <- sprintf("Mix_ThreeTypes_AThreshM_%1.2f_%1.2f_BThreshM_%1.2f_%1.2f_CThreshM_%1.2f_%1.2f_deltas_%1.2f_%1.2f_threshSlope_%d_Aalpha_%1.2f_%1.2f_Balpha_%1.2f_%1.2f_Calpha_%1.2f_%1.2f_quitP_%1.2f",
                      A_ThreshM[1], A_ThreshM[2], B_ThreshM[1], B_ThreshM[2], C_ThreshM[1], C_ThreshM[2],
                      deltas[1], deltas[2], threshSlope, 
                      A_alpha[1], A_alpha[2], B_alpha[1], B_alpha[2], C_alpha[1], C_alpha[2],
                      quitP[1])  # note quitp[1] = quitP[2]

####################
# Run simulation multiple times
####################
# Prep meta-lists for collection of group size simulations
groups_taskDist  <- list()
groups_taskCorr  <- list()
groups_taskStep  <- list()
groups_taskTally <- list()
groups_stim      <- list()
groups_entropy   <- list()

# Functions to set task update function based on mix
updateTask_ABC <- function(n_A, n_B, n_C, P_g, TaskMat, QuitProb) {
  X_g_A <- updateTaskPerformance(P_sub_g  = P_g[1:(n_A), ],
                                 TaskMat  = X_g[1:(n_A), ],
                                 QuitProb = quitP[1]) 
  X_g_B <- updateTaskPerformance(P_sub_g  = P_g[(n_A+1):(n_A+n_B), ],
                                 TaskMat  = X_g[(n_A+1):(n_A+n_B), ],
                                 QuitProb = quitP[2])
  X_g_C <- updateTaskPerformance(P_sub_g  = P_g[(n_A+1+n_B):nrow(P_g), ],
                                 TaskMat  = X_g[(n_A+1+n_B):nrow(X_g), ],
                                 QuitProb = quitP[3])
  rownames(X_g_B) <- paste0("v-", (n_A+1):(n_A+n_B))
  rownames(X_g_C) <- paste0("v-", (n_A+n_B+1):nrow(P_g))
  X_g <- rbind(X_g_A, X_g_B)
  X_g <- rbind(X_g, X_g_C)
  return(X_g)
}

updateTask_AB <- function(n_A, n_B, n_C, P_g, TaskMat, QuitProb) {
  X_g_A <- updateTaskPerformance(P_sub_g  = P_g[1:(n_A), ],
                                 TaskMat  = X_g[1:(n_A), ],
                                 QuitProb = quitP[1]) 
  X_g_B <- updateTaskPerformance(P_sub_g  = P_g[(n_A+1):nrow(P_g), ],
                                 TaskMat  = X_g[(n_A+1):nrow(X_g), ],
                                 QuitProb = quitP[2])
  rownames(X_g_B) <- paste0("v-", (n_A+1):nrow(P_g))
  X_g <- rbind(X_g_A, X_g_B)
  return(X_g)
}

updateTask_AC <- function(n_A, n_B, n_C, P_g, TaskMat, QuitProb) {
  X_g_A <- updateTaskPerformance(P_sub_g  = P_g[1:(n_A), ],
                                 TaskMat  = X_g[1:(n_A), ],
                                 QuitProb = quitP[1]) 
  X_g_C <- updateTaskPerformance(P_sub_g  = P_g[(n_C+1):nrow(P_g), ],
                                 TaskMat  = X_g[(n_C+1):nrow(X_g), ],
                                 QuitProb = quitP[2])
  rownames(X_g_C) <- paste0("v-", (n_A+1):nrow(P_g))
  X_g <- rbind(X_g_A, X_g_C)
  return(X_g)
}

updateTask_BC <- function(n_A, n_B, n_C, P_g, TaskMat, QuitProb) {
  X_g_B <- updateTaskPerformance(P_sub_g  = P_g[1:(n_B), ],
                                 TaskMat  = X_g[1:(n_B), ],
                                 QuitProb = quitP[2])
  X_g_C <- updateTaskPerformance(P_sub_g  = P_g[(n_B+1):nrow(P_g), ],
                                 TaskMat  = X_g[(n_B+1):nrow(P_g), ],
                                 QuitProb = quitP[3]) 
  rownames(X_g_C) <- paste0("v-", (n_B+1):nrow(P_g))
  X_g <- rbind(X_g_B, X_g_C)
  return(X_g)
}
    
updateTask_A <- function(n_A, n_B, n_C, P_g, TaskMat, QuitProb) {
  X_g <- updateTaskPerformance(P_sub_g    = P_g,
                               TaskMat    = X_g,
                               QuitProb   = quitP[1])
  return(X_g)
}

updateTask_B <- function(n_A, n_B, n_C, P_g, TaskMat, QuitProb) {
  X_g <- updateTaskPerformance(P_sub_g    = P_g,
                               TaskMat    = X_g,
                               QuitProb   = quitP[2])
  return(X_g)
}

updateTask_C <- function(n_A, n_B, n_C, P_g, TaskMat, QuitProb) {
  X_g <- updateTaskPerformance(P_sub_g    = P_g,
                               TaskMat    = X_g,
                               QuitProb   = quitP[3])
  return(X_g)
}



# Loop through group sizes
for (i in 1:length(Ns)) {
  
  # Set group size
  n <- Ns[i]
  
  # Loop through clonals mixes
  mix_taskDist <- list()
  mix_taskCorr <- list()
  
  for (mix_index in 1:nrow(mix_ratios)) {
    
    # Set mix
    n_A <- n * mix_ratios[mix_index, "A"]
    n_B <- n * mix_ratios[mix_index, "B"]
    n_C <- n * mix_ratios[mix_index, "C"]
    if (n_A %% 1 != 0 | n_B %% 1 != 0 | n_C %% 1 != 0) {
      print("ERROR: Mixes provided create fractions of individuals")
      break
    }
    
    # Set work efficiency by lines # updated 11/15/18
    input      <- c( rep(A_alpha, n_A), rep(B_alpha, n_B), rep(C_alpha, n_C) ) 
    alpha      <- matrix(input, ncol = m, byrow = T)
    
    # Choose task performance function
    if (n_A != 0 & n_B != 0 & n_C != 0) {
      mix <- "ABC"
      updateTask_mix <- updateTask_ABC
    } else if (n_A != 0 & n_B != 0){
      mix <- "AB"
      updateTask_mix <- updateTask_AB
    } else if (n_A != 0 & n_C != 0){
      mix <- "AC"
      updateTask_mix <- updateTask_AC
    } else if (n_B != 0 & n_C != 0){
      mix <- "BC"
      updateTask_mix <- updateTask_BC
    } else if (n_A != 0) {
      mix <- "A"
      updateTask_mix <- updateTask_A
    } else if (n_B != 0) {
      mix <- "B"
      updateTask_mix <- updateTask_B
    } else {
      mix <- "C"
      updateTask_mix <- updateTask_C
    }
    
    # Prep lists for collection of simulation outputs
    ens_taskDist  <- list()
    ens_taskCorr  <- list()
    ens_taskStep  <- list()
    ens_taskTally <- list()
    ens_entropy   <- list()
    ens_stim      <- list()
    
    # Run Simulations
    for (sim in 1:reps) {
      # Set initial probability matrix (P_g)
      P_g <- initiateProbMatrix(n = n, m = m)
      
      # Seed task (external) stimuli
      stimMat <- seedStimuls(InitialSVector = InitialStim, 
                             gens = gens)
      
      # Seed internal thresholds 
      if (n_A != 0 & n_B != 0 & n_C != 0) {
        threshMatA <- seedThresholds(n = n_A, 
                                     m = m, 
                                     ThresholdMeans = A_ThreshM, 
                                     ThresholdSDs = A_ThreshSD)
        rownames(threshMatA) <- paste0("A-", rownames(threshMatA))
        threshMatB <- seedThresholds(n = n_B, 
                                     m = m, 
                                     ThresholdMeans = B_ThreshM, 
                                     ThresholdSDs = B_ThreshSD)
        rownames(threshMatB) <- paste0("B-", rownames(threshMatB))
        threshMat <- rbind(threshMatA, threshMatB)
        threshMatC <- seedThresholds(n = n_C, 
                                     m = m, 
                                     ThresholdMeans = C_ThreshM, 
                                     ThresholdSDs = C_ThreshSD)
        rownames(threshMatC) <- paste0("C-", rownames(threshMatC))
        threshMat <- rbind(threshMat, threshMatC)
        rm(threshMatA, threshMatB, threshMatC)
      } else if (n_A != 0 & n_B != 0){
        threshMatA <- seedThresholds(n = n_A, 
                                     m = m, 
                                     ThresholdMeans = A_ThreshM, 
                                     ThresholdSDs = A_ThreshSD)
        rownames(threshMatA) <- paste0("A-", rownames(threshMatA))
        threshMatB <- seedThresholds(n = n_B, 
                                     m = m, 
                                     ThresholdMeans = B_ThreshM, 
                                     ThresholdSDs = B_ThreshSD)
        rownames(threshMatB) <- paste0("B-", rownames(threshMatB))
        threshMat <- rbind(threshMatA, threshMatB)
      } else if (n_A != 0 & n_C != 0){
        threshMatA <- seedThresholds(n = n_A, 
                                     m = m, 
                                     ThresholdMeans = A_ThreshM, 
                                     ThresholdSDs = A_ThreshSD)
        rownames(threshMatA) <- paste0("A-", rownames(threshMatA))
        threshMatC <- seedThresholds(n = n_C, 
                                     m = m, 
                                     ThresholdMeans = C_ThreshM, 
                                     ThresholdSDs = C_ThreshSD)
        rownames(threshMatC) <- paste0("C-", rownames(threshMatC))
        threshMat <- rbind(threshMatA, threshMatC)
      } else if (n_B != 0 & n_C != 0){
        threshMatB <- seedThresholds(n = n_B, 
                                     m = m, 
                                     ThresholdMeans = B_ThreshM, 
                                     ThresholdSDs = B_ThreshSD)
        rownames(threshMatB) <- paste0("B-", rownames(threshMatB))
        threshMatC <- seedThresholds(n = n_C, 
                                     m = m, 
                                     ThresholdMeans = C_ThreshM, 
                                     ThresholdSDs = C_ThreshSD)
        rownames(threshMatC) <- paste0("C-", rownames(threshMatC))
        threshMat <- rbind(threshMatB, threshMatC)
      } else if (n_A != 0) {
        threshMat <- seedThresholds(n = n, 
                                    m = m, 
                                    ThresholdMeans = A_ThreshM, 
                                    ThresholdSDs = A_ThreshSD)
        rownames(threshMat) <- paste0("A-", rownames(threshMat))
      } else if (n_B != 0) {
        threshMat <- seedThresholds(n = n, 
                                    m = m, 
                                    ThresholdMeans = B_ThreshM, 
                                    ThresholdSDs = B_ThreshSD)
        rownames(threshMat) <- paste0("B-", rownames(threshMat))
      } else {
        threshMat <- seedThresholds(n = n, 
                                    m = m, 
                                    ThresholdMeans = C_ThreshM, 
                                    ThresholdSDs = C_ThreshSD)
        rownames(threshMat) <- paste0("C-", rownames(threshMat))
      }

      
      # Start task performance
      X_g <- matrix(data = rep(0, length(P_g)), ncol = ncol(P_g))
      
      # Create cumulative task performance matrix
      X_tot <- X_g
      
      # Prep correlation step matrix
      X_prev <- matrix(data = rep(0, n * m), ncol = m)
      X_prevTot <- matrix(data = rep(0, n * m), ncol = m)
      taskCorr <- list()
      taskStep <- list()
      taskTally <- list()
      
      ####################
      # Simulate
      ####################
      # Run simulation
      for (t in 1:gens) {
        # Update stimuli
        stimMat <- update_stim(stim_matrix = stimMat, 
                               deltas = deltas, 
                               alpha = alpha, 
                               state_matrix = X_g, 
                               time_step = t)
        # Calculate task demand based on global stimuli
        P_g <- calcThresholdProbMat(TimeStep = t + 1, # first row is generation 0
                                    ThresholdMatrix = threshMat, 
                                    StimulusMatrix = stimMat, 
                                    nSlope = threshSlope) #!!! Change!!!
        # Update task performance
        X_g <- updateTask_mix(n_A, n_B, n_C, P_g, TaskMat, QuitProb)
        
        # Capture current task performance tally
        tally <- matrix(c(t, colSums(X_g)), ncol = ncol(X_g) + 1)
        colnames(tally) <- c("t", colnames(X_g))
        tally <- transform(tally, Inactive = n - sum(X_g), n = n, replicate = sim)
        taskTally[[t]] <- tally
        
        # Update total task performance profile
        X_tot <- X_tot + X_g
        
        # Create time step for correlation
        if (t %% corrStep == 0) {
          # Get tasks performance in correlation step
          X_step <- X_tot - X_prevTot
          # Add to ensemble list of task steps
          taskStep[[t / corrStep]] <- X_step
          # Calculate rank correlation if it is not the first step
          if(sum(X_prev) != 0) {
            # Normalize
            stepNorm <- X_step / rowSums(X_step)
            prevNorm <- X_prev / rowSums(X_prev)
            # Calculate ranks
            step_ranks <- calculateTaskRank(TaskStepMat = X_step)
            prev_ranks <- calculateTaskRank(TaskStepMat = X_prev)
            # Calculate Correlation
            rankCorr <- cor(prev_ranks, step_ranks, method = "spearman")
            # Put in list
            taskCorr[[(t / corrStep) - 1]] <- diag(rankCorr)
            names(taskCorr)[(t / corrStep) - 1] <- paste0("Gen", t)
          }
          # Update previous step total matrix
          X_prevTot <- X_tot
          # Update previous step total matrix
          X_prev <- X_step
        }
      }
      
      # Calculate Entropy
      # entropy <- mutualEntropy(TotalStateMat = X_tot)
      # entropy <- transform(entropy, n = n, replicate = sim)
      
      # Calculate total task distribution
      # totalTaskDist <- X_tot / rowSums(X_tot)
      totalTaskDist <- X_tot / gens
      totalTaskDist <- transform(totalTaskDist, Inactive = gens - rowSums(X_tot), 
                                 n = n, 
                                 replicate = sim, 
                                 Line = rownames(threshMat),
                                 Mix = mix)
      totalTaskDist$Line <- gsub("([A-Z]).*", "\\1", totalTaskDist$Line)
      
      # Create tasktally table
      # taskTally <- do.call("rbind", taskTally)
      
      # Create tasktally table
      # stimMat <- transform(stimMat, n = n, replicate = sim)
      
      # Add total task distributions, entropy values, and graphs to lists
      ens_taskDist[[sim]]  <- totalTaskDist
      ens_taskCorr[[sim]]  <- taskCorr
      # ens_entropy[[sim]]   <- entropy
      # ens_taskTally[[sim]] <- taskTally
      # ens_taskStep[[sim]]  <- taskStep
      # ens_stim[[sim]]      <- stimMat
      
      # Print simulation completed
      print(paste0("DONE: Mix = ", mix,  ", N = ", n, ", Simulation ", sim))
    }
    
    # Calculate mean correlation for each n
    runCorrs <- lapply(ens_taskCorr, function(x) {
      # Unlist
      runs <- do.call("rbind", x)
      # Calculate mean
      runMean <- matrix(data = rep(NA, m), ncol =  m)
      for (column in 1:m) {
        runMean[ , column] <- mean(runs[ , column], na.rm = TRUE)
      }
      colnames(runMean) <- colnames(runs)
      return(runMean)
    })
    runCorrs <- do.call("rbind", runCorrs)
    runCorrs <- transform(runCorrs, n = n, Mix = mix)
    
    # Add to list of lists
    mix_taskDist[[mix_index]]  <- ens_taskDist
    mix_taskCorr[[mix_index]]  <- runCorrs
    # groups_taskStep[[i]]  <- ens_taskStep
    # groups_taskTally[[i]] <- ens_taskTally
    # groups_stim[[i]]      <- ens_stim
    # groups_entropy[[i]]   <- ens_entropy
    
  }
  
  groups_taskDist[[i]]  <- mix_taskDist
  groups_taskCorr[[i]]  <- mix_taskCorr
}

# Bind and process
task_dist <- unlist(groups_taskDist, recursive = FALSE)
task_dist <- unlist(task_dist, recursive = FALSE)
task_dist <- do.call("rbind", task_dist)

task_corr <- unlist(groups_taskCorr, recursive = FALSE)
task_corr <- do.call("rbind", task_corr)


####################
# Save run
####################
save(task_dist, task_corr, file = paste0("output/Rdata/", file_name, ".Rdata"))



