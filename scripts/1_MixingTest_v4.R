################################################################################
#
# Modified version of 1_MixingTest.R, allowing for variation in threshSlope, quitP
#
################################################################################

rm(list = ls())
# params <- matrix(c(2,2,2,2,0.6,0.6,
#                    2,2,6,6,0.6,0.6, 
#                    2,6,6,2,0.6,0.6,
#                    2,2,1,1,0.6,0.6,
#                    2,1,1,2,0.6,0.6,
#                    2,2,2,2,0.6,1.0,
#                    2,2,2,2,0.6,0.2),
#                  nrow = 7, ncol = 6, byrow = TRUE)
# params <- matrix(c(2, 2, 6, 6, 0.6, 1.0,
#                    2, 6, 6, 2, 0.6, 1.0,
#                    2, 2, 6, 6, 1.0, 0.6,
#                    2,	6, 6, 2, 1.0,	0.6,
#                    2,	2, 1,	1, 0.6, 0.6,
#                    2,	1, 1,	2, 0.6, 1.0,
#                    2,	2, 1,	1, 0.6,	1.0,
#                    2,	1, 1,	2, 0.6,	1.0),
#                  nrow = 8, ncol = 6, byrow = TRUE)
# params <- matrix(c(2, 4, 1, 3, 0.6,	1.0,
#                    2,	4, 1,	3, 1.0,	0.6),
#                  nrow = 2, ncol = 6, byrow = TRUE)
# params <- matrix(c(6, 6, 2, 2, 0.6,	0.6, 10, 10, 15, 15,
#                    6, 6, 2, 2, 0.6,	0.6, 10, 15, 15, 10,
#                    6, 6, 2, 2, 0.6,	0.6, 15, 15, 10, 10),
#                    nrow = 3, ncol = 10, byrow = TRUE)

# params <- matrix(c(2, 2, 6, 6, 1,	1, 10, 10, 10, 10,
#                    2, 6, 6, 2, 1,	1, 10, 10, 10, 10,
#                    2, 2, 1, 1, 1,	1, 10, 10, 10, 10,
#                    2, 1, 1, 2, 1,	1, 10, 10, 10, 10),
#                  nrow = 4, ncol = 10, byrow = TRUE)

# params <- matrix(c(3, 3, 1, 1, 0.6,	0.6, 10, 10, 10, 10,
#                    3,	1, 1,	3, 0.6,	0.6, 10, 10, 10, 10),
#                  nrow = 2, ncol = 10, byrow = TRUE)

params <- matrix(c(2, 2, 6, 6, 0.6,	0.6, 10, 10, 10, 10),
                 nrow = 1, ncol = 10, byrow = TRUE)

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
  
  
  # Loop through group sizes
  for (i in 1:length(Ns)) {
    
    # Set group size
    n <- Ns[i]
    
    # Loop through clonals mixes
    mix_taskDist <- list()
    mix_taskCorr <- list()
    
    for (mix_index in 1:length(mixes)) {
      
      # Set mix
      mix <- mixes[mix_index]
      
      # Set stim levels (for different larvae)
      # if (mix == "A") {
      #   StimRates      <- c(0.6, 0.6)
      # } else if (mix == "AB") {
      #   StimRates      <- c(0.6, 0.6) * 1.3
      # } else if (mix == "B") {
      #   StimRates      <- c(0.6, 0.6) * 1.4
      # }
      
      # Set work efficiency by lines # updated 11/15/18
      if (mix == "A") {
        alpha      <- matrix(rep(A_alpha, n), ncol = m, byrow = T)
      } else if (mix == "AB") {
        input      <- c( rep(A_alpha, n/2), rep(B_alpha, n/2) ) 
        alpha      <- matrix(input, ncol = m, byrow = T)
      } else if (mix == "B") {
        alpha      <- matrix(rep(B_alpha, n), ncol = m, byrow = T)
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
        if (mix == "A") {
          threshMat <- seedThresholds(n = n, 
                                      m = m, 
                                      ThresholdMeans = A_ThreshM, 
                                      ThresholdSDs = A_ThreshSD)
          rownames(threshMat) <- paste0("A-", rownames(threshMat))
        } else if(mix == "B") {
          threshMat <- seedThresholds(n = n, 
                                      m = m, 
                                      ThresholdMeans = B_ThreshM, 
                                      ThresholdSDs = B_ThreshSD)
          rownames(threshMat) <- paste0("B-", rownames(threshMat))
        } else if(mix == "AB") {
          threshMatA <- seedThresholds(n = n / 2, 
                                       m = m, 
                                       ThresholdMeans = A_ThreshM, 
                                       ThresholdSDs = A_ThreshSD)
          rownames(threshMatA) <- paste0("A-", rownames(threshMatA))
          threshMatB <- seedThresholds(n = n / 2, 
                                       m = m, 
                                       ThresholdMeans = B_ThreshM, 
                                       ThresholdSDs = B_ThreshSD)
          rownames(threshMatB) <- paste0("B-", rownames(threshMatB))
          threshMat <- rbind(threshMatA, threshMatB)
          rm(threshMatA, threshMatB)
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
          if (mix == "A") {
            X_g <- updateTaskPerformance(P_sub_g    = P_g,
                                         TaskMat    = X_g,
                                         QuitProb   = quitP[1])  
          } else if(mix == "B") {
            X_g <- updateTaskPerformance(P_sub_g    = P_g,
                                         TaskMat    = X_g,
                                         QuitProb   = quitP[2])  
          } else if(mix == "AB") {
            X_g_A <- updateTaskPerformance(P_sub_g  = P_g[1:(n/2), ],
                                           TaskMat  = X_g[1:(n/2), ],
                                           QuitProb = quitP[1]) 
            X_g_B <- updateTaskPerformance(P_sub_g  = P_g[(n/2+1):nrow(P_g), ],
                                           TaskMat  = X_g[(n/2+1):nrow(P_g), ],
                                           QuitProb = quitP[2])
            rownames(X_g_B) <- paste0("v-", (n/2+1):nrow(P_g))
            X_g <- rbind(X_g_A, X_g_B)
            rm(X_g_A, X_g_B)
          }
          
          # Capture current task performance tally
          tally <- matrix(c(t, colSums(X_g)), ncol = ncol(X_g) + 1)
          colnames(tally) <- c("t", colnames(X_g))
          tally <- transform(tally, Inactive = n - sum(X_g), n = n, replicate = sim)
          taskTally[[t]] <- tally
          
          # Update total task performance profile
          X_tot <- X_tot + X_g
          X_tot
          
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
  
}
