################################################################################
#
# Model incorporating both thresholds and network dynamics
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")
source("scripts/__Util_TaskStateFunctions_OriginalFormat.R")

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
Ns             <- c(4, 16) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 2000 #number of generations to run simulation 
corrStep       <- 100 #number of time steps for calculation of correlation 
reps           <- 5 #number of replications per simulation (for ensemble) 
epsilon        <- 0 #weighting of social information relative to task demand information (must be between 0 and 1) 
AvgUpdateTime  <- 1 #average number of steps before a node does an update

# Threshold Parameters
mixes          <- c("A", "B", "AB")
A_ThreshM      <- c(10, 10) #population threshold means for clone line A !!Change!!
A_ThreshSD     <- A_ThreshM * 0 #population threshold standard deviations for clone line A !!Change!!
B_ThreshM      <- c(10, 10) #population threshold means for clone line B !!Change!!
B_ThreshSD     <- B_ThreshM * 0 #population threshold standard deviations for clone line B !!Change!!
InitialStim    <- c(0, 0) #intital vector of stimuli
StimRates      <- c(0.6, 0.6) #vector of stimuli increase rates  
threshSlope    <- 2 #exponent parameter for threshold curve shape  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active
phi            <- 0.01 #learning/forgetting of threshold
c              <- 0.3 #social reinforcer of threshold

# Social Network Parameters
p              <- 0.2 #probability of interacting with individual
q              <- 2 #probability of interacting with individual in same state relative to others




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
  
  # Loop through clonals mixes
  mix_taskDist <- list()
  
  for (mix_index in 1:length(mixes)) {
    
    # Set mix
    mix <- mixes[mix_index]
    
    # Set group size
    n <- Ns[i]
    
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
                             RateVector = StimRates, 
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
      taskCorr <- list()
      taskStep <- list()
      taskTally <- list()
      
      ####################
      # Simulate
      ####################
      # Run simulation
      for (t in 1:gens) {
        # Update stimuli
        for (j in 1:(ncol(stimMat)/2)) {
          # update stim
          stimMat[t + 1, j] <- globalStimUpdate(stimulus = stimMat[t, j],
                                                delta = stimMat[t, j + m], 
                                                alpha = alpha, 
                                                Ni = sum(X_g[ , j]), 
                                                n = n)
          # shift down delta (rate increases)
          stimMat[t + 1, j + m] <- stimMat[t, j + m]
        }
        
        # Calculate task demand based on global stimuli
        Tmat <- calcThresholdProbMat(TimeStep = t + 1, # first row is generation 0
                                     ThresholdMatrix = threshMat, 
                                     StimulusMatrix = stimMat, 
                                     nSlope = threshSlope)
        # Update task performance
        X_g <- updateTaskPerformance_og(P_sub_g    = Tmat,
                                        TaskMat    = X_g,
                                        QuitProb   = quitP)
        # Update social network
        g_adj <- temporalNetwork(X_sub_g = X_g,
                                 p = p, 
                                 bias = q)
        # Update threshold (reinforcement/forgetting)
        threshMat <- adjustThresholdsSocial(ThresholdMatrix = threshMat, 
                                            X_sub_g = X_g, 
                                            phi = phi,
                                            SocialNetwork = g_adj)
        
        # Capture current task performance tally
        tally <- matrix(c(t, colSums(X_g)), ncol = ncol(X_g) + 1)
        colnames(tally) <- c("t", colnames(X_g))
        tally <- transform(tally, Inactive = n - sum(X_g), n = n, replicate = sim)
        taskTally[[t]] <- tally
        
        # Create time step 
        if (t %% corrStep == 0) {
          # Get tasks performance in correlation step
          X_step <- X_tot - X_prev
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
          # Update previous step matrix
          X_prev <- X_tot
        }
        # Update total task performance profile
        X_tot <- X_tot + X_g
      }
      
      # Calculate Entropy
      entropy <- mutualEntropy(TotalStateMat = X_tot)
      entropy <- transform(entropy, n = n, replicate = sim)
      
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
      taskTally <- do.call("rbind", taskTally)
      
      # Create tasktally table
      stimMat <- transform(stimMat, n = n, replicate = sim)
      
      # Add total task distributions, entropy values, and graphs to lists
      ens_taskDist[[sim]]  <- totalTaskDist
      ens_entropy[[sim]]   <- entropy
      ens_taskCorr[[sim]]  <- taskCorr
      ens_taskTally[[sim]] <- taskTally
      ens_taskStep[[sim]]  <- taskStep
      ens_stim[[sim]]      <- stimMat
      
      # Print simulation completed
      print(paste0("DONE: Mix = ", mix,  ", N = ", n, ", Simulation ", sim))
    }
    
    # Calculate mean correlation for each run
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
    runCorrs <- transform(runCorrs, n = n)
    
    # Add to list of lists
    mix_taskDist[[mix_index]]  <- ens_taskDist
    # groups_taskCorr[[i]]  <- runCorrs
    # groups_taskStep[[i]]  <- ens_taskStep
    # groups_taskTally[[i]] <- ens_taskTally
    # groups_stim[[i]]      <- ens_stim
    # groups_entropy[[i]]   <- ens_entropy
    
  }
  
  groups_taskDist[[i]]  <- mix_taskDist
}


# trim out correlations for group size 1
if(1 %in% Ns) {
  drop <- length(Ns)
  groups_taskCorr <- groups_taskCorr[-(1:drop)]
}


####################
# Save run
####################
# rename
# ep_08nobreak_taskDist <- ens_taskDist
# ep_08nobreak_entropy <- ens_entropy
# filename <- "epsilon08nobreak"
# save(ep_08nobreak_taskDist, 
#      ep_08nobreak_entropy, 
#      file = paste0("output/Rdata/", filename,".Rdata"))



