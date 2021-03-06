################################################################################
#
# Model incorporating both thresholds and network dynamics
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
Ns             <- c(1, 2, 4, 6, 8, 12, 16) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 2000 #number of generations to run simulation 
corrStep       <- 100 #number of time steps for calculation of correlation 
reps           <- 2 #number of replications per simulation (for ensemble) !!Change!!
epsilon        <- 0.4 #weighting of social information relative to task demand information (must be between 0 and 1) 
AvgUpdateTime  <- 1 #average number of steps before a node does an update

# Threshold Parameters
ThreshM        <- c(10, 10) #population threshold means !!Change!!
ThreshSD       <- ThreshM * 0.1 #population threshold standard deviations !!Change!!
InitialStim    <- c(0, 0) #intital vector of stimuli
StimRates      <- c(0.8, 0.6) #vector of stimuli increase rates  
ExhaustM       <- 5 #population exhaustion threshold mean !!Change!!
ExhaustSD      <- ExhaustM * 0.1 #populution exhaustion threshold standard deviation !!Change!!
ExhaustRate    <- 0.6 #rate of exhaustion from performing task  
TaskExhaust   <- c(8, 1) #scales relative exhaustion rate of the two tasks  
threshSlope    <- 2 #exponent parameter for threshold curve shape 
alpha          <- m #efficiency of task performance

# Social Network Parameters
p              <- 0.2 #probability of interacting with individual in other states
q              <- 5 #probability of interacting with individual in same state relative to others
NeighborThresh <- FALSE #should stimulus be counted in the social information probabilities? 



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
    
    # Seed exhaustion stimuli
    exhaustStim <- matrix(rep(0, n), ncol = 1, dimnames = list(paste0("v-", 1:n), "ExhaustStim"))
    
    # Seed internal thresholds
    threshMat <- seedThresholds(n = n, 
                                m = m, 
                                ThresholdMeans = ThreshM, 
                                ThresholdSDs = ThreshSD)

    # Seed exhaustion thresholds
    exhaustThresh <- seedExhaustThresh(n = n, 
                                     ExhaustMean = ExhaustM, 
                                     ExhaustSD = ExhaustSD)
    
    
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
      # Update Exhaustion Stimuli
      exhaustStim <- updateExhaustStim(ExhaustStim = exhaustStim, 
                                       TaskExhaustVect = TaskExhaust, 
                                       ExhaustRate = ExhaustRate, 
                                       StateMat = X_g)
      # Update social network
      g_adj <- temporalNetwork(X_sub_g = X_g,
                               p = p, 
                               bias = q)
      # Calculate task demand based on social interactions
      Lmat <- calcSocialProbMat(SocialNetwork = g_adj, 
                                UseThreshold = NeighborThresh,
                                ThresholdMatrix = socialThreshMat, 
                                X_sub_g = X_g)
      # Calculate task demand based on global stimuli
      Tmat <- calcThresholdProbMat(TimeStep = t + 1, # first row is generation 0
                                   ThresholdMatrix = threshMat, 
                                   StimulusMatrix = stimMat, 
                                   nSlope = threshSlope)
      # Calculate rest demand based on exhaustion stimuli
      Emat <- calcExhaustDemand(ExhaustStim = exhaustStim, 
                                ExhaustThreshVector = exhaustThresh, 
                                nSlope = threshSlope)
      # Update Pg
      P_g <- updateProbMatrix(Lmatrix = Lmat, 
                              Tmatrix = Tmat, 
                              epsilon = epsilon,
                              Ematrix = Emat)
      # Update task performance
      X_g <- updateTaskPerformance(P_sub_g = P_g,
                                   TaskMat = X_g,
                                   UpdateTime = AvgUpdateTime)

      
      # Capture current task performance tally
      tally <- matrix(c(t, colSums(X_g)), ncol = ncol(X_g) + 1)
      colnames(tally) <- c("t", colnames(X_g))
      tally <- transform(tally, Inactive = n - sum(X_g), n = n, replicate = sim)
      taskTally[[t]] <- tally
      
      # Create time step 
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
      # Update total task performance profile
      X_tot <- X_tot + X_g
    }
    
    # Calculate Entropy
    entropy <- mutualEntropy(TotalStateMat = X_tot)
    entropy <- transform(entropy, n = n, replicate = sim)
    
    # Calculate total task distribution
    # totalTaskDist <- X_tot / rowSums(X_tot)
    totalTaskDist <- X_tot / gens
    totalTaskDist <- transform(totalTaskDist, Inactive = gens - rowSums(X_tot), n = n, replicate = sim)
    
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
    print(paste0("DONE: N = ", n, ", Simulation ", sim))
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
  groups_taskDist[[i]]  <- ens_taskDist
  groups_taskCorr[[i]]  <- runCorrs
  groups_taskStep[[i]]  <- ens_taskStep
  groups_taskTally[[i]] <- ens_taskTally
  groups_stim[[i]]      <- ens_stim
  groups_entropy[[i]]   <- ens_entropy
  
}

# trim out correlations for group size 1
if(1 %in% Ns) {
  groups_taskCorr <- groups_taskCorr[-1]
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



