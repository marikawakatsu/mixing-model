##################################################
#
# Global Stimulus Functions 
#
##################################################


####################
# Seed Stimuls
####################
seedStimuls <- function(InitialSVector, gens) {
  # Calculate number of blank spots to make
  repLength <- length(InitialSVector) * gens #intiial row does not count as gen
  # Build matrix
  stim <- matrix(data = c(InitialSVector, rep(NA, repLength)),
                 byrow = TRUE, 
                 nrow = (gens + 1))
  # Fix Names
  colnames(stim) <- paste0(rep("s", length(InitialSVector)), 1:length(InitialSVector))
  rownames(stim) <- paste0("Gen", 0:gens)
  # Return
  return(stim)
}

####################
# Stimulus Level
####################
# Update stimulus levels
update_stim <- function(stim_matrix, deltas, alpha, state_matrix, time_step) {
  # Get preliminary info
  n <- nrow(state_matrix)
  stim_values <- stim_matrix[time_step, ]
  # Calculate
  effective_act <- state_matrix / n
  effective_act <- effective_act * alpha
  active_amount <- colSums(effective_act)
  
  new_values <- stim_values + deltas - active_amount
  new_values[new_values < 0] <- 0
  # Insert
  stim_matrix[time_step + 1, ] <- new_values
  return(stim_matrix)
}

# Frequency dependent
globalStimUpdate <- function(stimulus, delta, alpha, Ni, n) {
  # Calculate
  s <- stimulus + delta - ( alpha * ( Ni / n ))
  # If negative, make zero
  if(s < 0.0001) {
    s <- 0
  }
  return(s)
}

# Density dependent (per capita)
globalStimUpdate_PerCap <- function(stimulus, delta, alpha, Ni, n, m) {
  # Calculate
  s <- stimulus + (alpha * (delta ) * (n / m)) - (Ni * alpha)
  # s <- stimulus + (alpha * (delta ) * (1 / (1 + quitP)) * (n / m)) - (Ni * alpha)
  # If negative, make zero
  if(s < 0.0001) {
    s <- 0
  }
  return(s)
}
