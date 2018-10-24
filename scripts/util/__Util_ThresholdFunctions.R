##################################################
#
# Threshold Probability of Performance 
#
##################################################


####################
# Seed task thresholds
####################
seedThresholds <- function(n, m, ThresholdMeans = NULL, ThresholdSDs = NULL) {
  # Loop through tasks and sample thresholds from normal dist
  threshMat <- lapply(1:length(ThresholdMeans), function(i) {
    threshList <- rtnorm(n = n, 
                         mean = ThresholdMeans[i], 
                         sd = ThresholdSDs[i], 
                         lower = 0)
    return(threshList)
  })
  threshMat <- do.call("cbind", threshMat)
  # Fix names
  colnames(threshMat) <- paste0("Thresh", 1:length(ThresholdMeans))
  rownames(threshMat) <- paste0("v-", 1:n)
  # Return
  return(threshMat)
}


####################
# Threshold function
####################
threshProb <- function(s, phi, nSlope) {
  T_vi <- (s^nSlope) / (s^nSlope + phi^nSlope)
}


####################
# Calculate Threshold
####################
calcThresholdProbMat <- function(TimeStep, ThresholdMatrix, StimulusMatrix, nSlope) {
  # select proper stimulus for this time step
  stimulusThisStep <- StimulusMatrix[TimeStep, ]
  # calculate threshold probabilities for one individual
  thresholdP <- lapply(1:nrow(ThresholdMatrix), function(i) {
    # select row for individual in threshold matrix
    indThresh <- ThresholdMatrix[i, ]
    # create task vector to be output and bound
    taskThresh <- rep(NA, length(indThresh))
    # loop through each task within individual
    for (j in 1:length(taskThresh)) {
      taskThresh[j] <- threshProb(s = stimulusThisStep[j], phi = indThresh[j], nSlope = nSlope)
    }
    return(taskThresh)
  })
  # bind and return
  thresholdP <- do.call("rbind", thresholdP)
  thresholdP[is.na(thresholdP)] <- 0 #fix NAs where thresh was 0 and stim was 0
  colnames(thresholdP) <- paste0("ThreshProb", 1:ncol(thresholdP))
  rownames(thresholdP) <- paste0("v-", 1:nrow(thresholdP))
  return(thresholdP)
}

####################
# Output Threshold Demands
####################
calcThresholdDetermMat <- function(TimeStep, ThresholdMatrix, StimulusMatrix) {
  # select proper stimulus for this time step
  stimulusThisStep <- StimulusMatrix[TimeStep, ]
  # calculate threshold probabilities for one individual
  thresholdP <- lapply(1:nrow(ThresholdMatrix), function(i) {
    # select row for individual in threshold matrix
    indThresh <- ThresholdMatrix[i, ]
    # create task vector to be output and bound
    taskThresh <- rep(0, length(indThresh))
    # loop through each task within individual
    for (j in 1:length(taskThresh)) {
      stim <- stimulusThisStep[j]
      thresh <- indThresh[j]
      if (stim > thresh) {
        taskThresh[j] <- 1
      }
    }
    return(taskThresh)
  })
  # bind and return
  thresholdP <- do.call("rbind", thresholdP)
  colnames(thresholdP) <- paste0("ThreshProb", 1:ncol(thresholdP))
  rownames(thresholdP) <- paste0("v-", 1:nrow(thresholdP))
  return(thresholdP)
}

####################
# Adjust thresholds due to social interactions
####################
adjust_thresh_social <- function(social_network, threshold_matrix, state_matrix, epsilon) {
  # Calculate "sum" of task states/probs of neighbors
  active_neighbors <- t(social_network) %*% state_matrix #active neighbors by category
  total_neighbors <- rowSums(active_neighbors) #total active neighbors
  # Calculate interacting individuals NOT performing that task
  not_sums <- total_neighbors - active_neighbors
  # Calculate net threshold effect
  net_effect <- not_sums - active_neighbors
  net_effect <- net_effect * epsilon
  # Adjust thresholds and return
  threshold_matrix <- threshold_matrix + net_effect
  # Catch minimum thresholds
  threshold_matrix[threshold_matrix < 0] <- 0
  # Return
  return(threshold_matrix)
}

####################
# Adjust thresholds due to social interactions--with maximum level
####################
adjust_thresholds_social_capped <- function(social_network, threshold_matrix, state_matrix, epsilon, threshold_max) {
  # Calculate "sum" of task states/probs of neighbors
  active_neighbors <- t(social_network) %*% state_matrix #active neighbors by category
  total_neighbors <- rowSums(active_neighbors) #total active neighbors
  # Calculate interacting individuals NOT performing that task
  not_sums <- total_neighbors - active_neighbors
  # Calculate net threshold effect
  net_effect <- not_sums - active_neighbors
  net_effect <- net_effect * epsilon
  # Adjust thresholds and return
  threshold_matrix <- threshold_matrix + net_effect
  # Catch minimum and maximum thresholds
  threshold_matrix[threshold_matrix < 0] <- 0
  threshold_matrix[threshold_matrix > threshold_max] <- threshold_max
  # Return
  return(threshold_matrix)
}

