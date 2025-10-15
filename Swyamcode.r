# --- SOCIAL STRUCTURE IN SEIR MODELS ---
#
# Repository link: https://github.com/mammadmammadov/Statistical-Programming-Practical-2
# Team Contribution: [Insert statement of proportional contribution by team members here]
#
# This R script implements a refined SEIR model to investigate the role of
# household structure and a heterogeneous contact network on epidemic dynamics.
# The code is designed to be efficient for large populations (n up to 10000)
# using vectorized operations where possible.
#
# States: S (Susceptible), E (Exposed), I (Infectious), R (Recovered/Immune)


# ==============================================================================
# 1. FUNCTION: create_households
# ==============================================================================
# PURPOSE:
# Creates a population structure by assigning n individuals to households,
# where household sizes are uniformly distributed between 1 and hmax.
#
# INPUTS:
# - n: Total population size (integer)
# - hmax: Maximum household size (integer, default 5)
#
# OUTPUT:
# - h: An n-vector of integers, where h[i] is the household ID of person i.

create_households <- function(n = 1000, hmax = 5) {
  
  # The strategy is to generate household sizes until their sum exceeds n, 
  # then adjust the last one to fit n exactly.
  
  # Step 1: Initialize a list of random sizes (1 to hmax) that sum to >= n.
  sizes <- c()
  while(sum(sizes) < n) {
    sizes <- c(sizes, sample(1:hmax, 1))
  }
  
  # Step 2: Adjust the last household size to ensure the total population is exactly n.
  sizes[length(sizes)] <- sizes[length(sizes)] - (sum(sizes) - n)
  
  # Step 3 (The "one line" construction): 
  # Use rep() to repeat each household ID (1 to length(sizes)) by its size,
  # and then sample() to randomly assign these IDs to the n individuals.
  h <- sample(rep(seq_along(sizes), times = sizes))
  
  return(h)
}


# ==============================================================================
# 2. FUNCTION: get.net
# ==============================================================================
# PURPOSE:
# Constructs the regular contact network based on individual sociability (beta)
# for people *not* in the same household. This function uses vectorization
# to avoid O(n^2) nested loops, which is required for efficiency.
#
# INPUTS:
# - beta: n-vector of sociability parameters (beta_i) for each person.
# - h: n-vector of household IDs.
# - nc: Average number of contacts per person (default 15).
#
# OUTPUT:
# - network: An adjacency list of length n, where network[[i]] is a vector
#            of indices of i's regular contacts.

get.net <- function(beta, h, nc = 15) {
  
  n <- length(beta)
  beta_mean <- mean(beta)
  
  # Denominator (d) for the probability expression: beta_bar^2 * (n - 1)
  d <- (beta_mean^2) * (n - 1)
  
  # Use an upper triangle index matrix for efficient generation of pairs (i, j) where i < j
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  i_indices <- idx[, 2] # Person i (column index)
  j_indices <- idx[, 1] # Person j (row index)
  
  # 1. Filter out household contacts (i and j are in different households)
  is_not_housemate <- h[i_indices] != h[j_indices]
  i_indices_valid <- i_indices[is_not_housemate]
  j_indices_valid <- j_indices[is_not_housemate]
  
  # 2. Calculate link probability for all valid pairs
  # P(link) = nc * beta_i * beta_j / d
  P_link <- (nc * beta[i_indices_valid] * beta[j_indices_valid]) / d
  
  # 3. Determine which links are created using vectorization
  # Use rbinom to generate 0 (no link) or 1 (link) for each pair
  is_linked <- rbinom(length(P_link), 1, P_link)
  
  # 4. Filter the established links
  linked_i <- i_indices_valid[is_linked == 1]
  linked_j <- j_indices_valid[is_linked == 1]
  
  # 5. Build the adjacency list (must be recorded in both directions)
  network <- vector("list", n)
  
  # Combine all outgoing links (i -> j) and incoming links (j -> i)
  all_i <- c(linked_i, linked_j)
  all_j <- c(linked_j, linked_i)
  
  # Use tapply to aggregate all contacts for each person ID (i.e., build the list efficiently)
  # This avoids the slow appending loop: network[[i]] <- c(network[[i]], j)
  network_data <- tapply(all_j, all_i, c, simplify = FALSE)
  
  # Merge the generated contacts back into the full 'network' list structure
  network[as.numeric(names(network_data))] <- network_data
  
  # Replace any NULL entries (people with no contacts) with an empty vector
  network[sapply(network, is.null)] <- list(integer(0))
  
  return(network)
}


# ==============================================================================
# 3. FUNCTION: nseir
# ==============================================================================
# PURPOSE:
# Simulates the SEIR model with household, network, and random mixing structure
# over nt days. Uses vectorization for efficient state transitions, particularly
# for the S -> E step.
#
# INPUTS:
# - beta: n-vector of sociability parameters.
# - h: n-vector of household IDs.
# - alink: Adjacency list (regular contact network) from get.net.
# - alpha: Vector of infection probabilities c(alpha_h, alpha_c, alpha_r).
# - delta: Daily probability of I -> R (default 0.2).
# - gamma: Daily probability of E -> I (default 0.4).
# - nc: Average contacts per person (default 15).
# - nt: Number of days to simulate (default 100).
# - pinf: Proportion of initial population to start in I state (default 0.005).
#
# OUTPUT:
# - A list with elements S, E, I, R (time series of counts) and t (time vector).

nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005){
  
  # Define states as integers for indexing
  S_STATE <- 1; E_STATE <- 2; I_STATE <- 3; R_STATE <- 4
  
  n <- length(beta)
  
  # Parse alpha probabilities
  alpha_h <- alpha[1] # Household transmission
  alpha_c <- alpha[2] # Network transmission
  alpha_r <- alpha[3] # Random mixing base rate
  
  # 1. INITIALIZATION
  
  # Status vector: start everyone as susceptible
  status <- rep(S_STATE, n)
  
  # Initial infected individuals
  initial_infected_count <- max(1, round(n * pinf))
  infected_indices <- sample(1:n, initial_infected_count)
  status[infected_indices] <- I_STATE
  
  # Initialize time series for counts
  S_t <- E_t <- I_t <- R_t <- numeric(nt + 1)
  
  # Record initial counts (Day 0)
  S_t[1] <- sum(status == S_STATE)
  E_t[1] <- sum(status == E_STATE)
  I_t[1] <- sum(status == I_STATE)
  R_t[1] <- sum(status == R_STATE)
  
  # Pre-calculate constants for random mixing (avoids repeated division)
  beta_mean <- mean(beta)
  # Denominator: d = beta_bar^2 * (n - 1)
  d <- (beta_mean^2) * (n - 1)
  # Pre-calculate the constant factor for random transmission probability: alpha_r * nc / d
  random_factor <- (alpha_r * nc) / d
  
  # Pre-calculate household membership to find housemates of an individual
  # This creates a list where house_members[[h[i]]] lists all members of i's household.
  house_members <- tapply(1:n, h, list, simplify = FALSE)
  
  # 2. SIMULATION LOOP
  
  for(day in 1:nt){
    
    # 2a. DETERMINE CURRENT STATE INDICES
    I_current <- which(status == I_STATE)
    E_current <- which(status == E_STATE)
    S_current <- which(status == S_STATE)
    
    # Check if the epidemic has died out
    if(length(I_current) == 0 && length(E_current) == 0){
      # If no exposed or infected, the simulation ends. Fill remaining days with current counts.
      S_t[day:nt + 1] <- S_t[day]
      E_t[day:nt + 1] <- E_t[day]
      I_t[day:nt + 1] <- 0
      R_t[day:nt + 1] <- R_t[day]
      break
    }
    
    # 2b. I -> R (RECOVERY)
    if(length(I_current) > 0){
      # Vectorized recovery: Check all current infected against delta probability
      recoveries <- I_current[runif(length(I_current)) < delta]
      if(length(recoveries) > 0){
        status[recoveries] <- R_STATE
      }
    }
    
    # 2c. E -> I (INCUBATION PERIOD ENDS)
    if(length(E_current) > 0){
      # Vectorized transition: Check all current exposed against gamma probability
      new_infections <- E_current[runif(length(E_current)) < gamma]
      if(length(new_infections) > 0){
        status[new_infections] <- I_STATE
      }
    }
    
    # 2d. S -> E (NEW EXPOSURES / INFECTION) - The most complex part
    
    # Initialize a logical vector to track which susceptibles get infected today.
    # We use this to ensure a person can only be infected once per day.
    newly_exposed <- logical(n) 
    
    # Only proceed if there are both susceptibles and infecteds
    if(length(S_current) > 0 && length(I_current) > 0){
      
      # --- HOUSEHOLD TRANSMISSION ---
      # Create a unique list of household IDs that contain at least one infected person.
      infected_household_ids <- unique(h[I_current])
      
      # Find all susceptible members in those infected households.
      # This vector contains susceptible indices that need to be checked for infection
      S_in_infected_households <- unlist(house_members[as.character(infected_household_ids)])
      S_to_check_hh <- intersect(S_current, S_in_infected_households)
      
      if(length(S_to_check_hh) > 0) {
        # Infection Probability: alpha_h (simple probability, irrespective of I count in HH)
        exposed_from_hh <- S_to_check_hh[runif(length(S_to_check_hh)) < alpha_h]
        newly_exposed[exposed_from_hh] <- TRUE
      }
      
      # --- REGULAR NETWORK TRANSMISSION ---
      
      # Create a list of all susceptible people who are contacts of *any* infected person.
      # Only consider those not already exposed from the household check.
      S_to_check_net <- S_current[!newly_exposed[S_current]]
      
      # Determine all contacts of all currently infected people (I_current)
      infected_contacts <- unlist(alink[I_current])
      
      # Intersect the contacts with the still-susceptible population
      S_in_infected_networks <- intersect(S_to_check_net, infected_contacts)
      
      if(length(S_in_infected_networks) > 0) {
        # Infection Probability: alpha_c (simple probability, irrespective of I count in net)
        exposed_from_net <- S_in_infected_networks[runif(length(S_in_infected_networks)) < alpha_c]
        newly_exposed[exposed_from_net] <- TRUE
      }
      
      # --- RANDOM MIXING TRANSMISSION ---
      
      # Only consider those not already exposed.
      S_to_check_random <- S_current[!newly_exposed[S_current]]
      
      if(length(S_to_check_random) > 0) {
        
        # Vectorized calculation: Sum of infection probabilities from all infected people (i)
        
        # Calculate P_infection from one infected person 'i' to ALL susceptibles 'j' in S_to_check_random
        # P_i_to_j = random_factor * beta_i * beta_j
        
        # Total probability of j being infected by ANY i is P_total = 1 - product(1 - P_i_to_j)
        
        # Calculate P_i_to_j matrix (I_current x S_to_check_random)
        # Random factor is already pre-calculated as (alpha_r * nc) / d
        
        # beta[I_current] is a vector of beta_i values for the infected
        # beta[S_to_check_random] is a vector of beta_j values for the susceptibles
        
        # The outer product calculates beta_i * beta_j for all pairs (i, j)
        beta_product_matrix <- outer(beta[I_current], beta[S_to_check_random])
        
        # P_i_to_j for all pairs (i, j)
        P_matrix <- random_factor * beta_product_matrix
        
        # Product of (1 - P_i_to_j) along the infected dimension (rows)
        # This gives the probability of NOT being infected by ANY of the infected people
        P_not_infected <- apply(1 - P_matrix, 2, prod)
        
        # Total probability of being infected is 1 - P_not_infected
        P_infected <- 1 - P_not_infected
        
        # Check against uniform random numbers to determine who gets infected
        exposed_from_random_idx <- runif(length(S_to_check_random)) < P_infected
        
        exposed_from_random <- S_to_check_random[exposed_from_random_idx]
        newly_exposed[exposed_from_random] <- TRUE
      }
      
      # Final S -> E transition for all individuals marked TRUE
      status[which(newly_exposed)] <- E_STATE
    }
    
    # 2e. RECORD COUNTS FOR THE DAY
    S_t[day + 1] <- sum(status == S_STATE)
    E_t[day + 1] <- sum(status == E_STATE)
    I_t[day + 1] <- sum(status == I_STATE)
    R_t[day + 1] <- sum(status == R_STATE)
  }
  
  # Return the time series
  return(list(S = S_t, E = E_t, I = I_t, R = R_t, t = 0:nt))
}


# ==============================================================================
# 4. FUNCTION: plot_simulation
# ==============================================================================
# PURPOSE:
# Plots the time series of the S, E, I, R states from the simulation results.
#
# INPUTS:
# - simulation_result: The list returned by the nseir function.
# - title: The title for the plot.
#
# OUTPUT:
# - A base R plot of the epidemic dynamics.

plot_simulation <- function(simulation_result, title) {
  
  # Determine the maximum population count for the y-axis limit
  y_max <- max(c(simulation_result$S, simulation_result$E, simulation_result$I, simulation_result$R))
  
  # Plot Susceptible (S) as the base layer
  plot(simulation_result$t, simulation_result$S, 
       type = "l", 
       col = "dodgerblue3",
       lwd = 2, 
       ylim = c(0, y_max), 
       xlab = "Days", 
       ylab = "Number of Individuals", 
       main = title)
  
  # Add Exposed (E), Infected (I), and Recovered (R) plots
  lines(simulation_result$t, simulation_result$E, col = "orange2", lwd = 2)
  lines(simulation_result$t, simulation_result$I, col = "red2", lwd = 2)
  lines(simulation_result$t, simulation_result$R, col = "green4", lwd = 2)
  
  # Add a legend
  legend("right", 
         legend = c("Susceptible", "Exposed", "Infected", "Recovered"), 
         col = c("dodgerblue3", "orange2", "red2", "green4"), 
         lty = 1,
         lwd = 2,
         cex = 0.8) # Smaller font for legend
}


# ==============================================================================
# 5. EXECUTION & SCENARIO COMPARISON
# ==============================================================================
# Code to run the model and answer the final part of the assignment.

# Set a consistent seed for reproducibility across all models
set.seed(42)

# --- Define General Parameters ---
population_size <- 1000 # Use 1000 for testing, but should scale to 10000+
max_household_size <- 5
average_contacts <- 15 # nc
simulation_days <- 100 # nt

# Fixed Disease Parameters
prob_recovery <- 0.2 # delta
prob_exposed_to_infected <- 0.4 # gamma
initial_infection_rate <- 0.005 # pinf

# --- Base Population Setup ---
# Sociability parameters, drawn from a vector of U(0,1) random variables (for Scenarios 1 & 2)
beta_base <- runif(population_size) 
# Household assignment (used for all scenarios)
h_base <- create_households(n = population_size, hmax = max_household_size) 

# --- Scenario 1: Full Model with Default Parameters (Heterogeneous Beta, All Mixing) ---
# alpha = (alpha_h=0.1, alpha_c=0.01, alpha_r=0.01)
alpha_model_1 <- c(0.1, 0.01, 0.01)
net_model_1 <- get.net(beta = beta_base, h = h_base, nc = average_contacts)

result_model_1 <- nseir(beta = beta_base, 
                        h = h_base, 
                        alink = net_model_1,
                        alpha = alpha_model_1, 
                        delta = prob_recovery,
                        gamma = prob_exposed_to_infected, 
                        nc = average_contacts, 
                        nt = simulation_days,
                        pinf = initial_infection_rate)


# --- Scenario 2: Only Random Mixing (Heterogeneous Beta, No Structure) ---
# Set alpha_h = 0, alpha_c = 0. Set alpha_r = 0.04 as specified.
alpha_model_2 <- c(0, 0, 0.04) 
# Note: The network structure (net) is irrelevant here since alpha_c = 0. 
# We use the same 'net' generation just to satisfy the function input, though it's unused in the S->E steps.
net_model_2 <- get.net(beta = beta_base, h = h_base, nc = average_contacts)

result_model_2 <- nseir(beta = beta_base, 
                        h = h_base, 
                        alink = net_model_2, # alink is unused when alpha_c=0
                        alpha = alpha_model_2, 
                        delta = prob_recovery,
                        gamma = prob_exposed_to_infected, 
                        nc = average_contacts, 
                        nt = simulation_days,
                        pinf = initial_infection_rate)


# --- Scenario 3: Full Model with Constant Beta (Homogeneous Beta, All Mixing) ---
# beta is replaced by its mean value for every element
beta_constant <- rep(mean(beta_base), population_size)
alpha_model_3 <- alpha_model_1

# A new network must be generated using the constant beta vector
net_model_3 <- get.net(beta = beta_constant, h = h_base, nc = average_contacts)

result_model_3 <- nseir(beta = beta_constant, 
                        h = h_base, 
                        alink = net_model_3,
                        alpha = alpha_model_3, 
                        delta = prob_recovery,
                        gamma = prob_exposed_to_infected, 
                        nc = average_contacts, 
                        nt = simulation_days,
                        pinf = initial_infection_rate)


# --- Scenario 4: Constant Beta and Only Random Mixing (Homogeneous Beta, No Structure) ---
# Combination of Scenarios 2 and 3 settings
alpha_model_4 <- alpha_model_2 # c(0, 0, 0.04)
beta_constant <- rep(mean(beta_base), population_size)

# Network is irrelevant since alpha_c=0
net_model_4 <- get.net(beta = beta_constant, h = h_base, nc = average_contacts)

result_model_4 <- nseir(beta = beta_constant, 
                        h = h_base, 
                        alink = net_model_4, # alink is unused when alpha_c=0
                        alpha = alpha_model_4, 
                        delta = prob_recovery,
                        gamma = prob_exposed_to_infected, 
                        nc = average_contacts, 
                        nt = simulation_days,
                        pinf = initial_infection_rate)


# --- Plotting All Models ---

# Setting up the plotting area to have 2 rows and 2 columns
par(mfrow=c(2,2), mar=c(4, 4, 2, 1) + 0.1) # Adjusted margins for better layout

# Plotting each model with appropriate titles
plot_simulation(result_model_1, "Scenario 1: Full Model (Heterogeneous Beta)")
plot_simulation(result_model_2, "Scenario 2: Only Random Mixing (Heterogeneous Beta)")
plot_simulation(result_model_3, "Scenario 3: Full Model (Constant Beta)")
plot_simulation(result_model_4, "Scenario 4: Only Random Mixing (Constant Beta)")


# --- Comment on Apparent Effect of Structure ---

# In Scenario 1 (Full Model, Heterogeneous Beta), the epidemic is generally smaller and 
# progresses slower than in Scenario 2 (Only Random Mixing). This suggests that the 
# household and network structure, combined with individual variability in sociability (beta), 
# effectively segregates the population. Infection is often confined within smaller, 
# tighter groups (households/networks) before spreading randomly, leading to a smaller overall epidemic size.

# When comparing Scenario 1 (Heterogeneous Beta, Full Model) with Scenario 3 (Constant Beta, Full Model),
# the heterogeneous beta leads to a slightly smaller peak and total infected count. This aligns with 
# the principle that higher variability often leads to smaller epidemics, as a few high-contact individuals 
# may get infected and recover quickly, effectively 'firewalling' the rest of the population.

# The most explosive epidemic is typically seen in Scenario 4 (Constant Beta, Random Mixing), 
# where everyone is essentially identical and mixes globally. The structural models (1 and 3) 
# delay the peak and reduce the maximum number of simultaneous infections compared to the
# purely random models (2 and 4), highlighting the importance of real-world social structure
# in mitigating spread and flattening the curve.
