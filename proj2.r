# SOCIAL STRUCTURE IN SEIR MODELS

# Repository link: https://github.com/mammadmammadov/Statistical-Programming-Practical-2

# In this project, we will simulate the spread of an infectious disease in a population using an SEIR model.

# The main states in the SEIR model are as follows: 
# S -> Susceptible (healthy but can get infected)
# E -> Exposed (infected but not yet infectious)
# I -> Infectious (can spread the disease)
# R -> Recovered (and immune)

# First, we will create a population with households. n people will be assigned to households of varying sizes, with a maximum household size of hmax.

# Second, we will create a contact network for the population based on individual sociability parameters (beta) and an average number of contacts (nc).

# Third, we will implement the SEIR model on the contact network. We will define parameters for infection probabilities (alpha), recovery probability (delta), and the probability of exposed individuals becoming infectious (gamma).
# When susceptible individuals come into contact with infectious individuals, they may become exposed based on the defined infection probabilities. This can happen through household contacts, regular contacts in the network, or random contacts.
# Then, we will simulate the spread of the disease over time to track the number of individuals in each state (S, E, I, R)

# Finally, we will run and plot the simulation under different scenarios, such as varying the infection probabilities and sociability parameters, to observe how these changes affect the spread of the disease.


start_time <- proc.time()


# *******************************************************************
# *******************************************************************


create_households <- function(n = 1000, hmax = 5) {
  
  # INPUTS:
  # n is the total population size
  # hmax is the maximum household size
  
  # OUTPUT:
  # a vector of length n assigning each individual to a household 
  # (for instance,if individual 1 is in household 3, then the first element of the vector will be 3)
  
  # PURPOSE:
  # this function creates households for a population of size n, with household sizes ranging from 1 to hmax
  
  
  household_sizes <- c()
  
  # we keep creating households until we have enough spaces for all n individuals
  while(sum(household_sizes) < n) {
    # we randomly sample a household size from 1 to hmax
    household_sizes <- c(household_sizes, sample(1:hmax, 1))
  }
  
  # if we create too many spaces, we remove the excess from the last household
  excess_space <- sum(household_sizes) - n
  household_sizes[length(household_sizes)] <- household_sizes[length(household_sizes)] - excess_space
  
  # now we create a vector that assigns each person to a household
  # we use the logic of rep() to repeat each household number according to its size
  # for example, if household_sizes = c(1, 3, 2), then we want to create a vector like c(1, 2, 2, 2, 3, 3) to accommodate  1st person in household 1, 2nd, 3rd, and 4th persons in household 2, and 5th and 6th persons in household 3
  # this also ensures that the total length of the vector is n, since sum(household_sizes) = n
  # finally, sample() is used to randomize the order of people in households
  h <- sample(rep(1:length(household_sizes), times = household_sizes))
  
  return(h)
}


# *******************************************************************
# *******************************************************************


# Now we will create the contact network for the population

get.net <- function(beta, h, nc = 15) {
  
  # INPUTS:
  # beta is a vector of sociability parameters
  # h is a vector assigning each individual to a household
  # nc is the average number of contacts
  
  # OUTPUT:
  # a list of length n, where each element is a vector of indices representing the contacts of that individual
  
  # PURPOSE:
  # this function creates a contact network for a population based on individual sociability parameters (beta) and an average number of contacts (nc)
  
  n <- length(beta)
  beta_mean <- mean(beta)
  
  # β_^2(n − 1) is assigned to a variable named "d" for easier use later on
  d <- (beta_mean^2) * (n - 1)
  
  # Use an upper triangle index matrix for efficient generation of pairs (i, j) where i < j
  # This avoids double checking
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  
  # Extract i and j indices from the index matrix
  i_indices <- idx[, 2] 
  j_indices <- idx[, 1]
  
  # make sure i and j are not in the same household if they are to be considered for a link
  is_not_housemate <- h[i_indices] != h[j_indices]
  i_indices <- i_indices[is_not_housemate]
  j_indices <- j_indices[is_not_housemate]
  
  # ncβiβj / β_^2(n − 1) gives the probability of a contact between individuals i and j
  probable_link  <- (nc * beta[i_indices] * beta[j_indices]) / d
  
  # we use binomial sampling to determine if a link is established
  # this gives a vector of 0s and 1s, where 1 indicates a link is established, and 0 indicates no link
  is_linked <- rbinom(length(probable_link), 1, probable_link)
  
  # extracting the pairs (i, j) where a link is established
  linked_i <- i_indices[is_linked == 1]
  linked_j <- j_indices[is_linked == 1]
  
  # initializing an empty list to hold the contact network
  network <- vector("list", n)
  
  # combining all outgoing links (i -> j) and incoming links (j -> i)
  links_from_i_to_j <- c(linked_i, linked_j)
  links_from_j_to_i <- c(linked_j, linked_i)
  
  # using tapply to create a list where each element corresponds to a person and contains their contacts
  network_data <- tapply(links_from_j_to_i, links_from_i_to_j, c, simplify = FALSE)
  
  # merging the generated contacts back into the full 'network' list structure
  network[as.numeric(names(network_data))] <- network_data
  
  # replacing any NULL entries (people with no contacts) with an empty vector
  network[sapply(network, is.null)] <- list(integer(0))
  
  return(network)
}


# ********************************************************************
# ********************************************************************


# Now, we will implement the SEIR model on a given contact network

nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005){
  
  # INPUTS:
  
  # beta is a vector of sociability parameters
  # h is a vector assigning each individual to a household
  # alink is an adjacency list representing the contact network
  # alpha is a vector of infection probabilities(within household, within regular network, random)
  # delta is a daily probability of infected -> recovered
  # gamma is a daily probability of exposed -> infected
  # nc is the average number of contacts
  # nt is the number of days to simulate
  # pinf is the proportion of the initial population to randomly start in the infected state
  
  # OUTPUT:
  # This function returns a list containing the time series of S, E, I, R counts and the time vector
  
  # PURPOSE:
  # This function simulates the spread of an infectious disease in a population using an SEIR model.
  # First, it initializes the state of each individual and sets up the parameters for the simulation.
  # Then, in each iteration, it updates the state of each individual based on the defined probabilities and contact network.
  # Finally, it returns the time series of the counts of individuals in each state (S, E, I, R) over the simulated period.
  
  # defining states as integers for easier handling
  s_state <- 1
  e_state <- 2
  i_state <- 3
  r_state <- 4
  
  # number of individuals in the population
  n <- length(beta)
  
  # alpha_h is the within household infection probability
  alpha_h <- alpha[1]
  
  # alpha_c is the within regular contact network infection probability
  alpha_c <- alpha[2]
  
  # alpha_r is the random infection probability
  alpha_r <- alpha[3]
  
  # this vector will hold the state of each individual, assuming everyone starts as susceptible
  status <- rep(s_state, n)
  
  # the number of initially infected individuals according to random sampling based on pinf
  initial_infected <- max(1, round(n * pinf))
  
  # randomly selecting indices to be infected
  infected_indices <- sample(1:n, initial_infected)
  
  # setting the selected individuals to infected state
  status[infected_indices] <- i_state
  
  # initializing vectors to hold the counts of S, E, I, R over time
  # we use nt+1 to include the initial state at time 0 (before the simulation starts)
  S_t <- E_t <- I_t <- R_t <- numeric(nt+1)
  
  # initial counts for susceptible and infected at the start of the simulation
  S_t[1] <- n - initial_infected
  I_t[1] <- initial_infected
  
  # calculating β_ 
  beta_mean <- mean(beta)
  
  # assigning it to a variable named "d" for easier use later on
  d <- (beta_mean^2) * (n - 1)
  
  # creating a list of household members for each individual
  household_members <- tapply(1:n ,h, list, simplify = FALSE)
  
  # simulating the disease spread over nt days
  for(day in 1:nt){
    
    # S_current, E_current, I_current are vectors of indices of individuals in each state
    S_current <- which(status == s_state)
    E_current <- which(status == e_state)
    I_current <- which(status == i_state)
    
    #  LOOKING AT STATE TRANSITIONS
    
    #  INFECTED TO RECOVERED
    if(length(I_current) > 0){
      
      # Randomly selecting individuals from I_current to transition to R state based on probability delta (δ)
      # runif(length(I_current)) < delta  gives a logical vector indicating which individuals will recover
      recoveries <- I_current[runif(length(I_current)) < delta]
      
      # updating the status of the recovered individuals from I to R
      if(length(recoveries) > 0){
        status[recoveries] <- r_state
      }
    }
    
    # EXPOSED TO INFECTED
    if(length(E_current) > 0){
      # Randomly selecting individuals from E_current to transition to I state based on probability gamma (γ)
      new_infections <- E_current[runif(length(E_current)) < gamma]
      
      # updating the status of the newly infected individuals from E to I
      if(length(new_infections) > 0){
        status[new_infections] <- i_state
      }
    }
    
    # SUSCEPTIBLE TO EXPOSED
    
    # initializing a logical vector to track new exposures
    new_exposed <- logical(n) 
    
    if(length(S_current) > 0 && length(I_current) > 0){
      
      
      # WITHIN HOUSEHOLD
      
      infected_household_indices <- unique(h[I_current])
      
      # identifying (potential) susceptibles in households with at least one infected individual
      susceptible_in_infected_households <- unlist(household_members[as.character(infected_household_indices)])
      
      # susceptibles to check for possible exposure
      susceptibles_to_check_household <- intersect(S_current, susceptible_in_infected_households)
      
      if(length(susceptibles_to_check_household) > 0){
        exposed_from_household <- susceptibles_to_check_household[runif(length(susceptibles_to_check_household)) < alpha_h]
        new_exposed[exposed_from_household] <- TRUE
      }
      
      # REGULAR CONTACT NETWORK
      
      # susceptibles to check for possible exposure in the contact network excluding those already exposed from household
      susceptibles_to_check_network <- S_current[!new_exposed[S_current]]
      
      infected_contacts <- unlist(alink[I_current])
      
      susceptibles_in_infected_contacts <- intersect(susceptibles_to_check_network, infected_contacts)
      
      if(length(susceptibles_in_infected_contacts) > 0){
        exposed_from_contacts <- susceptibles_in_infected_contacts[runif(length(susceptibles_in_infected_contacts)) < alpha_c]
        new_exposed[exposed_from_contacts] <- TRUE
      }
      
      # RANDOM MIXING
      
      susceptible_to_check_random <- S_current[!new_exposed[S_current]]
      
      if(length(susceptible_to_check_random) > 0){
        
        # we calculate beta_i * beta_j for all pairs of infected individuals i and susceptibles j
        beta_product_matrix <- outer(beta[I_current], beta[susceptible_to_check_random])
        
        # calculating the probability of random infection for each pair (i, j)
        p_matrix <- alpha_r * nc * beta_product_matrix / d
        
        # we need to calculate the probability that each susceptible j gets infected by at least one infected individual i
        # to find this, we first calculate the probability that j does not get infected by any i
        # then, we subtract this from 1 to get the probability that j gets infected by at least one i
        
        p_not_infected <- apply(1 - p_matrix, 2, prod)
        p_infected <- 1 - p_not_infected
        
        # now we randomly determine which susceptibles get infected based on the calculated probabilities
        exposed_from_random <- susceptible_to_check_random[runif(length(susceptible_to_check_random)) < p_infected]
        new_exposed[exposed_from_random] <- TRUE
      }
      
      status[which(new_exposed)] <- e_state
    }
    # updating counts for the day
    S_t[day + 1] <- sum(status == s_state)
    E_t[day + 1] <- sum(status == e_state)
    I_t[day + 1] <- sum(status == i_state)
    R_t[day + 1] <- sum(status == r_state)
  }
  
  # returning the time series of S, E, I, R counts, along with the time vector
  return(list(S = S_t, E = E_t, I = I_t, R = R_t, t = 0:nt))
  
}


# ************************************************************
# ************************************************************


# now, we will write a function to plot the results of the SEIR simulation

plot_simulation <- function(simulation_result, title) {
  
  # INPUTS:
  # simulation_result is the output of the nseir function
  # title is the title of the plot
  
  # OUTPUT:
  # a plot showing the number of individuals in each state (S, E, I, R) over time
  
  # PURPOSE:
  # this function takes the result of the SEIR simulation from the nseir() function and plots the number of individuals in each state (S, E, I, R) over time.
  
  
  plot(simulation_result$t, simulation_result$S, 
       type = "l", 
       col = "dodgerblue3",
       lwd = 2, 
       ylim = c(0, max(simulation_result$S, simulation_result$E, simulation_result$I, simulation_result$R)), xlab = "Days", ylab = "Number of Individuals", main = title)
  
  # adding Exposed, Infected, and Recovered plots
  lines(simulation_result$t, simulation_result$E, col = "orange2", lwd = 2)
  lines(simulation_result$t, simulation_result$I, col = "red2", lwd = 2)
  lines(simulation_result$t, simulation_result$R, col = "green4", lwd = 2)
  
  # adding a legend
  legend("right", 
         legend = c("Susceptible", "Exposed", "Infected", "Recovered"), 
         col = c("dodgerblue3", "orange2", "red2", "green4"), 
         lty = 1,
         lwd = 2)
}


# ************************************************************
# ************************************************************


# now we will run the full simulation and plot the results in four different scenarios

set.seed(42)

population_size <- 10000

max_household_size <- 5

average_contacts <- 15

simulation_days <- 100

infection_probabilities <- c(0.1, 0.01, 0.01)

# daily probability of infected -> recovered
prob_recovery <- 0.2 

# daily probability of exposed -> infected
prob_exposed_to_infected <- 0.4 

# proportion of the initial population to randomly start in the infected state
initial_infection_rate <- 0.005 

# sociability parameters taken from a vector of U(0,1) random variables
beta <- runif(population_size) 

# assigning households to individuals
h <- create_households(n = population_size, hmax = max_household_size) 

# creating the contact network (it is the same for the first and second models)
net <- get.net(beta = beta, h = h, nc = average_contacts)

# model 1 - default parameters
result_model_1 <- nseir(beta = beta, 
                        h = h, 
                        alink = net,
                        alpha = infection_probabilities, 
                        delta = prob_recovery,
                        gamma = prob_exposed_to_infected, 
                        nc = average_contacts, 
                        nt = simulation_days,
                        pinf = initial_infection_rate)

# model 2 - only random infection
infection_probabilities_modified <- c(0, 0, 0.04)
net <- get.net(beta = beta, h = h, nc = average_contacts)
result_model_2 <- nseir(beta = beta, 
                        h = h, 
                        alink = net,
                        alpha = infection_probabilities_modified, 
                        delta = prob_recovery,
                        gamma = prob_exposed_to_infected, 
                        nc = average_contacts, 
                        nt = simulation_days,
                        pinf = initial_infection_rate)

# model 3 - constant beta
beta_modified = rep(mean(beta), population_size)
net <- get.net(beta = beta_modified, h = h, nc = average_contacts)
result_model_3 <- nseir(beta = beta_modified, 
                        h = h, 
                        alink = net,
                        alpha = infection_probabilities, 
                        delta = prob_recovery,
                        gamma = prob_exposed_to_infected, 
                        nc = average_contacts, 
                        nt = simulation_days,
                        pinf = initial_infection_rate)

# model 4 - constant beta and random infection
infection_probabilities_modified <- c(0, 0, 0.04)
beta_modified = rep(mean(beta), population_size)
net <- get.net(beta = beta_modified, h = h, nc = average_contacts)
result_model_4 <- nseir(beta = beta_modified, 
                        h = h, 
                        alink = net,
                        alpha = infection_probabilities_modified, 
                        delta = prob_recovery,
                        gamma = prob_exposed_to_infected, 
                        nc = average_contacts, 
                        nt = simulation_days,
                        pinf = initial_infection_rate)

# plotting all models together in a 2x2 grid

# setting up the plotting area to have 2 rows and 2 columns
par(mfrow=c(2,2))

# plotting each model with appropriate titles
plot_simulation(result_model_1, "Default parameters")
plot_simulation(result_model_2, "Only random mixing")
plot_simulation(result_model_3, "Constant Beta")
plot_simulation(result_model_4, "Only random mixing and constant Beta")

# printing the time taken for the entire simulation
cat(sprintf("\nTotal execution time: %.2f seconds\n", (proc.time() - start_time)[3]))

# Analysis of the four SEIR plots:
# 1. Default parameters: All infection routes (household, regular network, random mixing) active.
#    > Produces a moderate infection peak and takes longer to reach equilibrium state.
# 2. Only random mixing: Infection spreads only throught random mixing.
#    > Leads to a slightly higher and earlier infection peak, with more people infected overall,
#      but the system also reaches equilibrium faster due to the rapid progression of the outbreak.
# 3. Constant Beta: All infection routes active, but each person has the same sociability (beta).
#    > Similar to the default model, showing a lower infection peak and slower return to equilibrium.
# 4. Only random mixing + Constant Beta: Infection spreads only throught random mixing, and everyone has identical sociability.
#    > Produces higher and earlier infection peak, followed by a rapid decline as the outbreak resolves quickly.
# Overall:
#   - Random mixing only increases infection intensity and accelerates epidemic resolution.
#   - Having interactions within the same groups (like households or regular networks) slow the spread and delay equilibrium, 
#     but fewer people get infected overall.
