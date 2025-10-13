# In this project, we will simulate the spread of an infectious disease in a population using an SEIR model.

# The main states in the SEIR model is as follows: 
# S -> Susceptible
# E -> Exposed (infected but not yet infectious)
# I -> Infectious
# R -> Recovered (and immune)

# The population will be structured into households and a contact network, and we will consider different modes of transmission.

# First, we will create a population with households and a contact network
# Then, we will simulate the spread of the disease over time such to track the number of individuals in each state (S, E, I, R) 
# Finally, we will plot the results of the simulation

#*******************************************************************

# INPUTS:
# n is the total population size
# hmax is the maximum household size

# OUTPUT:
# a vector of length n assigning each individual to a household 
# (for instance,if individual 1 is in household 3, then the first element of the vector will be 3)

# PURPOSE:
# this function creates households for a population of size n, with household sizes ranging from 1 to hmax

create_households <- function(n = 1000, hmax = 5) {
  
  
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


# Testing the create_households function

# population_size <- 1000
# max_household_size <- 5
# h <- create_households(n = population_size, hmax = max_household_size)
# print(head(h, 20))
# cat(sprintf("Total number of people: %d\n", length(h))) 
# cat(sprintf("Total number of households: %d\n", length(unique(h))))

# person_id <- 969
# household_id_of_person_1 <- h[1]
# members <- which(h == household_id_of_person_1)
# cat(sprintf("Person %d is in household %d\n", person_id, household_id_of_person_1))
# print(members)

# # this gives us the size of each household
# household_size_counts <- table(h)

# # Now, let's create a frequency table of those sizes.
# size_distribution <- table(household_size_counts)

# cat("\nHousehold Size Distribution\n")
# print(size_distribution)

# *******************************************************************
# *******************************************************************


# Now we will create the contact network for for the population

# INPUTS:
# beta is a vector of sociability parameters
# h is a vector assigning each individual to a household
# nc is the average number of contacts per 

# OUTPUT:
# a list of length n (adjacency list), where the ith element is a vector of the contacts of individual i

# PURPOSE:
# This function will create a random contact network based on individual sociability parameters

get.net <- function(beta, h, nc) {
  
  n <- length(beta)
  
  # initializing an empty adjacency list to store the contact network
  network <- vector("list", n)
  
  # β_ is the mean of the β_i where i = 1, ..., n
  beta_mean <- mean(beta)
  
  
  # β_^2(n − 1) is assigned to a variable named "d" for easier use later on
  d <- (beta_mean^2) * (n - 1)
  
  # looping over all pairs of individuals to determine if they are in a contact and not in the same household
  for (i in 1:(n-1)){
    for(j in (i+1): n){
      # checking if individuals i and j are in different households
      if(h[i] != h[j]) {
        
        # ncβiβj / β_^2(n − 1) gives the probability of a contact between individuals i and j
        probable_link <- nc * beta[i] * beta[j] / d
        
        # probable_link gives very small values (0.0X), so we use runif(1) < probable_link instead of runif(1) >= probable_link
        if(runif(1) < probable_link) {
          
          # we make sure to assign the contact in both directions
          network[[i]] <- c(network[[i]], j)
          network[[j]] <- c(network[[j]], i)
        }
      }
    }
  }
  return(network)
}

# ************************************************************
# ************************************************************

# Testing the get.net function

# set.seed(42)

# population_size <- 1000
# max_household_size <- 5 
# average_contacts <- 15
# h <- create_households(n = population_size, hmax = max_household_size)

# beta <- runif(population_size)

# network <- get.net(beta = beta, h = h, nc = average_contacts)

# print(network)


# contacts_of_1 <- network[[1]]
# household_id <- h[1]
# housemates <- which(h == 1)
# housemates <- housemates[housemates != 1] 
# has_contact_with_housemate <- any(contacts_of_1 %in% housemates)
# print(has_contact_with_housemate)

# ********************************************************************

# Now, we will implement the SEIR model on a given contact network

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


nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005){
  
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
      # looping over each infected individual to determine whom they infect
      for(i in I_current){
        
        # updating susceptibles not already exposed
        current_susceptibles <- S_current[!new_exposed[S_current]] 
        # no need to continue if no susceptibles left
        if(length(current_susceptibles) == 0){
          break
        }
        
        # WITHIN HOUSEHOLD TRANSMISSION
        
        # finding household members of infected individual i
        household_members <- which(h == h[i])
        
        # finding susceptibles in the same household
        household_susceptibles <- intersect(current_susceptibles, household_members)
        
        # infecting susceptibles in the household based on alpha_h
        if(length(household_susceptibles) > 0){
          infected_in_household <- household_susceptibles[runif(length(household_susceptibles)) < alpha_h]
          new_exposed[infected_in_household] <- TRUE
        }
        
        # WITHIN REGULAR CONTACT NETWORK TRANSMISSION
        
        # updating susceptibles not already exposed
        current_susceptibles <- S_current[!new_exposed[S_current]] 
        
        # no need to continue if no susceptibles left
        if(length(current_susceptibles) == 0){
          break 
        }
        
        # finding susceptibles in the contact network of infected individual i
        network_susceptibles <- intersect(current_susceptibles, alink[[i]])
        
        # infecting susceptibles in the contact network based on alpha_c
        if(length(network_susceptibles) > 0){
          infected_in_network <- network_susceptibles[runif(length(network_susceptibles)) < alpha_c]
          new_exposed[infected_in_network] <- TRUE
        }
        
        # RANDOM TRANSMISSION
        
        # updating susceptibles not already exposed
        current_susceptibles <- S_current[!new_exposed[S_current]] 
        
        # no need to continue if no susceptibles left
        if(length(current_susceptibles) == 0){
          break
        }
        # infecting random susceptibles based on probability alpha_r * nc * β_i * β_j / β_^2(n − 1)
        random_susceptibles <- current_susceptibles[runif(length(current_susceptibles)) < (alpha_r * nc * beta[i] * beta[current_susceptibles] / d)]
        new_exposed[random_susceptibles] <- TRUE
      }
    }
    
    # after checking all infected individuals, we apply the S -> E transitions
    status[which(new_exposed)] <- e_state
    
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
# Testing the nseir function

# set.seed(42)
# population_size <- 1000
# max_household_size <- 5
# average_contacts <- 15
# simulation_days <- 150

# infection_probabilities <- c(0.1, 0.01, 0.01)
# prob_recovery <- 0.2 # daily probability of infected -> recovered
# prob_exposed_to_infected <- 0.4 # daily probability of exposed -> infected
# initial_infection_rate <- 0.05 # proportion of the initial population to randomly start in the infected state

# beta <- runif(population_size) # sociability parameters
# h <- create_households(n = population_size, hmax = max_household_size) # assigning households to individuals
# net <- get.net(beta = beta, h = h, nc = average_contacts) # creating the contact network

# result <- nseir(beta = beta, h = h, alink = net, alpha = infection_probabilities, delta = prob_recovery,
#  gamma = prob_exposed_to_infected, nc = average_contacts, nt = simulation_days, pinf = initial_infection_rate)

# legend("right", legend = c("susceptible", "exposed", "infected", "recovered"), 
#        col = c("blue", "orange", "red", "green"), lty = 1)


# ************************************************************
# ************************************************************

# now, we will write a function to plot the results of the SEIR simulation

# INPUTS:
# simulation_result is the output of the nseir function
# title is the title of the plot

# OUTPUT:
# a plot showing the number of individuals in each state (S, E, I, R) over time

# PURPOSE:
# this function takes the result of the SEIR simulation from the nseir() function and plots the number of individuals in each state (S, E, I, R) over time.

plot_simulation <- function(simulation_result, title) {
  
  
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

population_size <- 1000

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

# sociability parameters, drawn from a vector of U(0,1) random variables
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