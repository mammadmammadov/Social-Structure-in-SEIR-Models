# In this project, we will simulate the spread of an infectious disease in a population using an SEIR model.
# SEIR model is as follows: 
# S -> Susceptible
# E -> Exposed (infected but not yet infectious)
# I -> Infectious
# R -> Recovered (and immune)

# First, we will create a population with households and a contact network
# Then, we will simulate the spread of the disease over time such that 
# Finally, we will plot the results of the simulation

#*******************************************************************

# First of all, we assign people to households

# n is the total population size
# hmax is the maximum household size
# the function returns a vector of length n, where the ith element is the household ID of individual i
create_households <- function(n = 1000, hmax = 5) {


    household_sizes <- c()

    while(sum(household_sizes) < n) {
        household_sizes <- c(household_sizes, sample(1:hmax, 1))
    }

    # if we create too many spaces, we remove the excess from the last household
    excess <- sum(household_sizes) - n
    household_sizes[length(household_sizes)] <- household_sizes[length(household_sizes)] - excess

    # now we create a vector that assigns each person to a household
    # for example, if household_sizes = c(1, 2, 3), then we want
    # household = c(1, 2, 2, 3, 3, 3)
    # sample() is used to randomize the order of people in households
    h <- sample(rep(1:length(household_sizes), times = household_sizes))

    return(h)
}

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


# Now we will create the contact network for for the population

# This function will create a random contact network based on individual sociability parameters

# beta is a n-length vector of sociability parameters for each individual
# h is a vector assigning each individual to a household
# nc is the average number of contacts per individual 

# the function returns a list of length n (adjacency list), where the ith element is a vector of the contacts of individual i

get.net <- function(beta, h, nc) {

    n <- length(beta)
    
    # initializing an empty adjacency list to store the contact network
    network <- vector("list", n)

    beta_mean <- mean(beta)
    # β_ is the mean of the β_i where i = 1, ..., n
    
    # β_^2(n − 1)
    d <- (beta_mean^2) * (n - 1)

    # looping over all pairs of individuals to determine if they are in a contact and not in the same household
    for (i in 1:(n-1)){
        for(j in (i+1): n){
            if(h[i] != h[j]) {

                # ncβiβj / β_^2(n − 1)
                probable_link <- nc * beta[i] * beta[j] / d

                # print(probable_link)
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

# This function simulates the SEIR model on a given contact network

# Input:

# beta: vector of sociability parameters
# h: vector assigning each individual to a household
# alink: adjacency list representing the contact network
# alpha: vector of infection probabilities(within household, within regular network, random)
# delta: daily probability of infected -> recovered
# gamma: daily probability of exposed -> infected
# nc: average number of contacts
# nt: number of days to simulate
# pinf: proportion of the initial population to randomly start in the infected state

# this function returns a list containing the time series of S, E, I, R counts and the time vector


nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005){
    
    # defining states as integers for easier handling
    s_state <- 1
    e_state <- 2
    i_state <- 3
    r_state <- 4
    n <- length(beta)

    # alpha_h is the within household infection probability
    alpha_h <- alpha[1]

    # alpha_c is the within regular contact network infection probability
    alpha_c <- alpha[2]

    # alpha_r is the random infection probability
    alpha_r <- alpha[3]

    # this vector will hold the state of each individual, assuming everyone starts as susceptible
    status <- rep(s_state, n)

    # the number of initially infected individuals
    initial_infected <- max(1, round(n * pinf))

    # randomly selecting indices to be infected
    infected_indices <- sample(1:n, initial_infected)

    # setting the selected individuals to infected state
    status[infected_indices] <- i_state

    # initializing vectors to hold the counts of S, E, I, R over time
    S_t <- E_t <- I_t <- R_t <- numeric(nt+1)

    # initial counts for succeptible and infected at the start of the simulation
    S_t[1] <- n - initial_infected
    I_t[1] <- initial_infected

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

            # Randomly selecting individuals from I_current to transition to R state based on probability delta
            # runif(length(I_current)) < delta gives a logical vector indicating which individuals will recover
            recoveries <- I_current[runif(length(I_current)) < delta]
            if(length(recoveries) > 0){
                status[recoveries] <- r_state
            }
        }

        # EXPOSED TO INFECTED
        if(length(E_current) > 0){
            new_infections <- E_current[runif(length(E_current)) < gamma]
            if(length(new_infections) > 0){
                status[new_infections] <- i_state
            }
        }

        # SUSCEPTIBLE TO EXPOSED
        new_exposed <- logical(n) # Initialize a logical vector to track new exposures

        if(length(S_current) > 0 && length(I_current) > 0){
            # looping over each infected individual to determine whom they infect
            for(i in I_current){
                current_susceptibles <- S_current[!new_exposed[S_current]] # Only consider susceptibles not already exposed
                
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
                current_susceptibles <- S_current[!new_exposed[S_current]] # Update susceptibles not already exposed
                
                # no need to continue if no susceptibles left
                if(length(current_susceptibles) == 0){
                    break 
                }

                network_susceptibles <- intersect(current_susceptibles, alink[[i]])
                # infecting susceptibles in the contact network based on alpha_c
                if(length(network_susceptibles) > 0){
                    infected_in_network <- network_susceptibles[runif(length(network_susceptibles)) < alpha_c]
                    new_exposed[infected_in_network] <- TRUE
                }

                # RANDOM TRANSMISSION
                current_susceptibles <- S_current[!new_exposed[S_current]] # Update susceptibles not already exposed
                # no need to continue if no susceptibles left
                if(length(current_susceptibles) == 0){
                    break
            }
                # infecting random susceptibles based on alpha_r * nc * β_i * β_j / β_^2(n − 1)
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

# This function will plot the results of the SEIR simulation

# simulation_result is the output of the nseir function

# this function does not return anything, it just creates a plot
plot_simulation <- function(simulation_result) {
    plot(simulation_result$t, simulation_result$S, type = "l", col = "blue", ylim = c(0, max(simulation_result$S, simulation_result$E, simulation_result$I, simulation_result$R)), 
         xlab = "Days", ylab = "Number of Individuals", main = "SEIR Model Simulation")
    lines(simulation_result$t, simulation_result$E, col = "orange", lwd = 2)
    lines(simulation_result$t, simulation_result$I, col = "red", lwd = 2)
    lines(simulation_result$t, simulation_result$R, col = "green", lwd = 2)
    legend("right", legend = c("susceptible", "exposed", "infected", "recovered"), 
           col = c("blue", "orange", "red", "green"), lty = 1)
}

# ************************************************************


# Testing the plot_simulation function

set.seed(42)

population_size <- 1000

max_household_size <- 5

average_contacts <- 15

simulation_days <- 150

infection_probabilities <- c(0.1, 0.02, 0.05)

prob_recovery <- 0.1 # daily probability of infected -> recovered

prob_exposed_to_infected <- 0.2 # daily probability of exposed -> infected

initial_infection_rate <- 0.005 # proportion of the initial population to randomly start in the infected state

beta <- runif(population_size) # sociability parameters

h <- create_households(n = population_size, hmax = max_household_size) # assigning households to individuals

net <- get.net(beta = beta, h = h, nc = average_contacts) #creating the contact network

result <- nseir(beta = beta, h = h, alink = net,
 alpha = infection_probabilities, delta = prob_recovery,
 gamma = prob_exposed_to_infected, nc = average_contacts, nt = simulation_days,
    pinf = initial_infection_rate)

plot_simulation(result)









