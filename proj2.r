# First of all, we assign people to households
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
    
    network <- vector("list", n)
    beta_mean <- mean(beta)
    # β_ is the mean of the β_i where i = 1, ..., n
    # β_^2(n − 1)

    d <- (beta_mean^2) * (n - 1)
    for (i in 1:(n-1)){
        for(j in (i+1): n){
            if(h[i] != h[j]) {
                # ncβiβj / β_^2(n − 1)
                probable_link <- nc * beta[i] * beta[j] / d
                # print(probable_link)
                if(runif(1) < probable_link) {
                    network[[i]] <- c(network[[i]], j)
                    network[[j]] <- c(network[[j]], i)
                }
            }
        }
    }
    return(network)
}

# Testing the get.net function

set.seed(42)

population_size <- 1000
max_household_size <- 5 
average_contacts <- 15
h <- create_households(n = population_size, hmax = max_household_size)

beta <- runif(population_size)

network <- get.net(beta = beta, h = h, nc = average_contacts)

print(network[977])


  
contacts_of_1 <- network[[1]]
household_id <- h[1]
housemates <- which(h == 1)
housemates <- housemates[housemates != 1] 
has_contact_with_housemate <- any(contacts_of_1 %in% housemates)
print(has_contact_with_housemate)

# ********************************************************************

# nseir(beta, h, alink, alpha = c(.1, .1, .1), gamma = c(.1, .1), days = 100, initial_infected = 10)
