# Version 1 of a complex COVID-19 ABM
# The key object generated with model outputs is called "report"

library(dplyr)
# library(R0)

run_in_par <- 1 # Set this variable to 1 to run in parallel, otherwise to 0 - requires the doParallel package to work
num_cores <- 6 # Specify how many cores to use
generate_html_output <- 1 # Set this to 1 if you would like to automatically generate an HTML file with key stats and graphs - it requires the rmarkdown package and the ac19.rmd file to be in your working directory

# How many times do you want to run the simulation? By default the first third will run with CTI and TaT =2, second third with CTI and TaT =2 and third third with CTI and TaT =8
num_runs <- 300

# Important parameters can be set here for ease of use. Most other parameters are set inside the run_model function
num_agents <- 10000 # Number of agents in the population
num_time <- 240 # Default time steps is days
initial_infections <- 20

start_time <- Sys.time()

# Initialise agents data frame
make_agents <- function(num_agents) {
  agents <- data.frame(
    ID = (c(1:num_agents)),
    gender = sample(0:1, num_agents, replace = T),
    exposed = 0,
    infectious = 0,
    symptomatic = 0,
    isolated = 0,
    tested_positive = 0,
    times_tested = 0,
    days_infectious = 0, # Will later contain days since became infectious
    days_exposed = 0, # Will later contain days since became infected
    time_for_infectiousness = 0,
    time_for_symptoms = 0,
    time_for_test = -20,
    time_for_stay_home = 0,
    time_for_icu = 0,
    time_for_cure = 0,
    time_for_end_isolation = 0,
    time_for_end_stay_home = 0,
    time_for_death = 0,
    will_need_icu = 0,
    will_die = 0,
    will_stay_home = 0,
    will_be_symptomatic = 0,
    stay_home = 0,
    age = sample(0:90, num_agents, replace = T), #  assigns agents an age from 0 to 90
    household = 0,
    block = 0,
    class = 0,
    workplace = 0,
    recovered = 0,
    taxi = 0,
    symptomatic = 0,
    icu = 0,
    dead = 0,
    dead_icu_full = 0,
    times_quarantined = 0,
    where_infected = 0, # 1=home, 2=block, 3=class, 4=work, 5=taxi
    birth_week = sample(1:52, num_agents, replace = T)
  )
  return(agents)
}

# Infect agents at start of simulation
start_infect <- function(agents_model, num_agents, initial_infections) {
  for (y in 1:initial_infections) {
    x <- sample(1:num_agents, size = 1, replace = T)
    agents_model$exposed[x] <- 1
  }
  return(agents_model)
}

# Put agents into households
make_households <- function(agents_model, household_size, num_agents) {
  x <- 0
  house_num <- 1
  while (x < num_agents) {
    size_next_house <- sample(2:((household_size * 2) + 2), size=1,  replace = T)
    if ((num_agents - x) <= ((household_size * 2) + 2)) {
      size_next_house <- (num_agents - x)
    }
    for (j in x:(x + size_next_house)) {
      agents_model$household[j] <- house_num
    }
    house_num <- house_num + 1
    x <- (x + size_next_house)
  }
  return(agents_model)
}

# Put households into blocks
make_blocks <- function(agents_model, households_per_block, num_agents) {
  for (x in 1:num_agents) {
    agents_model$block[x] <- round(agents_model$household[x] / households_per_block)
  }
  return(agents_model)
}

# Put agents into school classes
make_classes <- function(agents_model, class_size) {
  for (i in 1:12) { # For each of the 12 grades
    num_kids <- sum(agents_model$age == (i + 5))
    num_classes <- round(num_kids / class_size)
    if (num_kids < class_size) {
      num_classes <- 1
    }
    modified_class_size <- round(num_kids / num_classes)
    kids_left <- num_kids
    for (j in 1:num_classes) {
      if (j == num_classes) { # The size of the last or only class = kids_left
        m <- dplyr::filter(agents_model, age == (i + 5), class == 0)
        for (k in 1:kids_left) {
          agents_model$class[m$ID[k]] <- j
        }
      }
      if (j < num_classes) { # The size of all classes that are not the last or only class = modified_class_size
        m <- filter(agents_model, age == (i + 5), class == 0)
        for (k in 1:modified_class_size) {
          agents_model$class[m$ID[k]] <- j
        }
      }
      kids_left <- (kids_left - modified_class_size)
    }
  }
  return(agents_model)
}

# Put agents into workplaces
make_workplaces <- function(agents_model, employment_rate, num_workplaces) {
  temp1 <- dplyr::filter(agents_model, age < 61 && age > 18)
  work_age_agents <- length(temp1$ID)
  num_workers <- round(work_age_agents * employment_rate)
  temp2 <- dplyr::sample_n(temp1, num_workers)
  workplace_size <- round(num_workers / num_workplaces)
  workplace_num <- 1
  x <- 0
  while (x < num_workers) {
    size_next_workplace <- sample(2:((workplace_size * 2) + 2), size =1,  replace = T)
    if ((num_workers - x) <= ((workplace_size * 2) + 2)) {
      size_next_workplace <- (num_workers - x)
    }
    for (j in x:(x + size_next_workplace)) {
      agents_model$workplace[temp2$ID[j]] <- workplace_num
    }
    workplace_num <- workplace_num + 1
    x <- (x + size_next_workplace)
  }
  return(agents_model)
}

# Determines which working agents will regularly take taxis
make_taxis <- function(agents_model, regular_taxi_takers) {
  working_agents <- dplyr::filter(agents_model, workplace > 0)
  taxi_regulars <- dplyr::sample_frac(working_agents, size = regular_taxi_takers, replace = T)
  for (x in 1:length(taxi_regulars$ID)) {
    agents_model$taxi[taxi_regulars$ID[x]] <- 1
  }
  return(agents_model)
}

# Reads in days since infection and returns a test result accounting for the false negative percentage by day since infection
get_tested <- function(days_infected) {
  test_result <- 0
  k <- sample(1:100, size = 1, replace = T)
  # Below figures are based on this article in AIM https://www.acpjournals.org/doi/10.7326/M20-1495
  if (days_infected == 1) {
    if (k < 1) test_result <- 1
  } else if (days_infected == 2) {
    if (k < 1) test_result <- 1
  } else if (days_infected == 3) {
    if (k < 6) test_result <- 1
  } else if (days_infected == 4) {
    if (k < 33) test_result <- 1
  } else if (days_infected == 5) {
    if (k < 63) test_result <- 1
  } else if (days_infected == 6) {
    if (k < 76) test_result <- 1
  } else if (days_infected == 7) {
    if (k < 81) test_result <- 1
  } else if (days_infected == 8) {
    if (k < 81) test_result <- 1
  } else if (days_infected == 9) {
    if (k < 80) test_result <- 1
  } else if (days_infected == 10) {
    if (k < 78) test_result <- 1
  } else if (days_infected == 11) {
    if (k < 75) test_result <- 1
  } else if (days_infected == 12) {
    if (k < 71) test_result <- 1
  } else if (days_infected == 13) {
    if (k < 67) test_result <- 1
  } else if (days_infected == 14) {
    if (k < 63) test_result <- 1
  } else if (days_infected == 15) {
    if (k < 59) test_result <- 1
  } else if (days_infected == 16) {
    if (k < 55) test_result <- 1
  } else if (days_infected == 17) {
    if (k < 51) test_result <- 1
  } else if (days_infected == 18) {
    if (k < 47) test_result <- 1
  } else if (days_infected == 19) {
    if (k < 43) test_result <- 1
  } else if (days_infected == 20) {
    if (k < 39) test_result <- 1
  } else if (days_infected == 21) {
    if (k < 35) test_result <- 1
  } else if (days_infected == 22) {
    if (k < 31) test_result <- 1
  } else if (days_infected == 23) {
    if (k < 27) test_result <- 1
  } else if (days_infected == 24) {
    if (k < 19) test_result <- 1
  } else if (days_infected == 25) {
    if (k < 15) test_result <- 1
  } else if (days_infected == 26) {
    if (k < 11) test_result <- 1
  } else if (days_infected == 27) {
    if (k < 7) test_result <- 1
  } else if (days_infected == 28) {
    if (k < 3) test_result <- 1
  }
  return(test_result)
}

# This function first calibrates risk by age and then adjusts death risk according to death_risk variable and then returns an array of risk values by age
death_adjust <- function(death_rate) {
  # dr10s = death risk in agents in their 10s, etc. Current values are based on this Lancet study https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext
  dr0s <- 0
  dr10s <- 0.00695 / 100
  dr20s <- 0.0309 / 100
  dr30s <- 0.0844 / 100
  dr40s <- 0.161 / 100
  dr50s <- 0.595 / 100
  dr60s <- 1.93 / 100
  dr70s <- 4.28 / 100
  dr80s <- 7.8 / 100
  # The following adjusts the relative risk for different ages up or down based on the value given to the death-rate variable
  total_risk <- (dr0s + dr10s + dr20s + dr30s + dr40s + dr50s + dr60s + dr70s + dr80s) / 9
  dr0s <- dr0s * (death_rate / total_risk)
  dr10s <- dr10s * (death_rate / total_risk)
  if (dr10s > 1) dr10s <- 1
  dr20s <- dr20s * (death_rate / total_risk)
  if (dr20s > 1) dr20s <- 1
  dr30s <- dr30s * (death_rate / total_risk)
  if (dr30s > 1) dr30s <- 1
  dr40s <- dr40s * (death_rate / total_risk)
  if (dr40s > 1) dr40s <- 1
  dr50s <- dr50s * (death_rate / total_risk)
  if (dr50s > 1) dr50s <- 1
  dr60s <- dr60s * (death_rate / total_risk)
  if (dr60s > 1) dr60s <- 1
  dr70s <- dr70s * (death_rate / total_risk)
  if (dr70s > 1) dr70s <- 1
  dr80s <- dr80s * (death_rate / total_risk)
  if (dr80s > 1) dr80s <- 1
  death_risk_array <- c(dr10s, dr20s, dr30s, dr40s, dr50s, dr60s, dr70s, dr80s)
  return(death_risk_array)
}

# This function reads in an ICU agent's age and returns whether that agent will die based on the age-related death risk defined in this function
death_test <- function(age, death_risk_array) {
  will_die <- 0
  # Apply the age-based death risk to the agent
  if (age > 9 & age < 20) {
    will_die <- sample(0:1, size = 1, prob = c(1 - death_risk_array[1], death_risk_array[1]), replace = T)
  } else if (age > 19 & age < 30) {
    will_die <- sample(0:1, size = 1, prob = c(1 - death_risk_array[2], death_risk_array[2]), replace = T)
  } else if (age > 29 & age < 40) {
    will_die <- sample(0:1, size = 1, prob = c(1 - death_risk_array[3], death_risk_array[3]), replace = T)
  } else if (age > 39 & age < 50) {
    will_die <- sample(0:1, size = 1, prob = c(1 - death_risk_array[4], death_risk_array[4]), replace = T)
  } else if (age > 49 & age < 60) {
    will_die <- sample(0:1, size = 1, prob = c(1 - death_risk_array[5], death_risk_array[5]), replace = T)
  } else if (age > 59 & age < 70) {
    will_die <- sample(0:1, size = 1, prob = c(1 - death_risk_array[6], death_risk_array[6]), replace = T)
  } else if (age > 69 & age < 80) {
    will_die <- sample(0:1, size = 1, prob = c(1 - death_risk_array[7], death_risk_array[7]), replace = T)
  } else if (age > 79) {
    will_die <- sample(0:1, size = 1, prob = c(1 - death_risk_array[8], death_risk_array[8]), replace = T)
  }
  return(will_die)
}

# Calculates risk for each household and then exposes agents to that risk
house_transmit <- function(agents_model, house_risk, num_houses) {
  for (j in 1:num_houses) {
    temp <- dplyr::filter(agents_model, household == j & isolated == 0 & dead == 0 & icu == 0)
    num_infected <- sum(temp$infectious)
    if (num_infected > 0) {
      this_house_risk <- (num_infected * house_risk)
      num_in_house <- length(temp$ID)
      for (l in 1:num_in_house) {
        k <- sample(1:10000, size = 1, replace = T)
        if (k < this_house_risk & temp$recovered[l] == 0 & temp$days_exposed[l] == 0) {
          agents_model$exposed[temp$ID[l]] <- 1
          agents_model$where_infected[temp$ID[l]] <- 1
        }
      }
    }
  }
  return(agents_model)
}

# Calculates risk for a block ( or collection of households)  and then exposes agents to that risk
block_transmit <- function(agents_model, block_risk) {
  for (x in 1:max(agents_model$block)) {
    current_block <- dplyr::filter(agents_model, block == x & isolated == 0 & stay_home == 0 & dead == 0 & icu ==0)
    if (length(current_block$ID) > 0) {
      current_block_risk <- sqrt(block_risk * sum(current_block$infectious))
      for (y in 1:length(current_block$ID)) {
        k <- sample(1:10000, size = 1, replace = T)
        if (k < current_block_risk & current_block$infectious[y] == 0) {
          agents_model$exposed[current_block$ID[y]] <- 1
          agents_model$where_infected[current_block$ID[y]] <- 2
        }
      }
    }
  }
  return(agents_model)
}

# Calculates risk for each class and then exposes agents to that risk
class_transmit <- function(agents_model, class_risk) {
  for (k in 1:12) {
    temp1 <- dplyr::filter(agents_model, age == (k + 5 & stay_home == 0 & isolated == 0 & icu == 0 & dead == 0))
    num_classes <- max(temp1$class)
    for (j in 1:num_classes) {
      temp2 <- dplyr::filter(temp1, class == j)
      num_infected <- sum(temp2$infectious)
      if (num_infected > 0) {
        this_class_risk <- (num_infected * class_risk)
        num_in_class <- length(temp2$ID)
        for (l in 1:num_in_class) {
          k <- sample(1:10000, size = 1, replace = T)
          if (k < this_class_risk & temp2$recovered[l] == 0 & temp2$days_exposed[l] == 0) {
            agents_model$exposed[temp2$ID[l]] <- 1
            agents_model$where_infected[temp2$ID[l]] <- 3
          }
        }
      }
    }
  }
  return(agents_model)
}

# Calculates risk for each workplace and then exposes agents to that risk
work_transmit <- function(agents_model, work_risk) {
  temp1 <- dplyr::filter(agents_model, workplace > 0 & stay_home == 0 & isolated == 0 & dead == 0 & icu == 0)
  num_workplaces_mod <- max(agents_model$workplace)
  for (p in 1:num_workplaces_mod) {
    temp2 <- dplyr::filter(temp1, workplace == p)
    num_infected <- sum(temp2$infectious)
    if (num_infected > 0) {
      this_workplace_risk <- (num_infected * work_risk)
      num_in_workplace <- length(temp2$ID)
      for (l in 1:num_in_workplace) {
        k <- sample(1:10000, size = 1, replace = T)
        if (k < this_workplace_risk & temp2$recovered[l] == 0 & temp2$days_exposed[l] == 0) {
          agents_model$exposed[temp2$ID[l]] <- 1
          agents_model$where_infected[temp2$ID[l]] <- 4
        }
      }
    }
  }
  return(agents_model)
}

# Fills new taxis every timestep, calculates risk for each taxi, and then exposes agents to that risk
taxi_transmit <- function(agents_model, taxi_risk, taxi_capacity) {
  taxi_takers <- dplyr::filter(agents_model, taxi == 1 & isolated == 0 & stay_home == 0 & dead == 0 & icu == 0)
  if (length(taxi_takers$ID) > 24) {
    num_taxis <- floor(sum(taxi_takers$taxi) / taxi_capacity) - 1
    for (X in 1:num_taxis) {
      current_taxi <- dplyr::sample_n(taxi_takers, size = taxi_capacity, replace = F)
      if (sum(current_taxi$infectious) > 0) {
        current_taxi_risk <- taxi_risk * sum(current_taxi$infectious)
        at_risk <- dplyr::filter(current_taxi, infectious == 0 & exposed == 0 & recovered == 0 & dead == 0)
        if (length(at_risk$ID) > 0) {
          for (y in 1:length(at_risk$ID)) {
            k <- sample(1:10000, size = 1, replace = T)
            if (k < current_taxi_risk) {
              agents_model$exposed[at_risk$ID[y]] <- 1
              agents_model$where_infected[at_risk$ID[y]] <- 5
            }
          }
        }
      }
    }
  }
  return(agents_model)
}

# Let infectious disease progress and get cured or result in death
disease_progress <- function(agents_model, icu_capacity, time_step, TaT, tracing_efficacy, model_config, isolation_rate) {
  infectious_agents <- dplyr::filter(agents_model, infectious == 1 & dead == 0)
  for (x in 1:length(infectious_agents$ID)) {
    agents_model$days_infectious[infectious_agents$ID[x]] <- agents_model$days_infectious[infectious_agents$ID[x]] + 1
  }
  # End stay_home or end quarantine if it is time
  agents_to_release_2 <- dplyr::filter(agents_model, time_for_end_stay_home == time_step)
  if (length(agents_to_release_2$ID) > 0) {
    for (p in 1:length(agents_to_release_2$ID)) {
      agents_model$stay_home[agents_to_release_2$ID[p]] <- 0
    }
  }
  # End isolation if its time
  agents_to_release <- dplyr::filter(agents_model, time_for_end_isolation == time_step)
  if (length(agents_to_release$ID) > 0) {
    for (p in 1:length(agents_to_release$ID)) {
      agents_model$isolated[agents_to_release$ID[p]] <- 0
    }
  }
  # Some housekeeping regarding test frequency
  if (time_step == 20) {
    to_fix <- dplyr::filter(agents_model, time_for_test < 0)
    for (p in 1:length(to_fix$ID)) {
      agents_model$time_for_test[to_fix$ID[p]] <- 0
    }
  }
  # Let agents get test results if it is time for them to get the results. Note that time_for_test refers not to the day on which the test is done, but to the day of the result. time_for_test - TaT is the day of the test.
  agents_to_test <- dplyr::filter(agents_model, time_for_test == time_step & tested_positive == 0 & dead == 0)
  if (length(agents_to_test$ID) > 0) {
    for (j in 1:length(agents_to_test$ID)) {
      days_infected <- agents_to_test$days_exposed[j] + agents_to_test$days_infectious[j] - TaT
      if (agents_to_test$days_exposed[j] > 0) {
        agents_model$tested_positive[agents_to_test$ID[j]] <- get_tested(days_infected) # Only positive agents get sent to the get-tested function, negative agents test negative automatically and are counted in the next line
      }
      agents_model$times_tested[agents_to_test$ID[j]] <- agents_model$times_tested[agents_to_test$ID[j]] + 1
      if (model_config > 1 & agents_model$tested_positive[agents_to_test$ID[j]] == 1) { # This implements isolation and then contact tracing by assigning test result dates to contact agents
        agents_model$isolated[agents_to_test$ID[j]] <- 1
        agents_model$time_for_end_isolation[agents_to_test$ID[j]] <- time_step + 21
        # Trace household contacts
        all_house_contacts <- dplyr::filter(agents_model, household == agents_to_test$household[j] & tested_positive == 0 & time_for_test < (time_step - 14 & isolated == 0))
        house_contacts <- dplyr::sample_frac(all_house_contacts, size = tracing_efficacy, replace = T)
        if (length(house_contacts$ID) > 0) {
          for (x in 1:length(house_contacts$ID)) {
            agents_model$time_for_test[house_contacts$ID[x]] <- time_step + TaT + sample(1:2, size = 1, prob = c(0.5, 0.5), replace = T)
            agents_model$isolated[house_contacts$ID[x]] <- sample(0:1, size = 1, prob = c(1 - isolation_rate, isolation_rate), replace = T) # This and the following few lines describe whether agents will quarantine or just stay home
            if (agents_model$isolated[house_contacts$ID[x]] == 1) {
              agents_model$time_for_end_isolation[house_contacts$ID[x]] <- time_step + 14
            } else {
              agents_model$time_for_end_stay_home[house_contacts$ID[x]] <- time_step + 14
              agents_model$stay_home[house_contacts$ID[x]] <- 1
            }
          }
        }
        # Workplace tracing
        if (agents_to_test$workplace[j] > 0) {
          all_work_contacts <- dplyr::filter(agents_model, workplace == agents_to_test$workplace[j] & tested_positive == 0 & time_for_test < (time_step - 14 & isolated == 0))
          work_contacts <- dplyr::sample_frac(all_work_contacts, size = tracing_efficacy, replace = T)
          if (length(work_contacts$ID) > 0) {
            for (x in 1:length(work_contacts$ID)) {
              agents_model$time_for_test[work_contacts$ID[x]] <- time_step + TaT + sample(1:2, size = 1, prob = c(0.5, 0.5), replace = T)
              agents_model$isolated[work_contacts$ID[x]] <- sample(0:1, size = 1, prob = c(1 - isolation_rate, isolation_rate), replace = T) # This and the following few lines describe whether agents will quarantine or just stay home
              if (agents_model$isolated[work_contacts$ID[x]] == 1) {
                agents_model$time_for_end_isolation[work_contacts$ID[x]] <- time_step + 14
              } else {
                agents_model$time_for_end_stay_home[work_contacts$ID[x]] <- time_step + 14
                agents_model$stay_home[work_contacts$ID[x]] <- 1
              }
            }
          }
        }
        # School class tracing
        if (agents_to_test$class[j] > 0) {
          all_class_contacts <- dplyr::filter(agents_model, class == agents_to_test$class[j] & age == agents_to_test$age[j] & tested_positive == 0 & time_for_test < (time_step - 14 & isolated == 0))
          class_contacts <- dplyr::sample_frac(all_class_contacts, size = tracing_efficacy, replace = T)
          if (length(class_contacts$ID) > 0) {
            for (x in 1:length(class_contacts$ID)) {
              agents_model$time_for_test[class_contacts$ID[x]] <- time_step + TaT + sample(1:2, size = 1, prob = c(0.5, 0.5), replace = T)
              agents_model$isolated[class_contacts$ID[x]] <- sample(0:1, size = 1, prob = c(1 - isolation_rate, isolation_rate), replace = T) # This and the following few lines describe whether agents will quarantine or just stay home
              if (agents_model$isolated[class_contacts$ID[x]] == 1) {
                agents_model$time_for_end_isolation[class_contacts$ID[x]] <- time_step + 14
              } else {
                agents_model$time_for_end_stay_home[class_contacts$ID[x]] <- time_step + 14
                agents_model$time_for_end_stay_home[class_contacts$ID[x]] <- 1
              }
            }
          }
        }
      }
    }
  }
  # Let agents stay home if they are symptomatic and they are the type to stay home
  agents_to_stay_home <- dplyr::filter(agents_model, time_for_stay_home == time_step & dead == 0)
  if (length(agents_to_stay_home$ID) > 0) {
    for (j in 1:length(agents_to_stay_home$ID)) {
      agents_model$stay_home[agents_to_stay_home$ID[j]] <- 1
    }
  }
  # Let patients who need ICU go to ICU if it is time
  agents_to_icu <- dplyr::filter(agents_model, time_for_icu == time_step & dead == 0)
  if (length(agents_to_icu$ID) > 0) {
    for (j in 1:length(agents_to_icu$ID)) {
      if (sum(agents_model$icu) < icu_capacity) { # If there is room in ICU
        agents_model$icu[agents_to_icu$ID[j]] <- 1
      } else { # If ICU is full the agent dies
        agents_model$infectious[agents_to_icu$ID[j]] <- 0
        agents_model$dead[agents_to_icu$ID[j]] <- 1
        agents_model$dead_icu_full[agents_to_icu$ID[j]] <- 1
      }
    }
  }
  # Let patients be cured if it is there cure date
  agents_to_cure <- dplyr::filter(agents_model, time_for_cure == time_step & dead == 0)
  if (length(agents_to_cure$ID) > 0) {
    for (j in 1:length(agents_to_cure$ID)) {
      agents_model$infectious[agents_to_cure$ID[j]] <- 0
      agents_model$recovered[agents_to_cure$ID[j]] <- 1
      agents_model$icu[agents_to_cure$ID[j]] <- 0
      agents_model$isolated[agents_to_cure$ID[j]] <- 0
      agents_model$stay_home[agents_to_cure$ID[j]] <- 0
    }
  }
  # Let patients die if it is time for them to die
  agents_to_die <- dplyr::filter(agents_model, time_for_death == time_step & dead == 0)
  if (length(agents_to_die$ID) > 0) {
    for (j in 1:length(agents_to_die$ID)) {
      agents_model$infectious[agents_to_die$ID[j]] <- 0
      agents_model$dead[agents_to_die$ID[j]] <- 1
      agents_model$icu[agents_to_die$ID[j]] <- 0
      agents_model$isolated[agents_to_die$ID[j]] <- 0
    }
  }
  return(agents_model)
}

# Let exposed disease progress and become infectious - and on day 1 calculate natural disease course for the agent
exposed_progress <- function(agents_model, asymptomatic_rate, icu_rate, death_rate, stay_home_rate, time_step, TaT, day_first_infectious, v_day_first_infectious, days_infectious_asymptomatic, v_days_infectious_asymptomatic, symtoms_to_end_infectious, v_symtoms_to_end_infectious, symptoms_to_test, v_symptoms_to_test, symptoms_to_icu, v_symptoms_to_icu, death_risk_array) {
  exposed_agents <- dplyr::filter(agents_model, exposed == 1 & infectious == 0 & dead == 0)
  if (length(exposed_agents$ID) > 0) {
    for (x in 1:length(exposed_agents$ID)) {
      agents_model$days_exposed[exposed_agents$ID[x]] <- agents_model$days_exposed[exposed_agents$ID[x]] + 1
      # The below  checks whether it is time for an agent to become infectious - infectious agents are handled in another function
      if (agents_model$days_exposed[exposed_agents$ID[x]] > 2) {
        if (agents_model$time_for_infectiousness[exposed_agents$ID[x]] == time_step) {
          agents_model$exposed[exposed_agents$ID[x]] <- 0
          agents_model$infectious[exposed_agents$ID[x]] <- 1
        }
      }
    }
  }
  # Assign natural history of disease on day 1
  day_1_agents <- dplyr::filter(agents_model, days_exposed == 1 & dead == 0)
  if (length(day_1_agents$ID) > 0) {
    for (x in 1:length(day_1_agents$ID)) {
      # First we determine a few either or questions about the agent drawing from binomial distributions
      agents_model$will_be_symptomatic[day_1_agents$ID[x]] <- sample(0:1, size = 1, prob = c(asymptomatic_rate, 1 - asymptomatic_rate), replace = T) # Determines whether this agent will be symptomatic
      if (agents_model$will_be_symptomatic[day_1_agents$ID[x]] == 1) {
        agents_model$will_stay_home[day_1_agents$ID[x]] <- sample(0:1, size = 1, prob = c(1 - stay_home_rate, stay_home_rate), replace = T) # Determines whether this agent will stay home when symptomatic
      }
      if (agents_model$will_be_symptomatic[day_1_agents$ID[x]] == 1) {
        agents_model$will_need_icu[day_1_agents$ID[x]] <- sample(0:1, size = 1, prob = c(1 - icu_rate, icu_rate), replace = T) # Determines whether this agent will need to go to icu
      }
      if (agents_model$will_need_icu[day_1_agents$ID[x]] == 1) {
        agents_model$will_die[day_1_agents$ID[x]] <- death_test(day_1_agents$age[x], death_risk_array)
      }
      # The below determines on what dates the agent will undergo various disease transitions by drawing from a variety of distributions
      agents_model$time_for_infectiousness[day_1_agents$ID[x]] <- time_step - 1 + sample((day_first_infectious - v_day_first_infectious):(day_first_infectious + v_day_first_infectious), size = 1, replace = T) # The date the agent becomes infectious
      if (agents_model$will_be_symptomatic[day_1_agents$ID[x]] == 1) {
        agents_model$time_for_symptoms[day_1_agents$ID[x]] <- agents_model$time_for_infectiousness[day_1_agents$ID[x]] + sample((days_infectious_asymptomatic - v_days_infectious_asymptomatic):(days_infectious_asymptomatic + v_days_infectious_asymptomatic), size = 1, replace = T) # The date that agent becomes symptomatic
      }
      agents_model$time_for_test[day_1_agents$ID[x]] <- agents_model$time_for_symptoms[day_1_agents$ID[x]] + TaT + sample((symptoms_to_test - v_symptoms_to_test):(symptoms_to_test + v_symptoms_to_test), size = 1, replace = T)
      if (agents_model$will_stay_home[day_1_agents$ID[x]] == 1) {
        agents_model$time_for_stay_home[day_1_agents$ID[x]] <- (agents_model$time_for_symptoms[day_1_agents$ID[x]] + sample(1:5, size = 1, replace = T)) # The date that agent stays home
      }
      if (agents_model$will_need_icu[day_1_agents$ID[x]] == 1) {
        agents_model$time_for_icu[day_1_agents$ID[x]] <- agents_model$time_for_symptoms[day_1_agents$ID[x]] + sample(3:10, size = 1, replace = T)
      }
      if (agents_model$will_die[day_1_agents$ID[x]] == 0) { # If the agent will not die we need to provide a cure date
        agents_model$time_for_cure[day_1_agents$ID[x]] <- (time_step + sample(8:17, size = 1, replace = T))
        if (agents_model$time_for_icu[day_1_agents$ID[x]] > agents_model$time_for_cure[day_1_agents$ID[x]]) { # Just some house keeping to make sure agents do not go to hospital after getting cured
          agents_model$time_for_icu[day_1_agents$ID[x]] <- agents_model$time_for_cure[day_1_agents$ID[x]] - 1
        }
      } else { # If the agent will die we need to provide a death date
        agents_model$time_for_death[day_1_agents$ID[x]] <- agents_model$time_for_icu[day_1_agents$ID[x]] + sample(1:10, size = 1, replace = T)
        if (agents_model$time_for_icu[day_1_agents$ID[x]] > agents_model$time_for_death[day_1_agents$ID[x]]) { # Just some house keeping to make sure agents do not go to hospital after dying
          agents_model$time_for_icu[day_1_agents$ID[x]] <- agents_model$time_for_death[day_1_agents$ID[x]] - 1
        }
      }
    }
  }
  return(agents_model)
}

run_model <- function(current_run, num_runs, num_agents, num_time, initial_infections) {
  library(dplyr) # We need to reload the package here for it to work in parallel

  # Give values to key community variables
  household_size <- 4 # mean number of people in one household
  households_per_block <- 40 # mean number of households per block or geographic unit
  num_schools <- 1 # Number of schools in the community being simulated
  class_size <- 40 # mean number of pupils per class
  employment_rate <- 0.69 # Fraction of agents aged 19 to 65 who are employed
  num_workplaces <- 50 # Number of workplaces agents can be assigned to
  num_taxi_trips <- 100 # of taxi trips per time unit in the simulation
  regular_taxi_takers <- 0.8 # Fraction of employed persons who take taxis
  taxi_capacity <- 12 # Capacity of a single taxi
  icu_capacity <- 10
  TaT <- 2 # The default test turnaround time in days

  # Give values to key disease progression variables used in distributions lower down - currently using same values as SA modelling consortium
  day_first_infectious <- 4 # Mean day on which agents first become infectious
  v_day_first_infectious <- 2 # A variation factor for the above
  days_infectious_asymptomatic <- 2 # Mean number of days infectious prior to symptoms
  v_days_infectious_asymptomatic <- 1 # a variation factor for the above, in this case the previous value +1
  symtoms_to_end_infectious <- 5 # mean number of days infectious since became symptomatic
  v_symtoms_to_end_infectious <- 1 # a variation factor for the above
  symptoms_to_test <- 4 # mean days from first symtoms to test
  v_symptoms_to_test <- 3
  symptoms_to_icu <- 5
  v_symptoms_to_icu <- 1
  asymptomatic_rate <- 0.75 # fraction of infected people who are asymptomatic or not sick enough to seek care or stay home, which is slightly different from the strict understanding of asymptomatic
  # The following is not based on the SA modelling assumptions
  icu_rate <- 0.1 # Fraction of symptomatic agents requiring ICU (0.05 adjusted for the 0.75 asymptomatic rate above)
  death_rate <- 0.006 # This is the death rate in infected agents, providing ICU is available
  stay_home_rate <- 0.25 # the fraction of agents who stay home when experiencing symptoms
  tracing_efficacy <- 0.8 # The fraction of contacts who will be traced
  isolation_rate <- 0.8 # The fraction of traced agents who will quarantine - those who do not quarantine just stay at home and can still get infected from household contacts

  # Give values to transmission risk variables
  house_risk <- 180 # This is the risk per time step expressed as a fraction of 10 000 that any agent in a household would contract COVID-19 if one agent in the household has infectious TB
  block_risk <- 30 # This is the risk that any agent on a block would contract COVID-19 if one agent on the block has active COVID-19
  class_risk <- 30 # This is the risk that any agent in a school class would contract COVID-19 if one agent in the school class has infectious COVID-19
  taxi_risk <- 120 # This is the risk that any agent in a taxi would contract COVID-19  if one agent in that taxi has active COVID-19
  work_risk <- 20 # This is the risk that any agent in a workplace would contract COVID-19  if one agent at that workplace has infectious COVID-19

  death_rate <- death_rate * (1 / (icu_rate * (1 - asymptomatic_rate))) # We need this transformation since death-rate has to be scaled up since it only applies to the subset of agents needing ICU

  death_risk_array <- death_adjust(death_rate)

  model_config <- 1
  peak <- 0
  peak_time <- 0

  # Change model configuration if it is time to do so
  if (current_run > (num_runs / 3)) {
    model_config <- 2
  }
  if (current_run > ((num_runs / 3) * 2)) {
    model_config <- 3
    TaT <- 8
  }

  # Make agents and put them into settings by calling the above functions
  agents_model <- make_agents(num_agents)
  agents_model <- make_households(agents_model, household_size, num_agents)
  agents_model <- make_blocks(agents_model, households_per_block, num_agents)
  agents_model <- make_classes(agents_model, class_size)
  agents_model <- make_workplaces(agents_model, employment_rate, num_workplaces)
  agents_model <- make_taxis(agents_model, regular_taxi_takers)
  agents_model <- start_infect(agents_model, num_agents, initial_infections)

  # Calculating some variables before we start the loop
  num_houses <- max(agents_model$household) # This will have to be done inside the loop once the model allows households to change

  # The loop that takes us through the time steps and events in the model
  for (time_step in 1:num_time) {
    if (time_step == 1) {
      peak <- 1
      peak_time <- 1
    }

    # Run the transmission functions
    agents_model <- house_transmit(agents_model, house_risk, num_houses)
    agents_model <- block_transmit(agents_model, block_risk)
    agents_model <- class_transmit(agents_model, class_risk)
    agents_model <- work_transmit(agents_model, work_risk)
    agents_model <- taxi_transmit(agents_model, taxi_risk, taxi_capacity)

    # Run the disease progression functions
    agents_model <- disease_progress(agents_model, icu_capacity, time_step, TaT, tracing_efficacy, model_config, isolation_rate)
    agents_model <- exposed_progress(agents_model, asymptomatic_rate, icu_rate, death_rate, stay_home_rate, time_step, TaT, day_first_infectious, v_day_first_infectious, days_infectious_asymptomatic, v_days_infectious_asymptomatic, symtoms_to_end_infectious, v_symtoms_to_end_infectious, symptoms_to_test, v_symptoms_to_test, symptoms_to_icu, v_symptoms_to_icu, death_risk_array)

    if (peak < sum(agents_model$infectious) + sum(agents_model$exposed)) {
      peak <- sum(agents_model$infectious) + sum(agents_model$exposed)
      peak_time <- time_step
      total_infections_at_peak <- peak + sum(agents_model$recovered) + sum(agents_model$exposed)
    }

    if (time_step == num_time) {
      total_tests <- sum(agents_model$times_tested)
      total_positive_tests <- sum(agents_model$tested_positive)
      dead <- sum(agents_model$dead)
      total_infections <- sum(agents_model$recovered) + dead
      peak_infections <- peak
      dead_icu_full <- sum(agents_model$dead_icu_full)
      basic_reproductive_number <- 1 / ((num_agents - total_infections_at_peak) / num_agents)
      percent_detected <- (total_positive_tests / total_infections) * 100
      num_infected_home <- length(which(agents_model$where_infected == 1))
      num_infected_school <- length(which(agents_model$where_infected == 2))
      num_infected_block <- length(which(agents_model$where_infected == 3))
      num_infected_work <- length(which(agents_model$where_infected == 4))
      num_infected_taxi <- length(which(agents_model$where_infected == 5))
      num_quarantines <- sum(agents_model$times_quarantined)
      reporter <- c(current_run, dead, total_infections, peak_infections, peak_time, dead_icu_full, basic_reproductive_number, total_tests, percent_detected, model_config, num_infected_home, num_infected_block, num_infected_school, num_infected_work, num_infected_taxi, num_quarantines)
    }
  } # end timestep loop
  return(reporter)
} # End run_model function

# Run simulations either straight or in parallel
if (run_in_par == 0) {
  print("Running simulations...")
  for (current_run in 1:num_runs) {
    if (current_run == 1) {
      report <- rbind(run_model(current_run, num_runs, num_agents, num_time, initial_infections))
      cat("Simulation ", current_run, "completed. ")
    } else {
      report <- rbind(report, run_model(current_run, num_runs, num_agents, num_time, initial_infections))
      cat("Simulation ", current_run, "completed. ")
    }
  }
} else {
  library(doParallel)
  # Run simulations in parallel
  print("Running simulations in parallel...")
  current_run <- 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  report <- foreach(current_run = 1:num_runs, .combine = "rbind") %dopar% {
    run_model(current_run, num_runs, num_agents, num_time, initial_infections)
  }
  stopCluster(cl)
}

end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

report <- as.data.frame(report)
colnames(report)[1] <- "current_run"
colnames(report)[2] <- "dead"
colnames(report)[3] <- "total_infections"
colnames(report)[4] <- "peak_infections"
colnames(report)[5] <- "peak_time"
colnames(report)[6] <- "dead_icu_full"
colnames(report)[7] <- "basic_reproductive_number"
colnames(report)[8] <- "total_tests"
colnames(report)[9] <- "percent_detected"
colnames(report)[10] <- "model_config"
colnames(report)[11] <- "infected_home"
colnames(report)[12] <- "infected_block"
colnames(report)[13] <- "infected_school"
colnames(report)[14] <- "infected_work"
colnames(report)[15] <- "infected_taxi"
colnames(report)[16] <- "num_quarantines"

if (generate_html_output == 1) {
  library(rmarkdown)
  rmarkdown::render("ac19.rmd")
}

write.csv2(report, file = "C:/Users/Marcus/Dropbox/Programming/COVID_IBM_PAR_outputs.csv")
