# iterating sequence analysis

rm(list = ls())

#loading packages

library(tidyverse)
library(TraMineR)
library(fastcluster)
library(cluster)
library(WeightedCluster)
library(stringr)
library(NbClust)
library(modelr)

#read in traminer data
data(mvad)

#set up data from traminer website
mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school", 
                   "training")
mvad.labels <- c("employment", "further education", "higher education", 
                 "joblessness", "school", "training")
mvad.scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")


sa_fit <- function(x){
  #1. drawing a sample
  data_to_use <- x 
  
  # 2.creating sequence object
  mvad.seq <- seqdef(data_to_use, 17:86, alphabet = mvad.alphabet, states = mvad.scodes, 
                     labels = mvad.labels, xtstep = 6)

    # constructing distance matrix
  dist.om1 <- seqdist(mvad.seq, method = "OM", indel = 1, sm = "TRATE")
  
  #clustering
    clusterward1 <- agnes(dist.om1, diss = TRUE, method = "ward")
  
    #defining number of clusters to analyse
    k_clusters <- c(2:15)

    #3. constructing function for fit statistics
      medoids_fun <- function(z){
      results <- WeightedCluster::wcKMedoids(dist.om1, k = z, weights = NULL, initialclust = clusterward1)
      
      #exporting just the fit statistics
      results_k <- results$stats
      results_k <- as.data.frame(results_k)
      colnames(results_k) <- z
    
      #spitting out the results we want
    return(results_k)
  }

    # 5. running the medoids function for each number of clusters and saving
    results2 <- map(k_clusters, medoids_fun)
    #converting to dataframe
    results2 <- as.data.frame(results2)
    
    # attaching names and so on  
    stats_names <- as_tibble(c("PBC", "HG", "HGSD", "ASW", "ASWw", "CH", "R2", "CHsq", "R2sq", "HC"))
    stats_names <- stats_names %>% rename(test = value)
    results2$test <- stats_names$test
    
    # tidying up the data and converting to long
    results3 <- results2 %>%
      gather("k", "value", 1:14) %>%
      mutate(k = as.numeric(str_sub(k, 2))) %>%
      mutate(range = if_else(test == "CH" | test == "CHsq", "infinity", "one")) %>%
      mutate(best = if_else(test == "HC", "min", "max"))

      output <- list(results3, clusterward1, dist.om1) 
    
    #returning the tidy dataset
  return(output)
}

# mutate approach ----------------
by_year <- mvad %>% 
  group_by(birth_year) %>% 
  nest()

seq_obj <- function(x){
  data_to_use <- x 
  
  # 2.creating sequence object
  output <- seqdef(data_to_use, 17:86, alphabet = mvad.alphabet, states = mvad.scodes, 
                     labels = mvad.labels, xtstep = 6)
    return(output)
}

seq_dist <- function(x){
output <- seqdist(x, method = "OM", indel = 1, sm = "TRATE")
return(output)
}


sa_clust <- function(x){
  dist.om1 <- x
clusterward1 <- agnes(dist.om1, diss = TRUE, method = "ward")
return(clusterward1)
}

clust_medoids_fit <- function(x, y){
  dist.om1 <- as.matrix(x)
  clusterward1 <- y
  medoids_fun <- function(z){
    results <- WeightedCluster::wcKMedoids(dist.om1, k = z, weights = NULL, initialclust = clusterward1)
    
    #exporting just the fit statistics
    results_k <- results$stats
    results_k <- as.data.frame(results_k)
    colnames(results_k) <- z
    
    #spitting out the results we want
    return(results_k)
  }
  k_clusters <- c(2:15)
    # 5. running the medoids function for each number of clusters and saving
  results2 <- map(k_clusters, medoids_fun)
  #converting to dataframe
  results2 <- as.data.frame(results2)
  
  # attaching names and so on  
  stats_names <- as_tibble(c("PBC", "HG", "HGSD", "ASW", "ASWw", "CH", "R2", "CHsq", "R2sq", "HC"))
  stats_names <- stats_names %>% rename(test = value)
  results2$test <- stats_names$test
  
  # tidying up the data and converting to long
  results3 <- results2 %>%
    gather("k", "value", 1:14) %>%
    mutate(k = as.numeric(str_sub(k, 2))) %>%
    mutate(range = if_else(test == "CH" | test == "CHsq", "infinity", "one")) %>%
    mutate(best = if_else(test == "HC", "min", "max"))
  
  #returning the tidy dataset
  return(results3)
}

clust_medoids_cluster <- function(x, y, z){
  dist.om1 <- x
  clusterward1 <- y
  k_solution <- z
  medoids_fun <- function(z){
    results <- WeightedCluster::wcKMedoids(dist.om1, k = z, weights = NULL, initialclust = clusterward1)
      return(results)
  }
  results <- medoids_fun(k_solution)
  output <- results$clustering
  #returning the tidy dataset
  return(output)
}

wcKMedoids(dist.om1, k = 2, weights = NULL, initialclust = clusterward1)

si_plot <- function(x, y, z){
  mvad.seq <- x
  clustering <- y
  title <- z
    plot <- TraMineR::seqIplot(mvad.seq, group = clustering, border = NA, with.legend = "right", main = title)
    return(plot)
}


# making data frame ------------------------------------------------
years_list <- sample(1997:1999, size = 712, prob = c(1,3,7), replace = TRUE)




mvad$birth_year <- years_list


#feels like there should be a more striaghtforward implementation?
mvad2 <- modelr::bootstrap(mvad, 3) %>% 
  mutate(data = map(strap, as.data.frame)) %>% 
  select(data) %>% 
  unnest()


age_list <- sample(1:5, size = 2136, replace = TRUE)
mvad2$age <- age_list

mvad2 <- mutate(mvad2, cohort = birth_year - age)

# nest by year (allows us to fit the models for different years)
by_year <- mvad2 %>% 
  group_by(age, birth_year, cohort) %>% 
  nest()


# sequence analysis
data_cohort <- by_year %>% 
  group_by(age) %>% 
  mutate(sequence_object = map(data, seq_obj)) %>% 
  mutate(distance_matrix = map(sequence_object, seq_dist)) %>% 
  mutate(cluster_ward = map(distance_matrix, sa_clust)) %>% 
  mutate(fit = map2(distance_matrix, cluster_ward, clust_medoids_fit)) %>% 
  mutate(clustering = pmap(list(distance_matrix, cluster_ward, 3), clust_medoids_cluster))  %>% 
  mutate(best_k = map_dbl(fit, best_k_ch))
#  ungroup() %>% #had to ungroup this for some reason?
#  mutate(plot = pwalk(list(sequence_object, clustering, birth_year), si_plot))

#so this is all it takes to get the fit statistics!

data_cohort <- data_cohort %>% 
  mutate(best_k = map(fit, best_k)) %>% 
  select(birth_year, best_k)



# findings best cluster from fit statistics


# calculating turbulence and complexity index for sequences

data_cohort <- data_cohort %>% 
  mutate(turbulence = map2(sequence_object, FALSE, seqST)) %>% #second argument to turbulence function
  # is whether to normalise or not
  mutate(turbulence_normalised = map2(sequence_object, TRUE, seqST)) %>% 
  mutate(complexity_index = map2(sequence_object, FALSE, seqici))

data_cohort %>% 
  unnest(clustering, turbulence, turbulence_normalised, complexity_index) %>% 
  gather("measure", "statistic", 3:5) %>% 
  group_by(birth_year, clustering, measure) %>% 
  summarise(mean_statistic = mean(statistic)) %>% 
  ggplot(aes(x = birth_year, y = mean_statistic)) +
  geom_point(aes(colour = as.factor(clustering))) +
  facet_wrap(~ measure, scales = "free_y") +
  theme_bw()

boot_df <- data_cohort %>% 
  unnest(turbulence, turbulence_normalised, complexity_index) %>% 
  group_by(birth_year) %>% 
  modelr::bootstrap(100) %>% 
  mutate(data = map(strap, as.data.frame)) %>% #this might take a while on the actual data
  select(-strap) %>% 
  unnest() %>% 
  gather("measure", "value", 3:5) %>% 
  group_by(.id, birth_year, measure) %>% 
  summarise(mean_value = mean(value)) 


data_cohort %>% 
  unnest(turbulence, turbulence_normalised, complexity_index) %>% 
  gather("measure", "value", 2:4) %>% 
  group_by(birth_year, measure) %>% 
  summarise(mean_value = mean(value)) %>% 
  ggplot(aes(x = birth_year, y = mean_value, colour = measure)) +
  geom_point() +
  geom_point(data = boot_df, alpha = 0.05) +
  theme_bw() +
  facet_wrap(~ measure, scales = "free_y")
  

# so I want to calculate the mean as before, but also
# to bring with it the birth_year info

temp2 <- temp$strap[1]

as.data.frame(temp2)

# defining bootstrap funciton

seq_st_boot <- function(x, y) {
  x %>% 
    mutate(turbulence_normalised = map2(.data$sequence_object, y, seqST))
}


# mode fnction
modal <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


best_k <- function(x) {
  fit <- x
  best_k <- fit %>% 
    group_by(test) %>% 
    mutate(best_fit = ifelse(best == "max", k[value == max(value)], 0)) %>%  #need value == not just max(value)
    mutate(best_fit = ifelse(best == "min", k[value == min(value)], best_fit)) %>% 
    summarise(best_solution = first(best_fit))
  
  output <- best_k
  
  return(output)
}

best_k_ch <- function(x) {
  fit <- x
  best_k <- fit %>% 
    filter(test == "CH") %>% 
    mutate(best_fit = ifelse(best == "max", k[value == max(value)], 0)) %>%  #need value == not just max(value)
    mutate(best_fit = ifelse(best == "min", k[value == min(value)], best_fit)) %>% 
    summarise(best_solution = first(best_fit))
  
  output <- best_k[[1]]
  
  return(output)
}





# these give the best solution for the different tests
# now need to decide how to reconcile these - straight mode?
# probably need to be more theoretical than that

# need to incorporate duda hart into this as well
# perhaps easiest thing is to decide on the best indices first
# probably have to do some reading before i decide this


# to plot, presumably just add best fit column to dataset and then use that for number of clusters



# next steps:
# getting data in right APC structure to nest - yep
# extracting fit statistics from nested DF - yep
# determining which is the best solution?

sessionInfo()