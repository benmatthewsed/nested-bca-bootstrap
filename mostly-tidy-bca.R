#' ---
#' title: Bias-corrected and accelerated confidence intervals in the tidyverse
#' author: Ben Matthews
#' output: github_document
#' ---

library(boot)
library(tidyverse)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)


#' First I want to create my function to sit within `boot()`.

#' (See more info here - http://www.stat.ucla.edu/~rgould/252w02/bs2.pdf)

# bootstrapping subfunctin ------------------------------------------------


more_mean <- function(x, i){
  mean(x[i])
}

#' and then the function to wrap around `boot()`

my_boot <- 
  function(data, var, Rval){
    boot_dta <- boot(data[[var]], more_mean, R = Rval)
    return(boot_dta)
  }

#' This was by far the most annoying part. I didn't realise that `boot()` wanted
#' a set of values and *not* a vector. So [[]] is required to extract the
#' innards of the vector, rather than [] to present the vector itself. 
#' h/t http://r4ds.had.co.nz/vectors.html.

#' let's test this to see if it works.

test_boot <- my_boot(iris, "Sepal.Length", 500)

test_boot

#' looks good to me!


#' Now let's put together afunction to extract the confidence intervals from the boot object:

bca_cis <- function(boot_out){
  
  results <- boot.ci(boot_out, type = "bca")
  low_ci <- pluck(results, "bca", 4)
  high_ci <- pluck(results, "bca", 5)
  
  output <- 
    tribble(
      ~low_ci, ~high_ci,
      low_ci,   high_ci
    )
  
  return(output)
}

#' Let's try this on the object we just saw:

bca_cis(test_boot)

#' looks good!

#' package this all together and try it on an example on multiple groups


sep_length <- iris %>% 
  select(Sepal.Length, Species)

sep_length <- sep_length %>% 
  group_by(Species) %>% 
  nest()

sep_length <- sep_length %>% 
  mutate(boot_obj = pmap(list(data, "Sepal.Length", 500), my_boot)) %>% 
  mutate(results = map(boot_obj, bca_cis))

cis <- sep_length %>% 
  unnest(results) %>% 
  select(Species, 4:5)

cis

#' this gives us CIs in (probably unweildy) wide format. From here we can work
#' with them as we would anything else!


#' Let's calculate the mean in the normal way and then add the CIs
#' (NB it's probably possible to get the mean to come through with the CIS
#' and so avoid this step, but I'm not quite sure how :S)

results <- iris %>% 
  group_by(Species) %>% 
  summarise(mean = mean(Sepal.Length))

results <- left_join(results, cis, by = "Species")

results
