#=======================#
#==== 3.Base CKMR====####
#=======================#

#=====Emily Tryc=====#
#====May 27, 2026====#

# This uses the object "kinships' created from 1.Pair_Setup.R


## >Load libraries =======================================================####

library(tidyverse)
library(dplyr)
library(lubridate)
library(plyr)



  
  

## >JAGS model ============================================================####
# No survival estimation 
cat("model{
	
	surv ~ dunif(0.5, 0.95)

  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)

  # Adult numbers
  
	for(i in 1:years) {

	  Nadult[i] ~ dnorm(mu, 1/(sd^2)) T(0, 1e+09)

		TruePairs[i] ~ dbinom((RObase[i]*(surv^AgeDif[i]))/(Nadult[i]), Pair_viable_count[i])

	}


}", file = "HSPPOPsingleY.jags")









## >Define required objects ===============================================####
  
  ### Years to estimate ===================================================####
  Cohort_years <- kinships %>% 
    group_by(Cohort_2) %>% 
    dplyr::summarise(
      Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
    ) %>% 
    filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
    arrange(Cohort_2)
  
  
    #### Extract years ====================================================####
    #count to give to be i in JAGS
    years <- nrow(Cohort_years)
    
    
    #### Group together ===================================================####
    # by cohort year and age dif (differing probabilities)
    Cohort_years_cohort_age <- kinships %>%
      group_by(Cohort_2, AgeDif) %>%
      dplyr::summarise(
        Pair_viable_count_per_agedif = sum(HSPPOP_candidate > 0, na.rm = TRUE))
    

    
  

  ### Observations/year  ==================================================####
  observations <- kinships %>%
    left_join(Cohort_years_cohort_age, by = c("Cohort_2", "AgeDif"))
  
    #### Year index =======================================================####
    # to link kinship to years (just the number to match it, not the actual year)
    year_index <- factor(observations$Cohort_2, levels = Cohort_years$Cohort_2)
    observations <- observations %>%
      mutate(
        year_index = as.integer(year_index)
      )
  
  
  

    #### Sex info =========================================================####
    #Count the ones with known sex (U = 0)
    # Count females (F = 1)
    observations <- observations  %>% 
      mutate(
        sex_num = case_when(
          Sex_1 == "U" ~ NA_real_, # unknown to na
          Sex_1 == "M" ~ 1, # male is 1
          Sex_1 == "F" ~ 0 # female is 0
        )
      )
    observations <- observations %>%
      mutate(
        knownSex = as.numeric(!is.na(sex_num)),
        IsFemale = as.numeric(!is.na(sex_num) & (sex_num == 0))  # if sex_num encoded 0 = female
      )





    ### Combine data ======================================================####
    
      
      # Now collapse based on cohort 2, age difference, RObase, if sex is known, 
          # and female or male 
      
      group_vars_POP <- c("Cohort_2", "knownSex", "IsFemale", "RObase","AgeDif")
      
      
      collapsed_hsppop <- observations %>% 
        # filter(RObase == 1 | 2) %>% 
        group_by(year_index,across(all_of(group_vars_POP))) %>%
        dplyr::summarise(
          Pair_viable_count_per_agedif = n(),                   
          TruePairs = sum(true_pair, na.rm=TRUE),
          .groups = "drop"
        ) %>%
        arrange(Cohort_2) 
      
      
      collapsed_hsppop %>%  
        ggplot(aes(x = year_index, y = TruePairs)) +
        geom_line()
      
      
      collapsed_hsppop <- collapsed_hsppop %>% 
        filter(TruePairs != 0)
      
      data_POP = list( years = years,  # number of cohorts,
                       TruePairs = collapsed_hsppop$TruePairs,
                       AgeDif = collapsed_hsppop$AgeDif,
                       RObase = collapsed_hsppop$RObase,
                       nobs = nrow(collapsed_hsppop),
                       Pair_viable_count = collapsed_hsppop$Pair_viable_count_per_agedif,
                       year=collapsed_hsppop$year_index)
      
      

      
      
      
      
      
      

## >Run JAGS  =============================================================####
  ### Initial values ======================================================####
  inits_POP = function() {
    list(
      # surv = runif(1,0.01, 0.99),
      Nadult = (runif(years, min = 100, max = 2000)))
  }
  
  
  ### Parameter values ====================================================####
  params = c("Nadult")
  
  nburn <- 10000
  nchains <- 3
  niter <- 20000
  n.cores = 3
  
  ### Run  ================================================================####
  Out_CKMR_sim1_POP = jagsUI::jags(data_POP, 
                                 inits_POP,
                                 params, 
                                 "HSPPOPsingleY.jags", 
                                 n.burnin = nburn,
                                 n.chains = nchains,
                                 n.iter = niter,
                                 parallel = T,
                                 verbose = T)

  
  
  ### Plot and compare ===================================================####
  # Extract posterior summary as a data frame
  nhat_POP <- as.data.frame(Out_CKMR_sim1_POP$summary)
  
  # Keep only the Nadult rows
  nhat_POP <- nhat_POP[grep("^Nadult", rownames(nhat_POP)), ]
  
  # Add a Year index (1, 2, 3, ...)
  nhat_POP$Year <- 1:nrow(nhat_POP)
  
  # Select the useful columns
  nhat_df_POP <- nhat_POP %>%
    dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`)
  
  
  
  ggplot(nhat_df_POP, aes(x = Year, y = mean)) +
    geom_point(color = "blue") +
    geom_line(color = "blue") +
    labs(
      x = "Year (Cohort index)",
      y = expression(hat(N)[adult])) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, fill = "lightblue") +
    theme_bw() 
  

  
  
  
  
  
  
  
  
  
  