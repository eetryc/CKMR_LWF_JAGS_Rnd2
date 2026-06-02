#=========================#
#==== 1.Data setup ====####
#=========================#

#=====Emily Tryc=====#
#====May 27, 2026====#

## >Load Libraries =========================================================####
library(tidyverse)
library(plyr)


## >Load data ==============================================================####
original=read.csv("updated_LWF_Info.csv",h=T)

original_clean <- original %>% 
  drop_na(FinalAge)
original_clean <- original_clean %>% 
  mutate(FinalAge = as.numeric(FinalAge),
         SampleYear = as.numeric(SampleYear))


  ### Select columns ======================================================####
  
  LWF_final <- original_clean %>%
    mutate(Cohort = SampleYear - FinalAge) %>% 
    select(SampleID,Cohort,SampleYear,SampleDate,Sex) %>% 
    mutate(SampleDate = mdy(SampleDate),
           Month = month(SampleDate))
  
  
  
  
  ### Unfiltered pairs ====================================================####
  all_pair_comps <- LWF_final %>%
    rename_with(~ paste0(.x, "_1"), everything()) %>%      
    cross_join(LWF_final %>% rename_with(~ paste0(.x, "_2"), everything())) %>%
    filter(SampleID_1 != SampleID_2) %>% 
    #get rid of the ones that are the same comp, eg 1 and 2 vs 2 and 1 
    filter(SampleID_1 != SampleID_2) %>% 
    mutate(
      ID_min = pmin(SampleID_1, SampleID_2),
      ID_max = pmax(SampleID_1, SampleID_2)
    ) %>%
    filter(ID_min != ID_max) %>%
    distinct(ID_min, ID_max, .keep_all = TRUE) %>%
    select(-ID_min, -ID_max)
  
  
  
  
  
  ### POP filter ==========================================================####
  POpairs <- LWF_final %>%
    rename_with(~ paste0(.x, "_1"), everything()) %>% # add the columns fro each indiv
    cross_join(LWF_final %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
    mutate(conception = ymd(paste0(Cohort_2-1,"-10-01"))) %>% ###### Change as needed to fiddle with conception date
    # mutate to add the POP possibility (y/n) with the year and reproductive age filters
    mutate(HSPPOP_candidate = ifelse(Cohort_1 < Cohort_2 & 
                                       (Cohort_2 - Cohort_1 >= 5) & 
                                       (SampleYear_1 >= Cohort_2) & 
                                       (SampleDate_1 >= conception),"1","0")) %>% 
    #get rid of the ones that are the same comp, eg 1 and 2 vs 2 and 1 
    filter(SampleID_1 != SampleID_2) %>% 
    mutate(PairID = paste(SampleID_1, SampleID_2, sep = "_")) %>%
    mutate(older_first=ifelse(Cohort_1 < Cohort_2 |(Cohort_1 == Cohort_2 & SampleID_1 < SampleID_2), T,F) ) %>% 
    filter(older_first) %>%
    select(-older_first) %>% 
    mutate(AgeDif = Cohort_2 - Cohort_1) %>% 
    dplyr::rename(Individual_1 = SampleID_1,Individual_2 = SampleID_2) 
  
  
  
  # count numnber of yes or no to get a better idea of how they are changing
  POpairs %>% 
    dplyr::count(HSPPOP_candidate)
  
  
  POPonly <- POpairs %>% 
    filter(HSPPOP_candidate == 1) %>% 
    mutate(RObase=2)
  
  
  
  
  
  
  
  
  
  ### HSP filter ===========================================================####
  
  HSpairs <- LWF_final %>%
    rename_with(~ paste0(.x, "_1"), everything()) %>% # sort into col for each indiv
    cross_join(LWF_final %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
    mutate(conception = ymd(paste0(Cohort_2-1,"-10-01"))) %>% ###### Change as needed to fiddle with conception date
    mutate(HSPPOP_candidate = ifelse(Cohort_1 < Cohort_2 & (Cohort_1 != Cohort_2),"1","0")) %>% 
    #get rid of the ones that are the same comp, eg 1 and 2 vs 2 and 1 
    filter(SampleID_1 != SampleID_2) %>% 
    mutate(PairID = paste(SampleID_1, SampleID_2, sep = "_")) %>%
    mutate(older_first=ifelse(Cohort_1 < Cohort_2 |(Cohort_1 == Cohort_2 & SampleID_1 < SampleID_2), T,F) ) %>% 
    filter(older_first) %>%
    select(-older_first) %>% 
    mutate(AgeDif = Cohort_2 - Cohort_1) %>% 
    dplyr::rename(Individual_1 = SampleID_1,Individual_2 = SampleID_2) 
  
  # count numnber of yes or no to get a better idea of how they are changing
  HSpairs %>% 
    dplyr::count(HSPPOP_candidate)
  
  HSPonly <- HSpairs %>% 
    filter(HSPPOP_candidate == 1) %>% 
    mutate(RObase=4)
  
  
  
  
  
  
  
  ### Join tables =========================================================####
  kinships <- rbind(data.frame=POPonly,
                    data.frame=HSPonly) 

