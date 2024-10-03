#----------------------------------
## Libraries
#----------------------------------

library(tidyverse)
library(sf)
library(gmRi)
source("Code/enhance_r_funcs.R")


#------------------------------------
## Comparison of potential predators
#------------------------------------

potentials <- ecodata::species_groupings %>%
  filter(SOE.20 %in% c("Piscivore", "Benthivore", "Apex Predator")) %>% 
  distinct(COMNAME) # this is the list of fish predators, benthic predators, and apex predators. Qualitatively, it is similar (albeit more comprehensive) than the list from the food habits database

potentials <- as.vector(potentials$COMNAME)

preds <- read.csv("Data/FoodHabits_lobsterpredators.csv")
preds_upper <- str_trim(str_to_upper(preds$com_name))

comparison <- setdiff(potentials, preds_upper) # There are few species on this list that I feel like could be potential lobster predators, particularly in inshore waters. The ones that stand out to me are redfish, cunner, pollock, striped bass (which are poorly represented in the survey), summer flounder, tautog. Some of these are interesting because they are also species that are range expanding (summer flounder/fluke and tautog). 

# The other point that I think is worth considering is that larger lobster may be strong predators of small lobster. Maybe worth considering for the dSEM approach down the road? Not sure we want to include large lobster in an overlap analysis? Finally, in a similar way we know that cancer crabs have antagonistic interactions with lobster (particularly as smaller size classes). Not sure if we want to consider potential lobster competitors at this stage of the analysis. 

# Finally, I think its worth considering that ~75% of the predators that have been observed with lobster in their stomachs, evidence of lobster has only been observed in <10 individual stomachs. For instance, NOAA has sampled 62K whiting stomachs and only found lobster in 3 individuals!!!

# For now I will proceed with the food habits database characterization of lobster predators and use a size threshold of 20 cm (~8"). 

pred_subset <- preds$scientific_name

#----------------------------------
## Load & filter predators
#----------------------------------
ma <- read_rds("Data/TrawlSurvey_data/mass_weight_at_length.rds") %>%
  filter(scientific_name %in% pred_subset, 
         length_cm >= 20) %>%
  mutate(date = lubridate::date(date))
  
dfo <- read_rds("Data/TrawlSurvey_data/dfo_weight_at_length.rds")  %>%
  filter(scientific_name %in% pred_subset, 
         length_cm >= 20)

nefsc <- read_rds("Data/TrawlSurvey_data/nefsc_both_weight_at_length.rds") %>% 
  mutate(scientific_name = str_to_sentence(scientific_name)) %>%
  filter(scientific_name %in% pred_subset, 
         length_cm >= 20)

me <- read_rds("Data/TrawlSurvey_data/me_both_weight_at_length.rds") %>%
  filter(scientific_name %in% pred_subset, 
         length_cm >= 20) %>% 
  mutate(date = lubridate::date(date))

pred_df <- rbind(ma, dfo, nefsc, me) %>%
  mutate(season = str_to_sentence(season))

write_rds(pred_df, "Data/Derived/combined_and_filtered_predators.rds", compress = "gz")

#----------------------------------
## Load & filter lobster
#----------------------------------

ma_lob <- read_rds("Data/TrawlSurvey_data/mass_weight_at_length.rds") %>%
  filter(scientific_name == "Homarus americanus",
         length_cm < 8.2) %>%
  mutate(date = lubridate::date(date))

dfo_lob <- read_rds("Data/TrawlSurvey_data/dfo_weight_at_length.rds")  %>%
  filter(scientific_name == "Homarus americanus") %>%
  # mutate(length_cm = length_cm*10,
  #        weight_at_length = weight_at_length*1000) %>% # This is temporary fix!!!!
  filter(length_cm < 8.2)

nefsc_lob <- read_rds("Data/TrawlSurvey_data/nefsc_both_weight_at_length.rds") %>% 
  mutate(scientific_name = str_to_sentence(scientific_name)) %>%
  filter(scientific_name == "Homarus americanus",
         length_cm < 8.2)

me_lob <- read_rds("Data/TrawlSurvey_data/me_both_weight_at_length.rds") %>%
  filter(scientific_name == "Homarus americanus",
         length_cm < 8.2) %>%
  mutate(date = lubridate::date(date))

lob_df <- rbind(ma_lob, dfo_lob, nefsc_lob, me_lob) %>%
  mutate(season = str_to_sentence(season))

lob_df %>% 
  ungroup() %>%
  distinct(trawl_id, survey) %>%
  group_by(survey) %>%
  summarize(n = n())

lob_df %>% 
  group_by(survey) %>%
  summarize(total_observed = sum(number_at_length))

lob_df %>% 
  group_by(survey, year) %>% 
  summarize(total_observed = sum(number_at_length)) %>% 
  ggplot(aes(x = year, y = total_observed))+
  geom_line(aes(color = survey))+
  theme_classic()

lob_df %>% 
  group_by(survey, year) %>% 
  summarize(total_observed = sum(number_at_length)) %>% 
  ggplot(aes(x = year, y = total_observed))+
  geom_line(aes(color = survey))+
  facet_wrap(~survey, scales = "free")+
  theme_classic()


ggplot(lob_df, aes(x = weight_at_length))+
  geom_histogram()+
  facet_wrap(~survey, scales = "free")


write_rds(lob_df, "Data/Derived/combined_and_filtered_lobsters.rds", compress = "gz")

#---------------------------------------------
## Build helper file for covariate extraction
#---------------------------------------------


ma_meta <- read_rds("Data/TrawlSurvey_data/mass_weight_at_length.rds") %>%
  mutate(date = lubridate::date(date))

dfo_meta <- read_rds("Data/TrawlSurvey_data/dfo_weight_at_length.rds")

nefsc_meta <- read_rds("Data/TrawlSurvey_data/nefsc_both_weight_at_length.rds") %>% 
  mutate(scientific_name = str_to_sentence(scientific_name))

me_meta <- read_rds("Data/TrawlSurvey_data/me_both_weight_at_length.rds") %>%
  mutate(date = lubridate::date(date))

df_meta <- rbind(ma_meta, dfo_meta, nefsc_meta, me_meta) %>%
  mutate(season = str_to_sentence(season)) %>%
  ungroup() %>% 
  distinct(longitude, latitude, trawl_id, survey, season, year, date) %>%
  rename(ID = trawl_id, DATE = date, EST_YEAR = year, DECDEG_BEGLAT = latitude, DECDEG_BEGLON = longitude)

write_rds(df_meta, "Data/Derived/for_covariate_extraction.rds", compress = "gz")







# trawl_meta <- pred_df %>% 
#   ungroup() %>% 
#   distinct(longitude, latitude, trawl_id, survey, season, year, date) %>%
#   rename(ID = trawl_id, DATE = date, EST_YEAR = year, DECDEG_BEGLAT = latitude, DECDEG_BEGLON = longitude)
# 
# 
# tow_name_check<- c("ID", "DATE", "EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON")
# all(tow_name_check %in% names(trawl_meta))
# 
# trawl_meta <- trawl_meta %>%
#   drop_na(DECDEG_BEGLAT)
# 
# write_rds(trawl_meta, "Data/Derived/for_covariate_extraction.rds")



