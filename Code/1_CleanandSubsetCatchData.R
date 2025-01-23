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
potentials_df <- ecodata::species_groupings %>%
  filter(SOE.20 %in% c("Piscivore", "Benthivore", "Apex Predator"))

potentials <- potentials_df %>% 
  distinct(COMNAME) # this is the list of fish predators, benthic predators, and apex predators. Qualitatively, it is similar (albeit more comprehensive) than the list from the food habits database

potentials <- as.vector(potentials_df$COMNAME)

preds <- read.csv("Data/FoodHabits_lobsterpredators.csv")
preds_upper <- str_trim(str_to_upper(preds$com_name))

comparison <- setdiff(potentials, preds_upper) # There are few species on this list that I feel like could be potential lobster predators, particularly in inshore waters. The ones that stand out to me are redfish, cunner, pollock, striped bass (which are poorly represented in the survey), summer flounder, tautog. Some of these are interesting because they are also species that are range expanding (summer flounder/fluke and tautog). 

# The other point that I think is worth considering is that larger lobster may be strong predators of small lobster. Maybe worth considering for the dSEM approach down the road? Not sure we want to include large lobster in an overlap analysis? Finally, in a similar way we know that cancer crabs have antagonistic interactions with lobster (particularly as smaller size classes). Not sure if we want to consider potential lobster competitors at this stage of the analysis. 

# Finally, I think its worth considering that ~75% of the predators that have been observed with lobster in their stomachs, evidence of lobster has only been observed in <10 individual stomachs. For instance, NOAA has sampled 62K whiting stomachs and only found lobster in 3 individuals!!!

# For now I will proceed with the food habits database characterization of lobster predators and use a size threshold of 20 cm (~8"). 
pred_subset <- preds$scientific_name

#----------------------------------
## Load & filter unique tows
#----------------------------------
ma_tows <- read_rds("Data/TrawlSurvey_data/mass_weight_at_length.rds") |>
  ungroup() |>
  dplyr::select(longitude, latitude, trawl_id, season, year, survey, date) |>
  distinct() |>
  mutate(date = lubridate::date(date))

dfo_tows <- read_rds("Data/TrawlSurvey_data/dfo_weight_at_length.rds") |>
  ungroup() |>
  dplyr::select(longitude, latitude, trawl_id, season, year, survey, date) |>
  distinct() |>
  mutate(date = lubridate::date(date))

nefsc_tows <- read_rds("Data/TrawlSurvey_data/nefsc_both_weight_at_length.rds") |>
  ungroup() |>
  dplyr::select(longitude, latitude, trawl_id, season, year, survey, date) |>
  distinct() |>
  mutate(date = lubridate::date(date))

me_tows <- read_rds("Data/TrawlSurvey_data/me_both_weight_at_length.rds") |>
  ungroup() |>
  dplyr::select(longitude, latitude, trawl_id, season, year, survey, date) |>
  distinct() |>
  mutate(date = lubridate::date(date))

tows_df <- rbind(ma_tows, dfo_tows, nefsc_tows, me_tows) %>%
  mutate(season = str_to_sentence(season))

write_rds(tows_df, "Data/Derived/all_unique_tows.rds", compress = "gz")


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

# Quick plot?
potentials_df <- potentials_df |>
  mutate_if(is.character, to_sentence_case) |>
  mutate_if(is.factor, ~ to_sentence_case(as.character(.))) |>
  select(COMNAME, SCINAME)

pred_df <- pred_df |>
  left_join(potentials_df, by = c("scientific_name" = "SCINAME"))
head(pred_df)

# Species, survey, year, season, total biomass?
pred_summ_plot <- pred_df |>
  mutate(total_weight_at_length = number_at_length*weight_at_length) |>
  group_by(COMNAME, survey, year, season) |>
  summarize("total_biomass" = sum(total_weight_at_length))

ggplot(pred_summ_plot, aes(x = year, y = total_biomass, fill = COMNAME)) +
  geom_area() +
  facet_wrap(~season)

#----------------------------------
## Load, classify and explore lobster data
#----------------------------------
ma_lob <- read_rds("Data/TrawlSurvey_data/mass_weight_at_length.rds") %>%
  filter(scientific_name == "Homarus americanus") |>
  mutate(
    life_class = ifelse(length_cm < 8.2, "juvenile", "adult"),
    date = lubridate::date(date)
  )

dfo_lob <- read_rds("Data/TrawlSurvey_data/dfo_weight_at_length.rds") %>%
  filter(scientific_name == "Homarus americanus") %>%
  mutate(
    life_class = ifelse(length_cm < 8.2, "juvenile", "adult"),
    date = lubridate::date(date)
  )

nefsc_lob <- read_rds("Data/TrawlSurvey_data/nefsc_both_weight_at_length.rds") %>%
  mutate(scientific_name = str_to_sentence(scientific_name)) %>%
  filter(scientific_name == "Homarus americanus") |>
  mutate(
    life_class = ifelse(length_cm < 8.2, "juvenile", "adult"),
    date = lubridate::date(date)
  )

me_lob <- read_rds("Data/TrawlSurvey_data/me_both_weight_at_length.rds") %>%
  filter(scientific_name == "Homarus americanus") |>
   mutate(
    life_class = ifelse(length_cm < 8.2, "juvenile", "adult"),
    date = lubridate::date(date)
  )

lob_df <- rbind(ma_lob, dfo_lob, nefsc_lob, me_lob) %>%
  mutate(season = str_to_sentence(season))

# Some quick exploration
# Tows
lob_df %>% 
  ungroup() %>%
  distinct(trawl_id, survey) %>%
  group_by(survey) %>%
  summarize(n = n())

# Total observed
lob_df %>% 
  group_by(survey, life_class) %>%
  summarize(total_observed = sum(number_at_length))

lob_df %>% 
  group_by(survey, year, life_class) %>% 
  summarize(total_observed = sum(number_at_length)) %>% 
  ggplot(aes(x = year, y = total_observed))+
  geom_line(aes(color = survey))+
  theme_classic() +
  facet_wrap(~life_class)

lob_df %>% 
  group_by(survey, year, life_class) %>% 
  summarize(total_observed = sum(number_at_length)) %>% 
  ggplot(aes(x = year, y = total_observed))+
  geom_line(aes(color = survey))+
  facet_wrap(~survey + life_class, scales = "free", ncol = 2)+
  theme_classic()


ggplot(lob_df, aes(x = weight_at_length))+
  geom_histogram()+
  facet_wrap(~survey + life_class, scales = "free", ncol = 2) 


write_rds(lob_df, "Data/Derived/combined_lobsters_by_life_class.rds", compress = "gz")



