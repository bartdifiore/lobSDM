#####
## Quick trouble shooting -- lobster weirdness and covariate functions
#####

library(tidyverse)
library(gmRi)

# Get Box path to Box/Mills Lab/Projects/COCA19_Projections/data/combined/tidy_mod_data.rds
coca_path<- cs_path(box_group = "Mills Lab", subfolder = "Projects/COCA19_Projections/")
coca_dat<- readRDS(paste0(coca_path, "data/combined/tidy_mod_data.rds"))
summary(coca_dat)

# Get lobster
nefsc_code <- 301 # For reference, DFO code is 2550

lob_dat <- coca_dat |>
    filter(NMFS_SVSPP == nefsc_code)

lob_summ<- lob_dat |>
    group_by(EST_YEAR, SURVEY) |>
    summarize("Lob_Bio" = sum(BIOMASS), "Lob_Abund" = sum(ABUNDANCE))

ggplot() +
    geom_bar(data = lob_summ, aes(x = EST_YEAR, y = Lob_Abund, fill = SURVEY), stat = "identity", position = "dodge")

# Tows
load(paste0(coca_path, "data/dfo/raw/RV.GSCAT.Rdata"))
load(paste0(coca_path, "data/dfo/raw/RV.GSMISSIONS.Rdata"))

both <- GSMISSIONS |>
    left_join(GSCAT)

head(both)

tows <- unique(paste(both$MISSION, both$SETNO, both$YEAR, sep = "-"))
write.csv(tows, "~/Desktop/COCA19_DFO_Tows.csv")

### Covariate functions
# Moved over to "2_ExtractCovariates.R" function