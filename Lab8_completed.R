################################
###### In class exercises ######
################################

#### Load required packages ####

## to map proportion stocked fish
library(ggplot2)
library(scatterpie) # may need to install

library(tidyverse)
library(stringr)


#### Load required data ####
#load(url("https://github.com/jaredhomola/fw828_fall2025/raw/refs/heads/main/lab5.rdata"))

load("Lab8.Rdata")

## Inspect data objects
head(pedigree)
head(Lake_data)


### This rdata package contains...
## Datasets:
### pedigree: tibble with OffspringID, InferredMum, InferredDad, and Lake_code
### Lake_data: tibble with Lake_code, Lake_name, and Lat/Long

## Functions: all the functions in this lab are sourced from public GitHub repositories.

## The function used to create the breeding matrix is modified from https://github.com/ScribnerLab/SeaLampreyRapture
### brd.mat.safe() creates a binary breeding matrix. I modified this slightly to allow for large parent population sizes

## The functions used to populate and subsample the breeding matrix are found at https://github.com/nicksard/mater/
### brd.mat.fitness() populates breeding batrix based on user-defined fertilization rates
### mat.sub.sample() subsamples breeding batrix based on number of offspring sampled
### convert2ped() reshapes subsampled breeding matrix to pedigree format

## Functions used to calculate number of parents (Np) and number of breeders (Nb) are found at https://github.com/weiseell/NbdLamprey2
### Ns_calc() takes in a table of OffspringID, MotherID, and FatherID and calculates Np using the Chao estimator 
### Nb_SF() takes in a table of OffspringID, MotherID, and FatherID and calculates Nb using the sibship frequency method 

##################################
# Warm-up 
##################################


#### Parentage-based tagging  ####

## Create a new pbtHit column to identify whether each individual was stocked or naturally produced
pedigree <- pedigree |> 
  mutate(pbtHit = case_when(
    !str_detect(InferredDad, "[*#]") & !str_detect(InferredMum, "[*#]") ~ "stocked",
    TRUE ~ "naturally_produced"))


## Summarize proportion of individuals that were stocked
sum_dat <- pedigree |> 
  mutate(Lake_code = as.factor(Lake_code)) |> 
  group_by(Lake_code) %>% 
  count(pbtHit) %>% 
  pivot_wider(names_from = pbtHit, values_from = n) |> 
  mutate(proportion_stocked = stocked/sum(naturally_produced+stocked),
         prop0rtion_natural = naturally_produced/sum(naturally_produced+stocked))

## Join with lake data
loc_dat <- sum_dat %>% left_join(Lake_data)

#### Plot PBT results on map  ####

## Create polygon data of United States
states <- map_data("state")

## Check out object
head(states)

## Subset to Michigan Lower Peninsula
mi_map <- states |> 
  filter(region == "michigan",
         subregion == "south")

## plot proportion stocked on map 
loc_dat  |> 
  ggplot() + 
  geom_polygon(
    data = mi_map,
    aes(x = long, y = lat, group = group),
    fill = "gray", color = "darkgray") +
  geom_scatterpie(aes(x = Longitude, y = Latitude),
                  cols = c("naturally_produced", "stocked"),
                  pie_scale = 15)+
  theme_bw()


##################################
# Main exercise 
##################################

#### Create a breeding matrix  ####
### This activity follows the steps for breeding matrix creation from https://github.com/ScribnerLab/SeaLampreyRapture

set.seed(517)

## create breeding matrix
mat1 <- brd.mat.safe(moms = 100, dads = 100, lambda.low = 3, lambda.high = 3)
head(mat1[1:10,1:10])

## Minimum and maximum expected number of fertilized eggs per female
min.fert <- 25000
max.fert <- 100000

## Populate breeding matrix with numbers of offspring
mat.fitness <- brd.mat.fitness(mat1, min.fert = min.fert, max.fert = max.fert)
head(mat.fitness[1:10,1:10])

## Maximum expected number of offspring collected using specific gear
maxoff <- 500
minoff <- 10

## randomly selecting the number of offspring that will be sampled in a given sampling effort
n.off = round(runif(n = 1,min = minoff, max = maxoff))
n.off

# creating full pedigree with that breeding matrix, and sub-sampling each mate pair randomly
mat.sub <- mat.sub.sample(mat = mat.fitness, noff = n.off)
head(mat.sub)

# Remove zeroes
mat.sub <- mat.sub[mat.sub$off1 != 0,]

# converting to format usable by Np and Nb functions
ped2 <- convert2ped(mat.sub)
head(ped2)

#### Calculate Np and Nb
Np_est <- ped2 |> 
  as_tibble() |> 
  dplyr::rename(MotherID = mom, 
                FatherID = dad,
                OffspringID = off) |> 
  Ns_calc()
Np_est[[2]]$chao

Nb_est <- ped2 |> 
  dplyr::rename(MotherID = mom, 
                FatherID = dad,
                OffspringID = off) |> 
  Nb_SF()
Nb_est

