################################
###### In class exercises ######
################################

#### Load required packages ####

## to map proportion stocked fish
library(ggplot2)

library(tidyverse)
library(stringr)


#### Load required data ####
load(url("https://github.com/sdstraus/FW828_PedigreeBasedInference/raw/refs/heads/main/Lab8.rdata"))

## Inspect data objects
head(pedigree)
head(Lake_data)


### This rdata package contains...
## Datasets:
### pedigree: tibble with OffspringID, InferredMum, InferredDad, and Lake_code
### Lake_data: tibble with Lake_code, Lake_name, and Lat/Long

## Functions: 
## The function used to create the breeding matrix is modified from https://github.com/ScribnerLab/SeaLampreyRapture
### brd.mat.safe() creates a binary breeding matrix. I modified this slightly to allow for large parent population sizes

## The functions used to populate and subsample the breeding matrix are found at https://github.com/nicksard/mater/
### brd.mat.fitness() populates breeding matrix based on user-defined fertilization rates
### mat.sub.sample() subsamples breeding matrix based on number of offspring sampled
### convert2ped() reshapes subsampled breeding matrix to pedigree format

## Functions used to calculate number of parents (Np) and number of breeders (Nb) are found at https://github.com/weiseell/NbdLamprey2
### Ns_calc() takes in a table of OffspringID, MotherID, and FatherID and calculates Np using the Chao estimator 
### Nb_SF() takes in a table of OffspringID, MotherID, and FatherID and calculates Nb using the sibship frequency method 

## In addition to the above functions: 
### make_pie() creates pie-slice polygons for plotting hatchery contributions in ggplot
### simulate_breeding() combines many of the above functions to simulate multiple replicates of given conditions


##################################
# Part 1 
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
         proportion_natural = naturally_produced/sum(naturally_produced+stocked))

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

# create pies
pie_df <- loc_dat |> 
  rowwise() |>
  mutate(pie = list(make_pie(
    lon = Longitude,
    lat = Latitude,
    values = c(stocked = stocked, natural = naturally_produced)
  ))) |> 
  unnest(pie)

# plot map with pies
ggplot() +
  geom_polygon(
    data = mi_map,
    aes(x = long, y = lat, group = group),
    fill = "gray", color = "darkgray"
  ) +
  geom_polygon(
    data = pie_df,
    aes(x = x, y = y, group = interaction(Lake_code, part), fill = part),
    color = "black", linewidth = 0.2) +
  # coord_quickmap() +
  coord_fixed() +
  theme_bw() +
  labs(fill = "Origin") +
  theme(panel.grid = element_blank())

##################################
# Part 2
##################################

#### Create a breeding matrix  ####
### This activity follows the steps for breeding matrix creation from https://github.com/ScribnerLab/SeaLampreyRapture

set.seed(517)

## create breeding matrix
mat1 <- brd.mat.safe(moms = 10, dads = 10, lambda.low = 3, lambda.high = 3)
mat1
head(mat1[1:10,1:10])

## Minimum and maximum expected number of fertilized eggs per female
min.fert <- 10000
max.fert <- 50000

## Populate breeding matrix with numbers of offspring
mat.fitness <- brd.mat.fitness(mat1, min.fert = min.fert, max.fert = max.fert)
head(mat.fitness[1:10,1:10])

## Maximum expected number of offspring collected using specific gear
maxoff <- 5000
minoff <- 5000

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

#### plot ####

### simulate for n repetitions

test <- simulate_breeding(n_moms = 100, n_dads = 100, 
                  min_off = 10, max_off = 10, 
                  n_reps = 10, seed = 123)
test

multiple_sampling <- test[0, ]
offspring_sampled <- c(5, 10, 20)

for(samp in offspring_sampled){
  tmp <- simulate_breeding(n_moms = 100, n_dads = 100,
                           min_off = samp,
                           max_off = samp,
                           seed = 517, n_reps = 2)
  multiple_sampling <- rbind(multiple_sampling, tmp)}
multiple_sampling

## warnings ok for small offspring sample sizes

multiple_sampling |>
  mutate(nSampled = as.factor(nSampled)) |>
  ggplot(aes(x=nSampled, y = Ns))+
  geom_boxplot()+
  xlab("Number of offspring sampled")+
  ylab("Np")+
  theme_bw()

multiple_sampling |>
  mutate(nSampled = as.factor(nSampled)) |>
  ggplot(aes(x=nSampled, y = Nb))+
  geom_boxplot()+
  xlab("Number of offspring sampled")+
  theme_bw()


# sex_ratio <- test[0, ]
# sr_males <- c(50, 100, 200, 500)
# 
# for(n_males in sr_males){
#   tmp <- simulate_breeding(n_moms = 100, n_dads = n_males, 
#                            min_off = 500, 
#                            max_off = 500, 
#                            seed = 517, n_reps = 10)
#   sex_ratio <- rbind(sex_ratio, tmp)}
# 
# sex_ratio |> 
#   mutate(dad_ratio = as.factor(nDads/nMoms)) |> 
#   ggplot(aes(x = dad_ratio, y = Ns))+
#   geom_boxplot()+
#   ylab("Np")+
#   xlab("Number of Dads:Number of Moms")+
#   theme_bw()
# 
# sex_ratio |> 
#   mutate(dad_ratio = as.factor(nDads/nMoms)) |> 
#   ggplot(aes(x = dad_ratio, y = Nb))+
#   geom_boxplot()+
#   ylab("Nb")+
#   xlab("Number of Dads:Number of Moms")+
#   theme_bw()
#   
