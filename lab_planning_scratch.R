library(dplyr)
library(ggplot2)

# ---- Parameters ----
set.seed(123)                # for reproducibility
n_offspring <- 100            # total number of offspring
n_mothers <- 10               # number of distinct mothers
n_fathers <- 10              # number of distinct fathers
missing_prob <- 0.5          # proportion of missing MotherID or FatherID
sample_locations <- c("L1", "L2", "L3")  # sampling locations names

# ---- Generate parent IDs and assign hatcheries ----
mothers <- tibble(
  MotherID = paste0("M", 1:n_mothers))

fathers <- tibble(
  FatherID = paste0("F", 1:n_fathers))

# ---- Assign parents to offspring ----
pedigree <- tibble(
  OffspringID = paste0("O", 1:n_offspring),
  InferredMum = sample(mothers$MotherID, n_offspring, replace = TRUE),
  InferredDad = sample(fathers$FatherID, n_offspring, replace = TRUE),
  Lake_code = sample(sample_locations, n_offspring, replace = TRUE))

# ---- Randomly remove some parental info ----
pedigree <- pedigree %>%
  mutate(
    InferredMum = ifelse(runif(n_offspring) < missing_prob, "#", InferredMum),
    InferredDad = ifelse(InferredMum == "#", "*", InferredDad))

pedigree
# ---- Write to CSV ----

write.csv(pedigree, "pbt_demo.csv", row.names = FALSE)

library(dplyr)
library(tidyr)
library(stringr)

pedigree <- pedigree |> 
mutate(pbtHit = case_when(
  !str_detect(InferredDad, "[*#]") & !str_detect(InferredMum, "[*#]") ~ "stocked",
  TRUE ~ "naturally_produced"))


### Generate location data #####

Lake_data <- tibble(
  Lake_code = c("L1", "L2", "L3"),
  Lake_name = c("Houghton Lake", "Hubbard Lake", "Burt Lake"),
  Latitude = c(44.339260, 44.798656, 45.462890),
  Longitude = c(-84.730241, -83.553743, -84.267826))

# 


# --- Michigan basemap ---
states <- map_data("state")
head(states)
mi_map <- states |> 
  filter(region == "michigan",
         subregion == "south")

loc_dat <- pedigree %>%
  left_join(Lake_data)

sum_dat <- pedigree |> 
  mutate(Lake_code = as.factor(Lake_code)) |> 
  group_by(Lake_code) %>% 
  count(pbtHit) %>% 
  pivot_wider(names_from = pbtHit, values_from = n) |> 
  mutate(proportion_stocked = stocked/sum(naturally_produced+stocked),
         prop0rtion_natural = naturally_produced/sum(naturally_produced+stocked)) |> 
  left_join(Lake_data)

sum_dat  |> 
ggplot() + 
  geom_polygon(
    data = mi_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95", color = "gray60"
  ) +
  geom_scatterpie(aes(x = Longitude, y = Latitude),
    cols = c("naturally_produced", "stocked"),
    pie_scale = 15)+
  theme_bw()



##### calc metrics
ped2 |> 
  as_tibble() |> 
  dplyr::rename(MotherID = mom, 
                FatherID = dad,
                OffspringID = off) |> 
  Ns_calc()

ped2 |> 
  dplyr::rename(MotherID = mom, 
                FatherID = dad,
                OffspringID = off) |> 
  Nb_SF()


source("~/Documents/GitHub/SupCon/scripts/R/SF_functions/Ns_calc.R")
source("~/Documents/GitHub/SupCon/scripts/R/SF_functions/Nb_SF.R")

save(pedigree, Lake_data, brd.mat.safe, brd.mat.fitness, mat.sub.sample, convert2ped, 
     Ns_calc, Nb_SF, file = "Lab8.rdata")

brd.mat.safe <- function (moms = 100, dads = 100, lambda.low = 3, lambda.high = 3) 
{
  mat <- matrix(data = 0, nrow = dads, ncol = moms)
  lambda.mates <- sample(lambda.low:lambda.high, size = 1)
  
  # dad probabilities
  possible.dad.mates <- 1:dads
  df <- data.frame(mates = possible.dad.mates)
  suppressWarnings(df$prob <- (lambda.mates^df$mates) * (exp(-lambda.mates))/factorial(df$mates))
  dad.mate.probs <- df[df$prob != 0 & !is.na(df$prob), ]
  
  # mom probabilities
  possible.mom.mates <- 1:moms
  df <- data.frame(mates = possible.mom.mates)
  suppressWarnings(df$prob <- (lambda.mates^df$mates) * (exp(-lambda.mates))/factorial(df$mates))
  mom.mate.probs <- df[df$prob != 0 & !is.na(df$prob), ]
  
  mom.order <- paste(1:moms, "Female", sep = "_")
  dad.order <- paste(1:dads, "Male", sep = "_")
  parent.order <- sample(x = c(mom.order, dad.order), size = moms + dads, replace = FALSE)
  
  tmp <- data.frame(parent = parent.order, stringsAsFactors = FALSE)
  tmp$sex <- ifelse(grepl(pattern = "Female", x = tmp$parent), "Female", "Male")
  tmp$number <- as.numeric(gsub(pattern = "_.*", replacement = "", x = tmp$parent))
  tmp$mates_before <- NA
  
  for (i in 1:nrow(tmp)) {
    if (tmp$sex[i] == "Female") {
      my.col <- tmp$number[i]
      mates <- sample(dad.mate.probs$mates, size = 1, prob = dad.mate.probs$prob)
      tmp$mates_before[i] <- mates
      current.mates.list <- which(mat[, my.col] == 1)
      
      if (length(current.mates.list) == 0) {
        my.rows <- sample(x = possible.dad.mates, size = min(mates, length(possible.dad.mates)))
        mat[my.rows, my.col] <- 1
      } else {
        if (length(current.mates.list) < mates) {
          pdm2 <- possible.dad.mates[!(possible.dad.mates %in% current.mates.list)]
          mates2 <- mates - length(current.mates.list)
          my.rows <- sample(x = pdm2, size = min(mates2, length(pdm2)))
          mat[my.rows, my.col] <- 1
        }
        if (length(current.mates.list) > mates) {
          cml2 <- sample(x = current.mates.list, size = min(mates, length(current.mates.list)))
          rem.mates <- setdiff(current.mates.list, cml2)
          mat[rem.mates, my.col] <- 0
        }
      }
    } else {
      my.row <- tmp$number[i]
      mates <- sample(mom.mate.probs$mates, size = 1, prob = mom.mate.probs$prob)
      tmp$mates_before[i] <- mates
      current.mates.list <- which(mat[my.row, ] == 1)
      
      if (length(current.mates.list) == 0) {
        my.cols <- sample(x = possible.mom.mates, size = min(mates, length(possible.mom.mates)))
        mat[my.row, my.cols] <- 1
      } else {
        if (length(current.mates.list) < mates) {
          pmm2 <- possible.mom.mates[!(possible.mom.mates %in% current.mates.list)]
          mates2 <- mates - length(current.mates.list)
          my.cols <- sample(x = pmm2, size = min(mates2, length(pmm2)))
          mat[my.row, my.cols] <- 1
        }
        if (length(current.mates.list) > mates) {
          cml2 <- sample(x = current.mates.list, size = min(mates, length(current.mates.list)))
          rem.mates <- setdiff(current.mates.list, cml2)
          mat[my.row, rem.mates] <- 0
        }
      }
    }
  }
  
  return(mat)
}
