################################
###### In class exercises ######
################################

#### Load required packages ####
library(adegenet)
library(poppr)
library(pegas) ## May need to install
library(tidyverse)

#### Load required data ####
load(url("https://github.com/jaredhomola/fw828_fall2025/raw/refs/heads/main/lab5.rdata"))

## Inspect steelhead object
steelhead.genind
table(steelhead.genind@pop)

### This dataset contains microhaplotypes for six strains of steelhead that are stocked into the Great Lakes
### And one set of 45 steelhead captured in the Boardman River in NW Lower Michigan
### The Boardman fish require strain assignment and have a genind pop designation of "DownstreamBoardman"

## Filter for locus missingness
steelhead.maxLocMiss20.genind <- missingno(steelhead.genind,
                                           type = "loci",
                                           cutoff = 0.2)

## Filter for individual missingness
steelhead.maxLocMiss20.maxIndMiss20.genind <- missingno(steelhead.maxLocMiss20.genind,
                                                        type = "geno",
                                                        cutoff = 0.2)

## Filter for minor allele frequency > 0.01
maf <- minorAllele(steelhead.maxLocMiss20.maxIndMiss20.genind)

mafFilter <- maf %>% 
  as_tibble(rownames = "locus") %>% 
  filter(value < 0.01 | value == 1)

steelhead.maxLocMiss20.maxIndMiss20.minMaf01.genind <- steelhead.maxLocMiss20.maxIndMiss20.genind[loc = !(locNames(steelhead.maxLocMiss20.maxIndMiss20.genind) %in% mafFilter$locus)]

## Filter for HW proportions
hwe_full <- hw.test(steelhead.maxLocMiss20.maxIndMiss20.minMaf01.genind,
                    res.type = "full")

hwe_filtered <- as_tibble(hwe_full, rownames = "locus") %>% 
  mutate(q_value = p.adjust(`Pr(chi^2 >)`, method = "fdr")) %>% 
  filter(q_value < 0.001)

steelhead.maxLocMiss20.maxIndMiss20.minMaf01.hweFilt.genind <- steelhead.maxLocMiss20.maxIndMiss20.minMaf01.genind[loc = !(locNames(steelhead.maxLocMiss20.maxIndMiss20.minMaf01.genind) %in% hwe_filtered$locus)]

## Filter for LD
pairwise_ld <- pair.ia(steelhead.maxLocMiss20.maxIndMiss20.minMaf01.hweFilt.genind)
# Remove one locus from each pair with r2 > 0.2
high_ld_pairs <- which(pairwise_ld > 0.2, arr.ind = TRUE)
# Extract locus names4
ld_filtered <- high_ld_pairs %>% 
  as_tibble(rownames = "locusPair") %>% 
  separate(locusPair, into = c("locus1", "locus2")) %>% 
  distinct(locus1)

steelhead.maxLocMiss20.maxIndMiss20.minMaf01.hweFilt.ldFilt.genind <- steelhead.maxLocMiss20.maxIndMiss20.minMaf01.hweFilt.genind[loc = !(locNames(steelhead.maxLocMiss20.maxIndMiss20.minMaf01.hweFilt.genind) %in% ld_filtered$locus1)]

## Create a separate genind object for the DownstreamBoardman individuals
downstreamBoardman.genind <- steelhead.maxLocMiss20.maxIndMiss20.minMaf01.hweFilt.ldFilt.genind[pop(steelhead.maxLocMiss20.maxIndMiss20.minMaf01.hweFilt.ldFilt.genind) == "DownstreamBoardman"]

## Create a separate genind object for the reference individuals
reference.genind <- steelhead.maxLocMiss20.maxIndMiss20.minMaf01.hweFilt.ldFilt.genind[pop(steelhead.maxLocMiss20.maxIndMiss20.minMaf01.hweFilt.ldFilt.genind) != "DownstreamBoardman"]

## Construct initial DAPC model
dapc.ref <- dapc(reference.genind,
                 var.contrib = TRUE,
                 scale = FALSE,
                 n.pca = 8,
                 n.da = 5)

#optim.a.score(dapc.ref)

## Add test fish to perform the assignment test
pred.test <- predict.dapc(dapc.ref, newdata = downstreamBoardman.genind)
pred.test$posterior %>% 
  as_tibble(rownames = "SampleID")



#############################
#### Homework assignment ####
#############################

## The lakeTrout.genind contains data from seven lake trout strains stocked into Lake Michigan
## and one set of lake trout that have been tagged for a telemetry study.
## The telemetry researchers want to determine the strain of the fish they've tagged
## to determine whether strain is indicative of movement patterns. 

### Perform assignment testing using DAPC to assign each "tagged" fish to a strain.
## Remember to filter your data before performing the assignment test.

### Use tidyverse data wrangling functions to determine the 
### most probable strain assignment for every fish. 
### Report both the strain and the probability of assignment.