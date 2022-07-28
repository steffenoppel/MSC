## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculating seabird population growth rates
## R script to review MSC productivity and sensitivity analysis
## Steffen Oppel (steffen.oppel@rspb.org.uk), July 2022
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## code adapted from Carneiro et al. 2020. A framework for mapping the distribution of seabirds by integrating tracking, demography and phenology. Journal of Applied Ecology 57: 514- 525.
## available at https://github.com/anacarneiro/DensityMaps/blob/master/scripts/01_demography.R
## breeding success and survival estimates taken from https://github.com/anacarneiro/DensityMaps/blob/master/metadata_files/01_demography/metadata_demography_to_run_simulation.csv


#####################################################
################ LOADING PACKAGES ###################
#####################################################

library(popbio)
library(dplyr)
library(tidyverse)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select
library(readxl)



#####################################################
########## READ IN DEMOGRAPHY METADATA ##############
#####################################################

## GENERAL DIRECTORY
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\MSCreview\\MSC")  ## copy and paste here your working directory

seabirds<-read_excel("Seabird_data.xlsx", sheet="Species")
examples<-read_excel("Productivity_analysis_SO.xlsx", sheet="BIRDS")
demography<-read.csv("https://raw.githubusercontent.com/anacarneiro/DensityMaps/master/metadata_files/01_demography/metadata_demography_to_run_simulation.csv")

## summarise average breeding success and proportion of juvenile to adult survival
## this is a placeholder for using more appropriate species-specific data
## see further adjustments for seaducks, penguins, and boobies below (brood reduction etc.)
dem.fix<-demography %>% summarise(BS=mean(BreedingSuccess),Sj=mean(AnnualJuveSurvProb),Sa=mean(AnnualAdultSurvProb), AFB=mean(AgeFirstBreeding)) %>%
  mutate(prop.sj=Sj/Sa)

## EXTRACT AND COMBINE THE DATA NEEDED
seabirds<-seabirds %>% select("Family name","Genus","Scientific name" ,"Common name","Age at first breeding (mid of range, or avg value)","Fecundity (clutch size/frequency of breeding)","Adult mortality rate (lowest if range)") %>%
  rename(Family=`Family name`,Scientific=`Scientific name`,Species=`Common name`,AFB=`Age at first breeding (mid of range, or avg value)`,fec=`Fecundity (clutch size/frequency of breeding)`,mort=`Adult mortality rate (lowest if range)`) %>%
  mutate(mort=as.numeric(mort)) %>%
  mutate(mort=ifelse(is.na(mort),(1-dem.fix$Sa)*100,mort)) %>% ### insert generic average mortality for seabirds from Carneiro et al. 2020
  mutate(surv=(100-as.numeric(mort))/100) %>%
  mutate(fec=as.numeric(fec)) %>%
  mutate(fec=ifelse(is.na(fec),1,fec)) %>% ### insert generic 1 egg clutch for seabirds)
  mutate(fec=ifelse(Family=="Anatidae",fec*0.1386,fec*0.6)) ### multiply fecundity with an average seabird breeding success of 60% on predator free islands - taken from Table 2 in Caravaggi et al. 2019, but COULD BE IMPROVED by species-specific estimates. For ducks, breeding success is much lower.

examples<-examples %>% select("Common name","Productivity score") %>%
  rename(Species=`Common name`) %>%
  filter(!is.na(Species)) %>%
  left_join(seabirds,by="Species") %>%
  mutate(surv=ifelse(Species == "Surf scoter",0.77,surv)) %>%  ## replace inappropriate survival of surf scotes with data from Krementz et al. 1997
  mutate(AFB=ifelse(AFB=="UNK",4,AFB)) %>%  ## replace AFB with dem.fix$AFB average or just 4 years for now
  mutate(fec=ifelse(Species %in% c("Southern rockhopper penguin","Masked booby"),fec*0.5,fec)) ## boobies and penguins only rear one chick out of 2 eggs lais





##########################################################
############# CALCULATE LAMBDA FOR EACH SPECIES ##########
##########################################################


for (i in 1:nrow(examples)){
  
  #### a. DEFINE THE PARAMETERS FOR THIS SPECIES ----
  Sa <- as.numeric(examples$surv[i])            ### survival of adults
  Sj <- as.numeric(examples$surv[i]*dem.fix$prop.sj)       ### survival of juveniles (assumed to be fixed proportion of adult survival)
  fec <- as.numeric(examples$fec[i])                ### breeding success
  afb <- round(as.numeric(examples$AFB[i]), 0)    ### age at first breeding AS AN INTEGER (i.e. round to 0 decimal places)
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #### b. DEFINE THE MATRIX FOR THIS SPECIES ----
  ## Model follows Abraham et al. 2016 - see diagram in manuscript.
  ## Age classes are Adult, age 1, age 2, ..., age a.f.b - 1.
  ## Modelling a population of females (= half of the total species population)
  
  dim1 <- afb      ## number of age classes = a.f.b.
  dim2 <- dim1^2   ## n cells of matrix are a.f.b. * a.f.b.
  
  species.matrix <- (rep.int(0, dim2))                        ## construct Leslie matrix, filled with 0s
  species.matrix[dim2] <- expression(Sa)                ## last row  , last col   : adults in year t surviving to be adults in year t+1
  species.matrix[dim1] <- expression(fec*0.5)         ## row 1  , last col    : production of first year juveniles in year t+1 from adults in year t
  species.matrix[1 + dim1] <- expression(Sj)          ## row 2  , col 1 (first column): immatures of age a.f.b.-1 in year t surviving to become adults in year t+1
  if(afb>2){
  for (k in 2:(dim1-1)){
    species.matrix[k*dim1 + k] <- expression(Sj)              ## row k+1, col k    : fill in the off-diagonal Sj terms - juveniles/immatures surviving from year t to t+1
  }
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #### c. CALCULATING POPULATION GROWTH RATE ----
  
  #### Create Leslie matrix using the species' vital rates
  seabird.vr <- list(Sa=Sa, Sj=Sj,
                     fec=fec,
                     afb=afb) ## extract parameters for this simulation from the table of input params
  
  A <- matrix(sapply(species.matrix, eval, seabird.vr, NULL), nrow=sqrt(length(species.matrix)), byrow=TRUE) ## create Leslie matrix
  examples$lambda[i]<-lambda(A)
  
}


#####################################################
############# SUMMARISE OUTPUT ##########
#####################################################
hist(examples$lambda)
examples %>% arrange(lambda)
fwrite(examples,"Oppel_example_output_growth_rate.csv")
