# Chapter 3 preliminary analysis

# Occupancy and/or abundance models
# accounting for:
# 1. phenology
# 2. subtransect habitat use
# 3. clustering (for abundance)
# 4. detection probability (if possible)

library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(reshape)
library(reshape2)

# Simulate data
# Possibly use Reto's code here, or Matechou adaptation I made for phenology


# Organize monitoring data

# By week or by specific day?
# Week: could be specified in wide format
# Day: long format needed, but long might be necessary for subtransects


# Need untrimmed data with subsection counts
# Remember issues with duplicate species on same survey

rawdata <- read_csv("data/data.trim.csv")
dat <- rawdata %>% 
  rename(SiteID = SiteID.x) 


