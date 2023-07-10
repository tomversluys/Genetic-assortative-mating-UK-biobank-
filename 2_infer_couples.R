# this script infers couples (2)
###############################################################################

setwd("/Users/tomversluys/Documents")
rm(list=ls())
graphics.off()

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

# READ IN DATA
###############################################################################
# data <- fread("./newphd/data/ukbb/ukbb_old_subset.csv") 
data <- fread("<1_DATASET_HERE>") # load df output from previous script (1)

# begin couples identification
###############################################################################

data0 <- data %>%
  filter(Relation_to_people_in_house == 1 & 
           AcType > 0 & OwnRent > 0 &
           TimeatHome > 0 & Vehicles >= 0) %>%
  drop_na(HomeE, HomeN)

# create grouping variable to match pairs
###############################################################################
data1<-data0%>%
  mutate(geogroup=paste(HomeE, HomeN, sep=""))%>%
  mutate(CoupleID=paste(geogroup, Centre, AcType, OwnRent, TimeatHome, Vehicles, NoinHouse, sep="_")) %>%
  group_by(CoupleID) %>%
  mutate(CoupleID = cur_group_id()) %>%
  mutate(N_total=length(CoupleID)) %>%
  mutate(N_males_group = length(Sex[Sex==1]), N_females_group = length(Sex[Sex==0])) %>%
  group_by(geogroup) %>%
  mutate(density_proxy = length(geogroup)) %>%
  dplyr::select(-c(geogroup, density_proxy)) %>% data.table()

# limit to groups of two
###############################################################################
data2 <- data1 %>%
  group_by(CoupleID) %>%
  filter(n()==2) %>%
  as.data.table()

# create sexual orientation variable
###############################################################################
data3 <- data2 %>%
  group_by(CoupleID) %>%
  mutate(SexComposition = if_else(sum(Sex) %in% c(2, 0), "SameSex", "OppositeSex")) %>%
  data.table()

# recode ethnicity
###############################################################################
# wider to model
data4 <- data3[duplicated(CoupleID)]
data4 <- data4[order(CoupleID)]
names(data4) <- tolower(names(data4))
data5 <- data3[! ID %in% data4$id]
data5 <- data5[order(CoupleID)]
names(data5) <- toupper(names(data5))
combined <- cbind(data4, data5)

homo1 <- combined %>%
  filter(sexcomposition == "SameSex") %>%
  data.table
homo_original <- homo1 %>% mutate(Method = "Original")

hetero1 <- combined %>%
  filter(sexcomposition == "OppositeSex") %>%
  data.table
hetero_original <- hetero1 %>% mutate(Method = "Original")

# infer couples by new method (data returned to the ukbb)
###############################################################################
# check
couple_check <- fread("./newphd/data/ukbb/couples_derived.csv")

# WRANGLE
###############################################################################
# merge
derived_data2 <- left_join(couple_check, data, by = "ID") %>%
  drop_na(HomeE, HomeN) %>%
  group_by(CoupleID) %>%
  filter(n()==2) %>%
  data.table()

# make wide
###############################################################################
derived_data3 <- derived_data2[duplicated(CoupleID)]
derived_data3 <- derived_data3[order(CoupleID)]
names(derived_data3) <- tolower(names(derived_data3))
derived_data4 <- derived_data2[! ID %in% derived_data3$id]
derived_data4 <- derived_data4[order(CoupleID)]
names(derived_data4) <- toupper(names(derived_data4))
derived_combined <- cbind(derived_data3, derived_data4)
colnames(derived_combined)

# calculate birth distances between couples 
###############################################################################
derived_combined$birth_distance <- distHaversine(cbind(derived_combined$placebirth_long, derived_combined$placebirth_lat), 
                                                 cbind(derived_combined$PLACEBIRTH_LONG, derived_combined$PLACEBIRTH_LAT))
derived_combined <- derived_combined %>% 
  mutate(birth_distance = as.numeric(birth_distance)) 

sum(is.na(derived_combined$birth_distance)) # 1450
# split
homo2 <- derived_combined %>%
  filter(sexcomposition == "SameSex") %>%
  data.table
homo_new <- homo2 %>% mutate(Method = "New")

hetero2 <- derived_combined %>%
  filter(sexcomposition == "OppositeSex") %>%
  data.table
hetero_new <- hetero2 %>% mutate(Method = "New")

# save processed data
###############################################################################
fwrite(derived_combined, "<SAVE_TO_DIRECTORY>")
# fwrite(derived_combined, "./NATURE_PAPER/DATA/couples.csv")

# identify missing couples etc.
###############################################################################
new <- derived_combined %>% mutate(Method = "New") %>% 
  group_by(coupleid) %>% mutate(newid = id+ID) %>% group_by(newid) %>% filter(n() == 1) %>%
  select(id, ID, newid, coupleid, Method, sexcomposition) %>% data.table()

original <- combined %>% mutate(Method = "Original") %>% 
  group_by(coupleid) %>% mutate(newid = id+ID) %>% group_by(newid) %>% filter(n() == 1) %>%
  select(id, ID, newid, coupleid, Method, sexcomposition) %>% data.table()

# check false inclusions
original <- original %>% mutate(check = if_else(newid %in% new$newid, "true", "false")) %>%
  data.table()

original %>% group_by(sexcomposition, check) %>%
  tally()

# check false exclusions
new <- new %>% mutate(check = if_else(newid %in% original$newid, "included", "excluded")) %>%
  data.table()

new %>% group_by(sexcomposition, check) %>%
  tally()

