# this script processes the raw ukbb data prior to couple inference (1)
# data from three separate datasets is combined into a single df
###############################################################################
# setwd("/Users/tomversluys/Documents")
data_path <- "<YOUR_PATH_HERE>"
rm(list=ls())
graphics.off()

library(sp)
library(data.table)
library(dplyr)
library(geosphere)
library(rgdal)
library(dplyr)
library(tidyr)
library(ggplot2)
library(countrycode)
library(purrr)

#### LIST OF REQUIRED VARIABLES: 

# read in original ukbb data
###############################################################################
# biobank <- fread("./newphd/data/ukbb/ukbb_old_complete.csv")
biobank <- fread("<MAIN_UKBB_DATASET_HERE>") # load the main ukbb dataset

bio_subset <- biobank %>% select("eid", "31-0.0","34-0.0", "54-0.0", "670-0.0", "680-0.0",
                     "699-0.0", "709-0.0", "728-0.0", "20074-0.0", "20075-0.0", "50-0.0", "51-0.0", "21001-0.0", "21002-0.0",
                     "1687-0.0", "1697-0.0", "20022-0.0", "6138-0.0":"6138-0.5", "816-0.0", "21000-0.0", "20016-0.0", "2178-0.0", "2188-0.0",
                     "22006-0.0", "26411-0.0", "26418-0.0", "26428-0.0", "2405-0.0", "20499-0.0", "20500-0.0", "20122-0.0", "20126-0.0", "20127-0.0",  
                     "2129-0.0", "2149-0.0", "2159-0.0", "3669-0.0", "20116-0.0", "1239-0.0", "5959-0.0", "20403-0.0", "20414-0.0", "41227-0.0",
                     "6141-0.0", "22009-0.1":"22009-0.40", "22182-0.0", "20456-0.0", "129-0.0", "130-0.0", "845-0.0", 
                     "26410-0.0", "26412-0.0", "26413-0.0", "26414-0.0", "26415-0.0", "26416-0.0", "20118-0.0", "22003-0.0", 
                     "22004-0.0","22040-0.0", "53-0.0", "767-0.0", "47-0.0", "22034-0.0", "1647-0.0", "20115-0.0","738-0.0", "2744-0.0", 
                     "21022-0.0", "22027-0.0")
colnames(bio_subset)

colnames(bio_subset) <- c("EID", "Sex", "DOB", "Centre", "AcType", "OwnRent", "TimeatHome", "NoinHouse", "Vehicles",
                          "HomeE", "HomeN", "StandingHeight", "SeatedHeight", "BMI", "Weight", "BodSizeAge10", 
                          "HeightAge10", "BirthWeight", "Qualifications1", "Qualifications2", "Qualifications3", "Qualifications4",
                          "Qualifications5", "Qualifications6", "ManualJob", "Ethnicity", "FluidIntell", "OverallHealth", 
                          "LongstandingIllness", "GenEthGroup", "IncomeEng", "IncomeWal", "IncomeScot", "ChildFatherered",
                          "EverMentalhealthTreatment", "EverMentalhealthDistress", "Bipolar", "BIPorDepression", "Neur", "AnsSexHist", "NumberSexPart", "EverSSB",
                          "NumberSSBPart", "SmokingStatus", "CurrentSmoker", "PrevSmokeMostDatys", "AlchPerDay", "FreqAlch", "StatusBabyatBirth",
                          "R01",
                          "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
                          "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24",
                          "PC25", "PC26", "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", "PC33", "PC34", "PC35", "PC36",
                          "PC37", "PC38", "PC39", "PC40", "HLA", "Illicit/DrugAddiction", "PlaceBirthN", "PlaceBirthE", "AgeCompleteEducation",
                          "IndexMD", "EmploymentScore", "HealthScore", "EducationScore", "HouseScore", "CrimeScore", "HomeAreasPopDensity", "Heterozy", 
                          "HetPCcorrect", "WeeklyActivity", "Attendance_date", "length_working", "grip", "mins_activity", "birth_country_uk", "birth_country_non_uk", 
                          "household_income", "birth_weight_first_child", "recruitment_age", "missing_heterozygosity")


bio_subset <- bio_subset %>% select("EID", "Sex", "DOB", "Centre", "AcType", 
                                         "OwnRent", "TimeatHome", "NoinHouse", "Vehicles", "HomeE", 
                                         "HomeN","BMI",
                                         "Qualifications1", "Qualifications2":"Qualifications6", 
                                         "Ethnicity", "OverallHealth", 
                                         "ChildFatherered",
                                         "AnsSexHist", "NumberSexPart", "EverSSB", "NumberSSBPart", "SmokingStatus", 
                                         "R01", "PC1":"PC40", 
                                         "PlaceBirthN", "PlaceBirthE", "AgeCompleteEducation",
                                         "Heterozy", "HetPCcorrect",
                                         "Attendance_date", "birth_country_uk", 
                                         "birth_country_non_uk", "household_income", "birth_weight_first_child", "recruitment_age", "missing_heterozygosity")


bio_subset <- bio_subset %>% rename(Relation_to_people_in_house = R01, ID = EID) %>% data.table()

# add new fertility data
###############################################################################
# biobank <- fread("./newphd/data/ukbb/ukbb_old_complete.csv")
fertility_data <- fread("./newphd/data/ukbb/ukb670020.csv") # add the following variables "eid" "2774-0.0" "2734-0.0" "3829-0.0" "3839-0.0"

fertility_data <- fertility_data %>%
  select(c(1, 10, 6, 17, 21)) %>% data.table()
colnames(fertility_data) <- c("ID",
                              "Ever_mis_term_still",
                              "number_of_live_births", 
                              "number_of_stillbirths", 
                              "number_of_miscarriages")

bio_subset <- left_join(fertility_data, bio_subset, by = "ID")

# set non-male values for pregnancy variables to 0 
bio_subset1 <- bio_subset %>% 
  mutate(number_of_miscarriages = as.numeric(ifelse(Ever_mis_term_still == 0,
                                                    0, number_of_miscarriages))) %>%
  mutate(number_of_stillbirths = as.numeric(ifelse(Ever_mis_term_still == 0,
                                                   0, number_of_stillbirths))) %>%
  rowwise() %>%
  mutate(total_pregnancies = sum(number_of_live_births, number_of_stillbirths, 
                                 number_of_miscarriages, na.rm = T)) %>%
  data.table()

# add new fertility data
###############################################################################
age_data <- fread("./newphd/data/ukbb/ukbb_data_round_3.csv") # add the following variables: "eid"      "2754-0.0" "2754-1.0" "2754-2.0" "2754-3.0" "2764-0.0" "2764-1.0" "2764-2.0"
# "2764-3.0" "3872-0.0" "3872-1.0" "3872-2.0" "3872-3.0"
age_data <- age_data %>%
  select(c(1, 2, 6, 10)) %>% data.table()
colnames(age_data) <- c("ID",
                              "age_first_live_birth",
                              "age_last_live_birth", 
                              "age_prim_women_birth_first_child")

bio_subset <- left_join(bio_subset, age_data, by = "ID") %>% data.table()

# Convert UTM coordinates to long-lat
###############################################################################
# create function
wgs84 = "+init=epsg:4326"
bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 
+ellps=airy +datum=OSGB36 +units=m +no_defs'
ConvertCoordinates <- function(PlaceBirthE,PlaceBirthN) {
  out = cbind(PlaceBirthE,PlaceBirthN)
  mask = !is.na(PlaceBirthE)
  sp <-  sp::spTransform(sp::SpatialPoints(list(PlaceBirthE[mask],PlaceBirthN[mask]),proj4string=sp::CRS(bng)),sp::CRS(wgs84))
  out[mask,]=sp@coords
  out
}

###############################################################################
# note - the function below my be overridden by other packages, so run this script first
longlat_birth <- as.data.table(ConvertCoordinates(bio_subset$PlaceBirthE, bio_subset$PlaceBirthN)) %>%
  rename(PLaceBirth_long = PlaceBirthE, PLaceBirth_lat = PlaceBirthN) %>% data.table()

longlat_home <- as.data.table(ConvertCoordinates(bio_subset$HomeE, bio_subset$HomeN)) %>%
  rename(Home_long = PlaceBirthE, Home_lat = PlaceBirthN) %>% data.table()

bio_subset <- cbind(bio_subset, longlat_birth, longlat_home)

library(geosphere)
bio_subset$migration_distance <- distHaversine(cbind(bio_subset$PLaceBirth_long, bio_subset$PLaceBirth_lat),
                                                  cbind(bio_subset$Home_long, bio_subset$Home_lat))
bio_subset <- bio_subset %>%
  mutate(migration_distance = as.numeric(migration_distance))

# Recode qualifications to create "highest qualification variable"
# see https://link.springer.com/article/10.1007/s10519-019-09984-5 and Robinson (2017) on methods
###############################################################################

# recoding system
# 1	College or University degree (20)
# 2	A levels/AS levels or equivalent (13)
# 3	O levels/GCSEs or equivalent (10)
# 4	CSEs or equivalent (10)
# 5	NVQ or HND or HNC or equivalent (19)
# 6	Other professional qualifications eg: nursing, teaching (15)
# -7	None of the above (7)
# -3	Prefer not to answer (exclude)

bio_subset$NewQual <- 0
bio_subset$NewQual <- as.numeric(bio_subset$NewQual)
# Degrees (20 years)
bio_subset$NewQual[c(bio_subset$Qualifications1 %in% c(1) | bio_subset$Qualifications2 %in% c(1) | bio_subset$Qualifications3 %in% c(1) 
                     | bio_subset$Qualifications4 %in% c(1) | bio_subset$Qualifications5 %in% c(1) | bio_subset$Qualifications6 %in% c(1))] <- 20
# NVQ etc. (19 years)
bio_subset$NewQual[c(bio_subset$Qualifications1 %in% c(5) | bio_subset$Qualifications2 %in% c(5) | bio_subset$Qualifications3 %in% c(5) 
                     | bio_subset$Qualifications4 %in% c(5) | bio_subset$Qualifications5 %in% c(5) | bio_subset$Qualifications6 %in% c(5)) & bio_subset$NewQual == 0] <- 19
# Other professional Quals (15 years)
bio_subset$NewQual[c(bio_subset$Qualifications1 %in% c(6) | bio_subset$Qualifications2 %in% c(6) | bio_subset$Qualifications3 %in% c(6) 
                     | bio_subset$Qualifications4 %in% c(6) | bio_subset$Qualifications5 %in% c(6) | bio_subset$Qualifications6 %in% c(6)) & bio_subset$NewQual == 0] <- 15
# A Levels (13 years)
bio_subset$NewQual[c(bio_subset$Qualifications1 %in% c(2) | bio_subset$Qualifications2 %in% c(2) | bio_subset$Qualifications3 %in% c(2) 
                     | bio_subset$Qualifications4 %in% c(2) | bio_subset$Qualifications5 %in% c(2) | bio_subset$Qualifications6 %in% c(2)) & bio_subset$NewQual == 0] <- 13                                                                                           
# GCSEs etc. (10 years)
bio_subset$NewQual[c(bio_subset$Qualifications1 %in% c(3, 4) | bio_subset$Qualifications2 %in% c(3, 4) | bio_subset$Qualifications3 %in% c(3, 4) 
                     | bio_subset$Qualifications4 %in% c(3, 4) | bio_subset$Qualifications5 %in% c(3, 4) | bio_subset$Qualifications6 %in% c(3, 4)) & bio_subset$NewQual == 0] <- 10
# None of the above (7 years)
bio_subset$NewQual[c(bio_subset$Qualifications1 %in% c(-7) | bio_subset$Qualifications2 %in% c(-7) | bio_subset$Qualifications3 %in% c(-7) 
                     | bio_subset$Qualifications4 %in% c(-7) | bio_subset$Qualifications5 %in% c(-7) | bio_subset$Qualifications6 %in% c(-7)) & bio_subset$NewQual == 0] <- 7

# recode household income, setting values to midpoints of categories
bio_subset <- bio_subset %>% 
  mutate(household_income = as.numeric(household_income)) %>%
  mutate(household_income = case_when(household_income == 1 ~ 18000,
                                      household_income == 2 ~ 31000,
                                      household_income == 3 ~ 52000,
                                      household_income == 4 ~ 100000,
                                      household_income == 5 ~ 150000)) %>% data.table()

bio_subset <- bio_subset %>%
  mutate(old_dob = DOB, old_qual = NewQual, 
         old_bmi = BMI) %>%
  group_by(Sex, DOB) %>%
  mutate_at(c("BMI", "NewQual"),
            ~(scale(.) %>% as.vector)) %>%
  ungroup() %>% group_by(Sex) %>% mutate_at(c("DOB"), ~(scale(.) %>% as.vector)) %>%
  rename(BMIZScore = BMI, DOBZScore = DOB, NewQualZScore = NewQual) %>%
  data.frame()

derived_data2 <- bio_subset

# recode country
###############################################################################
derived_data2 <- derived_data2 %>% mutate(birth_country_uk = as.numeric(birth_country_uk), birth_country_non_uk = as.numeric(birth_country_non_uk)) %>%
  mutate(birth_country_uk = ifelse(birth_country_uk %in% c(-3, -2, -1, 6), NA, birth_country_uk)) %>%
  mutate(country = if_else(! is.na(birth_country_uk), birth_country_uk, birth_country_non_uk)) 

derived_data2 <- derived_data2 %>% mutate(birth_continent = case_when(between(country, 100, 152) ~ "Africa",
                                   between(country, 200, 251) ~ "Asia",
                                   between(country, -3, 6) | between(country, 300, 356)  ~ "Europe",
                                   between(country, 400, 427) ~ "North_America", # not USA
                                   between(country, 500, 506) ~ "Oceania",
                                   between(country, 600, 616) ~ "South_America")) %>%
  data.table()

# use external data for the remainer
###############################################################################

# recode country of birth and get coordinates
new_codes <- fread("./newphd/data/countries_codes_and_coordinates.csv")
codes <- fread("./newphd/data/coding.csv")

# first, get the general (non-ukbb) codes frame, remove quotes, rename, etc.
new_codes1 <- new_codes %>% mutate(across(everything(), ~ map_chr(.x, ~ gsub("\"", "", .x))))
new_codes1 <- new_codes1 %>% dplyr::select(1, 5, 6) %>% slice(2:257) %>%
  rename(country = countries_codes_and_coordinates, latitude_new_country = V5, longitude_new_country = V6) %>%
  mutate(country = as.factor(country), latitude_new_country = as.numeric(latitude_new_country),
         longitude_new_country = as.numeric(longitude_new_country))

# second, get old frame, recode, change caribbean, etc.
data <- derived_data2 %>%
  mutate(birth_country_non_uk = as.character(ifelse(is.na(birth_country_non_uk), "NA", birth_country_non_uk))) %>%
  mutate(PLaceBirth_long = ifelse(birth_country_non_uk == "109", -76.157227, PLaceBirth_long)) %>%
  mutate(PLaceBirth_lat = ifelse(birth_country_non_uk == "109", 15.326572, PLaceBirth_lat)) %>%
  dplyr::select(-country)

# # rename ukbb codes, etc.
codes1 <- codes %>%
  rename(birth_country_non_uk = coding) %>%
  mutate(meaning = if_else(meaning == "Tanzania", "Tanzania, United Republic of", meaning)) %>%
  mutate(meaning = if_else(meaning == "Burkina", "Burkina Faso", meaning)) %>%
  mutate(meaning = if_else(meaning == "USA", "United States", meaning)) %>%
  mutate(meaning = if_else(meaning == "Iran", "Iran, Islamic Republic of", meaning)) %>%
  mutate(meaning = if_else(meaning == "The Guianas", "French Guiana", meaning)) %>%
  mutate(meaning = if_else(meaning == "Myanmar (Burma)", "Myanmar", meaning)) %>%
  mutate(meaning = if_else(meaning == "Channel Islands", "United Kingdom", meaning)) %>%
  mutate(meaning = if_else(meaning == "Serbia/Montenegro", "Serbia", meaning)) %>%
  mutate(meaning = as.factor(meaning)) %>%
  mutate(birth_country_non_uk = as.factor(birth_country_non_uk)) %>%
  rename(country = meaning) %>%
  dplyr::select(birth_country_non_uk, country)

# bind ukbb data with ukbb codes to add country names
test1 <- merge(data, codes1, by = "birth_country_non_uk", all = TRUE) %>% data.frame()

# then bind that with general codes to add country coordinates
data <- merge(test1, new_codes1, by = "country", all = TRUE)

# check number of long/lat missing
sum(is.na(data$PLaceBirth_long)) #47,522
sum(is.na(data$PLaceBirth_lat)) #47,522

data <- data %>% mutate(PLaceBirth_long = if_else(is.na(PLaceBirth_long), longitude_new_country, PLaceBirth_long)) %>%
  mutate(PLaceBirth_lat = if_else(is.na(PLaceBirth_lat), latitude_new_country, PLaceBirth_lat)) %>%
  dplyr::select(-c(longitude_new_country, latitude_new_country)) %>%
  data.table()

# recheck
sum(is.na(data$PLaceBirth_long)) #9901
sum(is.na(data$PLaceBirth_lat))

# now check by country
data %>%
  filter(! country %in% new_codes1$country) %>%
  group_by(country) %>%
  summarise(sum(is.na(PLaceBirth_long)), sum(is.na(PLaceBirth_lat)), unique(country))
colnames(data)

d1 <- data %>% 
  mutate(V1 = as.numeric(ifelse(birth_country_non_uk == "NA", NA, paste(birth_country_non_uk)))) %>%
  mutate(birth_country_uk = as.numeric(birth_country_uk)) %>% dplyr::select(-birth_country_non_uk) %>%
  rename(birth_country_non_uk = V1)

# recode countries
d2 <- d1 %>% 
  mutate(birth_country_uk = ifelse(birth_country_uk %in% c(-3, -2, -1, 6), NA, birth_country_uk)) %>%
  mutate(country = if_else(! is.na(birth_country_uk), birth_country_uk, birth_country_non_uk)) 

d2 <- d2 %>% mutate(birth_continent = case_when(between(country, 100, 152) ~ "Africa",
                                                between(country, 200, 251) ~ "Asia",
                                                between(country, -3, 6) | between(country, 300, 356)  ~ "Europe",
                                                between(country, 400, 427) ~ "North_America", # not USA
                                                between(country, 500, 506) ~ "Oceania",
                                                between(country, 600, 616) ~ "South_America")) %>%
  data.table()

original_data <- d2 

# save processed data
fwrite(original_data, "<SAVE_TO_DIRECTORY>")
# fwrite(original_data, "./newphd/data/ukbb/ukbb_old_subset.csv")

                 