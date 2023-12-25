#libraries ####
library(svDialogs)
library(ggplot2)
library(tidyr)
library(plyr)
library(gganimate)
library(transformr)
library(dplyr)
library(tidyverse)
library(broom)
library(magick)
library(patchwork)
library(lubridate)
library(zoo)
library(ggsignif)
library(qpcR)
library(car) 
library(FSA)
library(ggrepel)
library(av)
library(viridis)
library(gifski)
library(tcltk)

#Load and sort data ####
input_directory <- tk_choose.dir(default = getwd(), caption = "Select the folder that contains all your plaque data (e.g. plaque_1")
setwd(input_directory)
dir.create(file.path("R_output/"))
output_directory <- paste0(file.path(input_directory, "R_output/"))
setwd(input_directory)
temp <- list.files(path = input_directory, pattern = ".csv", recursive = T)

#name conditions as appropriate
condition_1 <- "SMAD2KO"
condition_2 <- "SMAD3KO"
condition_3 <- "WT"

Plaque_Centres <- read.csv("Plaque Centres.csv")

#case_when inside mutate constructs a new column that is filled based on rules applied to existing data. In this case, addding the appropriate
#condition labels to each plaque
mydata_assembled <- data.table::rbindlist(lapply(temp, read.csv), use.names=FALSE)


scale <- 0.65   #Replace with whatever your scale is (in um/pixel)
interval <- 0.5 #set the interval between frames (in hours)

#Calculate per plaque metrics ####
data_uniqueID <- mydata_assembled %>% 
  unite("plaque_ID", Metadata_Date, Metadata_FoV_ID, remove=FALSE) %>% #creates a unique label for each cell
  unite("unique_ID", plaque_ID, TrackObjects_Label, remove=FALSE) %>%
  mutate(HPI=Metadata_Frame*interval) %>% #converts frames to hours
  mutate(TrackObjects_IntegratedDistance=TrackObjects_IntegratedDistance*scale) %>% #
  mutate(TrackObjects_Displacement=TrackObjects_Displacement*scale) %>% 
  mutate(TrackObjects_DistanceTraveled=TrackObjects_DistanceTraveled*scale) %>%
  mutate(speed= rollmean(TrackObjects_DistanceTraveled, 3, na.pad=TRUE)*interval) %>% #calculates a rolling mean and then uses interval to convert distance to speed
  group_by(unique_ID, HPI) %>% 
  mutate(duplicates = if_else((any(n()==1)), "NO", "YES")) %>%  #provides the output of a logical that can be used to filter and unique_IDs that were assigned to multiple cells in cellprofiler
  group_by(unique_ID) %>% 
  filter(!any(duplicates=="YES")) %>% 
  filter(n()>6) %>%   #remove any cells that were tracked for fewer than 6 frames (useful for getting rid of transient objects resulting from mis-segmentation)
  filter(HPI<=32)

number_of_cells <- n_distinct(data_uniqueID$unique_ID)

#calculates means for time points 
data_uniqueID_means_per_time <- data_uniqueID %>% 
  group_by(Metadata_Condition, plaque_ID, HPI) %>%
  summarise(mean_displacement=mean(TrackObjects_Displacement),
            mean_integrateddistance=mean(TrackObjects_IntegratedDistance),
            mean_distancetraveled=mean(TrackObjects_DistanceTraveled)
  )

#Calculate per cell metrics ####
# Function to calculate the straight-line distance between two points
calculate_distance <- function(x1, y1, x2, y2) {
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(distance)
}

data_uniqueID_per_cell_tail <- data_uniqueID %>% 
  group_by(unique_ID) %>% 
  arrange(HPI) %>% 
  slice_tail(n=1) %>% 
  left_join(Plaque_Centres, by=c("plaque_ID"),relationship = "many-to-many") %>% 
  dplyr::select(Metadata_Condition, plaque_ID,unique_ID, HPI,plaque_ID, TrackObjects_Lifetime, plaque_centre.x, plaque_centre.y, AreaShape_Center_X, AreaShape_Center_Y) %>% #capture radial position at final time epoint
  mutate(final_radial_distance = calculate_distance(plaque_centre.x, plaque_centre.y, AreaShape_Center_X, AreaShape_Center_Y))

data_uniqueID_per_cell_head <- data_uniqueID %>% 
  group_by(unique_ID) %>% 
  arrange(HPI) %>% 
  slice_head(n=1) %>% 
  left_join(Plaque_Centres, by=c("plaque_ID"),relationship = "many-to-many") %>% 
  dplyr::select(unique_ID,plaque_centre.x, plaque_centre.y, AreaShape_Center_X, AreaShape_Center_Y) %>% 
  mutate(starting_radial_distance = calculate_distance(plaque_centre.x, plaque_centre.y, AreaShape_Center_X, AreaShape_Center_Y)) %>% #capture radial position at first time epoint
  dplyr::select(unique_ID,starting_radial_distance)

data_uniqueID_per_cell<-data_uniqueID_per_cell_tail %>% 
  left_join(data_uniqueID_per_cell_head, by=c("unique_ID")) %>%
  filter(if_all(everything(), ~ !is.na(.))) %>% 
  mutate(radial_distance=(final_radial_distance-starting_radial_distance)*scale, #calculate total radial distance traveled per cell
         radial_velocity=radial_distance/(TrackObjects_Lifetime*interval)) #calculate average radial velocity over lifetime per cell

data_uniqueID_per_cell_summary <- data_uniqueID_per_cell %>% 
  group_by(Metadata_Condition) %>% 
  dplyr::summarise(mean_radial_distance=mean(radial_distance),
                   mean_radial_velocity=mean(radial_velocity))

data_uniqueID_CP_metrics <- data_uniqueID %>% 
  group_by(unique_ID) %>% 
  dplyr::filter(TrackObjects_IntegratedDistance==max(TrackObjects_IntegratedDistance)) %>% 
  dplyr::select(unique_ID, plaque_ID, Metadata_Condition,contains("TrackObjects"))  # Controls which plaque is displayed)

data_uniqueID_CP_metrics_summary <- data_uniqueID_CP_metrics %>% 
  group_by(Metadata_Condition) %>% 
  dplyr::summarise(mean_TrackObjects_IntegratedDistance=mean(TrackObjects_IntegratedDistance))

directionality_df<-data_uniqueID_per_cell %>% 
  left_join(data_uniqueID_CP_metrics, by="unique_ID") %>% 
  mutate(directionality=radial_distance/TrackObjects_IntegratedDistance)

