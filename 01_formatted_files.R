# new 01 script

# load data
library(tidyverse)
library(sizeSpectra)
library(readxl)

# downloaded formatted files from teams site on Dec 18 2024

### df_borneo complete: ####
# JPZ manually changed local copy of data
# sampling_year --> year
# sampling_month --> month
# mei_taille_units --> body_length_units


### df_france_irzetal_part1 and part2 ####
## fixed? ####
# site_france sheet renamed site_data
# need to email crew to update file on teams site


### df template lento-morin ####
# site_data does not have information for site "07"

### Loop to read in data  ####

name_target <- c("site",
                 "year",
                 "month",
                 "site_date",
                 "sampling_method",
                 "sample",
                 "sampling_area",
                 "organism_group",
                 "taxon", 
                 "body_mass",
                 "body_length",
                 "body_weight_units",
                 "body_length_units",
                 "count",
                 "multiplier")

file_paths <- list.files("data/formatted_files", 
                         pattern = "*.xlsx")

data_list <- list()
list_to_fix_names <- list()

tictoc::tic()
for(i in 1:length(file_paths)){
  # check colnames
  in_names <- names(
    read_excel(
      path = paste0("data/formatted_files/",
                    path = file_paths[i])))
  # read in size data
  if(identical(name_target, in_names)){
    dat_in <- read_excel(
      path = paste0("data/formatted_files/",
                    file_paths[i]),
      col_types = c("text",
                    "numeric", 
                    "numeric", 
                    "text",
                    "text",
                    "text", 
                    "numeric",
                    "text",
                    "text", 
                    "numeric",
                    "numeric", 
                    "text",
                    "text",
                    "numeric",
                    "numeric"))
    
    # read in site_data
    df_site <- read_excel(
      path = paste0("data/formatted_files/",
                    path = file_paths[i]),
      sheet = "site_data")
    df_site <- df_site %>%
      # add organism_groups
      select(site, organism_groups, geographical_latitude, geographical_longitude) %>%
      mutate(site = as.character(site),
             geographical_latitude = as.numeric(geographical_latitude), 
             geographical_longitude = as.numeric(geographical_longitude))
    
    # combine size and site
    dat_in <- left_join(dat_in, df_site, by = "site")
    
    # add dat_id (database file path) to dat_in
    dat_in$dat_id <- file_paths[i]
  } else { # this occurs when colnames doesn't match target
    {
      list_to_fix_names[[i]] <- file_paths[i]
      data_list[[i]] <- "Need to fix data"
      next
    }
  }
  # save dat_in to list
  data_list[[i]] <- dat_in
}
tictoc::toc() # takes ~ 3 minutes
# check files
list_to_fix_names

# code to remove bad hungary file for now
# data_list[[33]] <- NULL
# length(data_list)

# map(data_list,
#     \(df) df %>% pull(geographical_latitude) %>% is.double())

dat_df <- bind_rows(data_list)

dim(dat_df)
dat_df %>% select(dat_id, site_date) %>% unique() %>% dim()

dat_df <- dat_df %>%
  group_by(dat_id, site_date) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup()

range(dat_df$group_id)

# O'Gorman ####
# need to fix O'Gorman data once it is formatted 
# someone simulated this?? - DP 7/17/2024
# Looks to me like just duplicated body mass values for each individual?

# filter data ####
# remove NA body_weight_units 
dat_df <- dat_df %>%
  filter(!is.na(body_weight_units),
         !is.na(body_mass)) %>%
  filter(body_mass>0,
         count > 0)

dim(dat_df)


# body_weight_units ####
# convert everything to mg
# what are the body weight units?
dat_df$body_weight_units %>% unique() %>% sort()

dat_df %>%
  pull(organism_group) %>%
  unique() %>%
  sort()

dat_df |>
  select(organism_group, body_weight_units) |>
  distinct() %>%
  filter(organism_group == "Invertebrates" |
           organism_group == "Invertebrate" |
           organism_group == "invertebrates" )
dat_df |>
  filter(organism_group == "Invertebrates",
         body_weight_units == "g WW") |>
  select(dat_id) |>
  unique()

dat_df |>
  filter(body_weight_units == "g WW" |
           body_weight_units == "grams/wet" |
           body_weight_units == "mg wet weight" |
           body_weight_units == "mg WW" ) |>
  select(dat_id) |>
  unique()

# need to add code to update body_weight_units
dat_df <- dat_df %>%
  mutate(body_mass = case_when(
    # convert g to mg
    # convert wet weight to dry weight
    # Using an arbitrary value of 0.25 for now
    # assuming "dry_weight" is in mg
    body_weight_units == "dry_weight" ~ body_mass,
    body_weight_units == "g" ~ body_mass *1000,
    body_weight_units == "g WW" ~ body_mass *1000 *0.25,
    body_weight_units == "grams/wet" ~ body_mass *1000 *0.25,
    body_weight_units == "M.mg" ~ body_mass,
    body_weight_units == "mg" ~ body_mass,
    body_weight_units == "mg (dry_mass)" ~ body_mass,
    body_weight_units == "mg dry mass" ~ body_mass,
    body_weight_units == "mg DW" ~ body_mass,
    body_weight_units == "mg wet weight" ~ body_mass * 0.25,
    body_weight_units == "mg WW" ~ body_mass * 0.25,
    body_weight_units == "mg_DM" ~ body_mass,
    body_weight_units == "mg_dry_mass" ~ body_mass
  )) #%>% 
  # Filter out small body sizes
  #filter(body_mass > 0.0026) # not doing this anymore since MLE_bins code starts at peak anyways

# change weight units
dat_df <- dat_df %>%
  mutate(corrected_mass_units = "mg_dw")


# chnage multiplier in fish only ####

# multiplier * 1000

dat_df %>%
  pull(organism_groups) %>%
  unique() %>%
  sort()

dat_df <- dat_df %>%
  mutate(
    multiplier = case_when(
      organism_groups == "Fish" ~ multiplier * 1000, #multiply meters^2 by a large whole number
      organism_groups == "fish" ~ multiplier * 1000,
      .default = multiplier),
    organism_groups = case_when(
      organism_groups == "fish" ~ "Fish",
      organism_groups == "invertebrates" ~ "Invertebrates",
      organism_groups == "Invertebrates, Fish" ~ "Invertebrates + Fish",
      organism_groups == "Invertebrates + fish" ~ "Invertebrates + Fish",
      .default = organism_groups),
    organism_group = case_when(
      organism_group == "fish" ~ "Fish",
      organism_group == "Invertebrate" ~ "Invertebrates",
      organism_group == "invertebrates" ~ "Invertebrates",
      .default = organism_group
    ))


dat_df %>%
  filter(is.na(geographical_latitude)) %>%
  distinct(dat_id, site)

# remove sites without latitude data
dat_df <- dat_df %>%
  filter(!is.na(geographical_latitude))


# calculate ind_n ####
# ind_n is the individuals per unit area
# for fish only it is # / km^2
# when macroinvertebrates are present it is # /m^2

# ind_n per organism group? ####
dat_df <- dat_df %>%
  group_by(group_id, organism_group) %>% # need to add organism group here???
  mutate(ind_n = (count * multiplier) / n_distinct(sample))

dat_df %>%
  filter(is.na(ind_n)) %>%
  distinct(dat_id, site_date)

# dat_df %>%
#   filter(group_id == 15987) %>%
#   View()

# filter site_date with number and range ####
# No longer doing this here ####
# filter_vector <- dat_df %>%
#   group_by(group_id, site_date) %>%
#   summarise(n = n(), 
#             max_size = log10(max(body_mass)),
#             min_size = log10(min(body_mass))) %>%
#   mutate(size_range = (max_size - min_size)) %>%
#   filter(n > 100,
#          size_range >= 3) %>%
#   pull(site_date) %>%
#   unique()
# 
# filter_vector
# 
# dat_out <- dat_df |> 
#   filter(site_date %in% filter_vector) %>%
#   ungroup()

dat_out <- dat_df %>% ungroup()
dim(dat_df)
dim(dat_out)
n_distinct(dat_df$site_date)
n_distinct(dat_out$site_date)

dat_out %>%
  distinct(organism_group)
dat_out %>%
  distinct(organism_groups)

saveRDS(dat_out, "derived_data/formatted_files_stitched_filtered_Dec-2024.RDS")

