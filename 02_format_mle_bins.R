#rm(list=ls())

library(tidyverse)
library(furrr)
library(sizeSpectra)
source("new_GoF_functions.R")   # load in new functions
source("bin_data.R")

raw_orig <- readRDS("derived_data/formatted_files_stitched_filtered_Dec-2024.RDS")

raw_simp <- raw_orig |>
  ungroup()
range(raw_simp$group_id)

fish_and_inverts <- raw_simp |>
  filter(organism_groups == "Invertebrates + Fish") |>
  mutate(analysis_group = "Invertebrates + Fish")

raw_simp <- raw_simp |>
  mutate(analysis_group = organism_group)

dim(raw_simp) + dim(fish_and_inverts) 

raw_simp <- bind_rows(raw_simp, fish_and_inverts)
dim(raw_simp)

range(raw_simp$group_id)

raw_simp <- raw_simp |>
  group_by(group_id, analysis_group) |>
  mutate(analysis_id = cur_group_id()) |>
  ungroup()

# raw_simp |>
#   filter(dat_id == "df_NEON.xlsx") |>
#   select(site_date, group_id, analysis_id, analysis_group) %>%
#   distinct() %>%
#    View()

range(raw_simp$group_id)
range(raw_simp$analysis_id)

site_info <- raw_simp %>%
  select(site_date, 
         dat_id,
         group_id, 
         analysis_id, analysis_group,
         #organism_group, organism_groups,
         geographical_latitude) %>%
  distinct()

raw_simp_to_split <- select(raw_simp,
                   site_date,
                   analysis_id,
                   body_mass,
                   ind_n)


# smaller dataset to work out new functions
# dat_split <- raw_simp_to_split %>%
#   filter(analysis_id < 100) %>%
#   split(.$analysis_id)
# length(dat_split)

# full data set ~15k sites
dat_split <- raw_simp_to_split %>%
  split(raw_simp$analysis_id)
length(dat_split)

vecDiff = 20 # bigger value = more estimates, 
# but longer time and some huge CI's

possibly2 <- function (.f) {
  .f <- as_mapper(.f)
  function(...) {
    tryCatch(.f(...), error = function(e)  e$message)
  }
}

# new "bin_tibbles" fxn ####

tictoc::tic()
plan(cluster, workers = 10)

binned_tibbles <- dat_split |>
  future_map(possibly2(function(.data){
    bin_tibbles(.data)
    })
  )

plan(cluster, workers = 1)
tictoc::toc() 
# time for "bin_tibbles" fxn [only 99 sites] = 4 seconds
# all sites = 85 seconds


# mle_from tibbles function
tictoc::tic()
plan(cluster, workers = 10)

mle_res_no_cut <- binned_tibbles |>
  future_map(possibly2(function(.data){
    mle_from_no_cut_tibbles(.data)
  })
  )

plan(cluster, workers = 1)
tictoc::toc() 
# time for 99 sites = 5.85 seconds
# all sites = 137

# MLE results without cutting
mle_bin_no_cut_rows <- list()
for (i in 1:length(mle_res_no_cut)){
  if(is.character(mle_res_no_cut[[i]])){
    out <- data.frame(
      analysis_id = as.numeric(names(mle_res_no_cut[i])),
      min.x = NA,
      max.x = NA,
      n = NA,
      sum_bin_counts = NA,
      MLE.b = NA,
      MLE.b.l.ci = NA,
      MLE.b.u.ci = NA,
      p.val.k1 = NA,
      consistent.k1 = NA,
      p.val.k2 = NA,
      consistent.k2 = NA,
      error = mle_res_no_cut[[i]])
  } else{
    out <- mle_res_no_cut[[i]]
    out$error <- NA
  }
  mle_bin_no_cut_rows[[i]] <- out
}

no_cut_res_df <- bind_rows(mle_bin_no_cut_rows)
# mle_bin_res_df ####

dim(no_cut_res_df)

no_cut_res_df <- no_cut_res_df %>%
  left_join(site_info, join_by(analysis_id))
dim(no_cut_res_df)

names(no_cut_res_df)

# errors mle_bin_res ####
no_cut_res_df %>%
  group_by(error) %>%
  count()

#### vecDiff = 1 ####
# time = 165
# estimates = 14,197

# # A tibble: 5 × 2
# # Groups:   error [5]
# error                                          n
# <chr>                                      <int>
# 1 Need to make vecDiff larger - see R window   592
# 2 missing value where TRUE/FALSE needed         13
# not sure what this is 
# but it's only 13 estimates
# 3 sum(binCounts) >= minCounts is not TRUE       20
# not a high enough count even 
# after combining bins

# 4 too many open devices                       1564 
# pretty sure this is the same "vecDiff" error 
# but a result of 
# parallel computing - 
#  ,can't make a new plot when there 
# is more cores than graphics devices 

# 5 NA                                         14197 
# these are the estimates

#### vecDiff = 20 ####
# time = 1327 seconds (22 minutes)
# estimates = 15,170

# # A tibble: 5 × 2
# # Groups:   error [5]
# error                                          n
# <chr>                                      <int>
# 1 Need to make vecDiff larger - see R window   500
# 2 missing value where TRUE/FALSE needed         13
# 3 sum(binCounts) >= minCounts is not TRUE      192
# 4 too many open devices                        511
# 5 NA                                         15170


# MLE results without cutting
tictoc::tic()
plan(cluster, workers = 10)

mle_res_cut <- binned_tibbles |>
  future_map(possibly2(function(.data){
    mle_from_cut_tibbles(.data)
  })
  )

plan(cluster, workers = 1)
tictoc::toc() 
# all sites = 160 seconds

mle_bin_cut_rows <- list()
for (i in 1:length(mle_res_cut)){
  if(is.character(mle_res_cut[[i]])){
    out <- data.frame(
      analysis_id = as.numeric(names(mle_res_cut[i])),
      min.x = NA,
      max.x = NA,
      n = NA,
      sum_bin_counts = NA,
      MLE.b = NA,
      MLE.b.l.ci = NA,
      MLE.b.u.ci = NA,
      p.val.k1 = NA,
      consistent.k1 = NA,
      p.val.k2 = NA,
      consistent.k2 = NA,
      error = mle_res_cut[[i]])
  } else{
    out <- mle_res_cut[[i]]
    out$error <- NA
  }
  mle_bin_cut_rows[[i]] <- out
}

cut_res_df <- bind_rows(mle_bin_cut_rows)
# mle_bin_res_df ####

dim(cut_res_df)

cut_res_df <- cut_res_df %>%
  left_join(site_info, join_by(analysis_id))
dim(cut_res_df)

names(cut_res_df)

# errors mle_bin_res ####
cut_res_df %>%
  group_by(error) %>%
  count()

# save outputs ####
cut_res_df

# mle bin results dataframe
saveRDS(binned_tibbles,
        "derived_data/bins_and_biomass.RDS")
saveRDS(no_cut_res_df,
        "derived_data/mle_no_cut_results.RDS")
saveRDS(cut_res_df,
        "derived_data/mle_cut_results.RDS")



# biomass_df|>
#   ggplot(aes(x = abs(geographical_latitude),
#              y = biomass, 
#              color = analysis_group)) +
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~analysis_group) +
#   scale_y_log10()

### lbn_plot ####
# finish this with full data ####
# only got through testing and then got sidetracked. 

# mle_bin_res_df
# lines_test <- fit_one_list(dat_split[[1]])
# # df <- lines_test
# 
# lines_coef <- function(min.x,
#                        max.x,
#                        MLE.b,
#                        sum_bin_counts){ # bincount norm?
#   xmin = min.x
#   xmax = max.x
#   x_plb <- exp(seq(log(xmin),
#                    log(xmax), 
#                    length = 10000))
#   
#   y_plb <- dPLB(x_plb,
#                 b = MLE.b, 
#                 xmin = min(x_plb),
#                 xmax = max(x_plb)) * 
#     sum_bin_counts * x_plb
#   
#   coef_out <- coef(lm(y_plb ~ x_plb))
#   return(
#     #list(
#     fit = data.frame(fit_intercept = coef_out[1],
#                     fit_slope = coef_out[2],
#                     row.names = NULL)#,
#          # x_y_plb = data.frame(x_plb = x_plb,
#          #                      y_plb = y_plb))
#   )
# }
# 
# lines_coef(lines_test$min.x,
#            lines_test$max.x,
#            lines_test$MLE.b,
#            lines_test$sum_bin_counts)
# 
# sim_line_1 <- lines_coef(lines_test$min.x,
#                          lines_test$max.x,
#                          lines_test$MLE.b,
#                          lines_test$sum_bin_counts)
# 
# sim_line_2 <- lines_coef(lines_test$min.x,
#                          lines_test$max.x,
#                          lines_test$MLE.b - 1,
#                          lines_test$sum_bin_counts)
# sim_line_3 <- lines_coef(lines_test$min.x,
#                          lines_test$max.x,
#                          lines_test$MLE.b,
#                          lines_test$sum_bin_counts + 100)
# 
# sim_test <- bind_rows(sim_line_1, sim_line_2, sim_line_3)
# 
# sim_test %>%
#   rowwise() %>%
#   mutate(x_y = pmap(list(min.x = lines_test$min.x,
#                          max.x = lines_test$max.x,
#                          fit_intercept,
#                          fit_slope, 
#                          n = 100),
#                     .f = x_y_lines)) %>%
#   unnest(x_y) %>%
#   ggplot(aes(x = x, y = y, group = fit_intercept, color = fit_intercept)) +
#   geom_line() +
#   # scale_x_log10() +
#   # scale_y_log10() +
#   #xlim (0, 2000) +
#   NULL
# 
# lines_test_2 <- fit_one_list(dat_split[[2]])
# lines_test_3 <- fit_one_list(dat_split[[3]])
# lines_test_4 <- fit_one_list(dat_split[[4]])
# lines_test_5 <- fit_one_list(dat_split[[5]])
# 
# df_test <- bind_rows(lines_test_2,
#           lines_test_3,
#           lines_test_4,
#           lines_test_5,
#           lines_test)
# 
# coef_small <- df_test %>%
#   rowwise() %>%
#   mutate(coefs = pmap(list(min.x,
#                    max.x,
#                    MLE.b,
#                    sum_bin_counts),
#                     .f = lines_coef)) %>%
#   unnest(coefs)
# 
# x_y_lines <- function(min.x, max.x, fit_intercept, fit_slope, n = 1000){
#   x = seq(min.x, max.x, length.out = n)
#   y = fit_intercept + x * fit_slope
#   return(data.frame(x = x, y = y))
# }
# 
# # x_y_lines(coef_small[1, "min.x"],
# #           coef_small[1, "max.x"],
# #           coef_small[1, "fit_intercept"],
# #           coef_small[1, "fit_slope"])
# 
# coef_small %>%
#   sample_n(10) %>%
#   rowwise() %>%
#   mutate(x_y = pmap(list(min.x,
#                          max.x,
#                          fit_intercept,
#                          fit_slope, 
#                          n = 100),
#                     .f = x_y_lines)) %>%
#   unnest(x_y) %>%
#   ggplot(aes(x = x, y = y, group = analysis_id, color = analysis_id)) +
#   geom_line() +
#   # scale_x_log10() +
#   # scale_y_log10() +
#   #xlim (0, 2000) +
#   NULL
# 
# 
# coef_small <- mle_bin_res_df %>%
#   filter(!is.na(MLE.b)) %>%
#   rowwise() %>%
#   mutate(coefs = pmap(list(min.x,
#                            max.x,
#                            MLE.b,
#                            sum_bin_counts),
#                       .f = lines_coef)) %>%
#   unnest(coefs)
# 
# lines_df <- coef_small %>%
#     rowwise() %>%
#     mutate(x_y = pmap(list(min.x,
#                            max.x,
#                            fit_intercept,
#                            fit_slope,
#                            n = 100),
#                       .f = x_y_lines)) %>%
#     unnest(x_y) 
# 
# lines_df %>%
#   select(analysis_id, x, y) %>%
#   left_join(site_info %>%
#               select(analysis_id, analysis_group)) %>%
#   ggplot(aes(x = x, y = y, group = analysis_id, color = analysis_group)) +
#   geom_line() +
#   # scale_y_log10() +
#   # scale_x_log10() +
#   facet_wrap(~analysis_group, scales = "free") +
#   theme_bw()
# 
# # how can I extract the information to re-create the LBN plots? 
# # data in bin_tibble_list$binned
# ## boxes ##
# # how to make the boxes? geom_rect, left-right = binMin Max
# # top-bottom = high and low counts??
# ## lines ##
# # LBN uses dPLB() or something I think? "lines()" in base R plotting?
# 
# # seq x values
# # calculate y values 
# # B.PLB <- dPLB(x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB)) * sum(binTibble$Number) * x.PLB
# # fit lm(B.plb ~ x.PLB)
# # extract beta_0, beta_1
# # add to df output
# # use that data frame to make vector of x-y values and plot together
# 
# length(mle_res_id)

# save outputs ####
# mle bin results dataframe
saveRDS(mle_bin_res_df,
        "derived_data/mle_bin_cutoff_gof_result.RDS")
#mle_bin_res_df <- readRDS("derived_data/mle_bin_cutoff_gof_result.RDS")
saveRDS(mle_bin_no_cut_res_df,
        "derived_data/mle_bin_no_cutoff_gof_result.RDS")


# bin tibble lists
saveRDS(bin_tibble_list,
        "derived_data/list_bin_cutoff_tibbles.RDS")
# bin_tibble_list <- readRDS("derived_data/list_bin_cutoff_tibbles.RDS")
saveRDS(bin_tibble_no_cutoff_list,
        "derived_data/list_bin_no_cutoff_tibbles.RDS")
  
# biomass results 
saveRDS(biomass_df,
        "derived_data/sum_cutoff_biomass_df.RDS")
# biomass_df <- readRDS("derived_data/sum_biomass_df.RDS")
saveRDS(biomass_no_cut_df,
        "derived_data/sum_no_cutoff_biomass_df.RDS")
#


# #biomass_df sums make sense?
# biomass_df %>%
#   group_by(group_id) %>%
#   count() %>%
#   filter(n > 2)
# 
# biomass_df |>
#   filter(group_id == 15675)
# 
# biomass_no_cut_df %>%
#   group_by(group_id) %>%
#   count() %>%
#   filter(n > 2)
# 
# biomass_no_cut_df |>
#   filter(group_id == 15675)
# delete following? ####
# pretty sure these are old notes to myself. 

# # to do list july 2024 ----------------------------------------------------
# 
# # why are some collections not fitting?
# 
# mle_bin_res_df %>%
#   filter(is.na(MLE.b)) %>%
#   pull(group_id)
# 
# 
# 
# # possibly with error return ----------------------------------------------------
# 
# possibly2 <- function (.f) {
#   .f <- as_mapper(.f)
#   function(...) {
#     tryCatch(.f(...), error = function(e)  e$message)
#   }
# }
# 
# # make a vector of just NA groups
# na_groups <- mle_bin_res_df %>%
#   filter(is.na(MLE.b)) %>%
#   pull(group_id)
# 
# dat_na <- raw_simp %>%
#   filter(group_id %in% na_groups) 
# dim(raw_simp)
# dim(dat_na)
# dat_na <- dat_na %>%
#   split(dat_na$group_id)
# length(dat_na)
# 
# tictoc::tic()
# plan(cluster, workers = 10)
# 
# vecDiff = 2 # bigger = more estimates, but huge CI's
# # 2 = 57 seconds
# 
# na_bin <- dat_na |>
#   future_map(possibly2(function(.data){
#     fit_one_list(.data, vecDiff = vecDiff)  }
#   )
#   )
# 
# plan(cluster, workers = 1)
# tictoc::toc()
# 
# #na_bin[[1]] %>% class()
# 
# 
# na_bin_res_df <- bind_rows(lapply(na_bin, as.data.frame), .id = "id")
# # View(na_bin_res_df)
# # nrow(na_bin_res_df)
# # 
# # na_bin_res_df %>%
# #   filter(is.na(MLE.b)) %>%
# #   nrow()
# # 
# # names(na_bin_res_df)
# # head(na_bin_res_df)
# 
# na_bin_res_df <- na_bin_res_df %>%
#   rename(error = "X[[i]]")
# 
# na_bin_res_df %>%
#   group_by(error) %>%
#   count()
# 
# na_bin_res_df %>%
#   filter(str_detect(error, "sum"))
# 
# na_bin_res_df %>%
#   filter(str_detect(error, "sum"))
# 

# # looking at errors one by one --------------------------------------------
# 
# 
# 
# raw_simp_this_id <- dat_split[["12333"]]
# 
# # group_id before filtering out small rows
# dat_split[["65"]]
# # only one row of data
# dat_split[["84"]]
# # only three rows of data
# dat_split[["86"]]
# # 40 rows of data, error, need to make vecDiff larger
# dat_split[["90"]]
# # 16 rows of data, error, need to make vecDiff larger
# 
# dat_split[["9697"]]
# # make vecDiff larger
# 
# 
# dat_split[["10018"]]
# # error in combine_bins
# # only ~ 3 bins after peak
# dat_split[["12335"]]
# # error in combine_bins
# # only ~ 4 bins after peak
# 
# dat_split[["12333"]]
# # GoF_res_K2 NaNs produced
# 
# 
# # running functions one by one --------------------------------------------
# 
# 
# suppress.warnings = TRUE
# counts <- dplyr::select(raw_simp_this_id,
#                         x = body_mass,
#                         counts = ind_n)
# 
# # Prob has a peak. If first index is peak then still good.
# binned_with_peak <- bin_data(counts,
#                              binWidth = "2k")
# 
# index_peak <- which.max(binned_with_peak$binVals$binSumNorm)
# 
# binned <- binned_with_peak$binVals[index_peak:nrow(binned_with_peak$binVals), ]
# 
# # Note binned is just the tibble, shortened versino of
# # binned_with_peak$binVals, as we don't need $indiv for calcs.
# 
# num.bins <- nrow(binned)
# 
# # bin breaks are the minima plus the max of the final bin:
# binBreaks <- c(dplyr::pull(binned, binMin),
#                dplyr::pull(binned, binMax)[num.bins])
# 
# binCounts <- dplyr::pull(binned,
#                          binCount)
# 
# MLEbin.res <-  calcLike(negLL.PLB.binned,
#                         p = -1.5,
#                         w = binBreaks,
#                         d = binCounts,
#                         J = length(binCounts),   # = num.bins
#                         vecDiff = vecDiff,
#                         suppress.warnings = suppress.warnings)             # increase this if hit a bound
# 
# GoF_res_K1 <- GoF_PLB(bin_breaks = binBreaks,
#                       bin_counts = binCounts,
#                       b = MLEbin.res$MLE,
#                       K = 1)
# 
# GoF_res_K2 <- GoF_PLB(bin_breaks = binBreaks,
#                       bin_counts = binCounts,
#                       b = MLEbin.res$MLE,
#                       K = 2)
# 
# 
# 
# # running combined bins line by line --------------------------------------
# 
# combine_bins <- function(binBreaks,
#                          binCounts,
#                          minCounts = 5){
#   stopifnot(length(binBreaks) == length(binCounts) + 1)
#   stopifnot(sum(binCounts) >= minCounts)
#   
#   # Check each in turn, combine with subsequent one if binCounts[i] is < minCounts:
#   combined_breaks <- vector()
#   combined_counts <- vector()
#   
#   combined_counts_i <- 1      # counter for combined_counts vector
#   combined_breaks[1] <- binBreaks[1]
#   
#   i <- 1                      # counter for binCounts
#   while(i <= length(binCounts)){
#     if(binCounts[i] >= minCounts){
#       combined_counts[combined_counts_i] <- binCounts[i]
#       combined_breaks[combined_counts_i + 1] <- binBreaks[i + 1]
#       combined_counts_i <- combined_counts_i + 1
#       i <- i + 1
#     } else {
#       # Combine the next bins to get >=5 total count
#       cumul <- cumsum(binCounts[i:length(binCounts)])
#       if(max(cumul) >= minCounts){
#         first_to_five <- min(which(cumul >= minCounts))  # first index of remaining counts
#         #  to have cumulative count of
#         #  >= minCounts (default 5,
#         #  hence first_to_five),
#         #  will be index >=2
#         
#       } else {
#         first_to_five <- length(cumul)           # just take them all and fix afterwards
#       }
#       
#       combined_counts[combined_counts_i] <- cumul[first_to_five]
#       combined_breaks[combined_counts_i + 1] <- binBreaks[i + first_to_five]  # TODO CHECK THAT, may
#       # be off by 1
#       combined_counts_i <- combined_counts_i + 1
#       i <- i + first_to_five
#     }
#   }
#   
#   # But count in final bin may be <5, if so then combine with penultimate bin
#   # (which should be >=5 by definition)
#   M = length(combined_counts)
#   if(combined_counts[M] < minCounts){
#     combined_counts[M-1] <- combined_counts[M-1] + combined_counts[M]
#     combined_counts <- combined_counts[-M]
#     combined_breaks <- combined_breaks[-M]
#   }
#   
#   stopifnot(sum(combined_counts) == sum(binCounts))
#   
#   return(list(combined_breaks = combined_breaks,
#               combined_counts = combined_counts))
# }
# 
# 
# 
# 
# fit_one_list(dat_split[["1"]], vecDiff = 10, suppress.warnings = TRUE)

# old code from Perkins: --------------------------------------------------


# # set up dummy summary file for output to be saved
# s.out <- data.frame(
#   site_date=unique(raw_orig$site_date),
#   min.x = NA,
#   max.x = NA,
#   n = NA,
#   MLE.b =NA,
#   MLE.b.l.ci = NA,
#   MLE.b.u.ci =NA,
#   p.val.k1 = NA,
#   consistent.k1 = NA,
#   p.val.k2 = NA,
#   consistent.k2 = NA)
# 
# 
# 
# 
# # data set
# # res_all_group_id <- list()
# for(i in 1:length(group_id_vec)){
#   group_id_val <- group_id_vec[i]
#   
#   
#   res_this_id <- fit_one_group_id(raw_simp,
#                                   group_id_val)
#   
#   # ISD_bin_plot_nonoverlapping(binBreaks = res_this_id$binBreaks,
#   #                             binCounts = res_this_id$binCounts,
#   #                             b.MLE = res_this_id$MLEbin.res$MLE,
#   #                             b.confMin = res_this_id$MLEbin.res$conf[1],
#   #                             b.confMax = res_this_id$MLEbin.res$conf[2])
#   # 
#   # LBN_bin_plot(binValsTibble = res_this_id$binned,
#   #              b.MLE = res_this_id$MLEbin.res$MLE,
#   #              b.confMin = res_this_id$MLEbin.res$conf[1],
#   #              b.confMax = res_this_id$MLEbin.res$conf[2],
#   #              leg.text = "(a)",
#   #              log.xy = "x",
#   #              plot.binned.fitted = TRUE)
#   # 
#   # LBN_bin_plot(#binValsTibble = res_this_id$binned, 
#   #              binValsTibble <- binValsTibble[complete.cases(binValsTibble), ], #  removes missing bins to avoid plotting error
#   #              b.MLE = res_this_id$MLEbin.res$MLE,
#   #              b.confMin = res_this_id$MLEbin.res$conf[1],
#   #              b.confMax = res_this_id$MLEbin.res$conf[2],
#   #              leg.text = "(b)",
#   #              log.xy = "xy",
#   #              xLab = expression(paste("Body mass ", italic(x), "(mg)")),
#   #              plot.binned.fitted = TRUE)
#   # 
#   #print(res_this_id$MLEbin.res)
#   #print(res_this_id$GoF_K1)
#  #print(res_this_id$GoF_K2)
#   
#   # extract output to be saved
#   s.out$MLE.b[i] <- res_this_id$MLEbin.res$MLE
#   s.out$MLE.b.l.ci[i] <- res_this_id$MLEbin.res$conf[1]
#   s.out$MLE.b.u.ci[i] <- res_this_id$MLEbin.res$conf[2]
#   
#   site.df <- subset(raw_simp, group_id == group_id_vec[i])
# 
#   s.out$min.x[i] <- min(site.df$body_mass)
#   s.out$max.x[i] <- max(site.df$body_mass)
#   s.out$n[i] <- length(site.df$body_mass)
#   
#   s.out$p.val.k1[i] <- res_this_id$GoF_K1$Pvalue
#   s.out$consistent.k1[i] <- res_this_id$GoF_K1$consistent
#   s.out$p.val.k2[i] <- res_this_id$GoF_K2$Pvalue
#   s.out$consistent.k2[i] <- res_this_id$GoF_K2$consistent
# }
# 
# head(s.out)
# 
# 
