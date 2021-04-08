# devtools::install_github("andrewzm/FRK", "FRKTMB")
library("FRK") # FRKTMB branch (or >= v.2.0 if we have merged branches)
library("INLA")
library("mgcv")
library("LatticeKrig")
library("gstat")
library("ggplot2")
library("dplyr")
library("verification")
library("tidyr")
library("stringr")

## Names of packages used in the study
PACKAGES <- c("FRK", "INLA", "mgcv", "LKrig", "MRA", "gstat")

# ---- Data pre-processing ----

DIRECTORY <- "~/Dropbox/Basis_Function_Models_src/Code/Brazil_Malvinas_DEM/"
# DIRECTORY <- "~/Basis_Function_Models_src/Code/Brazil_Malvinas_DEM/"

## Load the data
load(paste0(DIRECTORY, "data/SST_sub_1000000.rda"))
load(paste0(DIRECTORY, "data/SST_sub_1000000_val.rda"))

## Rename the data
df_train <- SST_sub_1000000; rm(SST_sub_1000000)
df_test <- SST_sub_1000000_val; rm(SST_sub_1000000_val)

## Remove columns which will not be used
df_train$bias  <- df_test$bias  <- NULL
df_train$error <- df_test$error <- NULL

## Remove impossible locations, and remove repetitions
rm_invalid_locations <- function(df) {
  df %>%
    subset(lon <= 180 & lon >= -180 & lat <= 90 & lat >= -90) %>%
    distinct(lon, lat, .keep_all = TRUE) 
}

df_train <- rm_invalid_locations(df_train)
df_test <- rm_invalid_locations(df_test)

## Detrend the data following Zammit-Mangion and Rougier (2020).
## Here, we define the residuals from a linear model, which will be used as
## the response variable for the analysis.
df_train$lat2 <- df_train$lat^2
df_test$lat2 <- df_test$lat^2
df_train$Z <- residuals(lm(sst ~ 1 + lat + lat2, data = df_train))
df_test$Z <- residuals(lm(sst ~ 1 + lat + lat2, data = df_test))

## Keep a copy of the original training data, for plotting over the world
df_train_full <- df_train

## Subset the region corresponding to the Brazil-Malvinas Confluence zone:
between <- function(number,bounds) { number >= min(bounds) & number <= max(bounds)} 
subset_BM_confluence <- function(df, BM_box) {
  df %>% subset(
    between(lon, BM_box[, "lon"]) & between(lat, BM_box[, "lat"])
  )
}
BM_box <- cbind(lon = c(-60, -48), 
                lat = c(-50, -35))
df_train <- subset_BM_confluence(df_train, BM_box)
df_test  <- subset_BM_confluence(df_test, BM_box)


# ---- Data visualisation ----

nasa_palette <- c("#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7","#0f81f3",
                  "#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e","#aefc2a","#e9fc0d","#f6da0c","#f5a009",
                  "#f6780a","#f34a09","#f2210a","#f50008","#d90009","#a80109","#730005")

## Show the world map, but with a small box indicating the confluence zone. 
## Then, in the second panel, show the data that is in this region. Say that 
## we randomly split this data into a training and a test set. 

## A layer which will be added to most of the plots throughout the script:
map_layer <- geom_map(
  data = map_data("world"), map = map_data("world"),
  aes(group = group, map_id = region),
  fill = "white", colour = "black", size = 0.1
) 

## NB: restricting the colour scale to [-8, 8] prevents outliers dominating
# plot_of_SST_data_world <- ggplot(df_train[seq(1, nrow(df_train_full), 100), ]) +
plot_of_SST_data_world <- ggplot(df_train_full) +
  geom_point(aes(lon, lat, colour = pmin(pmax(Z, -8), 8)), pch = 46) +
  # scale_colour_distiller(palette = "Spectral", name = expression(degree*C)) +
  scale_colour_gradientn(colours = nasa_palette, name = expression(degree*C)) +
  xlab("Longitude (deg)") + ylab("Latitude (deg)") +
  xlim(c(-180, 180)) + ylim(c(-90, 90)) +
  map_layer + 
  geom_rect(aes(xmin = BM_box[1, "lon"], xmax = BM_box[2, "lon"],
                ymin = BM_box[1, "lat"], ymax = BM_box[2, "lat"]),
            colour = "red", alpha = 0) +
  theme_bw()  +
  coord_fixed(expand = FALSE, ratio = 1.3, ylim = c(-85, 90))

## Define a scale for the Brazil-Malvinas data
Brazil_Malv_colour_scale <- range(c(df_train$Z, df_test$Z))

## Reduce repetition
gg_basic <- ggplot() + 
  # scale_colour_distiller(palette = "Spectral", name = expression(degree*C),
  scale_colour_gradientn(colours = nasa_palette, name = expression(degree*C),
                         limits = Brazil_Malv_colour_scale) +
  xlab("Longitude (deg)") + ylab("Latitude (deg)") +
  map_layer + xlim(BM_box[, "lon"]) + ylim(BM_box[, "lat"]) + theme_bw() + 
  coord_fixed(expand = FALSE)

plot_of_SST_training_data <- gg_basic +
  geom_point(data = df_train, aes(lon, lat, colour = Z), pch = 46) 

plot_of_SST_testing_data <- gg_basic +
  geom_point(data = df_test, aes(lon, lat, colour = Z), pch = 46) 

ggpubr::ggarrange(plot_of_SST_data_world, 
                  plot_of_SST_training_data,
                  # plot_of_SST_testing_data, 
                  nrow = 1)
ggsave( 
  filename = "global_and_training_data2.png", device = "png", width = 8, height = 2.8,
  path = paste0(DIRECTORY, "img/")
)

## ---- Modelling functions ----

## Load the model fitting and prediction functions. 
## Each function fits the corresponding model using df_train, and predicts
## at locations pred.locs.
source(paste0(DIRECTORY, "modelling_functions/FRK.R"))
source(paste0(DIRECTORY, "modelling_functions/INLA.R"))
source(paste0(DIRECTORY, "modelling_functions/mgcv.R"))
source(paste0(DIRECTORY, "modelling_functions/LatticeKrig.R"))
source(paste0(DIRECTORY, "modelling_functions/gstat.R"))

## ---- Optimise arguments ----

source(paste0(DIRECTORY, "Diagnostic_fns.R"))

## See optimise_args.R

## ---- Run the models ----

## NB: Unfortunately, INLA cannot be broken down to a fit and predict function
## (see https://www.r-inla.org/faq#h.821k2r53fvx3). 
## This means we need to re-run the model fitting every time we want to predict
## over new locations.
## As a work-around, we could provide both sets of prediction locations at once;
## however, this would complicate things somewhat.



## Appends predictions and SE to pred.locs. 
## Requires the prediction function fn, prediction locations pred.locs, the 
## fitted object M (or training data in the case of INLA), and the method name.
pred_and_append_to_df <- function(fn, pred.locs, M, method_name) {
  tmp <- fn(pred.locs, M)
  pred.locs[, paste0("pred_", method_name)] <- tmp$pred
  pred.locs[, paste0("se_", method_name)] <- tmp$se
  return(pred.locs)
}

times <- list()
fitted_model_objects <- list()

times$FRK <- system.time({
  fitted_model_objects$FRK <- SST_FRK_fit(df_train, nres = 4)
  df_test <- pred_and_append_to_df(SST_FRK_pred, df_test, fitted_model_objects$FRK, "FRK")
})

times$INLA <- system.time(
  df_test <- pred_and_append_to_df(SST_INLA, df_test, df_train, "INLA")
)

times$mgcv <- system.time({
  fitted_model_objects$mgcv <- SST_mgcv_fit(df_train)
  df_test <- pred_and_append_to_df(SST_mgcv_pred, df_test, fitted_model_objects$mgcv, "mgcv")
})

times$LKrig <- system.time({
  fitted_model_objects$LKrig <- SST_LKrig_fit(df_train)
  df_test <- pred_and_append_to_df(SST_LKrig_pred, df_test, fitted_model_objects$LKrig, "LKrig")
})

times$gstat<- system.time(
  df_test <- pred_and_append_to_df(SST_gstat, df_test, df_train, "gstat")
)

df_test$se_gstat %>% is.na() %>% sum() # We will just omit these


#### MRA:
MRA_results <- read.csv(
  paste0(DIRECTORY, 
         # "results/MRA_res_v2/df_train_test_res_v2/dataset1_results_MRA.csv")
         "results/MRA_res_default_parameters/df_train_test_res/dataset1_results_MRA.csv")
  )
## Sanity check:
all(MRA_results$x == df_test$lon)
all(MRA_results$y == df_test$lat)
nrow(df_test) == nrow(MRA_results)

## There are some NAs and negative variances in MRA results:
MRA_results %>%
  summarise(prediction_NAs = sum(is.na(pred)),
            variance_NAs = sum(is.na(pred_var)),
            negative_variances = sum(pred_var < 0, na.rm = TRUE))

## This likely occured because the prediction locations coincided with knot 
## locations. The simplest solution is to jitter/offset these problematic 
## prediction locations, and then predict again. offsetting by 0.0002 degrees
## fixes all the problems. 
offset_predictions <- read.csv(
  paste0(DIRECTORY, 
         "results/MRA_res_v3/df_train_test_res_v3/dataset1_results_MRA_offset_0.0002.csv")
)

## Sanity check: Confirm that locations are the same
MRA_problematic_idx <- which(is.na(MRA_results$pred) | MRA_results$pred_var < 0)
MRA_results[MRA_problematic_idx, c("x", "y")] == (offset_predictions[, c("x", "y")] - 0.0002)
## We can confirm that the coordinates are the same, 
## with locations shifted up by 0.0002 (both longitude and latitude). 
## Now replace:
MRA_results[MRA_problematic_idx, ] <- offset_predictions

## Incorporate MRA into the study
df_test$pred_MRA <- MRA_results$pred
df_test$se_MRA <- sqrt(MRA_results$pred_var) 
MRA_time <-  read.csv(paste0(DIRECTORY, "results/MRA_res_default_parameters/df_train_test_res/MRA_parameters.csv"))
times$MRA <- c("elapsed" = MRA_time$totaltime) 
####

## save results for later use
write.csv(df_test, 
          file = paste0(DIRECTORY, "results/df_test.csv"),
          row.names = FALSE)
df_test <- read.csv(file = paste0(DIRECTORY, "results/df_test.csv"))


## save fitted model objects so we can predict over other locations if needed
save(fitted_model_objects, file = paste0(DIRECTORY, "fitted_model_objects.RData"))

## Save times: 
times <- sapply(times, function(x) unname(x["elapsed"]))
times <- times[PACKAGES] # Order times by PACKAGES just to be safe
times <- t(as.data.frame(times))
times <- times/60  # convert to minutes
write.csv(times,
          file = paste0(DIRECTORY, "results/times.csv"),
          row.names = FALSE)
times <- read.csv(file = paste0(DIRECTORY, "results/times.csv"))



## ---- Diagnostics (RMSPE, CRPS, IS90, COV90, time) ----

source(paste0(DIRECTORY, "Diagnostic_fns.R"))

df_test <- read.csv(paste0(DIRECTORY, "results/df_test.csv"))

## Function to create long form dataframe, useful for diagnostics and plotting
long_prediction_data_frame <- function(df) {
  data.frame(
    Method = rep(PACKAGES, each = nrow(df)),
    lon    = rep(df$lon, times = length(PACKAGES)), 
    lat    = rep(df$lat, times = length(PACKAGES)), 
    pred   = c(as.matrix(df[, paste0("pred_", PACKAGES)])),
    se     = c(as.matrix(df[, paste0("se_", PACKAGES)]))
  )
}

## Convert to long form:
df_test <- cbind(
  long_prediction_data_frame(df_test), 
  Z = df_test$Z
)


## Add time
tmp <- as.numeric(times)
names(tmp) <- names(times)
df_test$time <- as.numeric(tmp[as.character(df_test$Method)])

## Compute diagnostics
(
  diagnostics <- df_test %>% 
    group_by(Method) %>%
    drop_na() %>% # 11 of the gstat standard errors are NA
    summarise(
      RMSPE = RMSPE(Z, pred),
      COV90 = coverage90(Z, pred, se), 
      IS90 = IS90(Z, pred, se), 
      CRPS = verification::crps(Z, matrix(c(pred, se), ncol = 2))$CRPS, 
      time = time[1] 
    ) %>% 
    as.data.frame()
)

## Re-order the diagnostics to put gstat at the end:
idx <- c(which(diagnostics$Method != "gstat"), which(diagnostics$Method == "gstat"))
diagnostics <- diagnostics[idx, ]

write.csv(diagnostics, 
          file = paste0(DIRECTORY, "results/Diagnostics.csv"), 
          row.names = FALSE)
diagnostics <- read.csv(paste0(DIRECTORY, "results/Diagnostics.csv"))

xtable::xtable(
  diagnostics, 
  caption = "Root mean squared prediction error (RMSPE), empirical coverage 
  (COV90) and Interval Score (IS90) from a purported 90% prediction interval, 
  continuous ranked probability score (CRPS), and total run time (Time) 
  of predictions of Brazil-Malvinas confluence zone validation data using the 
  packages FRK, INLA, LatticeKrig, mgcv, and the MRA implemented with ...."
) %>%
  print(include.rownames = FALSE)



## ---- Generate predictions throughout the spatial domain D ----

## Run prediction scripts using fine-regular grid of pred locations
## throughout the spatial domain D.
grid_over_D <- expand.grid(
  lon = seq(-60, -48, length.out = 200), 
  lat = seq(-50, -35, length.out = 200)
)

grid_over_D <- pred_and_append_to_df(SST_FRK_pred, grid_over_D, fitted_model_objects$FRK, "FRK")
grid_over_D <- pred_and_append_to_df(SST_INLA, grid_over_D, df_train, "INLA")
grid_over_D <- pred_and_append_to_df(SST_mgcv_pred, grid_over_D, fitted_model_objects$mgcv, "mgcv")
grid_over_D <- pred_and_append_to_df(SST_LKrig_pred, grid_over_D, fitted_model_objects$LKrig, "LKrig")
system.time(grid_over_D <- pred_and_append_to_df(SST_gstat, grid_over_D, df_train, "gstat"))

MRA_grid_over_D <- read.csv(
  paste0(DIRECTORY, 
         "results/MRA_res_default_parameters/grid_over_D_res/dataset2_results_MRA.csv")
)

## NB: There are many NAs here:
sum(is.na(MRA_grid_over_D$pred))
MRA_grid_over_D %>% subset(is.na(pred))
## Weird; it is the first 201 locations, and then every 200 after that.
## It doesn't really matter for plotting, so we will just continue. 


grid_over_D$pred_MRA <- MRA_grid_over_D$pred
grid_over_D$se_MRA <- sqrt(MRA_grid_over_D$pred_var)

## save results for later use
write.csv(grid_over_D, 
          file = paste0(DIRECTORY, "results/grid_over_D.csv"), 
          row.names = FALSE)
grid_over_D <- read.csv(paste0(DIRECTORY, "results/grid_over_D.csv"))



# ---- Plot predictions and standard errors over D ----

## Convert to long form for plotting
df_long <- long_prediction_data_frame(grid_over_D)

## Re-order the levels of Method so that gstat appears last, as it does in the text
gstat_idx <- which(levels(df_long$Method) == "gstat")
df_long$Method <- factor(
  df_long$Method,
  levels = c(levels(df_long$Method)[-gstat_idx], "gstat")
)

## Replace LKrig with LatticeKrig (for plots)
levels(df_long$Method)[levels(df_long$Method)=="LKrig"] <- "LatticeKrig"

## Plot prediction
tmp1 <- ggplot(df_long) +
  # geom_tile(aes(x = lon, y = lat, fill = pmin(pmax(pred, -8), 8))) +
  geom_tile(aes(x = lon, y = lat, fill = pred)) +
  facet_wrap(vars(Method), nrow = 1) + 
  # scale_fill_distiller(palette = "Spectral") +
  scale_fill_gradientn(colours = nasa_palette) +
  labs(
    # fill = expression(widehat(p)[Y]["|"][bold(Z)]), 
    fill = "pred.", 
    x = "Longitude (deg)", y = "Latitude (deg)"
  ) +
  map_layer +
  theme_bw() +
  coord_fixed(expand = FALSE, xlim = c(-60, -48), ylim = c(-50, -35)) + 
  scale_x_continuous(breaks = c(-58, -54, -50))

## Plot standard error
tmp2 <- ggplot(df_long) +
  geom_raster(aes(x = lon, y = lat, fill = se)) +
  facet_wrap(vars(Method), nrow = 1) + 
  scale_fill_distiller(palette = "BrBG", direction = -1) +
  labs( fill = "s.e.", x = "Longitude (deg)", y = "Latitude (deg)") +
  map_layer +
  theme_bw() +
  coord_fixed(expand = FALSE, xlim = c(-60, -48), ylim = c(-50, -35)) + 
  scale_x_continuous(breaks = c(-58, -54, -50))

ggpubr::ggarrange(tmp1, tmp2, nrow = 2)

ggsave( 
  filename = "DEM_results.png", device = "png", width = 10, height = 5,
  path = paste0(DIRECTORY, "img/")
)


# ---- Coverage plot with/without fine-scale variance for varied number of basis functions ----

## NB: This section was interesting before I realised that I had forgotten to 
## include the measurement error variance when constructing the prediction 
## intervals! Now, it is not interesting, because both models achieve good 
## coverage irrespective of the number of basis functions (when including the 
## measurement error variance). Keep it in case it is useful in the future. 

# ## Here, we compare empirical coverage using FRK and LatticeKrig using different
# ## numbers of basis functions. The motivation is to assess the importance of 
# ## inclusion of the fine-scale variance term. Using FRK and LatticeKrig is good for 
# ## this task, as they implement similar models, but FRK includes a fine-scale 
# ## term while LatticeKrig does not. 
# source(paste0(DIRECTORY, "argument_selection/optimise_args_fns.R"))
# tmp <- test_arguments(
#   SST_FRK_fit, SST_FRK_pred, 
#   df_train, df_test, 
#   arguments = list(nres = 1:3)
# )
# FRK_coverage_study <- if(exists("FRK_coverage_study")) plyr::rbind.fill(FRK_coverage_study, tmp) else tmp
# FRK_coverage_study$Method <- "FRK"
# 
# ## Need to show coverage of LatticeKrig eventually becomes better
# ## Keep NC = 5 fixed (this makes the number of basis functions used for LKrig
# ## and FRK slightly more similar, so we can more directly compare the effect
# ## of including a fine-scale variation term)
# tmp <- test_arguments(
#   SST_LKrig_fit, SST_LKrig_pred, 
#   df_train, df_test, 
#   arguments = list(nlevel = 1:3, NC = 5)
# )
# LKrig_coverage_study <- if(exists("LKrig_coverage_study")) plyr::rbind.fill(LKrig_coverage_study, tmp) else tmp
# LKrig_coverage_study$Method <- "LatticeKrig"
# 
# coverage_study <- plyr::rbind.fill(FRK_coverage_study, LKrig_coverage_study)
# write.csv(coverage_study, 
#           file = paste0(DIRECTORY, "results/coverage_study.csv"), 
#           row.names = FALSE)
# coverage_study <- read.csv(paste0(DIRECTORY, "results/coverage_study.csv"))
# coverage_study$Method <- as.character(coverage_study$Method)
# 
# ## Plot
# coverage_study %>%
#   ggplot(., 
#          aes(x = NUM_BASIS_FUNCTIONS, y = COV90, group = Method, colour = Method)) + 
#   geom_line() + 
#   geom_point() + 
#   labs(x = "r", y = "Mean coverage from 90% prediction interval") + 
#   theme_bw()
# 
# ggsave( 
#   filename = "coverage_comparison.png", device = "png", width = 4.5, height =2.25,
#   path = paste0(DIRECTORY, "img/")
# )

