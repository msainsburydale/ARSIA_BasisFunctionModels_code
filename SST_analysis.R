devtools::install_github("andrewzm/FRK", "FRKTMB")
library("FRK")
library("INLA")
library("mgcv")
library("LatticeKrig")
library("gstat")
library("ggplot2")
library("dplyr")
library("verification")
library("tidyr")
library("stringr")
library("ggpubr") # ggarrange
library("xtable")


DIRECTORY <- "ENTER THE PATH TO THE DIRECTORY CONTAINING SST_analysis.R"
## e.g., 
## DIRECTORY <- "~/ARSIA_BasisFunctionModels_code/"

## Packages used in the study
PACKAGES <- c("FRK", "INLA", "mgcv", "LKrig", "MRA", "gstat")


# ---- Data pre-processing ----

## Load the data
load(paste0(DIRECTORY, "data/SST_sub_1000000.rda"))
load(paste0(DIRECTORY, "data/SST_sub_1000000_val.rda"))

## Rename the data
df_train <- SST_sub_1000000; rm(SST_sub_1000000)
df_test <- SST_sub_1000000_val; rm(SST_sub_1000000_val)

## Remove columns which will not be used
df_train$error <- df_test$error <- df_train$bias  <- df_test$bias  <- NULL

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
df_test$lat2  <- df_test$lat^2
df_train$Z    <- residuals(lm(sst ~ 1 + lat + lat2, data = df_train))
df_test$Z     <- residuals(lm(sst ~ 1 + lat + lat2, data = df_test))

## Keep a copy of the original training data for plotting over the world
df_train_full <- df_train

## Subset the region corresponding to the Brazil-Malvinas Confluence zone
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

## A layer which will be added to most of the plots throughout the script:
map_layer <- geom_map(
  data = map_data("world"), map = map_data("world"),
  aes(group = group, map_id = region),
  fill = "white", colour = "black", size = 0.1
) 

## A palette for the data and predictions
nasa_palette <- c(
  "#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7",
  "#0f81f3","#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e",
  "#aefc2a","#e9fc0d","#f6da0c","#f5a009","#f6780a","#f34a09","#f2210a",
  "#f50008","#d90009","#a80109","#730005"
)

## Figure 4: Left panel. Training data over the world map, with a small red box 
## indicating the confluence zone. 
plot_of_SST_data_world <- ggplot(df_train_full) +
  ## NB: restricting the colour scale to [-8, 8] prevents outliers dominating
  geom_point(aes(lon, lat, colour = pmin(pmax(Z, -8), 8)), pch = 46) +
  scale_colour_gradientn(colours = nasa_palette, name = expression(degree*C)) +
  xlab("Longitude (deg)") + ylab("Latitude (deg)") +
  xlim(c(-180, 180)) + ylim(c(-90, 90)) +
  map_layer + 
  geom_rect(aes(xmin = BM_box[1, "lon"], xmax = BM_box[2, "lon"],
                ymin = BM_box[1, "lat"], ymax = BM_box[2, "lat"]),
            colour = "red", alpha = 0) +
  theme_bw()  +
  coord_fixed(expand = FALSE, ratio = 1.3, ylim = c(-85, 90))

## Figure 4: Right panel. Training data in the confluence zone.
plot_of_SST_training_data <- ggplot() + 
  scale_colour_gradientn(colours = nasa_palette, name = expression(degree*C)) +
  xlab("Longitude (deg)") + ylab("Latitude (deg)") +
  map_layer + xlim(BM_box[, "lon"]) + ylim(BM_box[, "lat"]) + theme_bw() + 
  coord_fixed(expand = FALSE) +
  geom_point(data = df_train, aes(lon, lat, colour = Z), pch = 46) 

ggsave( 
  ggarrange(plot_of_SST_data_world, plot_of_SST_training_data, nrow = 1),
  filename = "global_and_training_data.pdf", device = "pdf", width = 8, height = 2.8,
  path = paste0(DIRECTORY, "img/")
)

# ---- Load model fitting and prediction functions ----

## Load the model fitting and prediction functions. 
## Each function fits the corresponding model using df_train, and predicts
## at locations pred.locs.
source(paste0(DIRECTORY, "modelling_functions/FRK.R"))
source(paste0(DIRECTORY, "modelling_functions/INLA.R"))
source(paste0(DIRECTORY, "modelling_functions/mgcv.R"))
source(paste0(DIRECTORY, "modelling_functions/LatticeKrig.R"))
source(paste0(DIRECTORY, "modelling_functions/gstat.R"))


# ---- Optimise arguments ----

## See optimise_args.R

# ---- Run the models ----

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
  fitted_model_objects$FRK <- SST_FRK_fit(df_train)
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

# ---- MRA ----

MRA_results <- read.csv(
  paste0(DIRECTORY, 
         "results/MRA/df_train_test_res/dataset1_results_MRA.csv")
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

## This likely occurred because the prediction locations coincided with knot 
## locations. The simplest solution is to jitter/offset these problematic 
## prediction locations, and then predict again. offsetting by 0.0002 degrees
## fixes all the problems. 
offset_predictions <- read.csv(
  paste0(DIRECTORY, 
         "results/MRA/df_train_test_res/dataset1_results_MRA_offset_0.0002.csv")
)

## Sanity check: Confirm that locations are the same
MRA_problematic_idx <- which(is.na(MRA_results$pred) | MRA_results$pred_var < 0)
max(MRA_results[MRA_problematic_idx, c("x", "y")] - (offset_predictions[, c("x", "y")] - 0.0002))
## We can confirm that the coordinates are the same (within machine precision), 
## with locations shifted up by 0.0002 (both longitude and latitude). 
## Now replace problematic locations:
MRA_results[MRA_problematic_idx, ] <- offset_predictions

## Add MRA to df_test and times
df_test$pred_MRA <- MRA_results$pred
df_test$se_MRA <- sqrt(MRA_results$pred_var) 
MRA_time <-  read.csv(paste0(DIRECTORY, "results/MRA/df_train_test_res/MRA_parameters.csv"))
times$MRA <- c("elapsed" = MRA_time$totaltime) 


# ---- Save the testing data with predictions and standard errors ----

## save results for later use
write.csv(df_test, 
          file = paste0(DIRECTORY, "results/df_test.csv"),
          row.names = FALSE)

## Save fitted model objects so we can predict over other locations if needed
# save(fitted_model_objects, file = paste0(DIRECTORY, "fitted_model_objects.RData"))

## Save times
times <- sapply(times, function(x) unname(x["elapsed"]))
times <- t(as.data.frame(times))
times <- times/60  # convert to minutes
write.csv(times, file = paste0(DIRECTORY, "results/times.csv"), row.names = FALSE)

# ---- Reload the saved data (useful if stepping through the script without running the models) ----

times   <- read.csv(file = paste0(DIRECTORY, "results/times.csv"))
df_test <- read.csv(file = paste0(DIRECTORY, "results/df_test.csv"))

# ---- Diagnostics (RMSPE, CRPS, IS90, COV90, time) ----

source(paste0(DIRECTORY, "Diagnostic_fns.R"))

## Function to create long form dataframe, useful for diagnostics
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
diagnostics
write.csv(diagnostics, 
          file = paste0(DIRECTORY, "results/Diagnostics.csv"), 
          row.names = FALSE)
diagnostics <- read.csv(paste0(DIRECTORY, "results/Diagnostics.csv"))

print(xtable(diagnostics), include.rownames = FALSE)



# ---- Generate predictions throughout the spatial domain D ----

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
grid_over_D <- pred_and_append_to_df(SST_gstat, grid_over_D, df_train, "gstat")

MRA_grid_over_D <- read.csv(
  paste0(DIRECTORY, 
         "results/MRA/grid_over_D_res/dataset2_results_MRA.csv")
)
grid_over_D$pred_MRA <- MRA_grid_over_D$pred
grid_over_D$se_MRA <- sqrt(MRA_grid_over_D$pred_var)

## save results for later use
write.csv(grid_over_D, 
          file = paste0(DIRECTORY, "results/grid_over_D.csv"), 
          row.names = FALSE)

# ---- Load predictions over D (useful if stepping through the script without running the models) ----

grid_over_D <- read.csv(paste0(DIRECTORY, "results/grid_over_D.csv"))


# ---- Plot predictions and standard errors over D ----

## Convert to long form for plotting
df_long <- long_prediction_data_frame(grid_over_D)

## Re-order levels of Method so that gstat appears last, as it does in the text
gstat_idx <- which(levels(df_long$Method) == "gstat")
df_long$Method <- factor(
  df_long$Method,
  levels = c(levels(df_long$Method)[-gstat_idx], "gstat")
)

## Replace LKrig with LatticeKrig (for plots)
levels(df_long$Method)[levels(df_long$Method)=="LKrig"] <- "LatticeKrig"

## Plot prediction
plot_predictions <- ggplot(df_long) +
  geom_tile(aes(x = lon, y = lat, fill = pred)) +
  facet_wrap(vars(Method), nrow = 1) + 
  scale_fill_gradientn(colours = nasa_palette) +
  labs(
    fill = "pred.", 
    x = "Longitude (deg)", y = "Latitude (deg)"
  ) +
  map_layer +
  theme_bw() +
  coord_fixed(expand = FALSE, xlim = c(-60, -48), ylim = c(-50, -35)) + 
  scale_x_continuous(breaks = c(-58, -54, -50))

## Plot standard error
plot_SE <- ggplot(df_long) +
  geom_raster(aes(x = lon, y = lat, fill = se)) +
  facet_wrap(vars(Method), nrow = 1) + 
  scale_fill_distiller(palette = "BrBG", direction = -1) +
  labs( fill = "s.e.", x = "Longitude (deg)", y = "Latitude (deg)") +
  map_layer +
  theme_bw() +
  coord_fixed(expand = FALSE, xlim = c(-60, -48), ylim = c(-50, -35)) + 
  scale_x_continuous(breaks = c(-58, -54, -50))

ggsave( 
  ggarrange(plot_predictions, plot_SE, nrow = 2),
  filename = "DEM_results.pdf", device = "pdf", width = 10, height = 5,
  path = paste0(DIRECTORY, "img/")
)


