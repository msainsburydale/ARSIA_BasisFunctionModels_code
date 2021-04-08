library("testarguments") # devtools::install_github("MattSainsbury-Dale/testarguments")

ARG_SELECTION_DIRECTORY <- paste0(DIRECTORY, "argument_selection/")

## We select arguments by splitting the training set into two; a training set, 
## and a validation set (i.e., we do NOT use the test set to choose function 
## arguments).
idx_train <- 1:nrow(df_train)
RNGversion("3.6.0")
set.seed(1)
idx_fit <- sample(idx_train, round(nrow(df_train)/2), replace = FALSE)
idx_val <- idx_train[-idx_fit]

df_fit <- df_train[idx_fit, ]
df_val <- df_train[idx_val, ]


# ---- LatticeKrig ----

## Wrapper to fit and predict in the same function, as required by test_arguments()
fun <- function(df_train, df_test, nlevel, NC) {
  M <- SST_LKrig_fit(df_train, nlevel = nlevel, NC = NC)
  SST_LKrig_pred(df_test, M)
}

tmp <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_Gaussian, 
  arguments = list(nlevel = 3:4, NC = 25)
)

LKrig_scores <- if(exists("LKrig_scores")) plyr::rbind.fill(LKrig_scores, tmp) else tmp

write.csv(LKrig_scores, 
          file = paste0(ARG_SELECTION_DIRECTORY, "LKrig_scores.csv"),
          row.names = FALSE)
LKrig_scores <- read.csv(file = paste0(ARG_SELECTION_DIRECTORY, "LKrig_scores.csv"))

plot_diagnostics(LKrig_scores, argument_names = c("NC", "nlevel")) 
## These are the two main arguments to optimize. 
## The curves suggest the diagnostics are still improving rather significantly, 
## so I will increase NC and nlevel further. 
plot_diagnostics(LKrig_scores, argument_names = c("NC")) 

ggsave( 
  filename = "LKrig_args.png", device = "png", width = 12, height = 6,
  path = ARG_SELECTION_DIRECTORY
)


# ---- mgcv ----

fun <- function(df_train, df_test, k) {
  M <- SST_mgcv_fit(df_train, k = k)
  SST_mgcv_pred(df_test, M)
}

tmp <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_Gaussian, 
  arguments = list(k = c(seq(100, 1000, by = 100),
                         seq(1200, 2000, by = 200),
                         seq(2250, 3000, by = 250)))
)

mgcv_scores <- if(exists("mgcv_scores")) plyr::rbind.fill(mgcv_scores, tmp) else tmp

write.csv(mgcv_scores, 
          file = paste0(ARG_SELECTION_DIRECTORY, "mgcv_scores.csv"),
          row.names = FALSE)
mgcv_scores <- read.csv(file = paste0(ARG_SELECTION_DIRECTORY, "mgcv_scores.csv"))

plot_diagnostics(mgcv_scores, argument_names = "k") + 
  scale_x_continuous(n.breaks = 5)

ggsave( 
  filename = "mgcv_args.png", device = "png", width = 12, height = 6,
  path = ARG_SELECTION_DIRECTORY
)


# ---- FRK ----

## Wrapper to fit and predict in the same function, as required by test_arguments()
fun <- function(df_train, df_test, nres, nBAUs) {
  M <- SST_FRK_fit(df_train, nres = nres, nBAUs = nBAUs)
  SST_FRK_pred(df_test, M)
}

tmp <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_Gaussian, 
  # arguments = list(nres = 4, nBAUs = seq(20000, 60000, by = 20000))
  arguments = list(nres = 1:2, nBAUs = 20000)
)

FRK_scores <- if(exists("FRK_scores")) plyr::rbind.fill(FRK_scores, tmp) else tmp

write.csv(FRK_scores, 
          file = paste0(ARG_SELECTION_DIRECTORY, "FRK_scores.csv"), 
          row.names = FALSE)
FRK_scores <- read.csv(file = paste0(ARG_SELECTION_DIRECTORY, "FRK_scores.csv"))

plot_diagnostics(FRK_scores, arg_names = c("nBAUs", "nres"))
ggsave( 
  filename = "FRK_args.png", device = "png", width = 12, height = 6,
  path = ARG_SELECTION_DIRECTORY
  )

plot_diagnostics(long_df, argument_names = c("nres"))

# ---- INLA ----

## Wrapper to fit and predict in the same function, as required by test_arguments()
fun <- function(df_train, df_test, max.edge.interior) {
  SST_INLA(df_train = df_train, pred_locs = df_test, max.edge.interior = max.edge.interior)
}

## NB: inla() crashes when I set sigma0 = 3 * sd(df_fit$Z)
tmp <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_Gaussian, 
  arguments = list(max.edge.interior = c(20))
  ## Run 1:
  # arguments = list(max.edge.interior = c(0.5, 1, 2, 3), 
  #                  sigma0 = c(0.5, 1, 2) * sd(df_fit$Z), 
  #                  range0 = c(0.5, 1, 2))
  ## Run 2:
  # arguments = list(max.edge.interior = c(0.3, 0.5), 
  #                  range0 = c(0.25, 0.5))
  ## Run 3:
  # arguments = list(max.edge.interior = c(0.2, 0.3, 0.5))
)

INLA_scores <- tmp
# INLA_scores <- if(exists("INLA_scores")) plyr::rbind.fill(INLA_scores, tmp) else tmp

write.csv(INLA_scores, 
          file = paste0(ARG_SELECTION_DIRECTORY, "INLA_scores_run2.csv"),
          row.names = FALSE)

# INLA_scores <- read.csv(file = paste0(ARG_SELECTION_DIRECTORY, "INLA_scores_run2.csv"))

plot_diagnostics(
  INLA_scores, 
  # argument_names = c("max.edge.interior", "sigma0", "range0")   # Run 1
  # argument_names = c("max.edge.interior", "range0") # Run 2
  argument_names = c("max.edge.interior") # Run 3
)
## Run 1) As expected, lower max.edge.interior seems to be the best.
## Run 2) Still a reasonable improvement by lowering max.edge.interior.
## Run 3) The time is significantly higher setting maz.edge.interior = 0.2;
#         there are still small improvements, 
ggsave( 
  filename = "INLA_args_maxedgelength.png", device = "png", width = 12, height = 6,
  path = ARG_SELECTION_DIRECTORY
)


plot_diagnostics(long_df, argument_names = c("sigma0")) 
# Run 1) Seems like sigma0 = 1 * sd(df_fit$Z) is the best. 
## Further, there is not much difference between sigma0 = 1 * sd(df_fit$Z) and
## sigma0 = 0.5 * sd(df_fit$Z). No need to explore further.

ggsave( 
  filename = "INLA_args_sigma0.png", device = "png", width = 12, height = 6,
  path = ARG_SELECTION_DIRECTORY
)

plot_diagnostics(long_df, argument_names = c("range0"))
## Run 1) range0 = 0.5 seems to be the best
## Run 2) no difference between range0 = 0.25 and range0 = 0.5 in terms of diagnostics;
## range0 = 0.25 is slightly faster, so set range0 = 0.25. I believe this is a 
## slightly less informative prior as well, which is good.

ggsave( 
  filename = "INLA_args_range0.png", device = "png", width = 12, height = 6,
  path = ARG_SELECTION_DIRECTORY
)