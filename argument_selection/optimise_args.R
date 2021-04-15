## testarguments is a package for testing and visualising the performance of a 
## prediction algorithm with different combinations of function arguments. 
## It is particularly useful if one suspects an interaction between argument 
## levels.
devtools::install_github("MattSainsbury-Dale/testarguments")
library("testarguments")

source(paste0(DIRECTORY, "Diagnostic_fns.R"))

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


# ---- FRK ----

fun <- function(df_train, df_test, nres, nBAUs) {
  M <- SST_FRK_fit(df_train, nres = nres, nBAUs = nBAUs)
  SST_FRK_pred(df_test, M)
}

FRK_scores <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_Gaussian, 
  arguments = list(nres = 3:4, nBAUs = seq(20000, 60000, by = 20000))
)

plot_diagnostics(FRK_scores, arg_names = c("nBAUs", "nres"))

# ---- INLA ----

## Wrapper to fit and predict in the same function, as required by test_arguments()
fun <- function(df_train, df_test, max.edge.interior) {
  SST_INLA(df_train = df_train, pred_locs = df_test, max.edge.interior = max.edge.interior)
}

## NB: inla() crashes when I set sigma0 = 3 * sd(df_fit$Z)
INLA_scores <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_Gaussian, 
  ## Run 1:
  # arguments = list(max.edge.interior = c(0.5, 1, 2, 3),
  #                  sigma0 = c(0.5, 1, 2) * sd(df_fit$Z),
  #                  range0 = c(0.5, 1, 2))
  ## Run 2:
  # arguments = list(max.edge.interior = c(0.3, 0.5), 
  #                  range0 = c(0.25, 0.5))
  ## Run 3:
  arguments = list(max.edge.interior = c(0.2, 0.3, 0.5))
)

plot_diagnostics(
  INLA_scores, 
  # arg_names = c("max.edge.interior", "sigma0", "range0")   # Run 1
  # arg_names = c("max.edge.interior", "range0") # Run 2
  arg_names = c("max.edge.interior") # Run 3
)

# ---- LatticeKrig ----

## Wrapper to fit and predict in the same function, as required by test_arguments()
fun <- function(df_train, df_test, nlevel, NC) {
  M <- SST_LKrig_fit(df_train, nlevel = nlevel, NC = NC)
  SST_LKrig_pred(df_test, M)
}

LKrig_scores <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_Gaussian, 
  arguments = list(nlevel = 2:4, NC = seq(5, 25, by = 5))
)

plot_diagnostics(LKrig_scores, arg_names = c("NC", "nlevel")) 


# ---- mgcv ----

fun <- function(df_train, df_test, k) {
  M <- SST_mgcv_fit(df_train, k = k)
  SST_mgcv_pred(df_test, M)
}

mgcv_scores <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_Gaussian, 
  arguments = list(k = c(seq(100, 1000, by = 100),
                         seq(1200, 2000, by = 200),
                         seq(2250, 3000, by = 250)))
)

plot_diagnostics(mgcv_scores, arg_names = "k") + scale_x_continuous(n.breaks = 5)