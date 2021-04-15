## NB: Unfortunately, INLA cannot be broken down to a fit and predict function
## (see https://www.r-inla.org/faq#h.821k2r53fvx3). 
## This means we need to re-run model fitting every time we want to predict
## over new locations. We could provide all sets of prediction locations at once;
## however, this would complicate things somewhat, and only one set of predictions
## are included in the run-time, so it is fine as is.
SST_INLA <- function(pred_locs, df_train, max.edge.interior = 0.3, 
                     sigma0 = 1 * sd(df_train$Z), range0 = 0.25) {
  
  ## max.edge.interior = 0.2 results in 25762 basis functions.
  ## max.edge.interior = 0.3 results in 11534 basis functions.
  ## max.edge.interior = 0.5 results in 4267 basis functions.
  
  ## Establish a boundary for the domain, D
  coords   <- as.matrix(df_train[, c("lon", "lat")])
  boundary <- inla.nonconvex.hull(as.matrix(coords))
  
  ## Triangulation of the domain
  ## NB: Here I set the max edge length of the exterior of the domain to be 50%
  ## larger than the max interior length, and the cutoff to be 10% of the 
  ## maximum length. Reducing the number of arguments to check over in this way
  ## makes it easier to select a good configuration. 
  max.edge.exterior <- 1.5 * max.edge.interior
  cutoff <- 0.1 * max.edge.interior
  mesh <- inla.mesh.2d(boundary = boundary,
                       max.edge = c(max.edge.interior, max.edge.exterior), 
                       cutoff = cutoff)
  
  cat("INLA using a triangular mesh with", mesh$n, "vertices (and hence basis functions)\n")
  
  ## Construct the SPDE on the mesh
  ## NB: I keep the probability associated with each prior (i.e., Psigma and 
  # Prange) fixed, and just optimise range0 and sigma0.
  spde <- inla.spde2.pcmatern(mesh = mesh, 
                              alpha = 2,
                              prior.range = c(range0, 0.01),
                              prior.sigma = c(sigma0, 0.01))
  
  ## Make an index which identifies spatial location 
  n_spatial <- mesh$n
  s_index   <- inla.spde.make.index(name = "spatial.field",
                                    n.spde = n_spatial)
  
  ## Fitting locations
  coords.fit <- df_train[ , c("lon", "lat")] %>% as.matrix()
  PHI <- inla.spde.make.A(mesh = mesh,
                          loc = coords.fit)
  
  ## First stack: Estimation
  n_data <- nrow(df_train)
  stack_est <- inla.stack(
    data    = list(z = df_train$Z),
    A       = list(PHI, 1),
    effects = list(s_index,
                   list(Intercept = rep(1, n_data))),
    tag = "est")
  
  df_pred <- data.frame(x = mesh$loc[,1],
                        y = mesh$loc[,2])
  n_pred <- nrow(df_pred)
  PHI_pred <- Diagonal(n = n_pred)
  
  ## Second stack: Prediction
  stack_pred <- inla.stack(
    data = list(z = NA), # NA means we want to predict it
    A = list(PHI_pred, 1),
    effects = list(s_index,
                   list(Intercept = rep(1, n_pred))),
    tag = "pred")
  
  stack <- inla.stack(stack_est, stack_pred)

  ## A matrix (need this for MC sampling)
  A_pred <- inla.stack.A(stack_pred) # can also be accessed with stack_pred$A
  
  ## Formula
  formula <- z ~ -1 + Intercept +
    f(spatial.field,
      model = spde,
      group = spatial.field.group)
  
  ## Fitting stage: almost all of the time is spent here
  inla.object <- inla(formula,
                     data = inla.stack.data(stack, spde = spde),
                     control.predictor = list(A = inla.stack.A(stack),
                                              compute = TRUE),
                     control.compute = list(config = TRUE)) # required for inla.posterior.sample()
  
  ## Prediction stage: this is comparatively fast
  ## Extract predictions and standard errors of basis function weights:
  index_pred <- inla.stack.index(stack, "pred")$data
  lp_mean    <- inla.object$summary.fitted.values$mean[index_pred]
  lp_sd      <- inla.object$summary.fitted.values$sd[index_pred]
  
  ## Locations we wish to predict over
  proj.grid <- inla.mesh.projector(
    mesh, 
    loc = as.matrix(pred_locs[, c("lon", "lat")])
    )

  ## Measurement error:
  ## Naive way (not correct):
  # 1 / inla.object$summary.hyperpar["Precision for the Gaussian observations", "mean"]
  
  ## The way recommended by https://stats.stackexchange.com/a/358897:
  m <- inla.object$internal.marginals.hyperpar[[1]]
  m.var <- inla.tmarginal(function(x) 1/exp(x), m)
  sigma2e <- inla.zmarginal(m.var)$mean
  
  ## Predictive variance: var(Z0 | Z) = var(Y0 | Z) + sigma2e
  var_Z0 <- (inla.mesh.project(proj.grid, lp_sd[1:mesh$n]))^2 + sigma2e
  
  return(data.frame(
    pred = inla.mesh.project(proj.grid, lp_mean[1:mesh$n]),
    se = sqrt(var_Z0)
  ))
}



