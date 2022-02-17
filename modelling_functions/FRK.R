SST_FRK_fit <- function(df_train, nres = 4, nBAUs = 60000) {
  
  ## Convert training data to SpatialPointsDataFrame
  zspdf <- df_train
  coordinates(zspdf) = ~ lon + lat   
  
  ## Construct Basis functions
  B <- auto_basis(data = zspdf, nres = nres)
  
  ## Assign NUM_BASIS_FUNCTIONS to the parent environment:
  NUM_BASIS_FUNCTIONS <- nbasis(B)
  assign("NUM_BASIS_FUNCTIONS", NUM_BASIS_FUNCTIONS, env = parent.frame()) 
  
  ## Construct a fine grid of BAUs (SpatialPixelsDataFrame)
  BAUs <- expand.grid(
    lon = seq(-60, -48, length.out = ceiling(sqrt(nBAUs))), 
    lat = seq(-50, -35, length.out = ceiling(sqrt(nBAUs)))
  ) 
  coordinates(BAUs) <- ~ lon + lat
  BAUs <- SpatialPixelsDataFrame(
    points = BAUs, 
    data = data.frame(fs = rep(1, length(BAUs)))
  )
  
  ## Construct and fit the SRE object
  M <- SRE(f = Z ~ 1, data = list(zspdf),  
           BAUs = BAUs, basis = B, K_type = "neighbour")
  
  M <- SRE.fit(M, method = "TMB")
  return(M)
  
  # ## Construct and fit the SRE object
  # ## NB: first need to fix FRK() not allowing users to set nres
  # return(FRK(f = Z ~ 1, data = list(zspdf), method = "TMB"))
}

SST_FRK_pred <- function(pred_locs, FRK_object) {
  
  ## Convert prediction locations to SpatialPointsDataFrame
  coordinates(pred_locs) = ~ lon + lat   
  pred <- predict(FRK_object, type = "link", 
                  newdata = pred_locs,
                  n_MC = 400, percentiles = NULL)
  
  ## Predictive variance: var(Z0 | Z) = var(Y0 | Z) + sigma2e
  var_Z0 <- (pred$newdata$RMSPE_Y)^2 + FRK_object@Ve[1, 1]
  
  return(data.frame(
    pred = pred$newdata$p_Y,
    se = sqrt(var_Z0)
  ))
}
