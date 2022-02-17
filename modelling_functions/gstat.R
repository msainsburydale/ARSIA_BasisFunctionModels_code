SST_gstat <- function(pred_locs, df_train) {
  
  ## Convert training data and prediction locations to SpatialPointsDataFrame
  coordinates(df_train) = ~ lon + lat   
  coordinates(pred_locs) = ~ lon + lat  
  
  ## Fit variogram 
  varg <- variogram(object = Z ~ 1, data = df_train)
  vfit <- fit.variogram(
    object = varg, 
    model = vgm(10, model = "Mat", range = 3, nugget = 0.1), 
    fit.sills = TRUE, 
    fit.ranges = TRUE, 
    fit.kappa = TRUE
  ) 
  # plot(varg, vfit)
  
 ## Predict
  system.time(
    pred <- krige(Z ~ 1, 
                  locations = df_train, 
                  newdata = pred_locs, 
                  model = vfit, 
                  nmax = 100)
  )

  return(data.frame(
    pred = pred$var1.pred,
    se = sqrt(pred$var1.var)
  ))
}
