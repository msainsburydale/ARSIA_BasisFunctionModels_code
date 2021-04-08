SST_mgcv_fit <- function(df_train, k = 2250) {

  f <- Z ~ te(lon, lat,            # inputs over which to smooth
              bs = "tp",           # type of bases
              k = k,             # knot count in each dimension
              d = 2)               # spatial basis dimension
  
  # return(gam(f, data = df_train))
  cat("mgcv: using bam()\n")
  return(bam(f, data = df_train))
}

SST_mgcv_pred <- function(pred_locs, mgcv_object) {
  
  tmp <- predict(mgcv_object, newdata = pred_locs, se.fit = TRUE)
  
  ## Predictive variance: var(Z0 | Z) = var(Y0 | Z) + sigma2e
  var_Z0 <- tmp$se.fit^2 + mgcv_object$sig2
  
  return(data.frame(
    pred = tmp$fit,
    se = sqrt(var_Z0)
  ))
}



