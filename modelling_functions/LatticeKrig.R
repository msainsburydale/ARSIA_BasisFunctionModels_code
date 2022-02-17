SST_LKrig_fit <- function(df_train, nlevel = 4, 
                          NC = 20, 
                          NC.buffer = NULL, a.wght = 7) {

  domain <- rbind(
    floor(c(min(df_train$lon), min(df_train$lat))), 
    ceiling(c(max(df_train$lon), max(df_train$lat)))
    )
  
  LKinfo <- LKrigSetup(
    domain, 
    nlevel = nlevel, 
    NC = NC,
    NC.buffer = NC.buffer, a.wght = a.wght,
    nu = 1.5
    )
  
  ## Assign NUM_BASIS_FUNCTIONS to the parent environment:
  NUM_BASIS_FUNCTIONS <- LKinfo$latticeInfo$m
  assign("NUM_BASIS_FUNCTIONS", NUM_BASIS_FUNCTIONS, env = parent.frame()) 
  
  LKrig_object <- LatticeKrig(
    dplyr::select(df_train, lon, lat), 
    df_train$Z, 
    LKinfo = LKinfo, 
    findAwght = FALSE
    )
  
  return(LKrig_object)
}

SST_LKrig_pred <- function(pred_locs, LKrig_object) {
   xnew <- dplyr::select(pred_locs, lon, lat)
   
   ## Predictive variance: var(Z0 | Z) = var(Y0 | Z) + sigma2e
   var_Z0 <- (as.vector(predictSE.LKrig(LKrig_object, xnew)))^2 + LKrig_object$sigma.MLE^2

   return(data.frame(
     pred = as.vector(predict.LKrig(LKrig_object, xnew)),
     se = sqrt(var_Z0)
   ))
}

