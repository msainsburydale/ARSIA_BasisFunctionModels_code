## Define the functions we will use for comparing across methods

## coverage from prediction interval with nominal coverage alpha
coverage90 <- function(Z, pred, se, alpha = 0.9) {
  z_alpha <- qnorm(1 - (1-alpha)/2)
  lower <- pred - z_alpha*se
  upper <- pred + z_alpha*se
  sum((Z < upper) & (Z > lower)) / length(Z)
}

## interval score from prediction interval with nominal coverage alpha
IS90 <- function(Z, pred, se, alpha = 0.9) {
  z_alpha <- qnorm(1 - (1-alpha)/2)
  pred_l <- pred - z_alpha*se
  pred_u <- pred + z_alpha*se
  ISs <- (pred_u - pred_l) + 2/(1-alpha) * (pred_l - Z) * (Z < pred_l) +
    2/(1-alpha) * (Z - pred_u) * (Z > pred_u)
  mean(ISs)
}

## Root-mean-squared prediction error
RMSPE <- function(z,pred) {
  Y <- (z - pred)^2
  sqrt(mean(Y))
}

## df should contain predictions in field named pred, prediction standard
## errors in se, and the validation data in Z. 
compute_diagnostics_Gaussian <- function(df) {
  
  if (any(is.na(df))) {
    warning("Removing rows with some NA entries.")
    # df <- df[rowSums(is.na(df)) == 0, ]
    df <- drop_na(df)
  }
  
  summarise(df,
            RMSPE = RMSPE(Z, pred),
            COV90 = coverage90(Z, pred, se), 
            IS90 = IS90(Z, pred, se), 
            CRPS = verification::crps(Z, matrix(c(pred, se), ncol = 2))$CRPS
  )
}
