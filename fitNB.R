#Fit NB function
fit_nb <- function(xs) {
  if(sum(xs > 0)) {
    fit <- MASS::fitdistr(xs, "negative binomial")
    fit$estimate
  } else {
    c(0, 0)
  }
}
