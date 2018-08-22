#Fit NB function
fit_nb <- function(xs) {
  if(sum(xs > 0)) {
    fit <- try(MASS::fitdistr(xs, "negative binomial"))
    if(class(fit) == "fitdistr") {
      fit$estimate
    } else {
      c(0, 0)
    }
  } else {
    c(0, 0)
  }
}
