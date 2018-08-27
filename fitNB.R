#Fit NB function
fit_nb <- function(xs) {
  suppressWarnings(
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
  )
}

# Test
# dist <- rnbinom(100, size = 0.1, mu = 2)
# fit <- fit_nb(dist)
