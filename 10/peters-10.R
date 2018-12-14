# Name: Christian Peters

# No. 1)
# ======

# load the data
load('nlkq_data.RData')

# solve the least squares problem
theta <- optim(c(0, 0), function(theta) {
  sum((theta[1] + theta[2]*X[,1] + theta[2]**2 * X[,2] - y)**2)
})$par

# print the results
print(paste0('theta_1 = ', theta[1], ', theta_2 = ', theta[2]))

# No. 3)
# ======

#' Extension of the BFGS optimization algorithm using multiple starting values.
#' 
#' @param fn    The function to optimize.
#' @param times The number of starting points.
#' @param lower The lower box constraint.
#' @param upper The upper box constraint.
#' 
#' @return The results of the best BFGS run.
multiStarts <- function(fn, times, lower, upper) {
  # check input parameters
  stopifnot(length(lower) == length(upper), times > 0)
  count <- 1
  repeat {
    # get a random starting point within the box constraints
    start <- runif(length(lower), min=lower, max=upper)
    # employ L-BFGS-B as the optimization method because it uses the box constraints
    current <- optim(start, fn, lower=lower, upper=upper, method="L-BFGS-B")
    if (!exists('best') || best$value > current$value) {
      best <- current
    }
    # return the best result
    if (count >= times) {
      return(best)
    }
    count <- count + 1
  }
}

testfun <- function(x) {
  sum(x**2 - 200 * cos(x))
}
x <- seq(-50, 50, 0.01)
plot(x, sapply(x, testfun), type='l')

simulation <- function (iterations) {
  successCount <- 0
  curIteration <- 1
  times <- 1
  repeat {
    res <- multiStarts(testfun, times, c(-50, -50), c(50, 50))
    if (isTRUE(all.equal(res$par, c(0, 0), tolerance=1e-5))) {
      successCount <- successCount + 1
    }
    if (curIteration >= iterations) {
      if (successCount / iterations < 0.95) {
        times <- times + 1
        curIteration <- 1
        successCount <- 0
        next
      } else {
        return(times)
      }
    }
    curIteration <- curIteration + 1
  }
}

print(simulation(100))