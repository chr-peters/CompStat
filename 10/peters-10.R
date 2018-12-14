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
# theta_1 = 0.0244228409444103, theta_2 = 1.66087324819251

# No. 2)
# ======

# mySimplex - Simplex-Verfahren
# Eingabe:
#    f           - Funktion, deren Optimum(Minimum) gefunden werden soll, soll
#                  aus dem R^K auf R abbilden
#    theta.start - numeric(k), Position des Start-Simplex
#    t           - numerisch, Groesse des Simplex
#    maxit       - maximale Anzahl an Iterationen (zweites Abbruchkriterium),
#                  natuerliche Zahl, default 1000
#
# Ausgabe:
#  Liste mit 2 benannten Elementen:
#    pars - numerische Matix, ausgewertete Parameter, einer je Zeile
#    vals - Funktionswerte an den Stelle 'par'

mySimplex <- function(f, theta.start, t = 1, maxit = 1000) {
  K <- length(theta.start)
  # Als erstes brauchen wir den Startsimplex:
  d1 <- t / K / sqrt(2) * (sqrt(K + 1) + K - 1)
  d2 <- t / K / sqrt(2) * (sqrt(K + 1) - 1)
  theta <- matrix(d2, ncol = K, nrow = K)
  diag(theta) <- d1
  # Achtung: Das + theta.start funktioniert nur, weil in jeder Spalte von
  # theta ein Punkt steht, und Matrizen spaltenweise aufgebaut sind.
  theta <- cbind(0, theta) + theta.start
  
  # Jetzt brauchen wir fuer jeden Punkt in theta den Funktionswert
  theta.vals <- apply(theta, 2, f)
  
  for (i in seq_len(maxit)) {
    # Bestimme neues Theta durch Spiegelung
    h <- which.max(theta.vals)
    theta.centroid <- rowMeans(theta[, -h])
    theta.new <- 2 * theta.centroid - theta[, h]
    # Simplex updaten
    theta[, h] <- theta.new
    theta.vals[h] <- f(theta.new)
  }
  best = which.min(theta.vals)
  return(list(pars = theta[, best], vals = theta.vals[best]))
}

print(mySimplex(function(x) sum(x^2), rep(10, 2)))

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

#' The test function for multiple starting values
testfun <- function(x) {
  sum(x**2 - 200 * cos(x))
}

# Plot the function for n=1
x <- seq(-50, 50, 0.01)
plot(x, sapply(x, testfun), type='l',
     main='Function to test the multiStarts optimization procedure')

# As we can see, this function lends itself well to test the multiStarts function
# because it has several local minima and only one global minimum.
# Without prior knowledge about the rough location of the global minimum, the only
# way to find it is to try multiple starting values and see which one produces
# the best results, which is exactly what the multiStarts functions does.
# We would expect that if the multiStarts function works properly, increasing
# the 'times' parameter will help finding the global minimum with a higher
# probability as we will see in the following simulation.

#' Simulation function to find an amount of starting points which reaches
#' a success rate of at least 95% in finding the global minimum of the
#' test function.
#' 
#' @param iterations Number of iterations to determine the success rate.
#' 
#' @return The number of starting points generated to achieve a >=95% success rate.
simulation <- function (iterations) {
  # count how often the global minimum was found
  successCount <- 0
  curIteration <- 1
  # current times value
  times <- 1
  
  repeat {
    # get the result for the current times value
    res <- multiStarts(testfun, times, c(-50, -50), c(50, 50))
    
    # check if the global minimum was found
    if (isTRUE(all.equal(res$par, c(0, 0), tolerance=1e-5))) {
      successCount <- successCount + 1
    }
    if (curIteration >= iterations) {
      # if the success rate is below 95%, increase times and start over
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

#print(simulation(100))

# So according to the simulation, roughly 120 starting points have to be
# generated in order to find the global optimum at least 95% of the time.