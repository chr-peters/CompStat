# Name: Christian Peters

library(numDeriv)
library(mvtnorm)
library(ggplot2)

# No. 2)
# ======

#' Implementation of a quasi-newton optimization algorithm using the DFP update.
#' 
#' @param f       The function to optimize.
#' @param start   The starting point.
#' @param tol     The minimum step size. Used as a stopping criterion.
#' @param maxIter Upper threshold for the maximum number of iterations.
#' 
#' @return A list with two entries: The result and the number of iterations.
dfpOptim <- function(f, start, tol=1e-4, maxIter=50) {
  # approximation of the inverse hessian matrix that is updated according
  # to the DFP update rule
  invHessian <- diag(length(start))
  last <- start
  curIter <- 1
  curGradient <- grad(f, start)
  
  repeat {
    # find the best stepsize
    ny <- optimize(function(x) {
        f(last - x * (invHessian %*% curGradient)[,1])
        # don't take steps bigger than a newton step
      }, interval=c(-1, 1))$minimum
    
    # take a step
    deltaTheta <- -ny * (invHessian %*% curGradient)[,1]
    new <- last + deltaTheta
    
    # check stopping conditions
    if (norm(new - last, type='2') < tol || curIter >= maxIter) {
      return(list(result=new, iterations=curIter))
    }
    
    # approximate the next gradient
    nextGradient <- grad(f, new)
    
    # update the inverse hessian approximation using the DFP update
    g <- nextGradient - curGradient
    
    dfpUpdate <- (deltaTheta %*% t(deltaTheta)) / (t(g) %*% deltaTheta)[1, 1] - 
      (invHessian %*% g %*% t(g) %*% invHessian) / (t(g) %*% invHessian %*% g)[1, 1]
    
    invHessian <- invHessian + dfpUpdate
    
    curGradient <- nextGradient
    last <- new
    curIter <- curIter + 1
  }
}

#' This function tests the DFP optimization algorithm by minimizing several
#' bivariate normal distributions with different correlation parameters and
#' examining the number of iterations that it takes to find the minimum.
testDfpOptim <- function() {
  # the different rho values to test
  rho <- seq(-0.95, 0.95, 0.05)
  # get the number of iterations it takes for each rho
  iterations <- sapply(rho, function(rho){
    sigma <- matrix(c(1, rho, rho, 1), nrow=2)
    mean <- c(0, 0)
    # start at a suitable point so that the algorithm converges
    start <- c(1, 1)
    res <- dfpOptim(function(x) -dmvnorm(x, mean=mean, sigma=sigma), start=start)
    if(!isTRUE(all.equal(mean, res$result, tolerance=1e-4))) {
      message(paste0('Error at rho = ', rho, ': The optimizer did not converge! ',
                     'Stopped at: [', res$result[1], ', ', res$result[2], ']'))
    }
    return(res$iterations)
  })
  print(
    ggplot(data.frame(rho=rho, iterations=iterations), aes(x=rho, y=iterations)) +
      ggtitle('Number of iterations of the DFP updating algorithm for different
              bivariate normal distributions') +
      geom_point() +
      scale_y_continuous(breaks=1:max(iterations)) +
      theme(plot.title = element_text(hjust = 0.5))
  )
}

testDfpOptim()

# No. 3)
# ======

singleSimulation <- function(method, b, n, k, eps, u) {
  # create the matrix X from random observations
  X <- matrix(rnorm(b*n), nrow=b, ncol=n)
  X[, n] <- 1
  y <- sample(c(1, 0), b, replace=TRUE)
  X <- X + k * y %*% t(rep(1, n))
  
  # the function to be minimized
  scoreFun <- function(theta) {
    score <- 0
    for (i in 1:b) {
      score <- score - (y[i] * t(X[i,]) %*% theta - log(1 + exp(t(X[i,]) %*% theta)))
    }
    return(score)
  }
  
  # find the starting value
  
  # first, use glm to find the optimum
  best <- glm(y ~ ., family=binomial(link='logit'), data=cbind(as.data.frame(X[,-n]), y=y))
  best <- unname(c(best$coefficients[-1], best$coefficients[1]))
  
  # now, get ny
  nySearch <- function(ny) {
    ((1 + eps) * scoreFun(best) - scoreFun(best + ny * u))**2
  }
  ny <- optimize(nySearch, c(-10, 10))$minimum
  
  # calculate the starting vector
  start <- best + ny * u
  
  # start the optimizer
  res <- optim(start, scoreFun, method=method)

  # return the relevant results
  return(list(iterations=unname(res$counts[2]), error=norm(res$par - best, type='2'),
              suboptimality=(scoreFun(res$par) - scoreFun(best))[1, 1]))
}

print(singleSimulation("CG", 50, 3, 0, 1, c(1, 0, 0)))