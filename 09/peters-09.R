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
      }, interval=c(0, 1))$minimum
    
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
  rho <- seq(-0.95, 0.95, 0.1)
  # get the number of iterations it takes for each rho
  iterations <- sapply(rho, function(rho){
    sigma <- matrix(c(1, rho, rho, 1), nrow=2)
    mean <- c(0, 0)
    # start at a suitable point so that the algorithm converges
    start <- c(0.75, sign(rho)*0.75)
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
      theme(plot.title = element_text(hjust = 0.5))
  )
}

testDfpOptim()

# Note No3: Quadratsummen, Korrelation

