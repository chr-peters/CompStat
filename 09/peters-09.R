# Name: Christian Peters

library(numDeriv)
library(mvtnorm)

# No. 2)
# ======

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
      }, interval=c(0, 1))$minimum
    
    print(ny)
    
    # take a step
    deltaTheta <- -ny * (invHessian %*% curGradient)[,1]
    new <- last + deltaTheta
    
    # check stopping conditions
    if (curIter >= maxIter) {
      return(new)
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

f <- function(x) -dmvnorm(x, c(0, 0), sigma=diag(2))

print(dfpOptim(f, c(1, 1), maxIter=50))

# Note No3: Quadratsummen, Korrelation

