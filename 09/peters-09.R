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

#' Helper function that executes a single simulation from the simulation study.
#' 
#' @param method The optimization algorithm to be used.
#' @param b      The number of observations.
#' @param n      The number of features.
#' @param k      Parameter that adjusts the condition number of the problem.
#' @param eps   Quality of the starting point.
#' @param u      Direction in which to deviate from the optimal starting point.
#' 
#' @return A list containing the number of function calls, the error as the euclidean
#'         distance to the optimum and the suboptimality of the solution.
singleSimulation <- function(method, b, n, k, eps, u) {
  # create the matrix X from random observations
  X <- matrix(rnorm(b*n), nrow=b, ncol=n)
  X[, n] <- 1
  y <- sample(c(1, 0), b, replace=TRUE)
  X <- X + k * y %*% t(rep(1, n))
  
  # restore the intercept
  X[, n] <- 1
  
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
  success <- TRUE
  tryCatch({
    best <- glm(y ~ ., family=binomial(link='logit'),
                data=cbind(as.data.frame(X[,-n]), y=y))
    best <- unname(c(best$coefficients[-1], best$coefficients[1]))
  }, warning = function(w){
    success <<- FALSE
  }, error = function(e){
    success <<- FALSE
  })
  
  if (!success) {
    # no optimum exists for that example
    stop('No optimum exists!')
  }
  
  # now, get ny
  nySearch <- function(ny) {
    ((1 + eps) * scoreFun(best) - scoreFun(best + ny * u))**2
  }
  ny <- optimize(nySearch, c(-10, 10), tol=1e-4)$minimum
  
  # calculate the starting vector
  start <- best + ny * u
  
  # start the optimizer
  res <- optim(start, scoreFun, method=method)

  # return the relevant results
  return(list(functionCalls=unname(res$counts[1]),
              error=norm(res$par - best, type='2'),
              suboptimality=(scoreFun(res$par) - scoreFun(best))[1, 1]))
}

#' This function is used to execute the simulation study.
#' 
#' @return A data.frame that contains the results.
simulationStudy <- function() {
  # create an empty data.frame to store the results
  result <- data.frame(method=character(), b=integer(), n=integer(), k=double(), 
                       eps=double(), functionCalls=integer(), error=double(),
                       suboptimality=double(), stringsAsFactors = FALSE)
  
  curIteration <- 1
  maxIterations <- 3072
  
  # now iterate over all the parameter combinations
  for (method in c("CG", "BFGS")) {
    for (b in c(100, 200, 500, 1000)) {
      for (n in c(5, 10, 15, 20)) {
        for (k in c(0, 0.25, 0.5, 1)) {
          # commenting this out to make it run in a reasonable amount of time
          eps <- 1
          #for (eps in c(1, 2, 3, 4)) {
            # same here
            #for (curDirection in 1:4) {
              # create a random direction
              u <- rnorm(n)
              u <- u / norm(u, type='2')
              
              # run a single simulation
              success <- TRUE
              tryCatch({
                simRes <- singleSimulation(method, b, n, k, eps, u)
              }, error = function(e){
                success <<- FALSE
              })
              
              # try the next example if no result could be obtained
              if (!success) {
                curIteration <- curIteration + 1
                next
              }
              
              # add the results to the data.frame
              result <- rbind(result, list(method=method, b=b, n=n, k=k, eps=eps,
                                           functionCalls=simRes$functionCalls,
                                           error=simRes$error,
                                           suboptimality=simRes$suboptimality),
                              stringsAsFactors = FALSE)
              
              # print progress
              percent <- as.integer(curIteration / maxIterations * 100)
              print(paste0('Iteration: ', curIteration,', ', percent, '% done.'))
              
              curIteration <- curIteration + 1
            #}
          #}
        }
      }
    }
  }
  return(result)
}

#result <- simulationStudy()

testResult <- function(result) {
  # do a t-test to see if the suboptimality of CG is better than that of BFGS
  bfgs <- result[result$method=='BFGS',]$suboptimality
  cg <- result[result$method=='CG',]$suboptimality
  
  # Zero Hypothesis: Mean of suboptimality of BFGS it lower than that of CG
  # Alternative: Suboptimality of BFGS has a greater mean than that of CG
  test <- t.test(bfgs, cg, alternative='greater')
  print(test)
}

#testResult(result)

# We can see that the CG algorithm produces significantly (p-value < 0.05)
# less suboptimal results than the BFGS algorithm according to the t-test.