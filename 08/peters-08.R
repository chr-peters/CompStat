# Name: Christian Peters

library(testthat)
library(ggplot2)
library(mvtnorm)

# No. 1)
# ======

#' Implementation of the quadratic interpolation optimization algorithm.
#' 
#' @param f     A function that is to be minimized.
#' @param left  The left endpoint of the search interval. It must be numeric.
#' @param right The right endpoint of the search interval. It must be numeric.
#' @param tol   Tolerance threshold. Used to adjust the precision of the optimizer.
#' 
#' @return The arguments of f for which f is minimal in the interval [left, right]
quadratic_interpolation <- function(f, left, right, tol=1e-4) {
  # type checks
  stopifnot(is.function(f), is.numeric(c(left, right)), left<=right, is.numeric(tol))
  
  # calculate initial values
  best <- (left+right)/2
  f_right <- f(right)
  f_left <- f(left)
  f_best <- f(best)
  new <- 1/2 * (f_right * (left**2 - best**2) + f_best*(right**2-left**2) +
                  f_left * (best**2 - right**2)) / (f_right*(left-best) + 
                                                    f_best*(right-left) + 
                                                    f_left * (best-right))
  # iterate until the desired precision is reached
  while (abs(right - left) > tol) {
    f_new <- f(new)
    if (f_new < f_best) {
      tmp <- new
      new <- best
      best <- tmp
      tmp <- f_new
      f_new <- f_best
      f_best <- tmp
    }
    if (new < best) {
      left <- new
      f_left <- f_new
    } else {
      right <- new
      f_right <- f_new
    }
    new_old <- new
    new <- 1/2 * (f_right * (left**2 - best**2) + f_best*(right**2 - left**2) +
                    f_left * (best**2 - right**2)) / (f_right*(left-best) + 
                                                      f_best*(right-left) + 
                                                      f_left * (best-right))
    
    # Also check if new is NA. This occurs when dividing by zero when computing
    # 'new' and can happen if 'left' or 'right' is the minimum.
    if (is.na(new) || isTRUE(all.equal(new, new_old)) || new <= left || new >= right) {
      return(best)
    }
  }
  return(best)
}

#' This function is used to test the quadratic interpolation algorithm.
test_quadratic_interpolation <- function() {
  test_that("Test of the quadratic interpolation algorithm.", {
    # simple quadratic function
    expect_equal(quadratic_interpolation(function(x) x**2, -1, 1), 0, tolerance=1e-4)
    # shifted 4'th degree polynomial
    expect_equal(quadratic_interpolation(function(x) (x-3)**4 + 1, 1, 4), 3, tolerance=1e-4)
    # find pi/2
    expect_equal(quadratic_interpolation(function(x) -sin(x), 1, 2), pi/2, tolerance=1e-4)
    # additional logarithm
    expect_equal(quadratic_interpolation(function(x) log((x+2)**2)+3, -10, 10), -2, tolerance=1e-4)
  })
}

test_quadratic_interpolation()

# No. 2)
# ======

# a)

plotMedianProblem <- function() {
  a <- c(15, 51, 33, 20, 66, 35, 72, 3, 34)
  x <- seq(0, 100, 1)
  print(
    ggplot(data.frame(x=x, y=sapply(x, function(beta) sum(abs(a-beta)))), aes(x=x, y=y)) +
      geom_point() +
      ggtitle("Median calculation by optimizing the sum of absolute errors") +
      scale_x_continuous("beta", breaks=2*0:50) +
      scale_y_continuous("Sum of the absolute errors") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}

plotMedianProblem()

# As we can see, the median is 34 because it minimizes the sum of the absolute errors.

# No. 3)
# ======

#' Implementation of the coordinate descent minimization algorithm.
#' 
#' @param f       The function to minimize.
#' @param start   The starting value.
#' @param tol     The minimal step difference. Used as a stopping condition.
#' @param maxIter The maximum number of iterations.
#' 
#' @return A local minimum of f.
coordinateDescent <- function(f, start, tol=1e-4, maxIter=20) {
  # dimensionality of the search space
  n <- length(start)
  # iteration counter
  curIter <- 1
  
  # initialize starting values
  last <- start
  last_stopCondition <- last
  f_last <- f(last)
  
  # internal value used to locate minima along coordinate axes
  eps <- 1e-6
  repeat {
    # determine k'th unit vector
    k <- 1 + curIter %% n
    e_k <- rep(0, n)
    e_k[k] <- 1
    
    # helper function that determines the endpoint of the search interval along
    # the current unit vector e_k where the sign of direction indicates the
    # search direction (either positive or negative)
    getIntervalEndpoint <- function(direction) {
      curIteration <- 1
      repeat {
        curEndpoint <- curIteration * sign(direction)
        f_new <- f(last + curEndpoint * e_k)
        f_test <- f(last + (curEndpoint - eps * sign(direction)) * e_k)
        if (f_new > f_test) {
          return(curEndpoint)
        }
        # increase the search distance exponentially
        curIteration <- curIteration * 2
      }
    }
    
    # minimize f along the vector
    if (f(last + eps * e_k) > f_last) {
      right <- 0
      left <- getIntervalEndpoint(-1)
    } else {
      left <- 0
      right <- getIntervalEndpoint(1)
    }
    
    # get the optimal step size along the current coordinate by using quadratic
    # interpolation
    ny <- quadratic_interpolation(function(ny) f(last + ny * e_k), left, right)
    last <- last + ny * e_k
    f_last <- f(last)
    
    # check for stop condition once all coordinate directions have been updated
    # once.
    if (curIter %% n == 0) {
      # stop if the last (complete) step was sufficiently small
      if (norm(last - last_stopCondition, type="2") < tol) {
        return(last)  
      }
      
      last_stopCondition <- last
    }
    # also stop if the maximum number of iterations was reached
    if (curIter >= maxIter) {
      return(last)
    }

    curIter <- curIter + 1
  }
}

# This function is used to test the coordinate descent algorithm.
test_coordinateDescent <- function() {
  test_that('Test of the coordinate descent algorithm.', {
    # simple one dimensional function
    expect_equal(coordinateDescent(function(x) (x-3)**2+1, 0), 3, tolerance=1e-4)
    # negative bivariate normal distribution with rho > 0
    expect_equal(coordinateDescent(function(x) -dmvnorm(x, c(1, 2),
                                                        matrix(c(1, 0.5, 0.5, 1),
                                                               nrow=2)),
                                   c(0, 0)),
                 c(1, 2), tolerance=1e-4)
  })
}

test_coordinateDescent()