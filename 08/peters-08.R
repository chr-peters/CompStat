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
    expect_equal(quadratic_interpolation(function(x) x**2, -1, 1), 0,
                 tolerance=1e-4)
    # shifted 4'th degree polynomial
    expect_equal(quadratic_interpolation(function(x) (x-3)**4 + 1, 1, 4), 3,
                 tolerance=1e-4)
    # find pi/2
    expect_equal(quadratic_interpolation(function(x) -sin(x), 1, 2), pi/2,
                 tolerance=1e-4)
    # additional logarithm
    expect_equal(quadratic_interpolation(function(x) log((x+2)**2)+3, -10, 10), -2,
                 tolerance=1e-4)
  })
}

test_quadratic_interpolation()

# No. 2)
# ======

# a)
plotMedianProblem <- function() {
  # init values
  a <- c(15, 51, 33, 20, 66, 35, 72, 3, 34)
  # init 'search space'
  x <- seq(0, 100, 1)
  print(
    # plot the optimization problem
    ggplot(data.frame(x=x, y=sapply(x, function(beta) sum(abs(a-beta)))),
           aes(x=x, y=y)) +
      geom_point() +
      ggtitle("Median calculation by optimizing the sum of the absolute errors") +
      scale_x_continuous("beta", breaks=2*0:50) +
      scale_y_continuous("Sum of the absolute errors") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}

plotMedianProblem()

# As we can see, the median is 34 because it minimizes the sum of the absolute errors.

# b)
# Iteration: 1
# left = min(a) = 3
# right = max(a) = 72
# best = 3 + phi * (72 - 3) = 45.642
# new = 3 + 72 - 45.642 = 29.358
# SWAP:
# best = 29.358
# new = 45.642
# ADJUST BORDERS:
# right = 45.642
#
# Iteration: 2
# new = 3 + 45.642 - 29.358 = 19.284
# SWAP:
# -
# ADJUST BORDERS:
# left = 19.284
#
# Iteration: 3
# new = 19.284 + 45.642 - 29.358 = 35.568
# SWAP:
# best = 35.568
# new = 29.358
# ADJUST BORDERS:
# left = 29.358
#
# Iteration: 4
# new = 29.358 + 45.642 - 35.568 = 39.432
# SWAP:
# -
# ADJUST BORDERS:
# right = 39.432
#
# Iteration: 5
# new = 29.358 + 39.432 - 35.568 = 33.222
# SWAP:
# best = 33.222
# new = 35.568
# ADJUST BORDERS
# right = 35.568
#
# Iteration: 6
# new = 29.358 + 35.568 - 33.222 = 31.704
# SWAP:
# -
# ADJUST BORDERS:
# left = 31.704
#
# Iteration: 7
# new = 31.704 + 35.568 - 33.222 = 34.05
# SWAP:
# best = 34.05
# new = 33.222
# ADJUST BORDERS
# left = 33.222
#
# Iteration: 8
# new = 33.222 + 35.568 - 34.05 = 34.74
# SWAP:
# -
# ADJUST BORDERS
# right = 34.74
#
# DONE! left = 33.222 <= MEDIAN <= right = 34.74
# => MEDIAN = 34

# c)
# Assumption: All the elements x_i are uniformly distributed in the interval
# [min(x), max(x)].
# Under this assumption, the GSS should terminate as soon as
# right-left <= (max(x)-min(x))/n, where n is the number of elements in x.
#
# In each step of GSS, the interval size is shrinked by the factor phi=0.618,
# which yields the following equation: (max - min)*phi^k <= (max-min)/n, where k
# is the number of iterations after which the GSS terminates.
#
# Solving this for k yields: k <= -ln(n)/ln(phi) = 2.078 * ln(n) which means
# that the number of iterations of the golden section median search is in O(ln(n)).
# So, on average GSS terminates after 2.078 * ln(n) iterations, assuming all
# the elements in x are uniformly distributed.

# d)
# When using GSS to compute the median, the algorithm terminates after
# O(ln(n)) iterations. It should be noted however, that in every single iteration
# the sum of the absolute errors has to be computed, which requires looping
# through the whole dataset once. Thus, the total complexity of the GSS
# median algorithm is O(n*ln(n)), just like an efficient sorting algorithm.
#
# In GSS, the number of iterations is only in O(ln(n)), if we can assume that
# the data is uniformly distributed. This means that in the worst case,
# GSS can be worse than O(n*ln(n)). Taking this into consideration, I would
# prefer using the sorting method, because it is more robust than the
# GSS method and always delivers its result in O(n*ln(n)), provided an
# efficient sorting algorithm is employed.

# e)
# In the median function of R, the first half of the data is partially sorted
# and then the mean of the two biggest elements of that half is retured.
# This is more intelligent than sorting the whole vector,
# because only half of the data has to be sorted.

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