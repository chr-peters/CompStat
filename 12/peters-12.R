# Name: Christian Peters

set.seed(1234)

# used for creating permutations
library(gtools)
library(ggplot2)

# No. 1)
# ======

sampleRNG <- function(n) {
  sample(2^31, n) / 2^31
}

#' Imlpementation of the gaptest.
#' 
#' @param generator  A random number generator.
#' @param alpha      Left boundary of the gap.
#' @param beta       Right boundary of the gap.
#' @param z          Maximum gap size.
#' @param numSamples Number of samples.
#' 
#' @return Results of the chi-squared test.
gapTest <- function(generator, alpha, beta, z, numSamples=1000) {
  # generate random numbers
  randomNumbers <- generator(numSamples)
  
  # get the gaps
  gaps <- integer(z + 1)
  
  curGap <- 0
  for (i in seq_along(randomNumbers)) {
    # test if the gap is broken
    if (alpha <= randomNumbers[i] && randomNumbers[i] <= beta) {
      # store the gap count at the corresponding index in the gaps vector
      # all gaps >= z are stored at the index z+1
      gapIndex <- min(curGap + 1, z + 1)
      gaps[gapIndex] <- gaps[gapIndex] + 1
      
      curGap <- 0
      next
    }
    
    # still outside the gap
    curGap <- curGap + 1
  }
  
  # calculate the probabilities for each gap
  probabilities <- double(z + 1)
  for (i in seq_along(gaps)) {
    probabilities[i] <- (beta - alpha) * (1 - beta + alpha)^(i - 1)
  }
  probabilities[z+1] <- (1 - beta + alpha)^z
  
  # do the chi-squared test
  test <- chisq.test(x=gaps, p=probabilities)
  
  return(test)
}

#' Implementation of the permutation test.
#' 
#' @param generator        A random number generator.
#' @param elementsPerGroup Number of elements for each group.
#' @param numGroups        Number of groups.
#' 
#' @return The results of the chi squared test.
permutationTest <- function(generator, elementsPerGroup, numGroups) {
  # organize the random numbers as a matrix
  randomNumbers <- matrix(generator(elementsPerGroup * numGroups),
                          nrow=elementsPerGroup, ncol=numGroups)
  
  # create vector to store the count of each permutation
  permutationCounts <- integer(factorial(elementsPerGroup))
  
  # now get the counts of each permutation
  indexPermutations <- permutations(elementsPerGroup, elementsPerGroup)
  
  for (i in seq_along(permutationCounts)) {
    # check each group if it corresponds to the current permutation
    permutationHits <- apply(randomNumbers, 2, function(group) {
      as.integer(!is.unsorted(group[indexPermutations[i, ]]))
    })
    
    # store the count
    permutationCounts[i] <- sum(permutationHits)
  }
  
  # do the chi-squared test
  test <- chisq.test(x=permutationCounts)
  
  return(test)
}

#' This function is used to test the RNG using both the gaptest as well as the
#' permutation test.
testGenerator <- function(runsPerTest) {
  # get the p-Values of the gaptest and plot their distribution using a histogram
  gapResults <-replicate(runsPerTest, gapTest(sampleRNG, 0.25, 0.75, 5)$p.value)
  hist(gapResults, xlab='p-value', main='Distribution of P-Values in the Gaptest',
       freq=FALSE)
  
  # do the same for the permutation test using multiple T parameters
  sequenceLength <- 10000
  for (t in 3:5) {
    permutationResults <- replicate(runsPerTest,
                                    permutationTest(
                                      sampleRNG, t, sequenceLength %/% t)$p.value)
    hist(permutationResults, xlab='p-value',
         main=paste0('Distribution of P-Values in the Permutationtest, T=', t),
         freq=FALSE)
  }
  
}

# testGenerator(100)

# As we can see, the p-values are uniformly distributed for both tests. This
# is a good indication that H0 holds true, since the p-values would not be
# uniformly distributed if this was not the case.
# These findings lead to the conclusion that the random number generator in
# fact generates uniformly distributed random numbers.

# No. 2)
# ======

#' This function draws n samples from a geometric distribution with prob being
#' the 'success' rate of the bernoulli experiment.
sampleGeom <- function(n, prob) {
  result <- integer(n)
  for (i in seq_along(result)) {
    # do the bernoulli experiments
    numTrials <- 0
    repeat {
      rand <- sampleRNG(1)
      if (rand <= prob) {
        break
      }
      numTrials <- numTrials + 1
    }
    result[i] <- numTrials
  }
  
  return(result)
}

#' This function is used to test the geometric random number generator.
testGeom <- function(prob, numSamples=10000) {
  # get the random numbers
  randomNumbers <- sampleGeom(numSamples, prob)
  
  # get the 'true' density
  xReference <- seq(0, max(randomNumbers), 1)
  yReference <- dgeom(xReference, prob)
  
  # plot the results
  hist(randomNumbers, xlab='x', ylab='f(x)',
       main=paste0('Comparison of Geometric Distribution Samplers with prob=', prob),
       freq=FALSE, right=FALSE)
  
  points(xReference, yReference, col='red', lwd=2)
  
  legend('topright', legend=c('sampleGeom', 'dgeom'), col=c('black', 'red'),
         lty=c(1, NA), pch=c(NA, 1))
}

testGeom(0.5)

#' This function draws n samples from an exponential distribution with parameter
#' rate.
sampleExp <- function(n, rate) {
  -log(1-sampleRNG(n))/rate
}

#' This function is used to test the exponential distribution random number generator.
testExp <- function(rate, numSamples=10000) {
  # get the random numbers
  randomNumbers <- sampleExp(numSamples, rate)
  
  # get the 'true' density
  xReference <- seq(0, max(randomNumbers), 0.01)
  yReference <- dexp(xReference, rate)
  
  # plot the results
  hist(randomNumbers, xlab='x', ylab='f(x)',
       main=paste0('Comparison of Exponential Distribution Samplers with rate=',
                   rate),
       freq=FALSE, right=FALSE)
  
  lines(xReference, yReference, col='red', lwd=2)
  
  legend('topright', legend=c('sampleExp', 'dexp'), col=c('black', 'red'), lty=1)
}

testExp(1)

#' This function draws n samples form a normal distribution with parameters
#' 'mean' and 'sd'. k represents the number of uniform random variables that
#' are added during the sampling procedure.
sampleNorm <- function(n, mean, sd, k=12) {
  replicate(n, (sum(sampleRNG(k)) - k/2) * sqrt(12/k) * sd + mean)
}

#' This function is used to test the normal distribution random number generator.
testNorm <- function(mean, sd, numSamples=10000) {
  # get the random numbers
  randomNumbers <- sampleNorm(numSamples, mean, sd)
  
  # get the 'true' density
  xReference <- seq(min(randomNumbers), max(randomNumbers), 0.01)
  yReference <- dnorm(xReference, mean, sd)
  
  # plot the results
  hist(randomNumbers, xlab='x', ylab='f(x)',
       main=paste0('Comparison of Normal Distribution Samplers with mean=', mean,
                   ' sd=', sd),
       freq=FALSE, right=FALSE)
  
  lines(xReference, yReference, col='red', lwd=2)
  
  legend('topright', legend=c('sampleNorm', 'dnorm'), col=c('black', 'red'), lty=1)
}

testNorm(0, 1)

#' This function checks the distribution of p-values of the shapiro test
#' for different values of k.
kPlots <- function(samplesPerTest=100, testsPerK=100, maxK=20) {
  # store the p-values here
  pValues <- double(testsPerK * maxK)
  # carry out the tests
  for (k in seq_len(maxK)) {
    for (i in seq_len(testsPerK)) {
      pValues[(k-1)*testsPerK+i] <- shapiro.test(
        sampleNorm(samplesPerTest, 0, 1, k))$p.value
    }
  }
  # plot the results using boxplots
  res <- data.frame(pValue=pValues, k=rep(seq_len(maxK), each=testsPerK))
  boxplot(pValue~k, data=res, xlab='k', ylab='p-value',
          main='Distribution of P-Values for Different Choices of k')
}

kPlots()

# As we can see in the plot, k=6 already yields quite satisfactory results.

# No. 3)
# ======

#' This function samples random numbers from a truncated normal distribution using
#' the rejection method.
#' 
#' @param n Number of samples to draw.
#' @param mean Mean of the normal distribution.
#' @param sd Standard deviation of the normal distribution.
#' @param a Lower interval border.
#' @param b Upper interval border.
#' 
#' @return A vector containing n random numbers of the truncated normal distribution.
sampleNormTruncated <- function(n, mean, sd, a, b) {
  res <- double(n)
  i <- 1
  # get the maximum of the normal distribution in the interval
  if (a <= mean && b >= mean) {
    kq <- dnorm(mean, mean=mean, sd=sd)
  } else if (a >= mean) {
    kq <- dnorm(a, mean=mean, sd=sd)
  } else {
    kq <- dnorm(b, mean=mean, sd=sd)
  }
  # try sampling random numbers until n numbers where generated
  while(i <= n) {
    x <- (b-a) * sampleRNG(1) + a
    u <- sampleRNG(1)
    if (u <= dnorm(x, mean=mean, sd=sd) / kq) {
      res[i] <- x
      i <- i + 1
    }
  }
  return(res)
}

#' This function is used to test the truncated normal distributed random number
#' generator.
testNormTruncated <- function(mean, sd, a, b, numSamples) {
  # get the random numbers
  randomNumbers <- sampleNormTruncated(numSamples, mean, sd, a, b)
  
  # get the 'true' density
  xReference <- seq(a-1, b+1, 0.01)
  yReference <- dnorm(xReference, mean, sd)
  
  # plot the results
  hist(randomNumbers, xlab='x', ylab='f(x)',
       main=paste0('Truncated Normal Distribution, mean=', mean, ', sd=', sd,
                   ', a=', a, ', b=', b),
       freq=FALSE, right=FALSE, xlim=c(a-1, b+1))
  
  lines(xReference, yReference, col='red', lwd=2)
  
  legend('topright', legend=c('sampleNormTruncated', 'dnorm'),
         col=c('black', 'red'), lty=1)
}

testNormTruncated(0, 0.5, -1, 1, 1000)