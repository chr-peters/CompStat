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
      if(is.unsorted(group[indexPermutations[i, ]])) {
        return(0)
      }
      return(1)
    })
    
    # store the count
    permutationCounts[i] <- sum(permutationHits)
  }
  
  # do the chi-squared test
  test <- chisq.test(x=permutationCounts)
  
  return(test)
}

testGenerator <- function(runsPerTest) {
  gapResults <-replicate(runsPerTest, gapTest(sampleRNG, 0.25, 0.75, 5)$p.value)
  hist(gapResults, xlab='p-value', main='Distribution of P-Values in the Gaptest',freq=FALSE)
  
  sequenceLength <- 10000
  for (t in 3:5) {
    permutationResults <- replicate(runsPerTest, permutationTest(sampleRNG, t, sequenceLength %/% t)$p.value)
    hist(permutationResults, xlab='p-value', main=paste0('Distribution of P-Values in the Permutationtest, T=', t), freq=FALSE)
  }
  
}

#testGenerator(200)

# No. 2)
# ======

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
  
  legend('topright', legend=c('sampleGeom', 'dgeom'), col=c('black', 'red'), lty=c(1, NA), pch=c(NA, 1))
}

testGeom(0.5)

sampleExp <- function(n, rate) {
  -log(1-sampleRNG(n))/rate
}

testExp <- function(rate, numSamples=10000) {
  # get the random numbers
  randomNumbers <- sampleExp(numSamples, rate)
  
  # get the 'true' density
  xReference <- seq(0, max(randomNumbers), 0.01)
  yReference <- dexp(xReference, rate)
  
  # plot the results
  hist(randomNumbers, xlab='x', ylab='f(x)',
       main=paste0('Comparison of Exponential Distribution Samplers with rate=', rate),
       freq=FALSE, right=FALSE)
  
  lines(xReference, yReference, col='red', lwd=2)
  
  legend('topright', legend=c('sampleExp', 'dexp'), col=c('black', 'red'), lty=1)
}

testExp(1)

sampleNorm <- function(n, mean, sd, k=12) {
  replicate(n, (sum(sampleRNG(k)) - k/2) * sqrt(12/k) * sd + mean)
}

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

# Test auf Normalverteilung:
# shapiro.test

# ein Boxplot pro k (p-Werte muessen gleichverteilt sein)