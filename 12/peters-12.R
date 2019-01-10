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

#print(gapTest(sampleRNG, 0.25, 0.75, 5))

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

#print(permutationTest(sampleRNG, 3, 1000))

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
  xReference <- seq(min(randomNumbers), max(randomNumbers))
  yReference <- dgeom(xReference, prob)
  
  # plot the results
  hist(randomNumbers, freq=FALSE, right=FALSE)
  
  lines(xReference, yReference)
}

testGeom(0.5)