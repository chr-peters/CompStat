# Name: Christian Peters

set.seed(1234)

# No. 1)
# ======

sampleRNG <- function(n) {
  sample(2^31, n) / 2^31
}

gapTest <- function(generator, alpha, beta, z, numSamples=1000) {
  # generate random numbers
  randomNumbers <- generator(numSamples)
  
  # get the gaps
  gaps <- double(z + 1)
  
  curGap <- 0
  for (i in seq_along(randomNumbers)) {
    # test if the gap is broken
    if (alpha <= randomNumbers[i] && randomNumbers[i] <= beta) {
      # Store the gap count at the corresponding index in the gaps vector.
      # All gaps >= z are stored at the index z+1
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

print(gapTest(sampleRNG, 0.25, 0.75, 5))