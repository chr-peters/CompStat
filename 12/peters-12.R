# Name: Christian Peters

set.seed(1234)

# No. 1)
# ======

sampleRNG <- function(n) {
  sample(2^31, n) / 2^31
}

# TODO: Adjust numSamples so that p * numSamples >= 5 for each bucket
gapTest <- function(generator, alpha, beta, z, numSamples=1000) {
  # generate random numbers
  randomNumbers <- generator(numSamples)
  
  # get the gaps
  gaps <- double(z + 1)
  
  curGap <- 0
  for (i in seq_along(randomNumbers)) {
    # test if the gap is broken
    if (alpha <= randomNumbers[i] && randomNumbers[i] <= beta) {
      # only store gaps <= z
      if (curGap <= z) {
        gaps[curGap+1] <- gaps[curGap+1] + 1
      }
      
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
  
  # TODO sum must be 1
  print(sum(probabilities))
  
  # do the chi-squared test
  test <- chisq.test(x=gaps, p=probabilities)
  
  print(test)
}

gapTest(sampleRNG, 0.25, 0.75, 5, numSamples=1000)