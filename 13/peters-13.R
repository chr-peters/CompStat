# Name: Christian Peters

library(microbenchmark)

# No. 3)
# ======

print(microbenchmark({
  x <<- matrix(rnorm(1e6), nrow=1000)
}))

print(microbenchmark({
  sums <<- rowSums(x)
}))

print(microbenchmark({
  sums[sums < -100] <- -100
  sums[sums > 100] <- 100
}))

print(microbenchmark({
  norm2 <<- sqrt(sum(sums**2))
}))
