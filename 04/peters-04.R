# Name: Christian Peters

# No 1)
#======
# i)
print(1.01 + 1.02 == 1.03)
# ii)
print(0.1*0.05/0.05 == 0.1)
# iii)
print(sqrt(2)**2 == 2)
# iv)
print(2 + 1e32 - 1e32 == 2)
# v)
print(exp(log(3)) == 3)
# vi)
print(choose(23, 2) == factorial(23) / (factorial(2)*factorial(21)))
# example 1: roots of a polynomial
# the polynomial (x-2)^2 = x^2 - 4x + 4 has 2 as the only root
print(polyroot(c(4, -4, 1))==2)
# example 2: sin and cos
# sin(x)**2 + cos(x)**2 = 1
print(sin(0.08)**2 + cos(0.08)**2 == 1)


# No 3)
#======
varTwoPass <- function(x) {
  return(sum((x-mean(x))**2)/(length(x)-1))
}

varTextBook <- function(x) {
  return((sum(x**2) - sum(x)**2/length(x))/(length(x)-1))
}

varAnalytical <- function(x) {
  n <- length(x)
  return((n*(n+1))/12.)
}

print(varAnalytical(1:1000) - var(1:1000))

varEval <- function() {
  maxExponent <- 26
  simCount <- 10
  errors_twoPass <- numeric(maxExponent*simCount)
  errors_textBook <- numeric(maxExponent*simCount)
  errors_R <- numeric(maxExponent*simCount)
  rowIndex <- 1
  for (exponent in 1:maxExponent) {
    for (curSim in 1:simCount) {
      permutation <- sample.int(2**exponent)
      errors_twoPass[rowIndex] <- abs(varTwoPass(permutation)
                                      - varAnalytical(permutation))
      rowIndex <- rowIndex + 1
    }
  }
  res <- data.frame()
}