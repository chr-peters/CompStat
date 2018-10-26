# No. 3:
#
# a) (9.8)_10 = (1001.11001100...)_2 -> (0.10100e(10000-1000))_2 = (0.10100e100)_2
#
# b) 1. Compare exponents -> e_w = (100)_2
#    2. Add mantissas:
#                        ( 0 . 1 0 1 0 0 )_2
#                      + ( 0 . 1 0 1 0 0 )_2
#                      ---------------------
#                      = ( 1 . 0 1 0 0 0 )_2
#    3. Normalize: (1.01000e100)_2 -> (0.101000e101)_2
#    4. Round: (0.101000e101)_2 -> (0.10100e101)_2
#
# c) (19.6)_10 = (10011.10011001...)_2 -> (0.10100e101)_2
#
# d) The rounding errors occured when converting the number (9.8)_10
#    to the inexact binary representation. When performing the addition,
#    no rounding errors happened, but due to the fact that the inputs were no
#    exact representations of the number (9.8)_10, the result was inexact as well.