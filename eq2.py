#######################
# Import extern files :

import numpy as N
from matplotlib import pyplot
from pylab import genfromtxt

#############################
# local function definition :


################
# main program :

mat = genfromtxt("lw_end.dat")
mat2 = genfromtxt("lw_endNL.dat")

x = mat[:,0]
N = mat[:,1]
N2 = mat2[:,1]


pyplot.plot(x, abs(N-N2), c='g', ls='-', label="water surface heights difference")
pyplot.legend()
pyplot.xlabel("x")
pyplot.ylabel("y")
pyplot.title("Difference between linear and non linear cases")

pyplot.savefig("lw_endNL.png")
