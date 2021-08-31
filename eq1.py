#######################
# Import extern files :

import numpy as Np
from matplotlib import pyplot
from pylab import genfromtxt

#############################
# local function definition :


################
# main program :

mat = genfromtxt("lw_end.dat")

x = mat[:,0]
N = mat[:,1]
u = mat[:,2]
D = mat[:,3]


pyplot.plot(x, N, c='g', ls='-', label="Wave form")
pyplot.plot(x, D, c='r', ls='-', label="Depth")
pyplot.legend()
pyplot.xlabel("x")
pyplot.ylabel("y=N(x)")
pyplot.title("Numerical approximation of long waves")

pyplot.savefig("lw_end.png")
