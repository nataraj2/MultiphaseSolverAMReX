#!/usr/bin/python

# Exact solution for the Riemann problem - Sod shock tube

from numpy import *
#from pylab import *
from matplotlib import rc, rcParams
from matplotlib.pyplot import *

data=loadtxt('total_no_cells.dat')
figure(1)
plot(data[:,0],data[:,1])

data=loadtxt('activity.txt')
figure(2)
plot(data[:,0],data[:,2],'k')
show()

