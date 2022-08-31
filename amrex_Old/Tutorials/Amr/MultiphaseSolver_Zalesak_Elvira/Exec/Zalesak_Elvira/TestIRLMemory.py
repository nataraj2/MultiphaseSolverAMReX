import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy import *

data1=loadtxt('testirl_1.txt')
data2=loadtxt('testirl_2.txt')

plt.plot(data1[:,0],data1[:,3],linewidth=2,color='k',label='Loop in sub1')
plt.plot(data2[:,0],data2[:,3],linewidth=2,color='b',label='No loop in sub1')

imgfilename="./Images/IRLMemoryUsage.png"

plt.xlabel('Time (sec)')
plt.ylabel('Memory usgae (MB)')
plt.gca().legend(loc=0)
plt.savefig(imgfilename)
plt.show()
