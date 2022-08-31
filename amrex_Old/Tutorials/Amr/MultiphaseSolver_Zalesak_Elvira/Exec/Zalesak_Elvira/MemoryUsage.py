import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy import *

data1=loadtxt('activity_noirl.txt')
data2=loadtxt('activity_irl.txt')

plt.plot(data1[:,0],data1[:,2],linewidth=2,color='k',label='No IRL')
plt.plot(data2[:,0],data2[:,2],linewidth=2,color='b',label='Single IRL call')

imgfilename="./Images/MemoryUsage.png"

plt.xlabel('Time (sec)')
plt.ylabel('Memory usgae (MB)')
plt.gca().legend(loc=0)
plt.savefig(imgfilename)
plt.show()
