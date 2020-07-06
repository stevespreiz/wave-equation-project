import numpy as np
import time
from matplotlib import pyplot as plt

f = open("out.txt")

x = np.array((f.readline().split("\t"))[2:-3],dtype = float)
u = np.array((f.readline().split("\t"))[2:-3],dtype = float)

for i in range(26):
    u = np.append(u,np.array((f.readline().split("\t"))[2:-3],dtype = float),axis = 0)
u = np.reshape(u, (26,11))

plt.plot(x,u[0,:])

for i in range(26):
    plt.plot(x,u[i,:])

plt.show()
