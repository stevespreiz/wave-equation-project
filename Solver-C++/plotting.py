import numpy as np
import time
from matplotlib import pyplot as plt

f = open("out.txt")

x = np.array((f.readline().split("\t"))[:-1],dtype = float)
u = np.array((f.readline().split("\t"))[:-1],dtype = float)

for i in range(34):
    u = np.append(u,np.array((f.readline().split("\t"))[:-1],dtype = float),axis = 0)
u = np.reshape(u, (35,11))

plt.plot(x,u[0,:])
plt.show()

for i in range(34):
    plt.plot(x,u[i,:])
    plt.show()
