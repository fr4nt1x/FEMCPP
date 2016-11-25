import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math

f = open('solution','r') 
x = []
y = []
z = []
for line in f.readlines():
    
    split = line.split()
    z.append( float(split[0]) )
    x.append(float(split[1]))
    y.append(float(split[2]))
fig1 = plt.figure()

ax1 = Axes3D(fig1)

ax1.plot_trisurf(x,y,z)
plt.show()
