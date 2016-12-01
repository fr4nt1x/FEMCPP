import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math

f = open('solution','r') 
x = []
y = []
z = []

sol = lambda x,y : 1+x*x + 2* y*y

for line in f.readlines():
    split = line.split()
    x.append(float(split[0]))
    y.append(float(split[1]))
    z.append(float(split[2]) )

error = []

for k,v in enumerate(x):
    error.append(abs(sol(v,y[k]) - z[k]))

print(max(error))
fig1 = plt.figure()
ax1 = Axes3D(fig1)

ax1.plot_trisurf(x,y,z)
plt.show()
