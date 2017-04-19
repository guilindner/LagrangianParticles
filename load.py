import h5py    # HDF5 support
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
fileName = "Particles1.h5"
f = h5py.File(fileName,  "r")

variables = ['x','y','z','vx','vy','vz']
for i in variables:
    globals()[i] = f['DNS']['BEAM'][i]
time = x.shape[0] 
posx = np.zeros(time)
posy = np.zeros(time)
posz = np.zeros(time)

fig = plt.figure()
ax = fig.gca(projection='3d')
for i in range(1):
    for j in range(0,time):
        posx[j] = x[j][i]
        posy[j] = y[j][i]
        posz[j] = z[j][i]
    ax.plot(posx,posy,posz,label='particle{}'.format(i))
ax.legend()    
plt.show()
