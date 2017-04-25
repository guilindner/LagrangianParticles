from scipy.stats import norm
import h5py    # HDF5 support
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import signal
import sys

fileName = "Particles1.h5"
f = h5py.File(fileName,  "r")

variables = ['x','y','z','vx','vy','vz']
for i in variables:
    globals()[i] = f['DNS']['BEAM'][i]
time = x.shape[0] 
posx = np.zeros(time)
posy = np.zeros(time)
posz = np.zeros(time)
velx = np.zeros(time)
data = np.zeros(time)
vel = np.zeros(time)
#mean, var, skew, kurt = norm.stats(moments='mvsk')

if sys.argv[1] == 'spectrum':
    for i in range(100):
        for j in range(0,time):
            data[j] += (vx[j][i]**2+vy[j][i]**2+vz[j][i]**2)**0.5
    data = data/100.
    fs = 20000
    f, Pxx= signal.welch(data, fs, window='hanning', nperseg=1024, 
    noverlap=16 , nfft=None, detrend='constant', return_onesided=True, 
    scaling='density', axis=-1)
    m = np.arange(100.,1000.)
    n = 100*m**(-2.)

    plt.figure()
    plt.loglog(f, np.sqrt(Pxx))
    plt.loglog(m,n)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Spectrum ')
    plt.grid(True, which="both")
    plt.show()
elif sys.argv[1] == 'trace': 
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
elif sys.argv[1] == 'pdf':
    print('pdf')
elif sys.argv[1] == 'fluc':
    for i in range(1):
        for j in range(0,time):
            vel[j] = (vx[j][i]**2+vy[j][i]**2+vz[j][i]**2)**0.5

else:
    print('Choose one: spectrum, trace or pdf')
#plt.show()
#fig, ax = plt.subplots(1,1)
#x = np.linspace(norm.ppf(0.01),norm.ppf(0.99,100))
#print(x)
#ax.plot(posx, norm.pdf(velx),'r-', lw=5, alpha=0.6, label='norm pdf')


