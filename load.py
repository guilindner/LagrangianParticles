import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.signal
from scipy.stats import norm, kurtosis

fileName = "Particles1.h5"
f = h5py.File(fileName,  "r")


variables = ['x','y','z','vx','vy','vz']
for i in variables:
    globals()[i] = f['DNS']['BEAM'][i]
time = x.shape[0] 
data = np.zeros(time)
posx = np.zeros(time)
posy = np.zeros(time)
posz = np.zeros(time)
if sys.argv[1] == 'spectrum':
    for i in range(100):
        for j in range(0,time):
            data[j] += (vx[j][i]**2+vy[j][i]**2+vz[j][i]**2)**0.5
    data = data/100.
    fs = 20000
    f, Pxx= scipy.signal.welch(data, fs, window='hanning', nperseg=1024, 
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
elif sys.argv[1] == 'spdf':
    print('single time statistics - pdf')
    
    #vel = (vx[0,:]**2+vy[0,:]**2+vz[0,:]**2)**0.5
    vel = vx[1,:]

    fig, ax = plt.subplots()
    ax.scatter(vel,norm.pdf(vel), c='b', label='time = 1')
    plt.legend()
    #ax.semilogy(stvx, norm.pdf(stvx,loc=0), label='norm pdf')
    ax.set_ylim(1e-4,1)
    ax.set_xlim(-6,6)
    ax.set_yscale('log')
    ax.set_xlabel('u_x [m/s]')
    ax.set_ylabel('PDF')
    ax.grid()
    plt.show()
elif sys.argv[1] == 'tpdf':
    print('two time statistics - pdf')
    #vel0 = (vx[0,:]**2+vy[0,:]**2+vz[0,:]**2)**0.5
    #vel = (vx[100,:]**2+vy[100,:]**2+vz[100,:]**2)**0.5
    #vel2 = (vx[500,:]**2+vy[500,:]**2+vz[500,:]**2)**0.5
    #vel3 = (vx[3000,:]**2+vy[3000,:]**2+vz[3000,:]**2)**0.5
    vel0 = vx[0,:]
    vel = vx[1,:]
    vel2 = vx[500,:]
    vel3 = vx[3000,:]   
    stvx = vel - vel0
    stvx2 = vel2 - vel0
    stvx3 = vel3 - vel0
    for i in range(10):
        print(stvx[i],stvx2[i],stvx3[i])
    fig, ax = plt.subplots()
    ax.scatter(stvx,norm.pdf(stvx), c='b', label='first')
    ax.scatter(stvx2,norm.pdf(stvx2), c='r', label='second')
    ax.scatter(stvx3,norm.pdf(stvx3), c='y', label='third')
    plt.legend()
    #ax.semilogy(stvx, norm.pdf(stvx,loc=0), label='norm pdf')
    #ax.set_ylim(1e-4,1)
    #ax.set_xlim(-6,6)
    #ax.set_yscale('log')
    ax.grid()
    plt.show()
elif sys.argv[1] == 'fluc':
    fig, ax = plt.subplots()
    ax.plot(vx[:,0], label='velocity x')
    ax.plot(vy[:,0], label='velocity y')
    ax.plot(vz[:,0], label='velocity z')
    ax.grid()
    ax.set_ylim(-6,6)
    ax.set_xlabel('time')
    ax.set_ylabel('velocity [m/s]')
    plt.legend()
    plt.show()

else:
    print('Choose one: spectrum, trace, pdf, fluc, ...')



