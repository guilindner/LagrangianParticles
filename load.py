import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.signal
from scipy.stats import norm, kurtosis, rv_continuous

# open file
fileName = "Particles1.h5"
f = h5py.File(fileName,  "r")


variables = ['x','y','z','vx','vy','vz','dxvx','dxvy','dxvz','Dvx']
for i in variables:
    globals()[i] = f['DNS']['BEAM'][i]
time = x.shape[0] 
data = np.zeros(time)
posx = np.zeros(time)
posy = np.zeros(time)
posz = np.zeros(time)
posx2 = np.zeros(time)
posy2 = np.zeros(time)
posz2 = np.zeros(time)
posx3 = np.zeros(time)
posy3 = np.zeros(time)
posz3 = np.zeros(time)

u_rms = np.mean(vx)
#u_fluc = vx[110,:] - u_rms
u_fluc = vx
#velx = np.mean(vx,axis=1)
velx = vx[:,0]
sigma_u = np.std(vx)
#acelx = np.mean(Dvx,axis=1)
acelx = Dvx[:,0]
print(velx.size)
corrx = np.corrcoef(vx)
deltaT = 0.004
Li = np.sum(np.mean(corrx,axis=0))/(3*np.var(vx))*deltaT
print('Integral scale:', Li)
Ti = Li/sigma_u
print('Integral time:', Ti)
A = 0.5 #0.5 - 1.0
epsilon = A*np.std(vx)**3/Li
print('Epsilon:',epsilon)
nu = 1.5e-5 #viscosity
Re_L = np.std(velx)*Li/nu
print('Integral Reynolds number:',Re_L)
eta = (nu**3/epsilon)**(1./4.) #kolmogorov lenght scale
print('Kolmogorov lenght scale:',eta)
tau_l = np.sqrt(2*np.mean(u_fluc**2)/np.mean(acelx**2))
print('tau_l ???', tau_l)
tau_eta = (nu/epsilon)**0.5
print('Kolmogorov time scale', tau_eta)
taylor = (15*nu/epsilon)**0.5 * sigma_u
print('Taylor microscale:',taylor)
Re_taylor = sigma_u*taylor/nu
print('Reynolds Taylor',Re_taylor)
Te = tau_l/u_rms #large-eddy turnover time
print('Large-eddy turnover time', Te)

Np = x[0].size # number of lagrangian tracers
print('Number of lagrangian tracers:',Np)

if sys.argv[1] == 'correlation':
    fig, ax = plt.subplots()
    plt.plot(corrx[0])
    ax.set_xlabel('time')
    ax.set_ylabel('correlation')
    ax.grid(which="both")
    plt.show()

if sys.argv[1] == 'time':
    Tl = np.zeros(3000)
    for i in range(3000):
        t = i
        ui = (vx[0,:])/np.std(vx[0,:])
        uit = (vx[0+t,:])/np.std(vx[0+t,:])
        Rii = np.mean(ui*uit)/np.var(vx[0,:])#(np.mean(ui**2)**0.5 * np.mean(uit**2)**0.5)
        Tl[i] = Rii
    print(Tl)
    print('area',(np.sum(Rii))/np.var(vx[0,:]))
    fig, ax = plt.subplots()
    #ax.set_xscale('log')
    ax.plot(Tl)
    plt.show()

if sys.argv[1] == 'time2':
    Tl = np.zeros(nech)
    for i in range(1,nech):
        for j in range(1, nup):
            Tl[i] += (velx[j]* velx[j+i])/nup
    print(Tl)
    print('area',(np.sum(Tl))/np.var(velx))
    fig, ax = plt.subplots()
    #ax.set_xscale('log')
    ax.plot(Tl)
    plt.show()

if sys.argv[1] == 'second':
    nech = 250
    nup = 1000
    x = np.arange(0,250)#/tau_eta
    s2u = np.zeros(nech)
    s2u2 = np.zeros(nech)
    for i in range(1,nech):
        for j in range(1,nup):
            s2u[i] += (vx[j+i,0]-vx[j,0])**2/nup
    fig, ax = plt.subplots()
    ax.plot(x,s2u)
    ax.set_xlabel('delta t')
    ax.set_ylabel('second order structure function')
    ax.set_xscale('log')
    ax.set_ylim(0,3.0)
    #ax.set_yscale('log')
    ax.grid(which="both")
    plt.show()

if sys.argv[1] == 'spectrum':
    for i in range(500):
        for j in range(0,time):
            data[j] += (vx[j][i]**2+vy[j][i]**2+vz[j][i]**2)**0.5
    data = data/100.
    fs = 20000
    f, Pxx= scipy.signal.welch(data, fs, window='hanning', nperseg=1024, 
    noverlap=16 , nfft=None, detrend='constant', return_onesided=True, 
    scaling='density', axis=-1)
    m = np.arange(100.,1000.)
    n = 300*m**(-2.)

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
    for j in range(0,time):
        posx[j] = x[j][55]
        posy[j] = y[j][55]
        posz[j] = z[j][55]
        posx2[j] = x[j][133]
        posy2[j] = y[j][133]
        posz2[j] = z[j][133]
        posx3[j] = x[j][433]
        posy3[j] = y[j][433]
        posz3[j] = z[j][433]
    ax.plot(posx,posy,posz,label='particle 1')
    ax.plot(posx2,posy2,posz2,label='particle 2')
    ax.plot(posx3,posy3,posz3,label='particle 3')
    ax.legend()    
    plt.show()
elif sys.argv[1] == 'spdf':
    print('single time statistics - pdf')
    
    #vel = (vx[0,:]**2+vy[0,:]**2+vz[0,:]**2)**0.5
    #vel = np.mean(vx,axis=0)
    vel = vx[1,:]
    vel = np.mean(vx,axis=0)
    fig, ax = plt.subplots()
    x = np.linspace(min(vel), max(vel), 100)
    y, bins = np.histogram(vel,bins=100)
    print(y.size)
    print(bins.size)
    ax.plot(bins[:-1], y)
    #plt.plot(x, mlab.normpdf(x, mean, sigma))
    #ax.scatter(vel,norm.logpdf(vel), c='b', label='time = 1')
    #ax.hist(vel, bins=100, normed=True, histtype = 'step')
    plt.legend()
    #ax.semilogy(stvx, norm.pdf(stvx,loc=0), label='norm pdf')
    #ax.set_ylim(1e-4,1)
    #ax.set_xlim(-6,6)
    #ax.set_yscale('log')
    ax.set_xlabel('u_x [m/s]')
    ax.set_ylabel('PDF')
    ax.grid(which="both")
    plt.show()
    
    
elif sys.argv[1] == 'together':
    
    fig, ax = plt.subplots()
    #data = np.mean(vx,axis=0)
    data = vx[0,:]
    (mu, sigma) = norm.fit(data)
    
    # the histogram of the data
    #n, bins, patches = ax.hist(data, 50, normed=1, facecolor='green', alpha=0.75) 
    # add a 'best fit' line
    x = np.linspace(-sigma*3,sigma*3,data.size)
    y = norm.logpdf(x, mu, sigma)
    ax.plot(x, y, 'r--', linewidth=2)
    ax.grid(which="both")
    #ax.set_yscale('log')
    
    plt.show()    
    
    
    
elif sys.argv[1] == 'tpdf':
    print('two time statistics - pdf')

    t = 0
    fig, ax = plt.subplots()
    for i in range(5):
        t = 10*(i+1)
        stvx = vx[t,:] - vx[0,:]
        x = (stvx)/np.sqrt(np.mean((stvx)**2))
        y = norm.pdf(stvx)/t
        ax.scatter(x,y)
    ax.set_xlabel('time')
    ax.grid(which="both")
    ax.set_xlim(-10,10)
    plt.show()
    ax.set_ylabel('log10 PDF vx(t)')

elif sys.argv[1] == 'tpdf2':
    print('two time statistics - pdf')

    t = 0
    fig, ax = plt.subplots()
    for i in range(5):
        t = 10*(i+1)
        stvx = vx[t,:] - vx[0,:]
        x = (stvx*t)/np.sqrt(np.mean((stvx*t)**2))
        y = norm.logpdf(stvx)
        ax.scatter(x,y)
    ax.set_xlabel('t*vx/<(t*vx)^2>^0.5')
    ax.set_ylabel('log10 PDF vx(t)')
    ax.grid(which="both")
    ax.set_xlim(-10,10)
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


elif sys.argv[1] == 'vpdf':

    fig, ax = plt.subplots()
    data = vx[10,:]
    (mu, sigma) = norm.fit(data)
    n, bins, patches = plt.hist(data, 100, normed=1)#, histtype='step')
    y = norm.pdf(bins, mu, sigma)
    ax.plot(bins, y, 'r--', linewidth=2)
    ax.set_xlim(-4,4)
    ax.set_ylabel('PDF')
    ax.set_xlabel('velocity [m/s]')
    ax.grid(which="both")
    plt.show()

elif sys.argv[1] == 'vpdf2':

    fig, ax = plt.subplots()
    data = np.sort(vx[100,:])
    fit = norm.pdf(data, np.mean(data), np.std(data))  
    ax.plot(data,fit,'-o')
    ax.set_xlim(-4,4)
    ax.set_ylabel('PDF')
    ax.set_xlabel('velocity [m/s]')
    ax.grid(which="both")
    plt.show()

elif sys.argv[1] == 'apdf':

    fig, ax = plt.subplots()
    data = Dvx[10,:]
    n, bins, patches = plt.hist(data, 100, normed=1, histtype='step')
    ax.set_xlim(-40,40)
    ax.grid(which="both")
    plt.show()

elif sys.argv[1] == 'apdf2':

    fig, ax = plt.subplots()
    data = np.sort(Dvx[100,:])
    fit = norm.pdf(data, np.mean(data), np.std(data))
    print(fit)
    ax.plot(data/np.std(data),fit,'-o')
    ax.set_xlim(-10,10)
    ax.set_ylabel('PDF')
    ax.set_xlabel('a/<a^2>^0.5')
    ax.grid(which="both")
    plt.show()

else:
    print('Choose one: spectrum, trace, pdf, fluc, ...')



