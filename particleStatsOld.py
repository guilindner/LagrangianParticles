import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.signal
from scipy.stats import norm, kurtosis, rv_continuous

# open file



class Domain:
    def __init__(self,diameter,height):
        self.diameter = diameter
        self.height = height
        self.volume = np.pi * (self.diameter/2)**2 * self.height
    def __repr__(self):
        return 'Domain: Diameter=%f, Height=%f' %(self.diameter,self.height)
    def draw(self):
        z = np.linspace(0, self.height, 20)
        theta = np.linspace(0, 2*np.pi, 20)
        theta_grid, z_grid=np.meshgrid(theta, z)
        x_grid = self.diameter*np.cos(theta_grid) #+ self.diameter/2
        y_grid = self.diameter*np.sin(theta_grid) #+ self.diameter/2
        return x_grid, y_grid, z_grid
        #fig = plt.figure()
        #ax = fig.add_subplot(projection='3d')
        #ax.plot_surface(x_grid, y_grid, z_grid, alpha=0.5)
        #plt.show()

class Particle():
    """reads a csv file with x, y and z coordinates"""
    def __init__(self,csvfile,deltaT):
        self.csvfile = csvfile
        self.deltaT = deltaT
        self.x,self.y,self.z = np.loadtxt(self.csvfile,delimiter=',',unpack=True)
        self.vx = np.gradient(self.x,deltaT)
        self.vy = np.gradient(self.y,deltaT)
        self.vz = np.gradient(self.z,deltaT)
        self.ax = np.gradient(self.vx,deltaT)
        self.ay = np.gradient(self.vy,deltaT)
        self.az = np.gradient(self.vz,deltaT)

    def divisions(self,domain,divisions):
        """ Divides the domain into equal parts and count the occurances"""
        size = domain.height/divisions
        counter = []
        for i in range(divisions):
            count = ((self.z >= i*size) & (self.z < (i+1)*size)).sum()
            counter.append(count)
        return counter

    def report(self):
        u_rms = np.mean(self.vx)
        #u_fluc = vx[110,:] - u_rms
        u_fluc = self.vx - u_rms
        #velx = np.mean(vx,axis=1)
        sigma_u = np.std(self.vx)
        #acelx = np.mean(Dvx,axis=1)

        print("Data size: "+len(self.vx))
        corrx = np.corrcoef(self.vx)

        Li = np.sum(np.mean(corrx))/(3*np.var(self.vx))*self.deltaT
        print('Integral scale:', Li)
        Ti = Li/sigma_u
        print('Integral time:', Ti)
        A = 0.5 #0.5 - 1.0
        epsilon = A*np.std(self.vx)**3/Li
        print('Epsilon:',epsilon)
        nu = 1.5e-5 #viscosity
        Re_L = np.std(self.vx)*Li/nu
        print('Integral Reynolds number:',Re_L)
        eta = (nu**3/epsilon)**(1./4.) #kolmogorov lenght scale
        print('Kolmogorov lenght scale:',eta)
        tau_l = np.sqrt(2*np.mean(u_fluc**2)/np.mean(self.ax**2))
        print('tau_l ???', tau_l)
        tau_eta = (nu/epsilon)**0.5
        print('Kolmogorov time scale', tau_eta)
        taylor = (15*nu/epsilon)**0.5 * sigma_u
        print('Taylor microscale:',taylor)
        Re_taylor = sigma_u*taylor/nu
        print('Reynolds Taylor',Re_taylor)
        Te = tau_l/u_rms #large-eddy turnover time
        print('Large-eddy turnover time', Te)

    def mean_STD(self,counter):
        """Mean of the Sojourn Time Distribution (STD)"""
        pass
    def variance():
        pass
    def standard_deviation():
        pass
    def skewness():
        pass
    def kurtosis():
        pass
    def spectrum():
        pass
    def trace(self,domain):
        x_grid, y_grid, z_grid = domain.draw()
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(x_grid, y_grid, z_grid, alpha=0.5)
        ax.plot(self.x,self.y,self.z,label='particle 1')
        ax.legend()
        plt.show()
    def __len__(self):
        return len(x)
    def __repr__(self):
        return 'Particle('+self.csvfile+','+str(self.deltaT)+')'


domain = Domain(1200,1200)
p0 = Particle('particle0.csv',0.004)
