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
    def __init__(self,csvfile):
        self.x,self.y,self.z = np.loadtxt(csvfile,delimiter=',',unpack=True)
        #self.x = x
        #self.y = y
        #self.z = z

    def divisions(self,domain,divisions):
        size = domain.height/divisions
        counter = []
        for i in range(divisions):
            count = ((self.z >= i*size) & (self.z < (i+1)*size)).sum()
            counter.append(count)
        return counter

    def mean_std(self,domain,divisions):
        size = domain.height/divisions

    def variance():
        pass
    def standard_deviation():
        pass
    def skewness():
        pass
    def kurtosis():
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
        return 'Particle object'


domain = Domain(1200,1200)
p0 = Particle('particle0.csv')
