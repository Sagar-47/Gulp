import os
from matplotlib import pyplot as plt
import numpy as np
import subprocess
from scipy.optimize import minimize
dataexp = np.loadtxt('abinit_phondis.dat').T
def fitfunction(pars):
    os.system('sed -e s/knn1/%f/g winkel_fit1.gin > lc_loop1.gin' % pars[0])
    os.system('sed -i.bak -e s/knn2/%f/g lc_loop1.gin' % pars[1])
    os.system('sed -i.bak -e s/knn3/%f/g lc_loop1.gin' % pars[2])
    os.system('sed -i.bak -e s/knn4/%f/g lc_loop1.gin' % pars[3])
    os.system('sed -i.bak -e s/knn5/%f/g lc_loop1.gin' % pars[4])
    os.system('sed -i.bak -e s/knn6/%f/g lc_loop1.gin' % pars[5])
    os.system('sed -i.bak -e s/knn7/%f/g lc_loop1.gin' % pars[6])
    os.system('sed -i.bak -e s/knn8/%f/g lc_loop1.gin' % pars[7])
    os.system('sed -i.bak -e s/knn9/%f/g lc_loop1.gin' % pars[8])
    os.system('gulp lc_loop1')
    data1 = np.loadtxt('phonGM.disp', skiprows=3, usecols=([1])).T
    data2 = np.loadtxt('phonMK.disp', skiprows=12, usecols=([1])).T
    data3 = np.loadtxt('phonKG.disp', skiprows=12, usecols=([1])).T
    data4 = np.concatenate((data1,data2,data3))
    #print data4.shape
    return sum((data4-dataexp)*(data4-dataexp))
x0 = [9.87,1.322,3.726,0.0436,1.57,7.80,2.81]
#k = minimize(fitfunction, x0, method='Nelder-Mead',options={'xtol': 1e-8, 'disp': True})
k = minimize(fitfunction, x0, method='Nelder-Mead',options={'xtol': 1e-11, 'disp': True})
print(k)
