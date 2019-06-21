import os
from matplotlib import pyplot as plt
import numpy as np
import subprocess
from scipy.optimize import minimize
dataexp = np.loadtxt('abinit_phondis.dat').T
#pars = np.zeros(8)
i=0
def fitfunction(pars):
    os.system('sed -e s/knn1/%f/g ZnS3.gin > lc_loop1.gin' % pars[0])
    os.system('sed -i.bak -e s/knn2/%f/g lc_loop1.gin' % pars[1])
    os.system('sed -i.bak -e s/knn3/%f/g lc_loop1.gin' % pars[2])
    os.system('sed -i.bak -e s/knn4/%f/g lc_loop1.gin' % pars[3])
    os.system('sed -i.bak -e s/knn5/%f/g lc_loop1.gin' % pars[4])
    os.system('sed -i.bak -e s/knn6/%f/g lc_loop1.gin' % pars[5])
    os.system('sed -i.bak -e s/knn7/%f/g lc_loop1.gin' % pars[6])
    os.system('sed -i.bak -e s/knn8/%f/g lc_loop1.gin' % pars[7])
    #os.system('sed -i.bak -e s/knn9/%f/g lc_loop1.gin' % pars[8])
    #os.system('sed -i.bak -e s/knn10/%f/g lc_loop1.gin' % pars[9])
    #os.system('sed -i.bak -e s/knn11/%f/g lc_loop1.gin' % pars[10])
    #os.system('sed -i.bak -e s/knn12/%f/g lc_loop1.gin' % pars[11])
    #os.system('sed -i.bak -e s/knn13/%f/g lc_loop1.gin' % pars[12])
    os.system('gulp < lc_loop1.gin > out')
    os.system( 'perl extract.pl' )
    data4 = np.loadtxt('disp.dat').T
    #data1 = np.loadtxt('phonGM.disp', skiprows=3, usecols=([1])).T
    #data2 = np.loadtxt('phonMK.disp', skiprows=12, usecols=([1])).T
    #data3 = np.loadtxt('phonKG.disp', skiprows=12, usecols=([1])).T
    #data4 = np.concatenate((data1,data2,data3))
    #print data4.shape
    return  np.sum((data4-dataexp)*(data4-dataexp), axis=None)
x0 = [4.4,48.0,1.0,28.0,1.0,28.0,27.9,2.8]
#print(fitfunction(x0))

k = minimize(fitfunction, x0, method='Nelder-Mead',options={'xtol': 1e-8, 'disp': True})
#k = minimize(fitfunction, x0, method='Nelder-Mead')
#,options={'xtol': 1e-11, 'disp': True}
print(k)
