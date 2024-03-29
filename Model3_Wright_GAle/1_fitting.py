import os
import numpy as np
import subprocess
from scipy.optimize import minimize
from scipy.optimize import basinhopping
dataexp = np.genfromtxt('abinit_phondis.dat').T

def fitfunction(pars):
	#pars[0] = abs(pars[0])
	#charge = -2-abs(pars[0])
	pars[1] = abs(pars[1])
	pars[3] = abs(pars[3])
	#pars[5] = abs(pars[5])
	#pars[2] = abs(pars[2])
	#pars[4] = abs(pars[4])
	#pars[6] = abs(pars[6])
	#pars[7] = abs(pars[7])
	#pars[8] = abs(pars[8])
	#pars[9] = abs(pars[9])
	os.system('sed -e s/knn1/%f/g ZnS1.gin > loop12.gin' % pars[0])
	os.system('sed -i.bak -e s/knn2/%f/g loop12.gin' % pars[1])
	os.system('sed -i.bak -e s/knn3/%f/g loop12.gin' % pars[2])
	os.system('sed -i.bak -e s/knn4/%f/g loop12.gin' % pars[3])
	os.system('sed -i.bak -e s/knn5/%f/g loop12.gin' % pars[4])
	#os.system('sed -i.bak -e s/knn6/%f/g loop12.gin' % pars[5])
	#os.system('sed -i.bak -e s/knn7/%f/g loop12.gin' % pars[6])
	#os.system('sed -i.bak -e s/knn8/%f/g loop12.gin' % pars[7])
	#os.system('sed -i.bak -e s/knn9/%f/g loop12.gin' % pars[8])
	#os.system('sed -i.bak -e s/knn10/%f/g loop12.gin' % pars[9])
	#os.system('sed -i.bak -e s/knn11/%f/g loop12.gin' % pars[10])
	#os.system('sed -i.bak -e s/knn12/%f/g loop12.gin' % pars[11])
	os.system('gulp < loop12.gin > out')
	os.system('perl extract.pl')
	os.system('python Space_Insertion.py')
	data4 = np.loadtxt('disp1.dat').T #problem with gulp, sometimes it generates 2 neighbouring values not separated by space for bad intial guesses 
	value = np.sum((data4-dataexp)*(data4-dataexp), axis=None)
	#print(pars)
	print(value)
	return value
#x0 = [10, 10,10, 1, 1000, 100, 100,10, 5, 100,5]
#x0 =  [1.0,-3.0,213.0,0.4,11413.0,0.1,664.3,10.5,1.0,129.0] 0.005,3
x0 = [672.28,0.39,1200.00,0.15,13.3]
#x0 = [13.3]
#fitfunction(x0)
#print(len(x0))
#xmin = [1., 1.]
#xmax = [11., 11.]
#bounds = [(low, high) for low, high in zip(xmin, xmax)]
#minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
#minimizer_kwargs = { "method": "Nelder-Mead" }
#k = basinhopping(fitfunction, x0, minimizer_kwargs=minimizer_kwargs,T=0.5, stepsize=5, niter=200, disp=True, niter_success=5)
k = minimize(fitfunction, x0, method='Nelder-Mead', options={'xtol': 1e-5, 'disp': True} )
#k = minimize(fitfunction, x0, method='Nelder-Mead')
#,options={'xtol': 1e-11, 'disp': True}
print(k)
