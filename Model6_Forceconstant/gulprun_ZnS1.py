import os
import numpy as np
import subprocess

""""dataexp = np.genfromtxt('abinit_phondis.dat').T"""
c=np.zeros(2)

def func(a) :
	c[0] = 4*8.617e-05*pow(a[1]*a[3],0.5)*pow(0.5*(a[0]+a[2]),12)
	c[1] = 4*8.617e-05*pow(a[1]*a[3],0.5)*pow(0.5*(a[0]+a[2]),6)
b = [knn1,knn2,knn3,knn4]
func(b)
os.system('sed -e s/pnn1/%f/g ZnS1.gin > loop12.gin' % c[0])
os.system('sed -i.bak -e s/pnn2/%f/g loop12.gin' % c[1])
b = [knn3,knn4,knn3,knn4]
func(b)
os.system('sed -i.bak -e s/pnn3/%f/g loop12.gin' % c[0])
os.system('sed -i.bak -e s/pnn4/%f/g loop12.gin' % c[1])
b = [knn1,knn2,knn1,knn2]
func(b)
os.system('sed -i.bak -e s/pnn5/%f/g loop12.gin' % c[0])
os.system('sed -i.bak -e s/pnn6/%f/g loop12.gin' % c[1])
"""os.system('gulp < loop12.gin > out')
os.system('perl extract.pl')
os.system('python Space_Insertion.py')
data4 = np.loadtxt('disp1.dat').T #problem with gulp, sometimes it generates 2 neighbouring values not separated by space for bad intial guesses 
value = np.sum((data4-dataexp)*(data4-dataexp), axis=None)
#print(pars)
print(value)"""

"""b = [0.02,17998.4,4.90,16.5]
func(b)
os.system('sed -e s/knn1/%f/g ZnS1.gin > loop12.gin' % c[0])
os.system('sed -i.bak -e s/knn2/%f/g loop12.gin' % c[1])
b = [4.90,16.5,4.90,16.5]
func(b)
os.system('sed -i.bak -e s/knn3/%f/g loop12.gin' % c[0])

os.system('sed -i.bak -e s/knn4/%f/g loop12.gin' % c[1])
b = [0.02,17998.4,0.02,17998.4]
func(b)
os.system('sed -i.bak -e s/knn5/%f/g loop12.gin' % c[0])

os.system('sed -i.bak -e s/knn6/%f/g loop12.gin' % c[1])
os.system('gulp < loop12.gin > out')"""
