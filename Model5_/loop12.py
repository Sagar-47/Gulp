import os
import numpy as np
import subprocess

""""dataexp = np.genfromtxt('abinit_phondis.dat').T"""
c=np.zeros(2)

def func(a) :
	c[0] = 4*8.617e-05*pow(a[1]*a[3],0.5)*pow(0.5*(a[0]+a[2]),12)
	c[1] = 4*8.617e-05*pow(a[1]*a[3],0.5)*pow(0.5*(a[0]+a[2]),6)
b = [2475.800000,225.500000,3.009900,1.732200]
func(b)
os.system('sed -e s/pnn1/%f/g ZnS1.gin > loop12.gin' % c[0])
os.system('sed -i.bak -e s/pnn2/%f/g loop12.gin' % c[1])
b = [3.009900,1.732200,3.009900,1.732200]
func(b)
os.system('sed -i.bak -e s/pnn3/%f/g loop12.gin' % c[0])
os.system('sed -i.bak -e s/pnn4/%f/g loop12.gin' % c[1])
b = [2475.800000,225.500000,2475.800000,225.500000]
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
os.system('sed -e s/2475.800000/%f/g ZnS1.gin > loop12.gin' % c[0])
os.system('sed -i.bak -e s/225.500000/%f/g loop12.gin' % c[1])
b = [4.90,16.5,4.90,16.5]
func(b)
os.system('sed -i.bak -e s/3.009900/%f/g loop12.gin' % c[0])

os.system('sed -i.bak -e s/1.732200/%f/g loop12.gin' % c[1])
b = [0.02,17998.4,0.02,17998.4]
func(b)
os.system('sed -i.bak -e s/0.000001/%f/g loop12.gin' % c[0])

os.system('sed -i.bak -e s/0.787340/%f/g loop12.gin' % c[1])
os.system('gulp < loop12.gin > out')"""
