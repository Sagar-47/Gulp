import numpy as np
import math
N = int(input("Total number of points:"))
n = int(input("Number of symmetry points:"))
theta1 = math.pi*float(input("angle between 1st and 2nd components(degrees) :"))/180
i = 0
a = np.zeros(3*n)
b = np.zeros((1,3),dtype= np.float)
d = np.zeros((1,3),dtype= np.float)
c = np.zeros(3*n,dtype= np.float )
f=0
#print(n)
while (i < 3*n):
	a[i] = float(input("enter symmetry 1st point:"))
	a[i+1] = float(input("enter symmetry 2nd point:"))
	a[i+2] = float(input("enter symmetry 3rd point:"))
	if(i>0) :
		mod = pow((((a[i]-a[i-3])**2) + ((a[i+1]-a[i-2])**2) + ((a[i+2]-a[i-1])**2) + (2*math.cos(theta1)*(a[i]-a[i-3])*(a[i+1]-a[i-2])) + (2*math.cos(theta1)*(a[i+1]-a[i-2])*(a[i+2]-a[i-1])) + (2*math.cos(theta1)*(a[i]-a[i-3])*(a[i+2]-a[i-1]))),1/2)
		#c = np.append(c, 100*mod)
		f=f+mod		
		c[i] = mod
		print(mod)			
	i=i+3
print(f)
i=0
while (i < 3*n):		
	if(i>0) :
		j=0
		e = math.ceil(N*c[i]/f)
		while (j<e) :
			d[0][0] = a[i-3] + j*((a[i] - a[i-3])/e)
			d[0][1] = a[i-2] + j*((a[i+1] - a[i-2])/e)
			d[0][2] = a[i-1] + j*((a[i+2] - a[i-1])/e)
			#print(d)
			b = np.append(b,d,axis=0)	
			j=j+1
		print(e)
	i=i+3
print(b)
print(a)
