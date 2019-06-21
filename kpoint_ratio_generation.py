import numpy as np
import math
n = int(input("number of symmetry points:"))
theta1 = math.pi*float(input("angle between 1st and 2nd components(degrees) :"))/180
theta2 = math.pi*float(input("angle between 2nd and 3rd components(degrees) :"))/180
theta3 = math.pi*float(input("angle between 3rd and 1st components(degrees) :"))/180
i = 0
k=-20
a = np.zeros(3*n)
#b = np.zeros((1,3),dtype= np.float)
c = np.zeros(3*n)
#d = np.zeros((1,3),dtype= np.float)
#print(n)
while (i < 3*n):
	a[i] = float(input("factor of 1st component (h) :"))
	a[i+1] = float(input("factor of 2nd component (k):"))
	a[i+2] = float(input("factor of 3rd component (l):"))
	if(i>0) :
		mod = pow((((a[i]-a[i-3])**2) + ((a[i+1]-a[i-2])**2) + ((a[i+2]-a[i-1])**2) + (2*math.cos(theta1)*(a[i]-a[i-3])*(a[i+1]-a[i-2])) + (2*math.cos(theta2)*(a[i+1]-a[i-2])*(a[i+2]-a[i-1])) + (2*math.cos(theta3)*(a[i]-a[i-3])*(a[i+2]-a[i-1]))),1/2)
		c = np.append(c, 100*mod)
	i=i+3
print(mod)
print(c)
