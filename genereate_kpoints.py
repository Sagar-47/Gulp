import numpy as np

n = int(input("enter number of symmetry points:"))
i = 0
a = np.zeros(3*n)
b = np.zeros((1,3),dtype= np.float)
c = np.zeros(3*n)
d = np.zeros((1,3),dtype= np.float)
print(n)
while (i < 3*n):
	a[i] = float(input("enter symmetry 1st point:"))
	a[i+1] = float(input("enter symmetry 2nd point:"))
	a[i+2] = float(input("enter symmetry 3rd point:"))
	j=0
	if(i>0) :
		c[i]= int(input("number of points :"))
		while (j<c[i]) :
			d[0][0] = a[i-3] + j*((a[i] - a[i-3])/20)
			d[0][1] = a[i-2] + j*((a[i+1] - a[i-2])/20)
			d[0][2] = a[i-1] + j*((a[i+2] - a[i-1])/20)
			b = np.append(b,d,axis=0)	
			j=j+1	
	i=i+3
print(b)
print(a)
