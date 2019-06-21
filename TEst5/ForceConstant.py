import os
import numpy as np
import subprocess

"""-0.06074
0.00356

-0.00187
0.00295
-0.00097

-0.00092
0.00047
-0.00041

0.00047
-0.00019
-0.00019
"""
n=int(input("No. of force constants to convert :"))
#b = [-0.06074,0.00356,-0.00187,0.00295,-0.00097,-0.00092,0.00047,-0.00041,0.00047,-0.00019,-0.00019]
j=0
a = np.zeros(n)
i=0
while(i<n):
	a[i] = float(input("force constant :"))
	a[i]=a[i]*97.33
	i=i+1
print(a)

