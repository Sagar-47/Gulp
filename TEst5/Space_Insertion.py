f=open("disp.dat",mode='r')
i=0
while(i<61) :
	np = f.readline()
	j=0
	c=0
	while(j<6):
		if (np.find('*') != -1):
			c=np.index("*")
			np = np[:c] + "99999.99" + np[(c+8):]
		c=np.index(".",c)+3
		np = np[:c] +" "+ np[c:]
		j=j+1
	if(i==0) :
		outf=open("disp1.dat",mode='w')
		outf.write(np)
	else :
		outf=open("disp1.dat",mode='a')
		outf.write(np)
	i=i+1

