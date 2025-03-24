#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:35:02 2017

@author: alexander
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

try:
    fname1=sys.argv[1]
    print(fname1)
except ValueError:
    print("No data name supplied")
    
try:
    fname2=sys.argv[2]
    print(fname2)
except ValueError:
    print("No data name supplied")

try:
    fname3=sys.argv[3]
    print(fname3)
except ValueError:
    print("No data name supplied")

if len(sys.argv) > 4:    
    try:
        beg=int(sys.argv[4])
        print(beg)
    except ValueError:
        print("Wrong input")
else:
    beg=30
    print(beg)
    
if len(sys.argv) > 5:    
    try:
        end=int(sys.argv[5])
        print(end)
    except ValueError:
        print("Wrong input")
else:
    end=50
    print(end)
    

    
    


loadit=np.loadtxt(fname1)
print(loadit.shape)


#fig3=plt.figure(3)
#ax3=fig3.gca(projection='3d')

#looky=45
#ax3.scatter(loadit[bsize*bsize*(looky-1):bsize*bsize*looky,2],loadit[bsize*bsize*(looky-1):bsize*bsize*looky,3],loadit[bsize*bsize*(looky-1):bsize*bsize*looky,5])
#plt.show()

#beg=50
#end=130
bsize=65
print(beg,end,end-beg,bsize*bsize*beg,bsize*bsize*end)
print(loadit[bsize*bsize*beg:bsize*bsize*end,:].shape)
avstruc=np.reshape(loadit[bsize*bsize*beg:bsize*bsize*end,:],(end-beg,bsize,bsize,9))



#print(avstruc.shape)
#print(avstruc)

avstruc=np.average(avstruc,axis=0)
#print(avstruc.shape)
avstruc=np.reshape(avstruc,(bsize*bsize,9))

#print(avstruc.shape)
#print(avstruc)
"""
fig4=plt.figure(4)
ax4=fig4.gca(projection='3d')
ax4.scatter(avstruc[:,2],avstruc[:,3],avstruc[:,5])
plt.show()
"""
#np.savetxt(fname2,avstruc)
f1=open(fname2,"w")
for i in range(0,avstruc.shape[0]):
	f1.write(" ".join(map(str, avstruc[i,:])))
	f1.write("\n")
	if ((i+1)%bsize == 0):
		f1.write("\n")
		f1.write("\n")
        
rbinsize=40
rmax=7.0
midxy=6.5
rmax2=rmax**2
rav=np.zeros((rbinsize,6))
rcount=np.zeros(rbinsize)
dR=rmax/rbinsize

for i in range(0,rbinsize):
    rav[i,0]=dR*i

for i in range(0,avstruc.shape[0]):
    rtemp2=(avstruc[i,2]-midxy)**2+(avstruc[i,3]-midxy)**2
    if (rtemp2 < rmax2):
        indR=int(np.sqrt(rtemp2)/dR)
        if indR < rbinsize:
            rav[indR,1:6]=rav[indR,1:6]+avstruc[i,4:9]
            rcount[indR]=rcount[indR]+1
                  
for i in range(0,rbinsize):
    #norm=dR**2*((i+1)**2-i**2)*np.pi
    if (rcount[i] > 0):
	rav[i,1:6]=rav[i,1:6]/rcount[i]
           
    
f3=open(fname3,"w")
for i in range(0,rbinsize):
	f3.write(" ".join(map(str, rav[i,:])))
	f3.write("\n")      


