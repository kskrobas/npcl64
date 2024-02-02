#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 09:40:21 2024

@author: fizyk
"""

from pandas import *

# Import libraries
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
 
'''
# Creating dataset
z = np.random.randint(100, size =(50))
x = np.random.randint(80, size =(50))
y = np.random.randint(60, size =(50))
 
# Creating figure
fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
 
# Creating plot
ax.scatter3D(x, y, z, color = "green")
plt.title("simple 3D scatter plot")
 
# show plot
plt.show()

'''

#  lattice parameter
lp=0.5
A=5
Vol=(2*A)**3

# build 3D lattice
v=np.arange(-A,A,lp)
x,y,z=np.meshgrid(v,v,v)
Ntot=x.size
aveD=Vol/Ntot

xrt=x.reshape(Ntot,)
yrt=y.reshape(Ntot,)
zrt=z.reshape(Ntot,)



#  integration step, radius range, integration range dr
step=0.002
rss=np.arange(lp/2,A-step,step)
dr=step*40

# constans
k0=4*np.pi*dr
k1=lp*0.25

wyniks=np.zeros((len(rss),),float)
# Creating figure
#fig,a = plt.subplots(figsize = (16, 9),nrows=1,ncols=1)
#ax=a
#ax = plt.axes(projection ="3d")

#ax.scatter3D(x,y,z)

# number of random collections
Mr=100
#normal distribution parameters
mu=0
sig=0.125


#------------------------------------------------------------------------------
print('start')

for i in range(0,Mr):
    '''
    xa=xrt+(0.5-np.random.rand(Ntot))*k1
    xb=yrt+(0.5-np.random.rand(Ntot))*k1
    xc=zrt+(0.5-np.random.rand(Ntot))*k1
    '''
    
    xa=xrt+(np.random.normal(mu,sig,Ntot))*k1
    xb=yrt+(np.random.normal(mu,sig,Ntot))*k1
    xc=zrt+(np.random.normal(mu,sig,Ntot))*k1
    
    
    rt=((xa**2+xb**2+xc**2)**0.5)
    
    for itr,rs in enumerate(rss):
        indm=np.where( (rt>rs) & (rt<rs+dr))
        ns=len(indm[0])      
        xs,ys,zs=xrt[indm[0]],yrt[indm[0]],zrt[indm[0]]
        #ax.scatter3D(xs,ys,zs,s=150)
        
        if ns!=0:
            #wsp=k0*rs**2/ns/aveD
            wyniks[itr]+=aveD*ns/(k0*rs**2)
            #print(itr,rs,f' {wsp:5.3f}')
        else:
            wyniks[itr]+=0  
                                    
print('stop')
#------------------------------------------------------------------------------       


#ax.set_box_aspect([1,1,1])

wyniks/=Mr
plt.plot(rss,wyniks,'-b')
plt.plot([0,A],[1,1],'--',color='gray')
plt.xlabel('r',fontsize=14)
plt.ylabel('g(r)',fontsize=14)
plt.title('g(r) distribution')



plt.show()
