#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 10:09:53 2024

@author: fizyk
"""

from pandas import *
import numpy as np
import matplotlib.pyplot as plt


plt.figure(figsize=(9,9))

#------------------------------------------------------------------------------
N=2**14
A=5.0
x=A*np.random.rand(N)+.025
y=A*np.random.rand(N)+.025
aveD=A**2/N

rt=(x**2+y**2)**0.5
dr=1/64
dr*=4
rss=np.arange(dr,A-dr,dr)
wyniks=np.ndarray((len(rss),),float)
#------------------------------------------------------------------------------

plt.subplot(2,2,1)



for itr,rs in enumerate(rss):
    indm=np.where( (rt>rs) & (rt<rs+dr))
    ns=len(indm[0])
    if ns!=0:
        locD=6.28*rs*dr/ns/4
        wsp=locD/aveD
        wyniks[itr]=1/wsp
        #print(' {wsp:5.3f}')
    else:
        wyniks[itr]=0
    
    xs,ys=x[indm],y[indm]
    plt.plot(xs,ys,'o',markersize=0.5)
    
    
plt.plot(x,y,'.k',markersize=0.25)
plt.axis('equal')    
plt.title('random/amorphic')


#------------------------------------------------------------------------------
plt.subplot(2,2,2)
plt.plot(rss,wyniks,'-b')
plt.plot([0,A],[1,1],'--',color='gray')
plt.ylim([0,2])


#------------------------------------------------------------------------------

step=1/300
rss=np.arange(0.1,A-step,step)

# value range ->  R, R+dr
dr=20*step


#  lattice parameter
lp=1/4

v=np.arange(-A,A,lp)
a,b=np.meshgrid(v,v)
Ntot=a.size
Area=(2*A)**2
aveD=Area/Ntot
wyniks=np.zeros((len(rss),),float)
art=a.reshape(Ntot,)
brt=b.reshape(Ntot,)

Mr=1000
k0=2*np.pi*dr
k1=1/16

mu=0
sig=0.02
#pcolors=['or','og','ob','om']
print('start')

for i in range(0,Mr):
    
    #xa=art+(0.5-np.random.rand(Ntot,))*k1
    #xb=brt+(0.5-np.random.rand(Ntot,))*k1
    
    xa=art+np.random.normal(mu,sig,Ntot)
    xb=brt+np.random.normal(mu,sig,Ntot)
    
    rt=((xa**2+xb**2)**0.5).reshape(Ntot,)

    
    for itr,rs in enumerate(rss):
        indm=np.where( (rt>rs) & (rt<rs+dr))
        ns=len(indm[0])

        if ns!=0:
            locD=k0*rs/ns
            
            wsp=locD/aveD
            wyniks[itr]+=1/wsp
            #print(itr,rs,f' {wsp:5.3f}')
        else:
            wyniks[itr]+=0         

    
    #aind.append(itr)
wyniks[0]=0    
wyniks/=Mr
print('stop')  

plt.axis('equal') 


plt.subplot(2,2,4)
plt.plot(rss,wyniks,'-b')
plt.plot([0,A],[1,1],'--',color='gray')

#------------------------------------------------------------------------------

dr=lp
rss=np.arange(1/20,A-dr,dr)

v=np.arange(-A,A,lp)
a,b=np.meshgrid(v,v)

xa=art+np.random.normal(mu,sig,Ntot)
xb=brt+np.random.normal(mu,sig,Ntot)
    
rt=((xa**2+xb**2)**0.5).reshape(Ntot,)
    
plt.subplot(2,2,3)
art=xa.reshape(Ntot,)
brt=xb.reshape(Ntot,)


for itr,rs in enumerate(rss):
    indm=np.where( (rt>rs) & (rt<rs+dr))
    ns=len(indm[0])      
    xs,ys=art[indm[0]],brt[indm[0]]
    plt.plot(xs,ys,'o',markersize=5)
    
    
plt.plot(art,brt,'ok',markersize=0.5)
plt.axis('equal')
plt.title('sc lattice + diffusion')



'''
clear
figure(1)
clf
N=0;
M=4;
v=N:M;
vs=N:M-1;

[x,y,z]=meshgrid(v,v,v);

[xs0,ys0,zs0]=meshgrid(vs,vs,v);
[xs1,ys1,zs1]=meshgrid(vs,v ,vs);
[xs2,ys2,zs2]=meshgrid(v ,vs ,vs);

%mp =[x(:),y(:),z(:)];


%%%%% sc
plot3(x(:),y(:),z(:),'or')
hold on


%%%%%  fcc
plot3(xs1(:)+0.5,ys1(:),zs1(:)+0.5,'ob')
plot3(xs2(:),ys2(:)+0.5,zs2(:)+0.5,'ob')
plot3(xs0(:)+0.5,ys0(:)+0.5,zs0(:),'ob')

hold off

cfcc=[x(:)      , y(:)       , z(:); ...
      xs0(:)+0.5 , ys0(:)+0.5 , zs0(:)    ; ...
      xs1(:)+0.5 , ys1(:)     , zs1(:)+0.5; ... 
      xs2(:)     , ys2(:)+0.5 , zs2(:)+0.5 ];




xlim([-1 M+1])
ylim([-1 M+1])
zlim([-1 M+1])
'''


