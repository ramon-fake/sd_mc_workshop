#!/usr/bin/env python3

# coding: utf-8

# In[1]:


#!/usr/bin/env python
# This import registers the 3D projection, but is otherwise unused.
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import math
from scipy.interpolate import griddata


# In[2]:

def rotrodrigues(vec, phi):
    s = (3,3)
    mat = np.zeros(s)
    a = math.cos((phi/2)*(math.pi/180))
    b = vec[0]*math.sin((phi/2)*(math.pi/180))
    c = vec[1]*math.sin((phi/2)*(math.pi/180))
    d = vec[2]*math.sin((phi/2)*(math.pi/180))
    mat[0][0] = a**2+b**2-c**2-d**2
    mat[0][1] = 2*(b*c-a*d)
    mat[0][2] = 2*(b*d+a*c)
    mat[1][0] = 2*(b*c+a*d)
    mat[1][1] = a**2+c**2-b**2-d**2
    mat[1][2] = 2*(c*d-a*b)
    mat[2][0] = 2*(b*d-a*c)
    mat[2][1] = 2*(c*d+a*b)
    mat[2][2] = a**2+d**2-b**2-c**2
    return mat


# read info from inpsd.dat
f = open('inpsd.dat', 'r')
cdata=f.read()

# Search for lattice repetition
centry=cdata.find('ncell')
clist=cdata[centry:].strip().split()
#print(clist[0:4])
nx=int(clist[1])
ny=int(clist[2])
nz=int(clist[3])

cdata=cdata[cdata.find('cell')+4:]
# Search for lattice
centry=cdata.find('cell')
clist=cdata[centry:].strip().split()
#f.close()
#x1=float(clist[1])
#y1=float(clist[2])
#z1=float(clist[3])
#x2=float(clist[4])
#y2=float(clist[5])
#z2=float(clist[6])
#x3=float(clist[7])
#y3=float(clist[8])
#z3=float(clist[9])

# Search for temperature
centry=cdata.find('Temp')
clist=cdata[centry:].strip().split()
temperature=float(clist[1])

# Search for magnetic field
#centry=cdata.find('hfield')
#clist=cdata[centry:].strip().split()
#field=float(clist[3])
#f.close()

#print("T=",temperature," B=",field)


# In[3]:


# load restart file
rod = rotrodrigues([math.sqrt(2)/2,-math.sqrt(2)/2,0],0.0)
momfile=glob.glob("restart.PdFeIr11.out")[0]
mom=np.loadtxt(momfile,usecols = (4,5,6),comments='#')
momtot=np.loadtxt(momfile,usecols = (3),comments='#')
mt=momtot[0]
#print(mom)
q=np.split(mom,len(mom))
pp = []
for i in range(len(mom)):
    qt=q[i].transpose(1,0)
    pt=np.dot(rod,qt)
    pp.append(pt.transpose(1,0))
mom=np.concatenate([pp[i] for i in range(len(pp))])
#print(mom)

coordfile=glob.glob("coord.right")[0]
coord=np.loadtxt(coordfile,usecols = (1,2,3),comments='#')

jobname=momfile.split('.')[1]
plotname='snapmap.'+jobname+'.png'
plotname2='scatmap.'+jobname+'.png'

# In[4]:


# define x,y coordinates
x=coord[:,0]
y=coord[:,1]

# define values
c=mom[:,2]+1.0e-12

### # find unique x,y entries
mu=np.unique(c)
### 
### # define grid.
### xi = np.linspace(xu.min(),xu.max(),2*(xu.size+1))
### yi = np.linspace(yu.min(),yu.max(),2*(yu.size+1))
### 
### # interpolate data to grid
### #si = griddata((x, y), s, (xi[None,:], yi[:,None]), method='cubic')
### ci = griddata((x, y), c, (xi[None,:], yi[:,None]), method='cubic')
### #sqci = griddata((x, y), sqc, (xi[None,:], yi[:,None]), method='cubic')
### 


# In[5]:



font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 18}

plt.rc('font', **font)
#plt.rc('figsize',(18,8))
###plt.xkcd()
fig,ax=plt.subplots(1,1,figsize=(6,9))
# contour the gridded data, plotting dots at the randomly spaced data points.
#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
#plt.imshow(ci.T,cmap=plt.cm.Wistia,extent=[yu.min(),yu.max(),xu.min(),xu.max()],vmin=0.0,vmax=1)
#plt.tripcolor(x,y,c,cmap=cm.gist_heat_r,shading='gouraud')
plt.tripcolor(x*0.384,y*0.384,c,cmap=cm.coolwarm,shading='gouraud',vmin=-1.0,vmax=1.0)
cbar = plt.colorbar(fraction=0.046, pad=0.04)
cbar.set_label('$M_z$ ($\mu_B$)')
plt.xlabel('$r_x$ (nm)')
plt.ylabel('$r_y$ (nm)')
ax.set_aspect('equal', 'box')
ax.autoscale_view()
#plt.xlim(0,265)
#plt.ylim(0,265)

#plt.title('$M_z$ for NX={:4d}, NY={:4d}'.format(nx,ny))
#plt.title('$S(q)$ for NX='+'{:4.2f}'.format(nx)+'K, B_z='+'{:4.2f}'.format(field)+'T')
fig.tight_layout()
plt.savefig(plotname)
#plt.show()

# In[10]:


### plt.figure(figsize=(16,8))
### # contour the gridded data, plotting dots at the randomly spaced data points.
### #CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
### #CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
### #plt.scatter(x,y,,ci.T,cmap=plt.cm.Wistia)
### plt.scatter(x,y,s=15,c=c,cmap=cm.coolwarm,marker='o')
### plt.colorbar
### plt.xlabel('$r_x$')
### plt.ylabel('$r_y$')
### 
### #plt.title('$S(q)^{1/2}$ for T='+'{:4.2f}'.format(temperature)+'K, B_z='+'{:4.2f}'.format(field)+'T')
### plt.title('$S(q)$ for NX={:4d}, NY={:4d}'.format(nx,ny))
### plt.savefig(plotname2)
### #plt.show()

