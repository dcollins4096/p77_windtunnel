
from dtools.starter1 import *

import dtools.davetools as dt
import tracer.trace_pc as t1
reload(t1)

import dtools.vis.pcolormesh_helper as pch
from tools import *

def analytic(rho1,rho2,L,z,dx):
    dx = (L-z)*(np.log(rho2)-np.log(rho1))
    return dx
def ann_eps(rho1,rho2,L,z):
    eps = (np.log(rho2)-np.log(rho1))
    return eps

dddx=[]
aaan=[]
eeep=[]
zzzz=[]
for zzz in [1,2,3,4,5,6,7,8]:
    N = 10
    rho1=1
    rho2=1.01
    xyz = [3,3,zzz]
    zzzz.append(zzz)
    length_units = 1.0 #m
    density_units = 1.204 #kg/m^3
    gladstone_dale = 2.3e-4 #m^3/kg
    cube = t1.get_cube_impulse(N, xyz, rho1=rho1,rho2=rho2)
    rays = t1.get_rays(N, 2, length_units)
    rays[2,:]+= 1e-3 #don't start at exactly the zone boundary, round off errors are strange)
    tracer = t1.tracer(cube,rays, length_units=length_units)
    tracer.march()
    tracer.get_dx()
    stuff = np.unique(tracer.Dx)
    ann = analytic(rho1,rho2,length_units, xyz[2]/N*length_units,length_units/N)
    eps = ann_eps(rho1,rho2,length_units, xyz[2]/N*length_units)
    mmm = stuff.max()
    dddx.append(mmm)
    aaan.append(ann)
    eeep.append(eps)
    #print('dx/ann dx',stuff.max()/eps)
    #print('dx/eps',stuff.max()/ann)
    #print('ann/eps',ann/eps)
    t1.image_dx(tracer,fname='%s/dx'%plot_dir)

plt.clf()
plt.plot(zzzz,dddx)
plt.plot(zzzz,aaan)
plt.savefig('%s/temp'%plot_dir)
    

