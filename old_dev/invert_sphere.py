

from dtools.starter1 import *

import dtools.davetools as dt
import tracer.trace_pc as t1
reload(t1)

import dtools.vis.pcolormesh_helper as pch
from tools import *

def ploot(tracer, rho2, fname):

    fig, axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]

    f1=tracer.cube.sum(axis=2)
    pl=ax0.imshow( f1.transpose())
    fig.colorbar(pl,ax=ax0)
    f2=rho2.real
    pl=ax1.imshow(f2.transpose())
    fig.colorbar(pl,ax=ax1)

    pch.simple_phase(f1.flatten(),f2.flatten(),ax=ax2,bins=[10,10])
    m = f1.min()
    x = f1.max()
    #ax2.plot([m,x],[m,x])


    fig.savefig(fname)

if 'cube_raw' not in dir():
    fptr = h5py.File('2_2/DD0026/data0026.cube.h5','r')
    cube_raw = fptr['Density'][()]
    fptr.close()
if 1:
    N = 128
    rho1=1
    rho2=1.1
    xyz = [3,3,zzz]
    length_units = 1.0 #m
    density_units = 1.204 #kg/m^3
    gladstone_dale = 2.3e-4 #m^3/kg
    #cube = t1.get_cube_impulse(N, xyz, rho1=rho1,rho2=rho2)
    cube = cube_raw*density_units*gladstone_dale+1
    rays = t1.get_rays(N, 2, length_units)
    tracer = t1.tracer(cube,rays, length_units=length_units)
    tracer.march()
    tracer.get_dx()

    rho2 = invert( tracer.Dx, tracer.Dy, mean=cube.sum())

    ploot(tracer, rho2, '%s/invert_turb'%plot_dir)
