
from dtools.starter1 import *

import dtools.davetools as dt
import tracer.trace_pc as t1
reload(t1)

import dtools.vis.pcolormesh_helper as pch
from tools import *

if 1:
    N = 10
    rho1=1
    rho2=1.01
    xyz = [3,3,5]
    length_units = 1.0 #m
    density_units = 1.204 #kg/m^3
    gladstone_dale = 2.3e-4 #m^3/kg
    cube = t1.get_cube_impulse(N, xyz, rho1=rho1,rho2=rho2)
    rays = t1.get_rays(N, 2, length_units)
    tracer = t1.tracer(cube,rays, length_units=length_units)
    tracer.march()
    tracer.get_dx()
    t1.image_dx(tracer,fname='%s/dx'%plot_dir)
    

