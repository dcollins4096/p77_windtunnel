


from dtools.starter1 import *

import dtools.davetools as dt
import tracer.trace_pc as t1
reload(t1)

import dtools.vis.pcolormesh_helper as pch
from tools import *


if 1:
    N = 128
    ldx = 2*np.pi/N
    #cube = t1.get_cubesin(N)
    cube = t1.get_cube_impulse(N,[3,3,4])
    proj = cube.sum(axis=2)
    #proj_ghost = np.zeros(nar(proj.shape))
    #proj_ghost[1:-1,1:-1]=proj
    proj_ghost = np.zeros(nar(proj.shape)+2)
    proj_ghost[1:-1,1:-1]=proj
    sl=slice(1,-1)
    proj_ghost[sl,sl]=np.log(proj)
    proj_ghost[0,:] =proj_ghost[-2,:]
    proj_ghost[-1,:]=proj_ghost[1,:]
    proj_ghost[:,0] =proj_ghost[:,-2]
    proj_ghost[:,-1]=proj_ghost[:,1]
    #p2(proj_ghost,'ghost')
    sl=slice(1,-1)
    #proj_ghost[sl,sl]=proj
    delta_x = (proj_ghost[2:,sl]-proj_ghost[:-2,sl])/(2*ldx)
    delta_y = (proj_ghost[sl,2:]-proj_ghost[sl,:-2])/(2*ldx)
    rho2 = invert(delta_x,delta_y, mean=proj.sum())
    ploot2( proj, rho2, "%s/invert_sin"%plot_dir)


    length_units=1
    rays = t1.get_rays(N, 2, length_units)
    tracer = t1.tracer(cube,rays, length_units=length_units)
    tracer.march()
    tracer.get_dx()
    ploot2( delta_x, tracer.Dx, '%s/dx_vs_dx'%plot_dir)

