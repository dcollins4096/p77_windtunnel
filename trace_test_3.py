

from dtools.starter1 import *
import tracer.trace_pc as t1
from dtools.math import brunt_tools as bt
reload(bt)
reload(t1)
plt.close('all')

N=128
dx = 1/N
n2 = 1.2
slope=0.6
off = 0.3
cube_xyz,cube = t1.get_cube1(N,ampl=.00,val=n2)
cube = bt.get_cubes2('2_2',26,do_rho_4=True)


pos = np.arange(dx*0.5,1,dx)
#ypos = np.zeros_like(xpos)+0.5
xpos,ypos=np.meshgrid(pos,pos)
xpos=xpos.flatten()
ypos=ypos.flatten()
zpos = np.zeros_like(xpos)+0.0
xyz = np.stack([xpos,ypos,zpos])

tracer = t1.tracer(cube,xyz)

tracer.march()

fig,axes=plt.subplots(2,2)
ax0=axes[0][0];ax1=axes[0][1]#;ax2=axes[0][2]
ax3=axes[1][0];ax4=axes[1][1]#; ax5=axes[1][2]
the_x = (cube_xyz[0]).sum(axis=1)/N
the_y = (cube_xyz[2]).sum(axis=1)/N
the_z = tracer.cube.sum(axis=1)
ax0.pcolormesh(the_x,the_y,the_z)
#ax0.plot(tracer.saver[0,:,:].transpose(), tracer.saver[2,:,:].transpose())
#ax1.scatter(tracer.saver[0,:,-1], tracer.saver[1,:,-1])
xf,yf=tracer.saver[0,:,-1], tracer.saver[1,:,-1]
import dtools.vis.pcolormesh_helper as pch
pch.simple_phase(xf.flatten(),yf.flatten(),bins=[64,64],ax=ax1)

fig.savefig('%s/sphere_trace'%plot_dir)
