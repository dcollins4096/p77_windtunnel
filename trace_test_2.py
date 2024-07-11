

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
#cube = bt.get_cubes2('2_2',26,do_rho_4=True)

cube_xyz*=0.1

length_units = 0.1 #m
density_units = 1.204 #kg/m^3
gladstone_dale = 2.3e-4 #m^3/kg
#gladstone_dale = 2.3e-2 #m^3/kg


for n in [70]:
    ypos = np.arange(dx*0.5,1,dx)
    xpos = np.zeros_like(ypos)+n/N

    #xpos,ypos=np.meshgrid(pos,pos)
    xpos=xpos.flatten()
    ypos=ypos.flatten()
    zpos = np.zeros_like(xpos)+0.0
    xyz = np.stack([xpos,ypos,zpos])

    #index = density_units*gladstone_dale*cube + 1
    index = cube
    xyz *= length_units

    tracer = t1.tracer(index,xyz,length_units=length_units)

    tracer.march()
    xf,yf,zf=tracer.saver[0,:,-1], tracer.saver[1,:,-1], tracer.saver[2,:,-1]
    x0,y0,z0=tracer.saver[0,:,0], tracer.saver[1,:,0], tracer.saver[2,:,0]
    xx,yy,zz=tracer.saver[0,:,:], tracer.saver[1,:,:],tracer.saver[2,:,:]

    fig,axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1]#;ax2=axes[0][2]
    ax3=axes[1][0];ax4=axes[1][1]#; ax5=axes[1][2]
    the_x = (cube_xyz[0]).sum(axis=2)/N
    the_y = (cube_xyz[1]).sum(axis=2)/N
    the_z = tracer.cube.sum(axis=2)
    ax0.pcolormesh(the_x,the_y,the_z)

    ax0.plot(xx.transpose(),yy.transpose())
    dx=.1/16
    L0=length_units
    ax0.set(xlim=[-dx,L0+dx],ylim=[-dx,L0+dx])
    #ax0.set(xlim=[-dx,L0+dx])#,ylim=[-dx,L0+dx])
    #ax0.set(xlim=[0.0115,0.012])


    the_x = (cube_xyz[1])[n,:,:]
    the_y = (cube_xyz[2])[n,:,:]
    the_z = tracer.cube[n,:,:]
    ax3.pcolormesh(the_x,the_y,the_z)
    ax3.plot(yy.transpose(),zz.transpose())
    ax3.set(xlim=[-dx,L0+dx],ylim=[-dx,L0+dx])

    the_x = (cube_xyz[0])[:,n,:]
    the_y = (cube_xyz[2])[:,n,:]
    the_z = tracer.cube[:,n,:]
    ax4.pcolormesh(the_x,the_y,the_z)
    ax4.plot(xx.transpose(),zz.transpose())
    ax4.set(xlim=[-dx,L0+dx],ylim=[-dx,L0+dx])


#ax0.plot(tracer.saver[0,:,:].transpose(), tracer.saver[2,:,:].transpose())
#ax1.scatter(tracer.saver[0,:,-1], tracer.saver[1,:,-1])
    import dtools.vis.pcolormesh_helper as pch
    #pch.simple_phase(xf.flatten(),yf.flatten(),bins=[64,64],ax=ax1)
    #ax1.scatter(xf.flatten(),yf.flatten())
    ax1.scatter((xf-x0).flatten(),(yf-y0).flatten())
    #ax1.set(xlim=[-dx,L0+dx],ylim=[-dx,L0+dx])

    fig.savefig('%s/sphere_trace_x%04d'%(plot_dir,n))
