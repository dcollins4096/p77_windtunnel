
from dtools.starter1 import *
import tracer.trace_pc as t1
reload(t1)
plt.close('all')

N=16
dx = 1/N
cube_xyz,cube = t1.get_cube1(N,ampl=.00,val=1.2)
cube_xyz,cube = t1.get_cube2(N,ampl=.00,val=1.2, slope=0.8)

xpos = np.arange(dx*0.5,1,dx)
ypos = np.zeros_like(xpos)+0.5
zpos = np.zeros_like(xpos)+0.0
xyz = np.stack([xpos,ypos,zpos])

tracer = t1.tracer(cube,xyz)
#tracer.make_periodic()
#tracer.make_constant()
if 1:
    tracer.march()

    fig,axes=plt.subplots(2,3)
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[0][2]
    ax3=axes[1][0];ax4=axes[1][1]
    ax0.imshow( tracer.cube.sum(axis=1))
    ax1.imshow( tracer.gx.sum(axis=1))
    ax2.imshow( tracer.gy.sum(axis=1))

    the_x = ((cube_xyz[0] + N/2)/N).sum(axis=1)/N
    the_y = ((cube_xyz[2] + N/2)/N).sum(axis=1)/N
    the_z = tracer.cube.sum(axis=1)
    ax3.pcolormesh(the_x,the_y,the_z)
    #ax3.imshow( tracer.cube.sum(axis=1),extent=[0,1,0,1])
    ok = (tracer.saver[0,:,:] >= 0 )*(tracer.saver[2,:,:]>=0)
    #got_here = np.where((tracer.saver[2,:,:]<0).any(axis=0))[0].min()
    #ax3.plot(tracer.saver[0,:,:got_here].transpose(), tracer.saver[2,:,:got_here].transpose())
    ax3.plot(tracer.saver[0,:,:].transpose(), tracer.saver[2,:,:].transpose())

    fig.savefig('%s/cube_gx'%plot_dir)
