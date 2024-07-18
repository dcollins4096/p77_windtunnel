
from dtools.starter1 import *
import tracer.trace_pc as t1
reload(t1)
plt.close('all')

N=10
dx = 1/N
n2 = 1.2
slope=0.3
off = 0.3
#cube_xyz,cube = t1.get_cube1(N,ampl=.00,val=n2)
cube_xyz,cube = t1.get_cube2(N,ampl=.00,val=n2, slope=slope,off=off)

xpos = np.arange(dx*0.5,1,dx/2)
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
    ax3=axes[1][0];ax4=axes[1][1]; ax5=axes[1][2]
    ax0.imshow( tracer.cube.sum(axis=1))
    ax1.imshow( tracer.gx.sum(axis=1))
    ax2.imshow( tracer.gy.sum(axis=1))

    #the_x = ((cube_xyz[0] + N/2)/N).sum(axis=1)/N
    #the_y = ((cube_xyz[2] + N/2)/N).sum(axis=1)/N
    the_x = (cube_xyz[0]).sum(axis=1)/N
    the_y = (cube_xyz[2]).sum(axis=1)/N
    the_z = tracer.cube.sum(axis=1)
    ax3.pcolormesh(the_x,the_y,the_z)
    #ax3.imshow( tracer.cube.sum(axis=1),extent=[0,1,0,1])
    ok = (tracer.saver[0,:,:] >= 0 )*(tracer.saver[2,:,:]>=0)
    #got_here = np.where((tracer.saver[2,:,:]<0).any(axis=0))[0].min()
    #ax3.plot(tracer.saver[0,:,:got_here].transpose(), tracer.saver[2,:,:got_here].transpose())
    ax3.plot(tracer.saver[0,:,:].transpose(), tracer.saver[2,:,:].transpose())


    #do I do snells law right?
    theta_in = np.arctan(slope)
    n1 = 1
    n2 = n2
    theta_out = np.arcsin( n1/n2*np.sin(theta_in))
    x_0 = tracer.saver[0,:,0]
    dZ = 1-(slope*x_0+off)
    theta_r = theta_in-theta_out
    Xnew = x_0 - dZ*np.tan(theta_r)
    ax4.scatter(tracer.saver[0,:,0], tracer.saver[0,:,-1],c='k')
    ax4.scatter(x_0,Xnew,c='r')
    ax3.plot(x_0,1-dZ,c='g')

    zmax=tracer.saver[2,:,:].max()
    for ni,xi in enumerate(x_0):
        ax5.plot([xi,xi],[0,1-dZ[ni]], c='r')
        ax5.plot([xi,Xnew[ni]],[1-dZ[ni],zmax],c='g')
    ax5.plot(tracer.saver[0,:,:].transpose(), tracer.saver[2,:,:].transpose(),c='b')

    fig.savefig('%s/cube_gx'%plot_dir)
