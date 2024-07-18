
from dtools.starter1 import *
import tracer.trace_pc as t1
reload(t1)
plt.close('all')

if 1:
    N=16
    dx = 1/N
    n2 = 1.2
    slope=0.3
    off = 0.3
    #cube_xyz,cube = t1.get_cube2(N,ampl=.00,val=n2, slope=slope,off=off)
    cube_xyz,cube = t1.get_cube1(N,ampl=.00,val=n2)

    length_units = 0.1 #m
    density_units = 1.204 #kg/m^3
    gladstone_dale = 2.3e-4 #m^3/kg

    pos = np.arange(dx*0.5,1,dx/3)
    xpos,ypos=np.meshgrid(pos,pos)
    xpos=xpos.flatten()
    ypos=ypos.flatten()
    zpos = np.zeros_like(xpos)+0.0
    xyz = np.stack([xpos,ypos,zpos])

    #index = density_units*gladstone_dale*cube + 1
    index=cube
    xyz *= length_units

    tracer = t1.tracer(index,xyz, length_units=length_units)

    print('march')
    tracer.march()
    print('marched')

if 1:
    Nrays=pos.size
    xf,yf,zf=tracer.saver[0,:,-1], tracer.saver[1,:,-1], tracer.saver[2,:,-1]
    x0,y0,z0=tracer.saver[0,:,0], tracer.saver[1,:,0], tracer.saver[2,:,0]
    xx,yy,zz=tracer.saver[0,:,:], tracer.saver[1,:,:],tracer.saver[2,:,:]
    for field in [x0,y0,xf,yf]:
        field.shape = Nrays,Nrays

    proj = cube.sum(axis=2)
    phat = np.fft.fftn(proj)

    dsize=1/N
    dx = np.zeros_like(proj)
    dy = np.zeros_like(proj)
    sl = slice(1,-1)
    dx[sl,sl]= (proj[2:,1:-1]-proj[:-2,1:-1])/dsize
    dy[sl,sl]= (proj[1:-1,2:]-proj[1:-1,:-2])/dsize
    Nrays = dx.shape[0]

    #dx = xf-x0
    #dy = yf-y0
    dxhat = np.fft.fftn(dx)
    dyhat = np.fft.fftn(dy)

    kvec = np.fft.fftfreq(Nrays)
    kx, ky = np.meshgrid(kvec,kvec,indexing='ij')
    k2 = kx**2+ky**2
    xvec = np.arange(0,1,1./Nrays)
    x0,y0=np.meshgrid(xvec,xvec,indexing='ij')

    dxhat2 = -1j*kx*phat
    dyhat2 = -1j*ky*phat
    #dx2 = np.fft.irfftn(dxhat2[:,:Nrays//2+1])
    dx2 = np.fft.ifftn(dxhat2)

    k2[0,0]=1 #to avoid 1/0

    D2 = -(1j*kx*dxhat+1j*ky*dyhat)/k2
    D2[0,0]=proj.mean()

    rho = np.fft.ifftn(D2)*dx.size
    print('Imag', np.abs(rho.imag).sum()/rho.size)
    print('Real', np.abs(rho.real).sum()/rho.size)
    rho = rho.real

if 1:
    fig,ax=plt.subplots(3,2,figsize=(8,12))
    ax0=ax[0][0];ax1=ax[0][1];ax2=ax[1][0];ax3=ax[1][1]
    ax4=ax[2][0];ax5=ax[2][1]


    #sl=slice(1,-1)
    sl=slice(None)
    the_x = (cube_xyz[0][sl,sl,sl]).sum(axis=2)/N
    the_y = (cube_xyz[1][sl,sl,sl]).sum(axis=2)/N
    the_z = tracer.cube.sum(axis=2)
    norm = mpl.colors.Normalize(vmin=the_z.min(),vmax=the_z.max())
    p=ax0.pcolormesh(the_x,the_y,the_z,norm=norm)
    fig.colorbar(p,ax=ax0)

    p=ax1.pcolormesh(x0,y0,rho)#,norm=norm)
    fig.colorbar(p,ax=ax1)

    #ax2.scatter(the_z.flatten(),rho.flatten())
    #p=ax2.imshow(dxhat.imag)
    ax2.imshow(dx)
    ax3.imshow(np.real(dx2))
    fig.colorbar(p,ax=ax2)
    #p=ax3.imshow(dxhat2.imag)
    fig.colorbar(p,ax=ax3)

    #rrr=dx.flatten()/dx2.real.flatten()
    rrr=-dx2.real.flatten()*14**2
    ax4.scatter( dx.flatten(), rrr)
    ax4.plot([-20,20],[-20,20])
    fig.savefig('%s/invert1'%(plot_dir))




