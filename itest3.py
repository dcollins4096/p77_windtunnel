

from dtools.starter1 import *
import dtools.davetools as dt
import tracer.trace_pc as t1
reload(t1)
plt.close('all')
import dtools.vis.pcolormesh_helper as pch

def pp(arr,ax,fig,**kwargs):
    p=ax.imshow(arr,**kwargs)
    fig.colorbar(p,ax=ax)

def p2(arr,fname):
    fig,ax=plt.subplots(1,1)
    pp(arr,ax,fig)
    fig.savefig('%s/%s'%(plot_dir,fname))
def p3(arr,fname,positive=False):
    fig,ax=plt.subplots(1,2,figsize=(8,4))
    aaa = np.abs(arr)

    if positive:
        norm = mpl.colors.Normalize(vmin = aaa[aaa>1e-7].min(),vmax=np.abs(arr).max())
    else:
        norm=None
    cmap=copy.copy(mpl.cm.get_cmap("jet"))
    cmap.set_under([0.5]*3)
    pp(arr.real,ax[0],fig,norm=norm,cmap=cmap)
    pp(arr.imag,ax[1],fig,norm=norm,cmap=cmap)
    fig.savefig('%s/%s'%(plot_dir,fname))
def stat(arr,label=''):
    print(label,np.abs(arr.real).sum(),np.abs(arr.imag).sum())



def invert(deltax,deltay, mean=0):
    dxhat = np.fft.fftn(deltax)
    dyhat = np.fft.fftn(deltay)
    k1 = np.fft.fftfreq(dxhat.shape[0])*N
    kx,ky = np.meshgrid(k1,k1,indexing='ij')
    k2 = kx**2+ky**2
    d2xhat = -1j*kx*dxhat-1j*ky*dyhat
    rhohat = np.zeros_like(dxhat)
    ok = k2>0
    rhohat[ok] = d2xhat[ok]/k2[ok]
    rhohat[k2==0] = mean
    print(rhohat[k2==0])
    rho = np.fft.ifftn(rhohat)
    return rho

if 1:
    N=128
    dx = 1/N
    n2 = 2
    slope=0.6
    off = 0.3
    cube_xyz,sphere = t1.get_cube1(N,ampl=.00,val=n2, center='off')
    if 'cube_raw' not in dir():
        fptr = h5py.File('2_2/DD0026/data0026.cube.h5','r')
        cube_raw = fptr['Density'][()]
        fptr.close()


    length_units = 0.1 #m
    density_units = 1.204 #kg/m^3
    gladstone_dale = 2.3e-4 #m^3/kg
    length_units = 1
    density_units =1
    gladstone_dale = 1

    cube_xyz*=length_units
    #cube=cube_raw
    cube_raw=sphere
    cube = density_units*gladstone_dale*cube_raw + 1
    cube = sphere

if 1:
    proj = cube.mean(axis=2)
    pdb.set_trace()
    projhat = np.fft.fftn(proj)

if 0:
    #perfect inversion.
    k = np.fft.fftfreq(N)*N
    kx, ky = np.meshgrid(k, k, indexing='ij')
    k2 = kx**2+ky**2

    epsx_hat = 1j*kx*projhat
    epsy_hat = 1j*ky*projhat
    epsx  = np.fft.ifftn(epsx_hat)
    epsy  = np.fft.ifftn(epsy_hat)
    rho = invert(epsx,epsy, mean=proj.mean()*N*N)
    p3(rho,'rho_back')
    stat(rho,'rho')
    stat(proj,'proj')
    stat(rho.real-proj,'zero')
    rho2=rho

if 1:
    L = length_units
    #harder inversion
    proj_ghost = np.zeros(nar(proj.shape)+2)
    proj_ghost[1:-1,1:-1]=proj
    sl=slice(1,-1)
    proj_ghost[sl,sl]=proj
    proj_ghost[0,:] =proj_ghost[-2,:]
    proj_ghost[-1,:]=proj_ghost[1,:]
    proj_ghost[:,0] =proj_ghost[:,-2]
    proj_ghost[:,-1]=proj_ghost[:,1]
    #p2(proj_ghost,'ghost')
    ldx = L/N
    delta_x = (proj_ghost[2:,sl]-proj_ghost[:-2,sl])/(2*ldx)
    delta_y = (proj_ghost[sl,2:]-proj_ghost[sl,:-2])/(2*ldx)
    rho2 = invert(delta_x,delta_y, mean=proj.sum())
    p3(rho2,'invert_harder')

if 1:
    #ray casting

    #set up rays
    pos = np.arange(dx*0.5,1,dx)
    xpos,ypos=np.meshgrid(pos,pos)
    xpos=xpos.flatten()
    ypos=ypos.flatten()
    zpos = np.zeros_like(xpos)+0.0
    xyz = np.stack([xpos,ypos,zpos])
    xyz *= length_units

    #cast rays
    tracer = t1.tracer(cube,xyz, length_units=length_units)
    print('march')
    tracer.march()
    print('marched')

    #catch rays
    xf,yf,zf=tracer.saver[0,:,-1], tracer.saver[1,:,-1], tracer.saver[2,:,-1]
    x0,y0,z0=tracer.saver[0,:,0], tracer.saver[1,:,0], tracer.saver[2,:,0]
    xx,yy,zz=tracer.saver[0,:,:], tracer.saver[1,:,:],tracer.saver[2,:,:]
    ddx = xf-x0
    ddy = yf-y0
    x0.shape = pos.size,pos.size
    y0.shape = pos.size,pos.size
    ddx.shape=pos.size,pos.size
    ddy.shape=pos.size,pos.size

    rho_live = invert(ddx,ddy, mean=proj.sum())
    fig,axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
    pp(proj,ax0,fig)
    ax0.set(title='proj')
    pp(delta_x.transpose(),ax1,fig)
    ax1.set(title='dproj/dx')
    pp(ddx,ax2,fig)
    ax2.set(title='ray shift, x')
    the_x = delta_x.transpose().flatten()
    the_y = ddx.flatten()
    pch.simple_phase(the_x,the_y,bins=[16,16],ax=ax3)
    #pp(delta_x.transpose()/ddx,ax3,fig)
    ax3.set(title='spectral')
    fig.tight_layout()
    fig.savefig('%s/dx'%plot_dir)



if 0:
    fig,axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
    ax0.hist(proj.flatten())
    ax1.hist(rho2.real.flatten())
    b = rho2.real.flatten()+0
    b.sort()
    cs= b.cumsum()
    ax2.plot(cs)
    fig.savefig('%s/hists'%plot_dir)


if 1:
    fig,axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
    ext = dt.extents()
    ext(proj)
    ext(rho2.real)
    norm = mpl.colors.Normalize(vmin=ext.minmax[0],vmax=ext.minmax[1])
    norm = None
    pp(proj,ax0,fig,norm=norm)
    ax0.set(title='original')
    pp(rho2.real,ax1,fig,norm=norm)
    ax1.set(title='reconstructed')
    pp(proj-rho2.real,ax2,fig)
    ax2.set(title='orig-recon')
    pp(rho2.imag,ax3,fig)
    ax3.set(title='recon imag')
    #pp(delta_x,ax0,fig)
    #pp(epsx.real,ax1,fig)
    fig.tight_layout()
    fig.savefig('%s/winner'%plot_dir)

