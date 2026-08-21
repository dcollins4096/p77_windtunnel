from dtools.starter1 import *

import dtools.vis.pcolormesh_helper as pch

def get_rays(N,dims,length):
    dx = length/N
    x = np.arange(0,length,dx)+0.5*dx
    x, y = np.meshgrid(x,x,indexing='ij')
    x=x.flatten()
    y=y.flatten()
    z = np.zeros_like(x)+dx*1e-3 #tiny offset to prevent zone edge confusion.
    xyz = np.stack([x,y,z])
    return xyz

def get_cube_impulse(N,xyz, rho1=1,rho2=1.01):
    cube = np.zeros([N,N,N])+rho1
    sl = slice(None)
    sss = [sl,sl,sl]
    for dim in range(len(xyz)):
        if xyz[dim] >= 0:
            sss[dim] = slice(xyz[dim],xyz[dim]+1)
    cube[tuple(sss)] = rho2
    return cube

def get_cubesin(N,k=[3,5],ampl=0,center=None):
    x,y,z = (np.mgrid[0:N,0:N,0:N]+0.5)/N
    sin = np.sin(2*np.pi*(k[0]*x+k[1]*y))
    cube = sin*1e-3+1

    return cube


def get_cubeslab(N,ampl=0,center=None):
    xyz = np.mgrid[0:N,0:N,0:N]/N
    if center is None:
        c = [0.5]*3
    else:
        c = [0.25, 0.35, 0.5]
    print(c)
    out = np.ones([N]*3)
    ok = xyz[1]*N==center
    out[ok]=xyz[0][ok]
    return xyz,out
def get_cube1(N,ampl=0,val=1.05,center=None):
    xyz = np.mgrid[0:N,0:N,0:N]/N
    if center is None:
        c = [0.5]*3
    else:
        c = [0.25, 0.35, 0.5]
    print(c)
    r2 = (xyz[0]-c[0])**2+(xyz[1]-c[1])**2+(xyz[2]-c[2])**2
    out = np.ones([N]*3)
    R = 0.25
    out[r2<R**2] = val
    if 1:
        np.random.seed(8675309)
        rando = np.random.random(out.size)*2*ampl-ampl
        rando.shape = out.shape
        out += rando
    return xyz,out
def get_cube2(N,ampl=0,val=1.05,slope=0.3,off=0.5):
    xyz = (np.mgrid[0:N,0:N,0:N])/N
    out = np.ones([N]*3)
    out[xyz[2] > xyz[0]*slope+off] = val
    if 1:
        np.random.seed(8675309)
        rando = np.random.random(out.size)*2*ampl-ampl
        rando.shape = out.shape
        out += rando
    return xyz,out
root='/data/cb1/Projects/P49_EE_BB/Simulations_128_downsample'
def get_cube_128(density,gladstone,sim='2_2'):

    frame = {'2_2':26,'1_1':30}[sim]
    fptr = h5py.File('%s/%s/DD%04d/data%04d.cube.h5'%(root,sim,frame,frame),'r')
    cube_raw = fptr['Density'][()]
    fptr.close()
    cube = cube_raw*density*gladstone+1
    return cube

def image_dx(tracer,fname):
    fig,axes=plt.subplots(1,3, figsize=(12,4))
    #fig,axes=plt.subplots(2,2)
    #ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
    ax0=axes[0];ax1=axes[1];ax2=axes[2]
    proj = tracer.cube.sum(axis=2)
    p=ax0.pcolormesh(tracer.xplane, tracer.yplane, proj)
    fig.colorbar(p,ax=ax0)
    
    ax0.set(title='proj')
    x0,y0,z0=tracer.saver[0,:,0], tracer.saver[1,:,0], tracer.saver[2,:,0]
    x0.shape = tracer.N[0], tracer.N[1]
    y0.shape = tracer.N[0], tracer.N[1]
    Z = tracer.Dx
    Z.shape = tracer.N[0], tracer.N[1]
    p=ax1.pcolormesh(x0,y0,Z)
    fig.colorbar(p,ax=ax1)
    Z = tracer.Dy
    Z.shape = tracer.N[0], tracer.N[1]
    p=ax2.pcolormesh(x0,y0,Z)
    fig.colorbar(p,ax=ax2)

    fig.tight_layout()
    fig.savefig(fname)


def ploot3(rho1, rho2, fname, unity=False, fit=False):

    fig, axes=plt.subplots(2,2,figsize=(8,8))
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]

    f1 = rho1
    pl=ax0.imshow( f1.transpose())
    fig.colorbar(pl,ax=ax0)
    f2=rho2.real
    pl=ax1.imshow(f2.transpose())
    fig.colorbar(pl,ax=ax1)

    pch.simple_phase(f1.flatten(),f2.flatten(),ax=ax2,bins=[10,10])
    if unity:
        m = f1.min()
        x = f1.max()
        ax2.plot([m,x],[m,x])
        #ax2.set_aspect('equal')
    if fit:
        the_x = f1.flatten()
        the_y = f2.flatten()
        pfit = np.polyfit(the_x,the_y,1)
        ax2.plot(the_x, pfit[0]*the_x+pfit[1],c='k')
        ax2.set_title('%0.1e x + %0.1e'%(pfit[0],pfit[1]))

    #ax3.hist((pfit[0]*f1+pfit[1]).flatten(),histtype='step',color='b')
    ax3.hist((f1).flatten(),histtype='step',color='b')
    ax3.hist(f2.flatten(),histtype='step',color='r')
    ok = np.abs(f2)>0
    r = (f1[ok]/f2[ok]).mean()
    if r < 1: 
        r = 1/r
    ax3.set_title('ratio = %0.1e'%r)


    fig.tight_layout()
    fig.savefig(fname)

def ploot2(rho1, rho2, fname, unity=False, fit=False):

    fig, axes=plt.subplots(1,3,figsize=(12,3))
    #ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
    ax0=axes[0];ax1=axes[1];ax2=axes[2]

    f1 = rho1
    pl=ax0.imshow( f1.transpose())
    fig.colorbar(pl,ax=ax0)
    f2=rho2.real
    pl=ax1.imshow(f2.transpose())
    fig.colorbar(pl,ax=ax1)

    pch.simple_phase(f1.flatten(),f2.flatten(),ax=ax2,bins=[10,10])
    if unity:
        m = f1.min()
        x = f1.max()
        ax2.plot([m,x],[m,x])
        #ax2.set_aspect('equal')
    if fit:
        the_x = f1.flatten()
        the_y = f2.flatten()
        pfit = np.polyfit(the_x,the_y,1)
        ax2.plot(the_x, pfit[0]*the_x+pfit[1],c='k')
        ax2.set_title('%0.1e x + %0.1e'%(pfit[0],pfit[1]))


    fig.tight_layout()
    fig.savefig(fname)

def sinner():
    #kvec = [5,6]
    kvec = [0,1]

    f=1#2*np.pi
    proj = np.sin(f*(kvec[0]*xx+kvec[1]*yy))
    kvec = [3,1]
    proj += np.sin(f*(kvec[0]*xx+kvec[1]*yy))

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
    k1 = np.fft.fftfreq(dxhat.shape[0])*dxhat.shape[0]
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
