from dtools.starter1 import *

import dtools.vis.pcolormesh_helper as pch

def ploot2(rho1, rho2, fname):

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
    m = f1.min()
    x = f1.max()
    ax2.plot([m,x],[m,x])
    #ax2.set_aspect('equal')


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
