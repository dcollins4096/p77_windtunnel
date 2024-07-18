

from dtools.starter1 import *
import dtools.davetools as dt
import tracer.trace_pc as t1
reload(t1)
plt.close('all')

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

if 1:
    L = 2*np.pi
    N = 64
    dx = L/N
    x = np.arange(0,L,dx)

    xx,yy= np.meshgrid(x,x,indexing='ij')

    #kvec = [5,6]
    kvec = [0,1]

    f=1#2*np.pi
    proj = np.sin(f*(kvec[0]*xx+kvec[1]*yy))
    proj += 2
    kvec = [3,1]
    #proj += np.sin(f*(kvec[0]*xx+kvec[1]*yy))

if 1:
    proj = np.ones([N,N])
    r2 = (xx-L/2)**2+(yy-L/2)**2
    dr = (0.25*L)**2
    proj[r2<dr]=1.2


if 1:
    k = np.fft.fftfreq(N)*N
    kx, ky = np.meshgrid(k, k, indexing='ij')
    k2 = kx**2+ky**2

if 1:
    projhat = np.fft.fftn(proj)


    p2(proj, 'P')
    p3(projhat, 'Phat')


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

if 0:
    #perfect inversion.
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
    #harder inversion
    proj_ghost = np.zeros(nar(proj.shape)+2)
    proj_ghost[1:-1,1:-1]=proj
    sl=slice(1,-1)
    proj_ghost[sl,sl]=proj
    proj_ghost[0,:] =proj_ghost[-2,:]
    proj_ghost[-1,:]=proj_ghost[1,:]
    proj_ghost[:,0] =proj_ghost[:,-2]
    proj_ghost[:,-1]=proj_ghost[:,1]
    p2(proj_ghost,'ghost')
    dx = L/N
    delta_x = (proj_ghost[2:,sl]-proj_ghost[:-2,sl])/(2*dx)
    delta_y = (proj_ghost[sl,2:]-proj_ghost[sl,:-2])/(2*dx)
    rho2 = invert(delta_x,delta_y, mean=proj.sum())
    p3(rho2,'invert_harder')



if 1:
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




if 0:
    #old code
    epsx_hat = -1j*kx*projhat
    epsy_hat = -1j*ky*projhat
    epsx  = np.fft.ifftn(epsx_hat)
    epsy  = np.fft.ifftn(epsy_hat)
    stat(epsx,'epsx')

    p3(epsx_hat, 'Epsx_hat')
    p3(epsx,'Epsx')
    d2x_hat = -kx**2*projhat
    d2y_hat = -ky**2*projhat
    stat(d2y_hat + ky**2*projhat,'should be zero')

    p3(d2x_hat,'d2x_hat')

    d2x = np.fft.ifftn(d2x_hat)
    p3(d2x,'d2x')

    d2_hat = np.zeros_like(d2x_hat)
    ok = k2>0
    d2_hat[ok] = (d2x_hat+d2y_hat)[ok]/k2[ok]
    d2rho = np.zeros_like(d2_hat)
    d2rho[ok] = -(-kx**2*projhat-ky**2*projhat)[ok]/k2[ok]
    stat(d2_hat-d2rho,'better be')
    stat(d2rho-projhat, 'hope')
    p3(d2rho,'d2rho')

    proj_return = np.fft.ifftn(d2rho)
    stat(proj_return)
    p3(proj_return,'ROTK')
    #d2rho  = -k2**2*projhat
    #stat(d2_hat,'d2')
    #stat(d2rho, 'drho')
