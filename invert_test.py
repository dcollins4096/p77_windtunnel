
from dtools.starter1 import *
import tracer.trace_pc as t1
reload(t1)
plt.close('all')


if 1:
    x = np.linspace(0,L-dx,N)
    N = 64
    L = 2*np.pi
    dx = L/N
    xx,yy=np.meshgrid(x,x,indexing='ij')
    kvec = [3,5]

    proj = np.sin(kvec[0]*xx+kvec[1]*yy)
    kvec = [-2,6]
    proj += np.sin(kvec[0]*xx+kvec[1]*yy)

if 1:
    kvec = np.fft.fftfreq(N)
    kx, ky = np.meshgrid(kvec,kvec,indexing='ij')
    k2 = kx**2+ky**2
    projhat=np.fft.fftn(proj)
    epsx = np.fft.ifftn(-1j*kx*projhat)
    epsy = np.fft.ifftn(-1j*ky*projhat)

    epsxhat = np.fft.fftn(epsx)
    epsyhat = np.fft.fftn(epsy)

    ddx_hat = -1j*kx*epsxhat
    zero = ddx_hat - kx**2*projhat
    print(np.abs(zero).sum())



if 1:

    fig,ax=plt.subplots(2,2,figsize=(8,8))

    #ax0=ax[0];ax1=ax[1];ax2=ax[2]
    ax0=ax[0][0];ax1=ax[0][1];ax2=ax[1][0];ax3=ax[1][1]

    #ax0.imshow(proj)
    #ax1.imshow(rho.real)
    p=ax0.imshow(projhat.real**kx**2)
    fig.colorbar(p,ax=ax0)
    wtf = projhat.imag**kx**2
    p=ax1.imshow(wtf)
    fig.colorbar(p,ax=ax1)
    p=ax2.imshow(ddx_hat.real)
    fig.colorbar(p,ax=ax2)
    p=ax3.imshow(ddx_hat.imag)
    fig.colorbar(p,ax=ax3)
    #ax3.imshow(zero)


    #ax2.imshow(deltax)
    #ax3.imshow( dproj_dx_hh.real)
    #print('imag',np.abs(dproj_dx_hh.imag).sum())

    fig.savefig('%s/poison'%plot_dir)



if 0:

    dsize=1/N
    deltax = np.zeros_like(proj)
    deltay = np.zeros_like(proj)
    sl = slice(1,-1)
    deltax[sl,sl]= (proj[2:,1:-1]-proj[:-2,1:-1])/dsize
    deltay[sl,sl]= (proj[1:-1,2:]-proj[1:-1,:-2])/dsize
    deltax = deltax[sl,sl]
    deltay = deltay[sl,sl]
    Nprim = deltax.shape[0]

if 0:


    dxhat = np.fft.fftn(deltax)
    dyhat = np.fft.fftn(deltay)
    Nprim = dxhat.shape[0]

    kvec = np.fft.fftfreq(Nprim)
    kx, ky = np.meshgrid(kvec,kvec,indexing='ij')
    k2 = kx**2+ky**2

if 0:
    #compare first derivatives.
    proj_hat = np.fft.fftn(proj)
    dproj_dx_hat = 1j*kx*proj_hat
    dproj_dx_hh = np.fft.ifftn(dproj_dx_hat)




if 0:

    D2 = -(1j*kx*dxhat+1j*ky*dyhat)/k2
    D2 = -1j*kx*dxhat
    D2[0,0]=proj.mean()

    rho = np.fft.ifftn(D2)*D2.size
    print('imag',np.abs(rho.imag).sum())

