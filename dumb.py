

from  dtools.starter1 import *


if 0:
    #sin.  works.
    N = 64
    L = 2*np.pi
    dx = L/N
    x = np.linspace(0,L-dx,N)
    #f = np.sin(2*np.pi*x)
    f = np.sin(x)
    fxN = (f[2:]-f[:-2])/(2*dx)
    k = np.fft.fftfreq(N)*N
    #k[N//2]=0
    fhat = np.fft.fft(f)
    fhatx = 1j*k*fhat
    fx = np.fft.ifft(fhatx)
    print('imag',np.abs(fx.imag).sum())

    fig,ax=plt.subplots(1,1)
    ax.plot(x,f,c='k')
    ax.plot(x,fx.real,c='g')
    ax.plot(x[1:-1],fxN,c='r')

    fig.savefig('%s/fft'%plot_dir)

if 1:
    #sin.  works.
    N = 64
    L = 2*np.pi
    dx = L/N
    x = np.linspace(0,L-dx,N)
    #f = np.sin(2*np.pi*x)
    #f = np.sin(x)
    f = np.zeros_like(x)
    f[ np.abs(x-np.pi)<1]=1e-3
    #f = np.abs(x-L/2)
    #f = np.sin(x)
    #f = np.sin(x)+np.sin(2*x)+np.sin(3*x)


    fxN = (f[2:]-f[:-2])/(2*dx)
    k = np.fft.fftfreq(N)*N
    fhat = np.fft.fft(f)
    fhatx = 1j*k*fhat
    fx = np.fft.ifft(fhatx)
    print('imag',np.abs(fx.imag).sum())

    fig,ax=plt.subplots(1,1)
    ax.plot(x,f,c='k')
    ax.plot(x,fx.real,c='g')
    ax.plot(x,fx.imag,c='b')
    ax.plot(x[1:-1],fxN,c='r')

    fig.savefig('%s/fft'%plot_dir)

if 0:
    N = 64
    L = 2*np.pi
    dx = L/N
    x = np.linspace(0,L-dx,N)
    xx,yy=np.meshgrid(x,x,indexing='ij')
    k = np.fft.fftfreq(N)*N
    kx,ky=np.meshgrid(k,k,indexing='ij')


    #kx[:,N//2+1]=0
    #kx[N//2+1,:]=0
    #f = np.sin(xx)
    f = np.zeros_like(xx)
    c0 = L/2
    r2 = (xx-c0)**2+(yy-c0)**2
    f[r2<0.25]=1


    fhat = np.fft.fftn(f)
    fhatx = 1j*kx*fhat
    #fx = np.fft.ifftn(fhatx)
    fx = np.fft.irfftn(fhatx[:,:N//2+1])

    fxN = np.zeros_like(f)*np.nan
    sl=slice(1,-1)
    fxN[sl,sl] = (f[2:,sl]-f[:-2,sl])/(2*dx)

    fig,ax=plt.subplots(2,2,figsize=(8,8))

    #ax0=ax[0];ax1=ax[1];ax2=ax[2]
    ax0=ax[0][0];ax1=ax[0][1];ax2=ax[1][0];ax3=ax[1][1]
    p=ax0.imshow(fxN)
    fig.colorbar(p,ax=ax0)
    p= ax1.imshow(fx.real)
    fig.colorbar(p,ax=ax1)
    ok = ~np.isnan(fxN)
    ax2.scatter( fx[ok].real.flatten(), fxN[ok].flatten())
    ax2.plot([-1,1],[-1,1])

    ax3.imshow(fx.imag)

    fig.savefig('%s/fft2'%plot_dir)

    



