

from dtools.starter1 import *
import dtools.davetools as dt
import tracer.trace_pc as t1
reload(t1)
plt.close('all')
import dtools.vis.pcolormesh_helper as pch
import tools

if 1:
    L = 2*np.pi
    L = 1
    N = 64
    dx = L/N
    x = np.arange(0,L,dx)

    xx,yy= np.meshgrid(x,x,indexing='ij')


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
    #p2(proj_ghost,'ghost')
    #ldx = L/N
    ldx = 2*np.pi/N
    print(ldx)
    delta_x = (proj_ghost[2:,sl]-proj_ghost[:-2,sl])/(2*ldx)
    delta_y = (proj_ghost[sl,2:]-proj_ghost[sl,:-2])/(2*ldx)
    rho2 = invert(delta_x,delta_y, mean=proj.sum())
    p3(rho2,'invert_harder')



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
    stat(proj-rho2.real)
    pp(proj-rho2.real,ax2,fig)
    ax2.set(title='orig-recon')
    pp(rho2.imag,ax3,fig)
    ax3.set(title='recon imag')
    #pp(delta_x,ax0,fig)
    #pp(epsx.real,ax1,fig)
    fig.tight_layout()
    fig.savefig('%s/winner'%plot_dir)

