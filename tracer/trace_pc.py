
from dtools.starter1 import *
from scipy.ndimage import gaussian_filter
from tools import *

from dtools.math import brunt_tools as bt
reload(bt)




def ddx(array, direction, dx):
    sl = slice(None)
    ii = slice(1,-1)
    ip1 = slice(2,None)
    im1 = slice(0,-2,None)
    rank = len(array.shape)
    iii =   [ii]*rank
    left =  [ii]*rank
    right = [ii]*rank
    left[direction]=ip1
    right[direction]=im1
    out = np.zeros_like(array)
    out[tuple(iii)] = (array[tuple(left)]-array[tuple(right)])/dx
    return out



class tracer():
    def __init__(self,density,rays,dim=2,length_units=1, density_units=1, gladstone_dale=1, skip_gd=False):
        self.density=density
        if not skip_gd:
            self.cube = self.density*density_units*gladstone_dale + 1
        else:
            self.cube=density
        self.rays=rays
        self.dim = dim
        self.length_units=length_units
        self.N = self.cube.shape
        self.xplane, self.yplane, z = get_rays(self.N[0], 2, length_units)
        self.xplane.shape = self.N[0],self.N[1]
        self.yplane.shape = self.N[0],self.N[1]
        self.density_units=density_units
        self.gladstone_dale=gladstone_dale
    def make_constant(self):
        shape = nar(self.cube1.shape)
        shape += 2
        self.cube = np.zeros(shape)
        self.cube[1:-1,1:-1,1:-1]=self.cube1
        self.cube[0,:,:]=self.cube[1,:,:]
        self.cube[:,0,:]=self.cube[:,1,:]
        self.cube[:,:,0]=self.cube[:,:,1]
        self.cube[-1,:,:]=self.cube[-2,:,:]
        self.cube[:,-1,:]=self.cube[:,-2,:]
        self.cube[:,:,-1]=self.cube[:,:,-2]
    def make_periodic(self):
        shape = nar(self.cube1.shape)
        shape += 2
        self.cube = np.zeros(shape)
        self.cube[1:-1,1:-1,1:-1]=self.cube1
        self.cube[0,:,:]=self.cube[-2,:,:]
        self.cube[:,0,:]=self.cube[:,-2,:]
        self.cube[:,:,0]=self.cube[:,:,-2]
        self.cube[-1,:,:]=self.cube[1,:,:]
        self.cube[:,-1,:]=self.cube[:,1,:]
        self.cube[:,:,-1]=self.cube[:,:,1]

    def get_dx_proj(self):
        N = self.N
        #z = (np.arange(N[2])+0.5)/N[2]
        z = (np.arange(N[2]))/N[2]
        z.shape = 1,1,z.size
        dz = 1/N[2]
        self.zproj = ((1-z)*np.log(self.cube)*dz).sum(axis=2)
        dx=4*np.pi/N[0] #domain is 2pi, Npixels, centered stencil.  What's the N for?
        self.dpdx = ddx(self.zproj,0,dx)
        self.dpdy = ddx(self.zproj,1,dx)

    def get_spectral(self):
        N = self.N[0]
        k = np.fft.fftfreq(N)*N
        kx, ky = np.meshgrid(k, k, indexing='ij')
        k2 = kx**2+ky**2

        self.get_dx_proj()
        hat = np.fft.fftn(self.zproj)
        self.spsx = np.fft.ifftn(1j*kx*hat)
        self.spsy = np.fft.ifftn(1j*ky*hat)

    def get_dx(self):
        xf,yf,zf=self.saver[0,:,-1], self.saver[1,:,-1], self.saver[2,:,-1]
        x0,y0,z0=self.saver[0,:,0], self.saver[1,:,0], self.saver[2,:,0]
        #xx,yy,zz=self.saver[0,:,:], self.saver[1,:,:],self.saver[2,:,:]
        #This 4pi is necessary, but not in the right place.
        self.Dx = (xf-x0)/(4*np.pi)
        self.Dy = (yf-y0)/(4*np.pi)
        shape = self.N[0:2]
        self.Dx.shape=shape
        self.Dy.shape=shape

    def image2(self, ftool,fname):
        fig,axes=plt.subplots(2,2, figsize=(8,8))
        ax0=axes[0][0];ax1=axes[0][1]
        ax2=axes[1][0];ax3=axes[1][1]

        f1 = self.density.mean(axis=2)
        f2 = self.recovered

        pl=ax0.imshow( f1.transpose())
        fig.colorbar(pl,ax=ax0)
        pl=ax1.imshow(f2.transpose())
        ax0.set_title('Original')
        fig.colorbar(pl,ax=ax1)
        ax1.set_title('Reconstructed')

        pch.simple_phase(f1.flatten(),f2.flatten(),ax=ax2,bins=[16,16])
        ax2.set_title('orig vs recons')

        bt.plot_brunt(ftool,method='full',ax=ax3)
        fig.tight_layout()
        fig.savefig(fname)





    def image(self, fname):
        fig,axes=plt.subplots(3,2, figsize=(8,12))
        ax0=axes[0][0];ax1=axes[0][1]
        ax2=axes[1][0];ax3=axes[1][1]
        ax4=axes[2][0];ax5=axes[2][1]

        f1 = self.density.mean(axis=2)
        f2 = self.recovered

        pl=ax0.imshow( f1.transpose())
        fig.colorbar(pl,ax=ax0)
        pl=ax1.imshow(f2.transpose())
        ax0.set_title('Original')
        fig.colorbar(pl,ax=ax1)
        ax1.set_title('Reconstructed')

        pch.simple_phase(f1.flatten(),f2.flatten(),ax=ax2,bins=[16,16])
        ax2.set_title('orig vs recons')

        pch.simple_phase( self.Dx.flatten(), self.Dy.flatten(), ax=ax3, bins=[self.N[0],self.N[1]])
        #ax3.pcolormesh(self.xplane, self.yplane, 
        xf,yf,zf=self.saver[0,:,-1], self.saver[1,:,-1], self.saver[2,:,-1]
        x0,y0,z0=self.saver[0,:,0], self.saver[1,:,0], self.saver[2,:,0]
        #pch.simple_phase( xf.flatten(), yf.flatten(), bins=[64,64],ax=ax3)

        x0.shape=self.N[0],self.N[1]
        y0.shape=self.N[0],self.N[1]
        ax4.pcolormesh(x0,y0,self.Dx, shading='nearest')










        fig.tight_layout()
        fig.savefig(fname)


    def invert(self):
        self.get_dx()
        self.get_dx_proj()
        total=self.zproj.sum()
        rho2 = invert( self.Dx, self.Dy, mean=total).real
        rho2 /= self.density_units*self.gladstone_dale
        rho2 *= 2
        self.recovered = rho2

    def march(self):

        N = nar(self.cube.shape)
        dx = 1/N*self.length_units
        dx.shape=dx.size,1

        Nray=self.rays.shape[1]
        self.eps = np.zeros([2,Nray])


        #generalize later
        self.dim1 = 0
        self.dim2 = 1

        #shorthand
        dim,dim1,dim2=self.dim,self.dim1,self.dim2
        cube=self.cube
        #cube = gaussian_filter(cube,1)

        gx = ddx(np.log(cube),0,dx[dim1])
        gy = ddx(np.log(cube),1,dx[dim2])
        nghost=0

        Nd = nar(gx.shape)
        self.gx=gx
        self.gy=gy

        self.rays[2,:]+=dx[dim]*nghost
        zstep_array=(np.arange(Nd[0])+nghost)*dx[dim]
        ray_shape = self.rays.shape
        self.saver = np.zeros([ray_shape[dim1],ray_shape[dim2],zstep_array.size+1])-1
        self.saver[:,:,0]=self.rays

        


        for nz,zstep in enumerate(zstep_array):
            next_z = zstep+dx[dim]

            marchers = slice(None)
            any_marching = True
            eject=False
                
            xyz = self.rays[:,marchers]
            ijk = (xyz//dx-nghost).astype('int')
            ijk = np.minimum(ijk,Nd[0]-nghost-1)
            ijk = np.maximum(ijk,0)
            this_gx = gx[tuple(ijk)]
            this_gy = gy[tuple(ijk)]

            dz = np.ones_like(this_gx)*dx[dim]
            self.eps[0] += dz*this_gx#+0.5*dz**2*this_dgx_dz
            self.eps[1] += dz*this_gy#+0.5*dz**2*this_dgy_dz
            ok = np.abs(self.eps[1])>0

            this_shiftx = dz*self.eps[0,:]
            this_shifty = dz*self.eps[1,:]

            shift = np.stack([this_shiftx, this_shifty, dz])
            self.rays[:,marchers] = self.rays[:,marchers] + shift


            self.saver[:,:,nz+1]=self.rays
            if eject:
                print('EJECT')
                break






