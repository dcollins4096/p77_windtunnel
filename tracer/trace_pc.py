
from dtools.starter1 import *
from scipy.ndimage import gaussian_filter
from tools import *

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






def ddx(array, direction, dx):
    sl = slice(None)
    ii = slice(1,-1)
    ip1 = slice(2,None)
    im1 = slice(0,-2,None)
    iii = [ii,ii,ii]
    left =  [ii,ii,ii]
    right = [ii,ii,ii]
    left[direction]=ip1
    right[direction]=im1
    out = np.zeros_like(array)
    out[tuple(iii)] = (array[tuple(left)]-array[tuple(right)])/dx
    return out



class tracer():
    def __init__(self,cube,rays,dim=2,length_units=1):
        self.cube1=cube
        self.cube=cube
        self.rays=rays
        self.dim = dim
        self.length_units=length_units
        self.N = cube.shape
        self.xplane, self.yplane, z = get_rays(self.N[0], 2, length_units)
        self.xplane.shape = self.N[0],self.N[1]
        self.yplane.shape = self.N[0],self.N[1]
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

    def get_dx(self):
        xf,yf,zf=self.saver[0,:,-1], self.saver[1,:,-1], self.saver[2,:,-1]
        x0,y0,z0=self.saver[0,:,0], self.saver[1,:,0], self.saver[2,:,0]
        #xx,yy,zz=self.saver[0,:,:], self.saver[1,:,:],self.saver[2,:,:]
        self.Dx = xf-x0
        self.Dy = yf-y0

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
            ijk = np.minimum(ijk,Nd[0]-nghost)
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






