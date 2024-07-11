
from dtools.starter1 import *

def get_cube1(N,ampl=0,val=1.05):
    xyz = np.mgrid[0:N,0:N,0:N]/N
    c = [0.5]*3
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

class tracer():
    def __init__(self,cube,rays,dim=2,length_units=1):
        self.cube1=cube
        self.cube=cube
        self.rays=rays
        self.dim = dim
        self.length_units=length_units
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

    def march(self):

        N = nar(self.cube.shape)
        print(N)
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


        #stelcils for derivatives
        im1 = [slice(1,-1),slice(1,-1),slice(1,-1)]
        ip1 = [slice(1,-1),slice(1,-1),slice(1,-1)]
        jm1 = [slice(1,-1),slice(1,-1),slice(1,-1)]
        jp1 = [slice(1,-1),slice(1,-1),slice(1,-1)]
        km1 = [slice(1,-1),slice(1,-1),slice(1,-1)]
        kp1 = [slice(1,-1),slice(1,-1),slice(1,-1)]

        im1[dim1] = slice(0,-2)
        ip1[dim1] = slice(2,None)
        jm1[dim2] = slice(0,-2)
        jp1[dim2] = slice(2,None)
        km1[dim]  = slice(0,-2)
        kp1[dim]  = slice(2,None)
        ip1[dim] = slice(None)
        im1[dim] = slice(None)
        jp1[dim] = slice(None)
        jm1[dim] = slice(None)
        kp1[dim1]= slice(None)
        km1[dim1]= slice(None)
        kp1[dim2]= slice(None)
        km1[dim2]= slice(None)

        #take derivatives
        gx = (np.log(cube[tuple(ip1)])-np.log(cube[tuple(im1)]))/dx[dim1]
        gy = (np.log(cube[tuple(jp1)])-np.log(cube[tuple(jm1)]))/dx[dim1]
        dgx_dz = (gx[tuple(kp1)] - gx[tuple(km1)])/dx[dim]
        dgy_dz = (gy[tuple(kp1)] - gy[tuple(km1)])/dx[dim]
        gx = gx[:,:,1:-1]
        gy = gy[:,:,1:-1]
        from scipy.ndimage import gaussian_filter
        #gx = gaussian_filter(gx,1)
        #gy = gaussian_filter(gy,1)

        Nd = nar(gx.shape)
        Dx = np.sign(gx)
        Dy = np.sign(gy)
        self.gx=gx
        self.gy=gy

        zstep_array=np.arange(Nd[0])*dx[dim]
        print(dx)
        ray_shape = self.rays.shape
        self.saver = np.zeros([ray_shape[dim1],ray_shape[dim2],zstep_array.size+1])-1
        self.saver[:,:,0]=self.rays

        for nz,zstep in enumerate(zstep_array):
            next_z = zstep+dx[dim]

            marchers = slice(None)
            any_marching = True
            eject=False
            while any_marching:

                
                ijk = (self.rays[:,marchers]//dx-1).astype('int')
                ijk = np.minimum(ijk,Nd[0]-1)
                ijk = np.maximum(ijk,0)
                this_Dx = Dx[tuple(ijk)]
                this_Dy = Dy[tuple(ijk)]
                this_gx = gx[tuple(ijk)]
                this_gy = gy[tuple(ijk)]
                this_dgx_dz = dgx_dz[tuple(ijk)]
                this_dgy_dz = dgy_dz[tuple(ijk)]

                dz = np.ones_like(this_gx)*dx[dim]
                self.eps[0] += dz*this_gx#+0.5*dz**2*this_dgx_dz
                self.eps[1] += dz*this_gy#+0.5*dz**2*this_dgy_dz
                this_shiftx = dz*self.eps[0,:]
                this_shifty = dz*self.eps[1,:]

                shift = np.stack([this_shiftx, this_shifty, dz])
                self.rays[:,marchers] = self.rays[:,marchers] + shift
                #marchers = self.rays[dim] < next_z
                #any_marching = marchers.any()
                any_marching=False
                if any_marching and False:

                    eject=True
                    print("EJECT")
                    pdb.set_trace()
                    break


            self.saver[:,:,nz+1]=self.rays
            if eject:
                print('EJECT')
                break






