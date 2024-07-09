
from dtools.starter1 import *

def get_cube1(N,ampl=0,val=1.05):
    xyz = np.mgrid[0:N,0:N,0:N]
    c = [N//2]*3
    r2 = (xyz[0]-c[0])**2+(xyz[1]-c[1])**2+(xyz[2]-c[2])**2
    out = np.ones([N]*3)
    R = N//4
    out[r2<R**2] = val
    if 1:
        np.random.seed(8675309)
        rando = np.random.random(out.size)*2*ampl-ampl
        rando.shape = out.shape
        out += rando
    return xyz,out
def get_cube2(N,ampl=0,val=1.05,slope=0.3):
    xyz = np.mgrid[0:N,0:N,0:N]-N//2
    c = [N//2]*3
    r2 = (xyz[0]-c[0])**2+(xyz[1]-c[1])**2+(xyz[2]-c[2])**2
    out = np.ones([N]*3)
    R = N//4
    out[xyz[2] > xyz[0]*slope] = val
    if 1:
        np.random.seed(8675309)
        rando = np.random.random(out.size)*2*ampl-ampl
        rando.shape = out.shape
        out += rando
    return xyz,out

class tracer():
    def __init__(self,cube,rays,dim=2):
        self.cube1=cube
        self.cube=cube
        self.rays=rays
        self.dim = dim
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
        N = nar(self.cube1.shape)
        dx = 1/N
        dx.shape=dx.size,1

        Nray=self.rays.shape[1]
        self.eps = np.zeros([2,Nray])


        #generalize later
        self.dim1 = 0
        self.dim2 = 1

        #shorthand
        dim,dim1,dim2=self.dim,self.dim1,self.dim2
        cube=self.cube


        ip1 = [slice(1,-1),slice(1,-1),slice(1,-1)]
        jp1 = [slice(1,-1),slice(1,-1),slice(1,-1)]
        im1 = [slice(1,-1),slice(1,-1),slice(1,-1)]
        jm1 = [slice(1,-1),slice(1,-1),slice(1,-1)]

        ip1[dim1] = slice(2,None)
        jp1[dim2] = slice(2,None)
        im1[dim1] = slice(0,-2)
        jm1[dim2] = slice(0,-2)
        #ip1[dim] = slice(None)
        #im1[dim] = slice(None)
        #jp1[dim] = slice(None)
        #jm1[dim] = slice(None)
        gx = (np.log(cube[tuple(ip1)])-np.log(cube[tuple(im1)]))/dx[dim1]
        gy = (np.log(cube[tuple(jp1)])-np.log(cube[tuple(jm1)]))/dx[dim1]
        Nd = nar(gx.shape)
        Dx = np.sign(gx)
        Dy = np.sign(gy)
        self.gx=gx
        self.gy=gy

        zstep_array=np.linspace(0,1,N[dim]+1)
        ray_shape = self.rays.shape
        self.saver = np.zeros([ray_shape[dim1],ray_shape[dim2],zstep_array.size+1])-1
        self.saver[:,:,0]=self.rays

        for nz,zstep in enumerate(zstep_array):
            print('nz',nz)
            next_z = zstep+dx[dim]

            marchers = slice(None)
            any_marching = True
            eject=False
            while any_marching:

                ijk = (self.rays[:,marchers]//dx).astype('int')
                ijk = np.minimum(ijk,Nd[0]-1)
                ijk = np.maximum(ijk,0)
                this_Dx = Dx[tuple(ijk)]
                this_Dy = Dy[tuple(ijk)]
                this_gx = gx[tuple(ijk)]
                this_gy = gy[tuple(ijk)]

                if 0:
                    #print('ok',self.rays[dim1]*16)
                    to_edge=(dx[dim1]*(ijk[dim1]+0.5*(1-this_Dx)))

                    deltax = np.abs(self.rays[dim1][marchers] - dx[dim1]*(ijk[dim1]+0.5*(1+this_Dx)))
                    deltay = np.abs(self.rays[dim2][marchers] - dx[dim2]*(ijk[dim2]+0.5*(1+this_Dy)))

                    zcx = np.zeros_like(deltax)+dx[dim1]
                    ok = np.abs(this_Dx)>0
                    zcx[ok] = np.sqrt(deltax[ok]/np.abs(this_gx[ok]))
                    zcy = np.zeros_like(deltay)+dx[dim2]
                    ok = np.abs(this_Dy)>0
                    zcy[ok] = np.sqrt(deltax[ok]/np.abs(this_gx[ok]))
                    dz = np.zeros_like(deltax)+dx[dim1]
                    dz = np.minimum(zcx, dz)
                    dz = np.minimum(zcy, dz)
                dz = np.ones_like(this_gx)*dx[dim]
                self.eps[0] += dz*this_gx
                self.eps[1] += dz*this_gy
                print(ijk[0])
                this_shiftx = dz*self.eps[0,:]
                this_shifty = dz*self.eps[1,:]

                shift = np.stack([this_shiftx, this_shifty, dz])
                self.rays[:,marchers] = self.rays[:,marchers] + shift
                marchers = self.rays[dim] < next_z
                any_marching = marchers.any()
                if any_marching:

                    eject=True
                    print("EJECT")
                    break


            self.saver[:,:,nz+1]=self.rays
            if eject:
                print('EJECT')
                break






