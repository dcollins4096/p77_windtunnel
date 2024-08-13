

from dtools.starter1 import *

import dtools.davetools as dt
import tracer.trace_pc as t1
reload(t1)

import dtools.vis.pcolormesh_helper as pch
from tools import *
plt.close('all')

def z_weight(cube):
    N = cube.shape
    z = (np.arange(N[2])+0.5)/N[2]
    z.shape = 1,1,z.size
    z=1
    return (z*cube).sum(axis=2)


for zzz in [4]:
    N = 128
    rho1=1
    rho2=1.1
    xyz = [3,3,zzz]
    length_units = 1.0 #m
    density_units = 1.204 #kg/m^3
    gladstone_dale = 2.3e-4 #m^3/kg
    #cube = t1.get_cube_impulse(N, xyz, rho1=rho1,rho2=rho2)
    cube = t1.get_cube_128(density_units, gladstone_dale)
    #cube = t1.get_cubesin(N)
    rays = t1.get_rays(N, 2, length_units)
    tracer = t1.tracer(cube,rays, length_units=length_units)
    tracer.march()
    tracer.get_dx()
    tracer.get_spectral()

    tracer.get_dx_proj()
    total=tracer.zproj.sum()

    rho_spect = invert( tracer.spsx, tracer.spsy, mean = total)
    ploot2( tracer.zproj, rho_spect.real,'%s/invert_impulse_test_spectral'%plot_dir)

    ploot2( tracer.spsx.real, tracer.dpdx, '%s/spectral_vs_centered'%plot_dir)
    if 0:
        ok = np.abs(tracer.dpdx.real) > 0
        vert = tracer.spsx.real[ok]/tracer.dpdx[ok]
        plt.clf()
        plt.hist(vert.flatten())
        plt.savefig('%s/spec_to_real'%plot_dir)
    rho3 = invert( tracer.dpdx, tracer.dpdy, mean = total)
    ploot2( tracer.zproj, rho3.real,'%s/invert_impulse_test1'%plot_dir)
    fig,ax=plt.subplots(1,1)
    import dtools.math.equal_probability_binner as epb
    ok = np.abs(tracer.zproj)>0
    xxxx= (rho3.real[ok]/tracer.zproj[ok]).flatten()
    if 0:
        ax.hist(xxxx)
        xxxx.sort()
        yyyy=np.arange(xxxx.size)/xxxx.size
        tax=ax.twinx()
        tax.plot(xxxx,yyyy,c='k')
    stuff=epb.equal_prob(xxxx,16,ax)

    fig.savefig('%s/recon_over_real'%plot_dir)
    #print(tracer.Dx[tracer.Dx>0])
    #print(tracer.dpdx[tracer.dpdx>0])
    rho2 = invert( tracer.Dx, tracer.Dy, mean=total)

    ploot2(tracer.dpdx, tracer.Dx, '%s/Dx_vs_dx_impulse'%plot_dir)
    ploot2(tracer.zproj, rho2, '%s/invert_impulse'%plot_dir)
    one = (tracer.cube.mean(axis=2)-1)/(2*(1/N))
    ploot2(one,rho2,'%s/how_far_off'%plot_dir)
    fig2,ax47=plt.subplots(1,2)
    xxx = one.flatten()
    yyy = rho2.flatten()
    ok = np.abs(xxx)>0
    rat = yyy[ok]/xxx[ok]
    print(rat.mean())
    

    
