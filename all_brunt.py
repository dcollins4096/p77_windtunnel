
from dtools.starter1 import *
from downsample import volavg

import simulation
reload(simulation)
import simulation_info.all_sims
import dtools.math.brunt_tools as bt
reload(bt)
import tracer.trace_pc as t1
reload(t1)
import p77_windtunnel.tools as p77tools
reload(p77tools)

import bucket
def plot_all_brunt(sim_list, projax=0, clobber=False):
    if 'ftool' not in bucket.things:
        bucket.things['ftool']={}
    if 'tracer' not in bucket.things:
        bucket.things['tracer']={}
    for nsim,sim in enumerate(sim_list):
        this_sim=simulation.corral[sim]
        this_sim.load()
        frame = this_sim.ann_frames[-1]

        if sim not in bucket.things['ftool'] or clobber:
            rho = this_sim.load_small_rho(frame)
            density_units = 1.204 #kg/m^3
            gladstone_dale = 2.3e-4 #m^3/kg
            #index = rho*density_units*gladstone_dale +1
            N = rho.shape[0]
            rays = t1.get_rays(N, 2, 1)
            tracer = t1.tracer(rho,rays, length_units=1, density_units=density_units,\
                              gladstone_dale=gladstone_dale)
            tracer.march()
            tracer.invert()
            tracer.image('%s/shadowgraph_%s'%(plot_dir,sim))
            rho2 = tracer.recovered

            if 0:
                z = (np.arange(N))/N
                z.shape = 1,1,z.size
                dz = 1/N
                local_zproj = ((1-z)*np.log(index)*dz).sum(axis=2)
                local_zproj /= density_units*gladstone_dale
                local_zproj *= 2
                #rho2 /= gladstone_dale
                #rho2 = (tracer.cube.mean(axis=2))-1
                #rho3 = tracer.zproj
                #rho2 = local_zproj
                #rho3 = ((index-1).mean(axis=2))*128/2
                #rho1 = (rho*density_units*gladstone_dale).mean(axis=2)/2
                rho1 = rho.mean(axis=2)
                #rho3 = rho.mean(axis=2)
                p77tools.ploot3(rho1,rho2,'%s/ab_%s'%(plot_dir,sim), fit=True,unity=True)

            #pdb.set_trace()
            ftool = bt.fft_tool(rho,rho2)
            ftool.do3()
            ftool.do2()
            #print(ftool.ps2.da/ftool.ps3.da)
            bucket.things['ftool'][sim]=ftool
            bucket.things['tracer'][sim]=tracer


    ncol = 3
    nrow = np.ceil(len(sim_list)/3).astype('int')
    nrow = max([nrow,1])
    fig,axes = plt.subplots(nrow,ncol,figsize=(12,12))

    if nrow==1:
        axes=nar([axes])
    for nsim,sim in enumerate(sim_list):
        this_sim=simulation.corral[sim]
        this_sim.load()
        frame = this_sim.ann_frames[-1]
        ftool=bucket.things['ftool'][sim]
        nc = nsim%ncol
        nr = nsim//ncol
        ax = axes[nr][nc]
        bt.plot_brunt(ftool,method='full',ax=ax)
    for ax in axes.flatten():
        ax.set(xticks=[],yticks=[])
    fig.subplots_adjust(wspace=0,hspace=0,left=0,right=1,top=1,bottom=0)
    #fig.tight_layout()
    fig.savefig('%s/all_brunt'%(plot_dir))

def image2(sim_list, projax=0):

    fig,axes = plt.subplots(1,2)
    ax0=axes[0];ax1=axes[1]

    for nsim,sim in enumerate(sim_list):
        ftool=bucket.things['ftool'][sim]
        tracer = bucket.things['tracer'][sim]
        tracer.image2(ftool, '%s/win_image_%s'%(plot_dir,sim))
        #ftool.sigmas_std()

def plot_std(sim_list, projax=0):

    fig,axes = plt.subplots(1,2)
    ax0=axes[0];ax1=axes[1]

    for nsim,sim in enumerate(sim_list):
        this_sim=simulation.corral[sim]
        this_sim.load()
        ftool=bucket.things['ftool'][sim]
        #ftool.sigmas_std()
        bt.sigmas_std(ftool)
        ax0.scatter(ftool.std_x3d, ftool.std_k2dk, c=[this_sim.color], marker=this_sim.marker, s=this_sim.marker_size*100)
        ax0.set(xlabel=r'$\sigma_{3d}$', ylabel='inferred')
        ax1.scatter(this_sim.Ms_mean,1-ftool.std_k2dk/ftool.std_x3d , c=[this_sim.color], marker=this_sim.marker, s=this_sim.marker_size*100)
        ax1.set(xlabel=r'$M_s$', ylabel = r'error')

    fig.tight_layout()
    fig.savefig('%s/all_std'%(plot_dir))



def plot_sigmas(sim_list, projax=0):

    fig,axes = plt.subplots(2,2)

    for nsim,sim in enumerate(sim_list):
        this_sim=simulation.corral[sim]
        this_sim.load()
        ftool=bucket.things['ftool'][sim]
        axes[0][0].scatter(this_sim.Ms_mean,1-ftool.sigma_x3d/ftool.sigma_k3d, c=[this_sim.color], marker=this_sim.marker,s=this_sim.marker_size*100)
        axes[0][0].set(xlabel='Ms',ylabel=r'$1-\sigma_{3x}/\sigma_{3k}$', title='3x vs 3k')
        axes[0][1].scatter(this_sim.Ms_mean,1-ftool.sigma_x2d/ftool.sigma_k2d, c=[this_sim.color], marker=this_sim.marker,s=this_sim.marker_size*100)
        axes[0][1].set(xlabel='Ms',ylabel=r'$1-\sigma_{2x}/\sigma_{2k}$', title='2x vs 2k')
        axes[1][0].scatter(this_sim.Ma_mean, 1-ftool.sigma_k3d/ftool.sigma_k2dk, c=[this_sim.color], marker=this_sim.marker,s=this_sim.marker_size*100)
        axes[1][0].set(xlabel=r'$M_a$',ylabel=r'$1-\sigma_{3k}/\sigma_{k2k}$', title='3k vs k2k')
        axes[1][1].scatter( this_sim.Ma_mean, ftool.ratio_1,c=[this_sim.color], marker=this_sim.marker,s=this_sim.marker_size*100)
        axes[1][1].set(xlabel=r'$M_a$',ylabel='ratio',title='goal')


    fig.tight_layout()
    fig.savefig('%s/all_sigmas'%(plot_dir))


