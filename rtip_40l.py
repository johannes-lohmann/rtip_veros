import matplotlib.pyplot as pl
import numpy as np
import xarray as xr
from eofs.xarray import Eof
from matplotlib.collections import LineCollection
import matplotlib as mpl
from matplotlib import gridspec
from scipy.interpolate import UnivariateSpline

pl.rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'size'   : 17})

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

mpl.rcParams['axes.xmargin'] = 0.03
mpl.rcParams['axes.ymargin'] = 0.03

mpl.rcParams['axes.unicode_minus'] = False

def MAIN():

        ot_hyst = xr.open_dataarray("overturning_s40l360q.nc")

        t_amoc_hyst = []; amoc_hyst = []; amoc_hyst360q = [ot_hyst.isel(time=0, depth=slice(0,32), lat=slice(24,40)).max()/1000000.]

        for i in range(len(ot_hyst.time)):
                t0 = ot_hyst.time[i]/360.
                t_amoc_hyst.append(t0)
                amoc_hyst.append(ot_hyst.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                if (t0-5500)%300==0:
                        amoc_vals = np.asarray([ot_hyst.isel(time=i-j, depth=slice(0,32), lat=slice(24,40)).max()/1000000. for j in range(5)])
                        amoc_hyst360q.append(np.mean(amoc_vals))

        F_hyst = np.linspace(0.85694,2.0567,29)
        F_hyst = np.concatenate([F_hyst, np.linspace(2.0138,0.7712,30)])
        
        
        warming(F_hyst, amoc_hyst360q)
        cooling(F_hyst, amoc_hyst360q)
        initial_cond(F_hyst, amoc_hyst360q)
        test_eof()
        edge_state()
        clustering()
        test_random_sequence()
        
        pl.show()

def edge_state():

        ot_c = xr.open_dataarray("overturning_rtipSc0.nc")
        t_amoc_c, amoc_c = amoc_timeseries(ot_c)
        t_amoc_c, aabw_c = aabw_timeseries(ot_c)
        ot_c01 = xr.open_dataarray("overturning_rtipSc01.nc")
        t_amoc_c01, amoc_c01 = amoc_timeseries(ot_c01)
        t_amoc_c01, aabw_c01 = aabw_timeseries(ot_c01)
        ot_c1 = xr.open_dataarray("overturning_rtipSc1.nc")
        t_amoc_c1, amoc_c1 = amoc_timeseries(ot_c1)
        t_amoc_c1, aabw_c1 = aabw_timeseries(ot_c1)
        t_amoc_c1 = t_amoc_c1[210:]
        amoc_c1 = amoc_c1[210:]
        aabw_c1 = aabw_c1[210:]
        ot_c006 = xr.open_dataarray("overturning_rtipSc006.nc")
        salt_temp006 = xr.open_dataset("salt_temp_rtipSc006.nc")
        t_amoc_c006, amoc_c006 = amoc_timeseries(ot_c006)
        t_amoc_c006, aabw_c006 = aabw_timeseries(ot_c006)

        fig=pl.figure()
        pl.subplot(211)
        pl.plot(t_amoc_c, amoc_c)
        pl.plot(t_amoc_c01, amoc_c01)
        pl.plot(t_amoc_c1, amoc_c1)
        pl.plot(t_amoc_c006, amoc_c006)
        pl.subplot(212)
        pl.plot(t_amoc_c, aabw_c)
        pl.plot(t_amoc_c01, aabw_c01)
        pl.plot(t_amoc_c1, aabw_c1)
        pl.plot(t_amoc_c006, aabw_c006)

        fig=pl.figure()
        pl.plot(amoc_c, aabw_c)
        pl.plot(amoc_c01, aabw_c01)
        pl.plot(amoc_c1, aabw_c1)
        pl.plot(amoc_c006, aabw_c006)

        fig=pl.figure()
        pl.subplot(311)
        pl.plot(t_amoc_c006, amoc_c006)
        pl.subplot(312)
        (salt_temp006['sst_SA']-salt_temp006['sst_NA']).plot()
        pl.subplot(313)
        (salt_temp006['sss_NA']-salt_temp006['sss_SA']).plot()

def initial_cond(F_hyst, amoc_hyst360q):

        ### 600 year spinup; 5y averages
        ot_init0 = xr.open_dataarray("overturning_rtipSc0.nc")
        salt_temp_init0 =  xr.open_dataset("salt_temp_rtipSc0.nc")

        ### 100 year branch of spinup; 2y averages
        ot_init = xr.open_dataarray("overturning_rtipSc0i.nc")
        salt_temp_init =  xr.open_dataset("salt_temp_rtipSc0i.nc")

        ### 10 year branch of 100y spinup; 2-month averages
        ot_initZ = xr.open_dataarray("overturning_rtipzSp.nc")
        salt_temp_initZ =  xr.open_dataset("salt_temp_rtipzSp.nc")
        t_amoc_initZ, amoc_initZ = amoc_timeseries(ot_initZ)

        N = 46
        n_lines = N+1

        print(salt_temp_init)

        t_amoc_init1, amoc_init1 = amoc_timeseries(ot_init0)
        k = len(amoc_init1)        
        t_amoc_init0, amoc_init0 = amoc_timeseries(ot_init)
        amoc_init = amoc_init1[:120] + amoc_init0
        t_amoc_init = t_amoc_init1[:120] + t_amoc_init0
        amoc_init= amoc_init[119:]
        t_amoc_init= t_amoc_init[119:]

        t_aabw_init1, aabw_init1 = aabw_timeseries(ot_init0)
        t_aabw_init0, aabw_init0 = aabw_timeseries(ot_init)
        aabw_init = aabw_init1 + aabw_init0
        t_aabw_init = t_aabw_init1 + t_aabw_init0
        aabw_init= aabw_init[k-1:]
        t_aabw_init= t_aabw_init[k-1:]

        '''
        ###include original initial condition
        sst_pc1 = salt_temp_init['sst_pc1'].values
        sst_pc2 = salt_temp_init['sst_pc2'].values
        sst_pc3 = salt_temp_init['sst_pc3'].values
        sst_pc4 = salt_temp_init['sst_pc4'].values
        sst_pc5 = salt_temp_init['sst_pc5'].values

        temp = salt_temp_init['temp'].values

        fig=pl.figure()
        #pl.plot(amoc_init[1:], sst_pc2, color='black')
        pl.plot(amoc_init[1:], temp, color='black')
        #pl.plot(amoc_init[1:], aabw_init[1:], color='black')
        ax = pl.axes()
        ax.set_prop_cycle([pl.cm.viridis(i) for i in np.linspace(0, 1, N)])
        for i in range(N):
                #pl.plot(amoc_init[1+i], sst_pc2[i], 'o')
                pl.plot(amoc_init[1+i], temp[i], 'o')
                #pl.plot(amoc_init[1+i], aabw_init[1+i], 'o')

        fig=pl.figure()
        ax=pl.subplot(911)
        pl.plot(t_amoc_init, amoc_init)
        ax.set_prop_cycle([pl.cm.viridis(i) for i in np.linspace(0, 1, n_lines)])
        for i in range(N+1):
                pl.plot(t_amoc_init[i], amoc_init[i], 'o')
        pl.subplot(912)
        (salt_temp_init['sst_SA']-salt_temp_init['sst_NA']).plot()
        #salt_temp_init['sst_SA'].plot()
        pl.subplot(913)
        (salt_temp_init['sss_NA']-salt_temp_init['sss_SA']).plot()
        #salt_temp_init['sss_SA'].plot()
        pl.subplot(914)
        salt_temp_init['seaice'].plot()
        pl.subplot(915)
        salt_temp_init['sst_pc1'].plot()
        pl.subplot(916)
        salt_temp_init['sst_pc2'].plot()
        pl.subplot(917)
        salt_temp_init['sst_pc3'].plot()
        pl.subplot(918)
        salt_temp_init['sst_pc4'].plot()
        pl.subplot(919)
        #salt_temp_init['sst_pc5'].plot()
        pl.plot(t_aabw_init, aabw_init)


        ot0 = xr.open_dataarray("overturning_rtipSc007.nc")
        t_amoc0, amoc0, aabw0 = amoc_aabw_timeseries(ot0)
        t_amoc = [t_amoc0]; amoc = [amoc0]; aabw = [aabw0]
        if amoc0[-1]<4.:
                tip=[1]
        else:
                tip=[0]


        for i in range(N):
                ot = xr.open_dataarray("overturning_rtipSc007i%s.nc"%i)
                t_amoc0, amoc0, aabw0 = amoc_aabw_timeseries(ot)
                t_amoc.append(t_amoc0)
                amoc.append(amoc0); aabw.append(aabw0)
                if amoc0[-1]<4.:
                        tip.append(1)
                else:
                        tip.append(0)

        fig=pl.figure()
        ax = pl.axes()
        ax.set_prop_cycle([pl.cm.viridis(i) for i in np.linspace(0, 1, n_lines)])
        for i in range(N+1):
                #fig=pl.figure()
                pl.plot(t_amoc[i], amoc[i])

        #fig=pl.figure()
        #ax = pl.axes()
        #ax.set_prop_cycle([pl.cm.plasma(i) for i in np.linspace(0, 1, n_lines)])
        #for i in range(N+1):
        #        pl.plot(t_amoc[i], aabw[i])


        fig=pl.figure()
        for i in range(N):
                amoc0 = amoc[i]; aabw0 = aabw[i]
                if (i==6 or i==12 or i==18 or i==23 or i==27 or i==32 or i==41 or i==42):
                        pl.plot(amoc0[14], aabw0[14], '<', color='gray') 
                elif tip[i]==0:
                        pl.plot(amoc0[14], aabw0[14], 's', color='crimson') 
                else:
                        pl.plot(amoc0[14], aabw0[14], 'o', color='black') 
                        

        print tip

        def make_fig(t_amoc, amoc, i, width=6.5):
                r0 = np.array(len(amoc)*[i])
                pl.plot(np.asarray(t_amoc)-t_amoc[0], r0, alpha=0.0)
                points2 = np.array([np.asarray(t_amoc)-9100., r0]).transpose().reshape(-1,1,2)#-t_amoc[0]
                segs2 = np.concatenate([points2[:-1],points2[1:]],axis=1)
                lc2 = LineCollection(segs2, cmap=pl.get_cmap('Spectral'),  norm = mpl.colors.Normalize(vmin=2.,vmax=7.4)
)
                lc2.set_array(np.asarray(amoc))
                lc2.set_linewidth(width)
                pl.gca().add_collection(lc2)
                if i==45:
                        ax2 = fig.colorbar(lc2)
                        ax2.set_label('AMOC max. (Sv)')

        fig=pl.figure(figsize=(12,5))
        pl.subplots_adjust(left=0.1, bottom=0.14, right=0.97, top=0.98, wspace=0.02, hspace=0.0)

        spec = gridspec.GridSpec(ncols=2, nrows=1,
                         width_ratios=[1, 5])

        ax0 = fig.add_subplot(spec[0])
        pl.plot(amoc_init[:N], t_amoc_init[:N], 'x-', color='black')
        pl.ylabel('Branch-off time (years)'); pl.xlabel('AMOC max. (Sv)')

        ax1 = fig.add_subplot(spec[1])

        for i in range(N):
                t_amoc1 = t_amoc[i]; amoc1 = amoc[i]
                t_amoc0 = np.concatenate((t_amoc_init1[100:119], t_amoc_init[0:i+1], t_amoc1))
                amoc0 = np.concatenate((amoc_init1[100:119], amoc_init[0:i+1], amoc1))
                make_fig(t_amoc0, amoc0, i)
                pl.plot([(i+1)*2.]*50, np.linspace(i-0.5, i+0.5, 50), color='black')

                print t_amoc0

                #make_fig(t_amoc[i], amoc[i], i)
        ax1.axes.get_yaxis().set_visible(False)
        ax1.spines['left'].set_visible(False)

        pl.xlabel('Simulation time (years)');pl.xlim(-100.,700.)
        '''

        M = 60

        edge_reals = [6, 11, 14, 15, 17, 19, 20, 23, 24, 29, 35, 40, 46, 54, 57]
        tip = []; X=[]; t_amoc=[]; Y = []; amoc = []; aabw = []
        for i in range(M):
                ot = xr.open_dataarray("overturning_rtipz%s.nc"%i)
                st = xr.open_dataset("salt_temp_rtipz%s.nc"%i)
                t_amoc0, amoc0, aabw0 = amoc_aabw_timeseries(ot)
                Y.append(st['rho_sub_NA'].values-st['rho_sub_SA'].values)
                t_amoc.append(t_amoc0)
                amoc.append(amoc0)
                X.append(st['rho_sub_SA'].values)
                aabw.append(aabw0)
                if amoc0[-1]<4.:
                        tip.append(1)
                else:
                        tip.append(0)

        print(st)

        fig=pl.figure(figsize=(7.5,6))
        pl.subplots_adjust(left=0.15, bottom=0.12, right=0.97, top=0.98, wspace=0.0, hspace=0.2)
        ax=pl.subplot(411)
        pl.plot(t_amoc_initZ, amoc_initZ); pl.ylabel('AMOC max.')
        ax.set_prop_cycle([pl.cm.plasma(i) for i in np.linspace(0, 1, M)])
        for i in range(M):
                pl.plot(t_amoc_initZ[i], amoc_initZ[i], 'o')
        pl.subplot(412)
        (salt_temp_initZ['sst_SA']-salt_temp_initZ['sst_NA']).plot()
        pl.subplot(413)
        (salt_temp_initZ['sss_NA']-salt_temp_initZ['sss_SA']).plot()
        pl.subplot(414)
        salt_temp_initZ['seaice'].plot()


        #warm state final
        ot_c1 = xr.open_dataarray("overturning_rtipSc016.nc")
        st_c1 = xr.open_dataset("salt_temp_rtipSc016.nc")
        t_amoc_c1, amoc_c1, aabw_c1 = amoc_aabw_timeseries(ot_c1)
        XW = st_c1['rho_sub_SA'].values; XW = XW[300:]
        YW = st_c1['rho_sub_NA'].values-st_c1['rho_sub_SA'].values; YW=YW[300:]

        #warm state initial
        XI = salt_temp_init0['rho_sub_SA'].values
        YI = salt_temp_init0['rho_sub_NA'].values-salt_temp_init0['rho_sub_SA'].values

        #edge state
        ot_c006 = xr.open_dataarray("overturning_rtipSc006.nc")
        salt_temp006 = xr.open_dataset("salt_temp_rtipSc006.nc")
        t_amoc_c006, amoc_c006, aabw_c006 = amoc_aabw_timeseries(ot_c006)
        XE = salt_temp006['rho_sub_SA'].values; XE=XE[200:]
        YE = salt_temp006['rho_sub_NA'].values-salt_temp006['rho_sub_SA'].values; YE=YE[200:]

        fig=pl.figure(figsize=(8,5))
        pl.subplots_adjust(left=0.08, bottom=0.14, right=0.97, top=0.98, wspace=0.0, hspace=0.2)
        ax = pl.axes()
        pl.plot(t_amoc_init1, amoc_init1, color='black')
        pl.plot(t_amoc_init, amoc_init, color='black')
        ax.set_prop_cycle([pl.cm.plasma(i) for i in np.linspace(0, 1, M)])
        for i in range(M):
                if i in edge_reals:
                        pl.plot(t_amoc[i][1:], amoc[i][1:],'--', linewidth=0.7)
                else:
                        pl.plot(t_amoc[i][1:], amoc[i][1:], linewidth=0.7)
        pl.ylabel('AMOC max. (Sv)'); pl.xlabel('Simulation time (years)')


        fig=pl.figure()
        ax = pl.axes()
        ax.set_prop_cycle([pl.cm.plasma(i) for i in np.linspace(0, 1, M)])
        for i in range(M):
                if i in edge_reals:
                        pl.plot(t_amoc[i], X[i],'--')
                else:
                        pl.plot(t_amoc[i], X[i])
        pl.plot(t_amoc_init1, XI, color='black')
        pl.plot(t_amoc_c1[300:], XW, color='crimson')
        pl.plot(t_amoc_c006[200:], XE, color='gray')

        fig=pl.figure()
        ax = pl.axes()
        ax.set_prop_cycle([pl.cm.plasma(i) for i in np.linspace(0, 1, M)])
        for i in range(M):
                pl.plot(t_amoc[i], Y[i])
        pl.plot(t_amoc_init1, YI, color='black')
        pl.plot(t_amoc_c1[300:], YW, color='crimson')
        pl.plot(t_amoc_c006[200:], YE, color='gray')
        
        fig=pl.figure()
        ax = pl.axes()
        ax.set_prop_cycle([pl.cm.plasma(i) for i in np.linspace(0, 1, M)])
        pl.plot(XI, YI, color='black')
        pl.plot(XW, YW, color='crimson')
        pl.plot(XE, YE, color='gray')
        
        for i in range(M):
                X0 = X[i]; Y0 = Y[i]
                pl.plot(X0[1], Y0[1], 'o', markersize=10, markerfacecolor='None', color='orange')
                if i in edge_reals:
                        tip[i] = 0.5
                        pl.plot(X0[14], Y0[14], '<', color='gray') 
                elif tip[i]==0:
                        pl.plot(X0[14], Y0[14], 's', color='crimson') 
                else:
                        pl.plot(X0[14], Y0[14], 'o', color='black') 
                pl.plot(X0[-1], Y0[-1], 'o', markersize=10, color='green')

        fig=pl.figure()
        for i in range(M):
                X0 = X[i]; Y0 = Y[i]
                if i in edge_reals:
                        tip[i] = 0.5
                        pl.plot(X0[14], Y0[14], '<', color='gray') 
                elif tip[i]==0:
                        pl.plot(X0[14], Y0[14], 's', color='crimson') 
                else:
                        pl.plot(X0[14], Y0[14], 'o', color='black') 

        print(tip)

        def make_fig(t_amoc, amoc, i, width=6.5):
                r0 = np.array(len(amoc)*[i])
                pl.plot(np.asarray(t_amoc)-t_amoc[0], r0, alpha=0.0)
                points2 = np.array([np.asarray(t_amoc)-9100., r0]).transpose().reshape(-1,1,2)
                segs2 = np.concatenate([points2[:-1],points2[1:]],axis=1)
                lc2 = LineCollection(segs2, cmap=pl.get_cmap('Spectral'),  norm = mpl.colors.Normalize(vmin=2.,vmax=7.4)
)
                lc2.set_array(np.asarray(amoc))
                lc2.set_linewidth(width)
                pl.gca().add_collection(lc2)
                if i==59:
                        ax2 = fig.colorbar(lc2)
                        ax2.set_label('AMOC max. (Sv)')

        fig=pl.figure(figsize=(12,5))
        pl.subplots_adjust(left=0.1, bottom=0.14, right=0.97, top=0.98, wspace=0.02, hspace=0.0)

        spec = gridspec.GridSpec(ncols=2, nrows=1,
                         width_ratios=[1, 5])

        ax0 = fig.add_subplot(spec[0])
        pl.plot(amoc_initZ, t_amoc_initZ, 'x-', color='black')
        pl.ylabel('Branch-off time (years)'); pl.xlabel('AMOC max. (Sv)')

        ax1 = fig.add_subplot(spec[1])
        cut_idcs = [23*[5], 30*[8], 7*[10]]; cut_idcs = [item for sublist in cut_idcs for item in sublist]
        amoc_sdev = []; t_amoc_sdev = []
        for i in range(M):
                t_amoc1 = t_amoc[i]; amoc1 = amoc[i]

                t_amoc0 = np.concatenate((t_amoc_init1[100:119], t_amoc_init[0:cut_idcs[i]], t_amoc1))
                amoc0 = np.concatenate((amoc_init1[100:119], amoc_init[0:cut_idcs[i]], amoc1))

                amoc_sdev.append(amoc1[1:]); t_amoc_sdev.append(t_amoc1[1:])
                
                make_fig(t_amoc0, amoc0, i)
                pl.plot([(i+1)/6.]*50, np.linspace(i-0.5, i+0.5, 50), color='black')

        ax1.axes.get_yaxis().set_visible(False)
        ax1.spines['left'].set_visible(False)
        pl.xlim(-100,500.)

        pl.xlabel('Simulation time (years)')

        amoc_sdev_new = []
        fig=pl.figure()
        for i in range(M):

                if i<23:
                        amoc_test = amoc_sdev[i]

                elif 22<i<53:
                        amoc_test = np.concatenate(([np.nan],amoc_sdev[i][:-1]))

                else:
                        amoc_test = np.concatenate(([np.nan,np.nan],amoc_sdev[i][:-2]))

                print(len(amoc_test), amoc_test[0:3], t_amoc_sdev[i][0])
                pl.plot(t_amoc_sdev[0],amoc_test)
                amoc_sdev_new.append(amoc_test)
        
        sdev = np.asarray([np.nanstd(np.asarray(amoc_sdev_new)[:,i]) for i in range(len(amoc_sdev_new[0]))])

        for i in range(len(sdev)):
                print(t_amoc_sdev[0][i]-9106., sdev[i])

        fig=pl.figure(figsize=(12,6.5))
        pl.subplot(211)
        pl.plot(np.asarray(t_amoc_sdev[0])-9106., sdev)
        pl.plot(100*[76.],np.linspace(min(sdev),max(sdev),100), color='black')
        pl.plot(100*[86.],np.linspace(min(sdev),max(sdev),100), color='black')
        pl.plot(100*[160.],np.linspace(min(sdev),max(sdev),100), color='black')
        pl.ylabel('Ensemble Stand. Dev.'); pl.xlim(-100,500)
        pl.subplot(212)
        pl.plot(np.asarray(t_amoc_sdev[0])-9106., sdev)
        pl.plot(100*[76.],np.linspace(min(sdev),max(sdev),100), color='black')
        pl.plot(100*[86.],np.linspace(min(sdev),max(sdev),100), color='black')
        pl.plot(100*[160.],np.linspace(min(sdev),max(sdev),100), color='black')
        pl.yscale('log')
        pl.ylabel('Ensemble Stand. Dev.'); pl.xlabel('Simulation time (years)')
        pl.xlim(-100,500)


def test_random_sequence():

        ### 70 years zoom.
        #tip = [0, 0, 1, 0, 0, 1, 0.5, 0, 1, 0, 0, 0.5, 1, 0, 0.5, 0.5, 1, 0.5, 1, 0.5, 0.5, 1, 1, 0.5, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 1, 0, 0, 0.5, 1, 0, 0, 1, 0.5, 0, 0, 0, 0, 1, 0.5, 0, 0, 1, 1, 0, 0, 1, 0.5, 0, 0, 0.5, 1, 1]

        ### instant zoom.
        #tip = [0, 0.5, 0, 0, 1, 0.5, 0, 0.5, 1, 0.5, 0.5, 0, 0.5, 0, 0, 1, 0.5, 1, 0.5, 1, 0, 1, 1, 1, 1, 0.5, 0.5, 0, 0, 0.5, 0.5, 0.5, 0, 0, 0.5, 1, 1, 0, 0.5, 0.5, 1, 1, 0.5, 0.5, 1, 0.5, 1, 0, 1, 0.5, 0.5, 0, 1, 0.5, 0.5, 1, 1, 0.5, 1, 0.5]

        ### instant zoom 2
        tip = [0.5, 0, 0.5, 0, 0.5, 0, 0.5, 0, 1, 1, 0, 0, 0.5, 0, 0.5, 0.5, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0.5, 1, 0.5, 0.5, 1, 0, 0,  0.5, 1, 1, 0.5, 1, 0.5, 0, 0, 0.5, 1, 0]

        M = len(tip)
        N = 1000000

        def calc_p(tip):
                p = np.zeros(3)
                for i in xrange(M):
                        if tip[i]==0:
                                p[0]+=1
                        elif tip[i]==0.5:
                                p[1]+=1
                        else:
                                p[2]+=1
                p = p/M
                return p

        def calc_joint_p(tip):
                p_joint = np.zeros((3,3))
                for i in xrange(1,M):
                        if (tip[i]==0 and tip[i-1]==0): p_joint[0,0]+=1
                        elif (tip[i]==0.5 and tip[i-1]==0.5): p_joint[1,1]+=1
                        elif (tip[i]==1 and tip[i-1]==1): p_joint[2,2]+=1
                        elif (tip[i]==0 and tip[i-1]==0.5): p_joint[0,1]+=1
                        elif (tip[i]==0 and tip[i-1]==1): p_joint[0,2]+=1
                        elif (tip[i]==0.5 and tip[i-1]==0): p_joint[1,0]+=1
                        elif (tip[i]==0.5 and tip[i-1]==1): p_joint[1,2]+=1
                        elif (tip[i]==1 and tip[i-1]==0): p_joint[2,0]+=1
                        elif (tip[i]==1 and tip[i-1]==0.5): p_joint[2,1]+=1
                idcs = np.where(p_joint==0)
                p_joint[idcs] = 0.00000001

                p_joint = p_joint/M
                p_joint = p_joint/np.sum(p_joint)
                return p_joint

        def realization(p):
                real=np.empty(M)
                rands = np.random.uniform(size=M)
                for i in xrange(M):
                        if rands[i]<p[0]:
                                real[i] = 0
                        elif p[0]<rands[i]<p[0]+p[1]:
                                real[i] = 0.5
                        elif p[0]+p[1]<rands[i]:
                                real[i] = 1
                return real

        def calc_entr(p):
                return -np.sum([p[i]*np.log(p[i]) for i in range(3)])

        def calc_cond_entr(p, p_joint):
                return -np.sum([[p_joint[l,j]*np.log(p_joint[l,j]/p[j]) for l in range(3)] for j in range(3)])

        p_data = calc_p(tip)
        p_joint_data = calc_joint_p(tip)

        print(p_data)
        print(p_joint_data)

        entr_data = calc_entr(p_data)
        joint_entr_data = calc_cond_entr(p_data, p_joint_data)
        print('Entropy', entr_data)
        print('Conditional Entropy', joint_entr_data)

        joints = np.empty(N)
        for k in range(N):
                real1 = realization(p_data)
                p1 = calc_p(real1)
                p_joint1 = calc_joint_p(real1)
                bla = calc_cond_entr(p1, p_joint1)
                joints[k] = bla

        data_statistic = joint_entr_data
        statistics = joints
        statistics_sort = np.sort(statistics)
        p_index = np.argmin(np.abs(statistics_sort - data_statistic))
        pvalue_precise = 1-float(p_index)/(len(statistics))

        print('p-value: ',pvalue_precise)

        bin_no=100
        hist, bins = np.histogram(statistics, bins=bin_no, density=True)
        binwidth = bins[1] - bins[0]
        width = 0.85 * binwidth
        center = (bins[:-1] + bins[1:]) / 2
        
        hist_p, center_p, pvalue = calc_pvalue(statistics,hist,center,binwidth,data_statistic)
        fig2=pl.figure()
        pl.title('p = %s'%pvalue_precise)
        ax1 = pl.subplot(111)
        pl.bar(center, hist, align='center', width=width, alpha=0.7,label='simulation')
        pl.bar(center_p, hist_p, align='center', width=width, alpha=0.7,color='gold')
        pl.plot(50*[data_statistic],np.linspace(0,max(hist),50),'--', color='crimson',linewidth=2.0,label='data')
        pl.xlabel('Conditional Entropy')
        pl.ylabel('P')


def clustering():
        
        start_idcs = [23*[3], 30*[2], 7*[1]]; start_idcs = [item for sublist in start_idcs for item in sublist]

        M = 60
        amoc = []; aabw = []
        edge_reals = [6, 11, 14, 15, 17, 19, 20, 23, 24, 29, 35, 40, 46, 54, 57]

        t_amoc=[];amoc=[]; eq=[]; aabw=[]
        for i in range(M):
                ot = xr.open_dataarray("overturning_rtipz%s.nc"%i)
                t_amoc0, amoc0, eq0, aabw0 = amoc_aabw_timeseries_ext(ot)
                amoc.append(amoc0);eq.append(eq0);aabw.append(aabw0);t_amoc.append(t_amoc0)

        variables = ['salt_sub_NA', 'salt_sub_SA', 'temp_sub_NA', 'temp_sub_SA', 'sst_NA', 'sst_SA', 'sss_NA', 'sss_SA', 'rho_sub_NA', 'rho_sub_SA', 'rho_NA', 'rho_SA', 'seaice', 'AMOC', 'AABW', 'eq', 'delta_rho', 'delta_rho_sub', 'delta_sst', 'delta_temp_sub', 'delta_salt_sub']

        variables=['rho_SA','rho_sub_NA','salt_sub_SA']

        #warm state final
        ot_c1 = xr.open_dataarray("overturning_rtipSc016.nc")
        st_c1 = xr.open_dataset("salt_temp_rtipSc016.nc")
        t_amoc_c1, amoc_c1, eq_c1, aabw_c1 = amoc_aabw_timeseries_ext(ot_c1)

        #edge state
        ot_c006 = xr.open_dataarray("overturning_rtipSc006.nc")
        salt_temp006 = xr.open_dataset("salt_temp_rtipSc006.nc")
        t_amoc_c006, amoc_c006, eq_c006, aabw_c006 = amoc_aabw_timeseries_ext(ot_c006)

        #warm state initial: spinup; 5y averages
        ot_init0 = xr.open_dataarray("overturning_rtipSc0.nc")
        salt_temp_init0 =  xr.open_dataset("salt_temp_rtipSc0.nc")
        t_amoc_init0, amoc_init0, eq_init0, aabw_init0 = amoc_aabw_timeseries_ext(ot_init0)


        for j in range(len(variables)):
                fig=pl.figure()
                pl.subplots_adjust(left=0.1, bottom=0.14, right=0.97, top=0.98, wspace=0.02, hspace=0.0)
                X=np.empty(M); Y=np.empty(M)
                for i in range(M):
                        st = xr.open_dataset("salt_temp_rtipz%s.nc"%i)
                        t_amoc0=t_amoc[i]; amoc0=amoc[i]
                        if variables[j]=='eq':
                                X0 = eq[i]
                                XE = eq_c006; XE=XE[200:]
                                XW = eq_c1; XW=XW[300:]
                                XI = eq_init0
                        elif variables[j]=='AMOC':
                                X0 = amoc[i]
                                XE = amoc_c006; XE=XE[200:]
                                XW = amoc_c1; XW=XW[300:]
                                XI = amoc_init0
                        elif variables[j]=='AABW':
                                X0 = aabw[i]
                                XE = aabw_c006; XE=XE[200:]
                                XW = aabw_c1; XW=XW[300:]
                                XI = aabw_init0
                        elif variables[j]=='delta_rho':
                                X0 = st['rho_NA'].values-st['rho_SA'].values
                                XE = salt_temp006['rho_NA'].values-salt_temp006['rho_SA'].values; XE=XE[200:]
                                XW = st_c1['rho_NA'].values-st_c1['rho_SA'].values; XW=XW[300:]
                                XI = salt_temp_init0['rho_NA'].values-salt_temp_init0['rho_SA'].values

                        elif variables[j]=='delta_sst':
                                X0 = st['sst_NA'].values-st['sst_SA'].values
                                XE = salt_temp006['sst_NA'].values-salt_temp006['sst_SA'].values; XE=XE[200:]
                                XW = st_c1['sst_NA'].values-st_c1['sst_SA'].values; XW=XW[300:]
                                XI = salt_temp_init0['sst_NA'].values-salt_temp_init0['sst_SA'].values

                        elif variables[j]=='delta_rho_sub':
                                X0 = st['rho_sub_NA'].values-st['rho_sub_SA'].values
                                XE = salt_temp006['rho_sub_NA'].values-salt_temp006['rho_sub_SA'].values; XE=XE[200:]
                                XW = st_c1['rho_sub_NA'].values-st_c1['rho_sub_SA'].values; XW=XW[300:]
                                XI = salt_temp_init0['rho_sub_NA'].values-salt_temp_init0['rho_sub_SA'].values

                        elif variables[j]=='delta_temp_sub':
                                X0 = st['temp_sub_NA'].values-st['temp_sub_SA'].values
                                XE = salt_temp006['temp_sub_NA'].values-salt_temp006['temp_sub_SA'].values; XE=XE[200:]
                                XW = st_c1['temp_sub_NA'].values-st_c1['temp_sub_SA'].values; XW=XW[300:]
                                XI = salt_temp_init0['temp_sub_NA'].values-salt_temp_init0['temp_sub_SA'].values

                        elif variables[j]=='delta_salt_sub':
                                X0 = st['salt_sub_NA'].values-st['salt_sub_SA'].values
                                XE = salt_temp006['salt_sub_NA'].values-salt_temp006['salt_sub_SA'].values; XE=XE[200:]
                                XW = st_c1['salt_sub_NA'].values-st_c1['salt_sub_SA'].values; XW=XW[300:]
                                XI = salt_temp_init0['salt_sub_NA'].values-salt_temp_init0['salt_sub_SA'].values
                        else:
                                X0 = st[variables[j]].values
                                XE = salt_temp006[variables[j]].values; XE=XE[200:]
                                XW = st_c1[variables[j]].values; XW=XW[300:]
                                XI = salt_temp_init0[variables[j]].values

                        if j==0:
                                print(i, t_amoc0[0], t_amoc0[14], t_amoc0[start_idcs[i]])


                        Y0 = eq[i]

                        Y[i] = Y0[14]

                        if i==0:
                                YE = eq_c006; YE=YE[200:]
                                YW = eq_c1; YW=YW[300:]
                                YI = eq_init0#; YI=YI[150:-20]
                                pl.plot(XI, YI, ',', color='blue')
                                pl.plot(XW, YW, ',', color='crimson')
                                pl.plot(XE, YE, ',', color='gray')

                        pl.plot(X0[start_idcs[i]], Y0[start_idcs[i]], 'x', color='blue')
                        pl.plot(X0[start_idcs[i]:start_idcs[i]+14], Y0[start_idcs[i]:start_idcs[i]+14], color='black', linewidth=0.5, alpha=0.25)
                        if i in edge_reals:
                                tip = 0.5
                                pl.plot(X0[start_idcs[i]+13], Y0[start_idcs[i]+13], '<', color='gray') 
                                pl.plot(X0[-1], Y0[-1], 'x', color='gray')
                        elif amoc0[-1]<4.:
                                tip = 1
                                pl.plot(X0[start_idcs[i]+13], Y0[start_idcs[i]+13], 'o', color='black')
                                pl.plot(X0[-1], Y0[-1], 'x', color='black')

                        else:
                                tip = 0
                                pl.plot(X0[start_idcs[i]+13], Y0[start_idcs[i]+13], 's', color='crimson') 
                                pl.plot(X0[-1], Y0[-1], 'x', color='crimson')

                        
                pl.xlabel(variables[j])



def cooling(F_hyst, amoc_hyst360q):

        ot_c = xr.open_dataarray("overturning_rtipSc0.nc")
        ot_c1 = xr.open_dataarray("overturning_rtipSc1.nc")
        ot_c05 = xr.open_dataarray("overturning_rtipSc05.nc")
        ot_c03 = xr.open_dataarray("overturning_rtipSc03.nc")
        ot_c025 = xr.open_dataarray("overturning_rtipSc025.nc")
        ot_c02 = xr.open_dataarray("overturning_rtipSc02.nc")
        ot_c019 = xr.open_dataarray("overturning_rtipSc019.nc")
        ot_c018 = xr.open_dataarray("overturning_rtipSc018.nc")
        ot_c017 = xr.open_dataarray("overturning_rtipSc017.nc")
        ot_c016 = xr.open_dataarray("overturning_rtipSc016.nc")
        ot_c015 = xr.open_dataarray("overturning_rtipSc015.nc")
        ot_c014 = xr.open_dataarray("overturning_rtipSc014.nc")
        ot_c013 = xr.open_dataarray("overturning_rtipSc013.nc")
        ot_c012 = xr.open_dataarray("overturning_rtipSc012.nc")
        ot_c011 = xr.open_dataarray("overturning_rtipSc011.nc")
        ot_c01 = xr.open_dataarray("overturning_rtipSc01.nc")
        ot_c009 = xr.open_dataarray("overturning_rtipSc009.nc")
        ot_c008 = xr.open_dataarray("overturning_rtipSc008.nc")
        ot_c007 = xr.open_dataarray("overturning_rtipSc007.nc")
        ot_c006 = xr.open_dataarray("overturning_rtipSc006.nc")
        ot_c0055 = xr.open_dataarray("overturning_rtipSc0055.nc")
        ot_c005 = xr.open_dataarray("overturning_rtipSc005.nc")
        ot_c004 = xr.open_dataarray("overturning_rtipSc004.nc")
        ot_c003 = xr.open_dataarray("overturning_rtipSc003.nc")
        ot_c002 = xr.open_dataarray("overturning_rtipSc002.nc")
        ot_c0015 = xr.open_dataarray("overturning_rtipSc0015.nc")
        ot_c001 = xr.open_dataarray("overturning_rtipSc001.nc")
        ot_c0005 = xr.open_dataarray("overturning_rtipSc0005.nc")
        ot_c0001 = xr.open_dataarray("overturning_rtipSc0001.nc")

        t_amoc_c, amoc_c = amoc_timeseries(ot_c)
        t_amoc_c1, amoc_c1 = amoc_timeseries(ot_c1)
        t_amoc_c05, amoc_c05 = amoc_timeseries(ot_c05)
        t_amoc_c03, amoc_c03 = amoc_timeseries(ot_c03)
        t_amoc_c025, amoc_c025 = amoc_timeseries(ot_c025)
        t_amoc_c02, amoc_c02 = amoc_timeseries(ot_c02)
        t_amoc_c019, amoc_c019 = amoc_timeseries(ot_c019)
        t_amoc_c018, amoc_c018 = amoc_timeseries(ot_c018)
        t_amoc_c017, amoc_c017 = amoc_timeseries(ot_c017)
        t_amoc_c016, amoc_c016 = amoc_timeseries(ot_c016)
        t_amoc_c015, amoc_c015 = amoc_timeseries(ot_c015)
        t_amoc_c014, amoc_c014 = amoc_timeseries(ot_c014)
        t_amoc_c013, amoc_c013 = amoc_timeseries(ot_c013)
        t_amoc_c012, amoc_c012 = amoc_timeseries(ot_c012)
        t_amoc_c011, amoc_c011 = amoc_timeseries(ot_c011)
        t_amoc_c01, amoc_c01 = amoc_timeseries(ot_c01)
        t_amoc_c009, amoc_c009 = amoc_timeseries(ot_c009)
        t_amoc_c008, amoc_c008 = amoc_timeseries(ot_c008)
        t_amoc_c007, amoc_c007 = amoc_timeseries(ot_c007)
        t_amoc_c006, amoc_c006 = amoc_timeseries(ot_c006)
        t_amoc_c0055, amoc_c0055 = amoc_timeseries(ot_c0055)
        t_amoc_c005, amoc_c005 = amoc_timeseries(ot_c005)
        t_amoc_c004, amoc_c004 = amoc_timeseries(ot_c004)
        t_amoc_c003, amoc_c003 = amoc_timeseries(ot_c003)
        t_amoc_c002, amoc_c002 = amoc_timeseries(ot_c002)
        t_amoc_c0015, amoc_c0015 = amoc_timeseries(ot_c0015)
        t_amoc_c001, amoc_c001 = amoc_timeseries(ot_c001)
        t_amoc_c0005, amoc_c0005 = amoc_timeseries(ot_c0005)
        t_amoc_c0001, amoc_c0001 = amoc_timeseries(ot_c0001)



        Fc1 = np.linspace(0.045,0.051,200)
        Fc1 = np.concatenate([Fc1,200*[0.051]])

        Fc05 = np.linspace(0.045,0.051,100)
        Fc05 = np.concatenate([Fc05,40*[0.051]])

        Fc03 = np.linspace(0.045,0.051,60)
        Fc03 = np.concatenate([Fc03,80*[0.051]])

        Fc025 = np.linspace(0.045,0.051,50)
        Fc03 = np.concatenate([Fc03,90*[0.051]])

        Fc02 = np.linspace(0.045,0.051,40)
        Fc02 = np.concatenate([Fc02,60*[0.051]])

        Fc019 = np.linspace(0.045,0.051,38)
        Fc019 = np.concatenate([Fc019,102*[0.051]])

        Fc018 = np.linspace(0.045,0.051,36)
        Fc018 = np.concatenate([Fc018,104*[0.051]])

        Fc017 = np.linspace(0.045,0.051,34)
        Fc017 = np.concatenate([Fc017,106*[0.051]])

        Fc016 = np.linspace(0.045,0.051,32)
        Fc016 = np.concatenate([Fc016,108*[0.051]])

        Fc015 = np.linspace(0.045,0.051,30)
        Fc015 = np.concatenate([Fc015,70*[0.051]])

        Fc014 = np.linspace(0.045,0.051,28)
        Fc014 = np.concatenate([Fc014,112*[0.051]])

        Fc013 = np.linspace(0.045,0.051,26)
        Fc013 = np.concatenate([Fc013,74*[0.051]])

        Fc012 = np.linspace(0.045,0.051,24)
        Fc012 = np.concatenate([Fc012,76*[0.051]])

        Fc011 = np.linspace(0.045,0.051,22)
        Fc011 = np.concatenate([Fc011,78*[0.051]])

        Fc01 = np.linspace(0.045,0.051,20)
        Fc01 = np.concatenate([Fc01,120*[0.051]])

        Fc009 = np.linspace(0.045,0.051,18)
        Fc009 = np.concatenate([Fc009,82*[0.051]])

        Fc008 = np.linspace(0.045,0.051,16)
        Fc008 = np.concatenate([Fc008,84*[0.051]])

        Fc007 = np.linspace(0.045,0.051,14)
        Fc007 = np.concatenate([Fc007,86*[0.051]])

        Fc006 = np.linspace(0.045,0.051,12)
        Fc006 = np.concatenate([Fc006,88*[0.051]])

        Fc0055 = np.linspace(0.045,0.051,11)
        Fc0055 = np.concatenate([Fc0055,69*[0.051]])

        Fc005 = np.linspace(0.045,0.051,10)
        Fc005 = np.concatenate([Fc005,70*[0.051]])

        Fc004 = np.linspace(0.045,0.051,8)
        Fc004 = np.concatenate([Fc004,112*[0.051]])

        Fc003 = np.linspace(0.045,0.051,6)
        Fc003 = np.concatenate([Fc003,194*[0.051]])

        Fc002 = np.linspace(0.045,0.051,4)
        Fc002 = np.concatenate([Fc002,196*[0.051]])

        Fc0015 = np.linspace(0.045,0.051,3)
        Fc0015 = np.concatenate([Fc0015,137*[0.051]])

        Fc001 = np.linspace(0.045,0.051,2)
        Fc001 = np.concatenate([Fc001,138*[0.051]])

        Fc0005 = np.linspace(0.045,0.051,1)
        Fc0005 = np.concatenate([Fc0005,139*[0.051]])
        
        fact = 28.5646


        fig=pl.figure()
        pl.subplot(221)
        pl.plot(t_amoc_c, amoc_c, label='Spinup')
        pl.plot(t_amoc_c1, amoc_c1, label='1000y ramp')
        pl.plot(t_amoc_c05, amoc_c05, label='500y')
        pl.plot(t_amoc_c03, amoc_c03, label='300y')
        pl.plot(t_amoc_c025, amoc_c025, label='250y')
        pl.plot(t_amoc_c02, amoc_c02, label='200y')
        pl.plot(t_amoc_c019, amoc_c019, label='190y')
        pl.plot(t_amoc_c018, amoc_c018, label='180y')
        pl.plot(t_amoc_c017, amoc_c017, label='170y')
        pl.plot(t_amoc_c016, amoc_c016, label='160y')
        pl.legend(loc='best')

        pl.subplot(222)
        pl.plot(t_amoc_c015, amoc_c015, label='150y')
        pl.plot(t_amoc_c014, amoc_c014, label='140y')
        pl.plot(t_amoc_c013, amoc_c013, label='130y')
        pl.plot(t_amoc_c012, amoc_c012, label='120y')
        pl.plot(t_amoc_c011, amoc_c011, label='110y')
        pl.legend(loc='best')

        pl.subplot(223)

        pl.plot(t_amoc_c01, amoc_c01, label='100y')
        pl.plot(t_amoc_c009, amoc_c009, label='90y')
        pl.plot(t_amoc_c008, amoc_c008, label='80y')
        pl.plot(t_amoc_c007, amoc_c007, label='70y')
        pl.plot(t_amoc_c006, amoc_c006, label='60y')
        pl.legend(loc='best')

        pl.subplot(224)
        pl.plot(t_amoc_c0055, amoc_c0055, label='55y')
        pl.plot(t_amoc_c005, amoc_c005, label='50y')
        pl.plot(t_amoc_c004, amoc_c004, label='40y')
        pl.plot(t_amoc_c003, amoc_c003, label='30y')
        pl.plot(t_amoc_c002, amoc_c002, label='20y')
        pl.plot(t_amoc_c0015, amoc_c0015, label='15y')
        pl.plot(t_amoc_c001, amoc_c001, label='10y')
        pl.plot(t_amoc_c0005, amoc_c0005, label='5y')
        pl.plot(t_amoc_c0001, amoc_c0001, label='1y')
        pl.legend(loc='best')

        fig=pl.figure(figsize=(10,12))
        pl.subplots_adjust(left=0.08, bottom=0.07, right=0.97, top=0.98, wspace=0.0, hspace=0.06)
        ax=pl.subplot(511)
        pl.plot(t_amoc_c015, amoc_c015, label='150y')
        pl.legend(loc='best'); pl.xlim(9000,16200);pl.ylim(2.0,7.2)
        ax.spines['bottom'].set_visible(False); ax.get_xaxis().set_visible(False)
        ax.spines["left"].set_position(("outward", 3));pl.ylabel('Max. AMOC (Sv)')
        ax=pl.subplot(512)
        pl.plot(t_amoc_c01, amoc_c01, label='100y')
        pl.legend(loc='best'); pl.xlim(9000,16200);pl.ylim(2.0,7.2)
        ax.spines['bottom'].set_visible(False); ax.get_xaxis().set_visible(False)
        ax.spines["left"].set_position(("outward", 3));pl.ylabel('Max. AMOC (Sv)')
        ax=pl.subplot(513)
        pl.plot(t_amoc_c006, amoc_c006, label='60y')
        pl.legend(loc='best'); pl.xlim(9000,16200);pl.ylim(2.0,7.2)
        ax.spines['bottom'].set_visible(False); ax.get_xaxis().set_visible(False)
        ax.spines["left"].set_position(("outward", 3));pl.ylabel('Max. AMOC (Sv)')
        ax=pl.subplot(514)
        pl.plot(t_amoc_c0055, amoc_c0055, label='55y')
        pl.legend(loc='best'); pl.xlim(9000,16200);pl.ylim(2.0,7.2)
        ax.spines['bottom'].set_visible(False); ax.get_xaxis().set_visible(False)
        ax.spines["left"].set_position(("outward", 3));pl.ylabel('Max. AMOC (Sv)')
        ax=pl.subplot(515)
        pl.plot(t_amoc_c001, amoc_c001, label='10y')
        pl.legend(loc='best'); pl.xlim(9000,16200);pl.ylim(2.0,7.2)
        ax.spines["bottom"].set_position(("outward", 3))
        ax.spines["left"].set_position(("outward", 3))
        pl.xlabel('Simulation time (years)');pl.ylabel('Max. AMOC (Sv)')

        n_lines = 26
        fig=pl.figure(figsize=(9,7))
        pl.subplots_adjust(left=0.08, bottom=0.12, right=0.97, top=0.98, wspace=0.0, hspace=0.06)
        pl.plot(t_amoc_c, amoc_c, label='Spinup')
        ax = pl.axes()
        ax.set_prop_cycle(color=[pl.cm.viridis(i) for i in np.linspace(0, 1, n_lines)])

        pl.plot(t_amoc_c03, amoc_c03, label='300y')
        pl.plot(t_amoc_c025, amoc_c025, label='250y')
        pl.plot(t_amoc_c02, amoc_c02, label='200y')
        pl.plot(t_amoc_c019, amoc_c019, label='190y')
        pl.plot(t_amoc_c018, amoc_c018, label='180y')
        pl.plot(t_amoc_c017, amoc_c017, label='170y')
        pl.plot(t_amoc_c016, amoc_c016, label='160y')
        pl.plot(t_amoc_c015, amoc_c015, label='150y')
        pl.plot(t_amoc_c014, amoc_c014, label='140y')
        pl.plot(t_amoc_c013, amoc_c013, label='130y')
        pl.plot(t_amoc_c012, amoc_c012, label='120y')
        pl.plot(t_amoc_c011, amoc_c011, label='110y')
        pl.plot(t_amoc_c01, amoc_c01, label='100y')
        pl.plot(t_amoc_c009, amoc_c009, label='90y')
        pl.plot(t_amoc_c008, amoc_c008, label='80y')
        pl.plot(t_amoc_c007, amoc_c007, label='70y')
        pl.plot(t_amoc_c006, amoc_c006, label='60y')
        pl.plot(t_amoc_c0055, amoc_c0055, label='55y')
        pl.plot(t_amoc_c005, amoc_c005, label='50y')
        pl.plot(t_amoc_c004, amoc_c004, label='40y')
        pl.plot(t_amoc_c003, amoc_c003, label='30y')
        pl.plot(t_amoc_c002, amoc_c002, label='20y')
        pl.plot(t_amoc_c0015, amoc_c0015, label='15y')
        pl.plot(t_amoc_c001, amoc_c001, label='10y')
        pl.plot(t_amoc_c0005, amoc_c0005, label='5y')
        pl.plot(t_amoc_c0001, amoc_c0001, label='1y')
        ax.spines["bottom"].set_position(("outward", 3))
        ax.spines["left"].set_position(("outward", 3))
        pl.xlabel('Simulation time (years)');pl.ylabel('Max. AMOC (Sv)')
        
        to_file = np.asarray([t_amoc_c008, amoc_c0001[:140], amoc_c0005[:140], amoc_c001[:140], amoc_c0015[:140], amoc_c002[:140], amoc_c003[:140], amoc_c004[:140], amoc_c005[:140], amoc_c0055[:140], amoc_c006[:140], amoc_c007[:140], amoc_c008[:140], amoc_c009[:140], amoc_c01[:140], amoc_c011[:140], amoc_c012[:140], amoc_c013[:140], amoc_c014[:140], amoc_c015[:140], amoc_c016[:140], amoc_c017[:140], amoc_c018[:140], amoc_c019[:140], amoc_c02[:140], amoc_c025[:140], amoc_c03[:140]])#.transpose()
        
        #np.savetxt('veros_timeseries_rates.txt', to_file)
        np.savetxt('veros_timeseries_spinup.txt', np.asarray([t_amoc_c[:121], amoc_c[:121]]))

        fig=pl.figure()
        pl.plot(t_amoc_c, amoc_c, label='Spinup')
        pl.plot(t_amoc_c1, amoc_c1, label='1000y ramp')
        pl.plot(t_amoc_c016, amoc_c016, label='160y')
        pl.plot(t_amoc_c0001, amoc_c0001, label='1y')
        pl.legend(loc='best')
        

        fig=pl.figure(figsize=(0.5,7))
        pl.subplots_adjust(left=0.03, bottom=0.03, right=0.3, top=0.98, wspace=0.0, hspace=0.06)
        ax=pl.subplot(111)
        vals = [0,5,10,15,20,30,40,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300]

        cmap = mpl.colors.ListedColormap([pl.cm.viridis(i) for i in np.linspace(0, 1, n_lines)][::-1])
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, boundaries=vals+[350], ticks=vals, orientation='vertical')

        cb.set_label('Ramping duration (years)')

        fig=pl.figure(figsize=(13,4))
        pl.subplots_adjust(left=0.07, bottom=0.18, right=0.96, top=0.96, wspace=0.2, hspace=0.06)
        ax=pl.subplot(121)
        ax.spines["bottom"].set_position(("outward", 4))
        ax.spines["left"].set_position(("outward", 4))
        pl.plot(F_hyst[:29], amoc_hyst360q[:29], 'o-', color='black')
        pl.plot(F_hyst[29:len(amoc_hyst360q)], amoc_hyst360q[29:], 'o-', color='black')
        pl.plot(100*[1.285], np.linspace(min(amoc_hyst360q), max(amoc_hyst360q), 100), color='gray')
        pl.plot(100*[1.456], np.linspace(min(amoc_hyst360q), max(amoc_hyst360q), 100), color='gray')

        points = np.array([fact*Fc016, amoc_c016[:len(Fc016)]]).transpose().reshape(-1,1,2)
        segs = np.concatenate([points[:-1],points[1:]],axis=1)
        lc = LineCollection(segs, cmap=pl.get_cmap('Spectral'))#,norm=plt.Normalize(250, 1500))
        lc.set_array(np.asarray(t_amoc_c016[:len(Fc016)])-9100.)
        lc.set_linewidth(3)
        pl.gca().add_collection(lc)

        axcb = fig.colorbar(lc)
        axcb.set_label('time (years)')

        pl.xlim(0.9,1.7)
        pl.xlabel('Freshwater forcing (PSU)')
        pl.ylabel('Max. AMOC (Sv)')

        ax=pl.subplot(122)
        ax.spines["bottom"].set_position(("outward", 4))
        ax.spines["left"].set_position(("outward", 4))
        pl.plot(F_hyst[:29], amoc_hyst360q[:29], 'o-', color='black')
        pl.plot(F_hyst[29:len(amoc_hyst360q)], amoc_hyst360q[29:], 'o-', color='black')
        pl.plot(100*[1.285], np.linspace(min(amoc_hyst360q), max(amoc_hyst360q), 100), color='gray')
        pl.plot(100*[1.456], np.linspace(min(amoc_hyst360q), max(amoc_hyst360q), 100), color='gray')

        points1 = np.array([fact*Fc014, amoc_c014]).transpose().reshape(-1,1,2)
        #print(points1)
        segs1 = np.concatenate([points1[:-1],points1[1:]],axis=1)
        #print(segs1)
        lc1 = LineCollection(segs1, cmap=pl.get_cmap('Spectral'))
        lc1.set_array(np.asarray(t_amoc_c014)-9100.)
        lc1.set_linewidth(3)
        pl.gca().add_collection(lc1)

        axcb = fig.colorbar(lc1)
        axcb.set_label('time (years)')

        pl.xlim(0.9,1.7)
        pl.xlabel('Freshwater forcing (PSU)')
        pl.ylabel('Max. AMOC (Sv)')

        def make_fig(t_amoc, amoc, ramp_dur, width, bar=False):
                t_amoc = np.concatenate((t_amoc_c[100:119], t_amoc))
                amoc = np.concatenate((amoc_c[100:119], amoc))
                r0 = np.array(len(amoc)*[ramp_dur])
                if bar==True:
                        pl.plot(np.asarray(t_amoc)-9100., r0, alpha=0.0)
                points2 = np.array([np.asarray(t_amoc)-9100., r0]).transpose().reshape(-1,1,2)
                segs2 = np.concatenate([points2[:-1],points2[1:]],axis=1)
                lc2 = LineCollection(segs2, cmap=pl.get_cmap('Spectral'),  norm = mpl.colors.Normalize(vmin=2.,vmax=7.4)
)
                lc2.set_array(np.asarray(amoc))
                lc2.set_linewidth(width)
                pl.gca().add_collection(lc2)
                if bar==True:
                        ax2 = fig.colorbar(lc2)
                        ax2.set_label('AMOC')
                

        fig=pl.figure(figsize=(12,5))
        pl.subplots_adjust(left=0.1, bottom=0.14, right=0.97, top=0.98, wspace=0.02, hspace=0.0)

        spec = gridspec.GridSpec(ncols=2, nrows=1,
                         width_ratios=[1, 5])

        ax0 = fig.add_subplot(spec[0])
        ax1 = fig.add_subplot(spec[1])


        w3=50
        w2=10.
        w1=5.

        make_fig(t_amoc_c03, amoc_c03, 300., width=w3)
        make_fig(t_amoc_c025, amoc_c025, 250., width=w3)
        make_fig(t_amoc_c02, amoc_c02, 200., width=w3)
        make_fig(t_amoc_c019, amoc_c019, 190., width=w2)
        make_fig(t_amoc_c018, amoc_c018, 180., width=w2)
        make_fig(t_amoc_c017, amoc_c017, 170., width=w2)
        make_fig(t_amoc_c016[:140], amoc_c016[:140], 160., width=w2)
        make_fig(t_amoc_c015[:140], amoc_c015[:140], 150., width=w2)
        make_fig(t_amoc_c014, amoc_c014, 140., width=w2)
        make_fig(t_amoc_c013, amoc_c013, 130., width=w2)
        make_fig(t_amoc_c012, amoc_c012, 120., width=w2)
        make_fig(t_amoc_c011, amoc_c011, 110., width=w2)
        make_fig(t_amoc_c01[:140], amoc_c01[:140], 100., width=w2)
        make_fig(t_amoc_c009, amoc_c009, 90., width=w2)
        make_fig(t_amoc_c008, amoc_c008, 80., width=w2)
        make_fig(t_amoc_c007, amoc_c007, 70., width=w2)
        make_fig(t_amoc_c006[:140], amoc_c006[:140], 60., width=w2)
        make_fig(t_amoc_c005, amoc_c005, 50., width=w2)
        make_fig(t_amoc_c0055[:140], amoc_c0055[:140], 55., width=w1)
        make_fig(t_amoc_c004, amoc_c004, 40., width=w2)
        make_fig(t_amoc_c003[:140], amoc_c003[:140], 30., width=w2)
        make_fig(t_amoc_c002[:140], amoc_c002[:140], 20., width=w2)
        make_fig(t_amoc_c0015, amoc_c0015, 15., width=w1)
        make_fig(t_amoc_c001[:140], amoc_c001[:140], 10., width=w1)
        make_fig(t_amoc_c0005, amoc_c0005, 5., width=w1)
        make_fig(t_amoc_c0001, amoc_c0001, 0., width=w1, bar=True)
        
        pl.xlabel('Simulation time (years)');pl.ylabel('Ramping duration (years)')
        pl.xlim(-100.,700.)

        '''
        N = len(ot_c.time)
        fig=pl.figure()
        for i in range(16):
                pl.subplot(4,4,i+1)
                ot_c.isel(time=int(i*N/16.)).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        '''

        '''
        ### Lorenz map from 1k run (after ramping, 200 data points)
        maxima = []; t_max=[]
        for i in range(1,199):
                if ((amoc_c1[-200+i-1]<amoc_c1[-200+i]) and (amoc_c1[-200+i]>amoc_c1[-200+i+1])):
                        maxima.append(amoc_c1[-200+i]); t_max.append(t_amoc_c1[-200+i])

        fig=pl.figure()
        pl.plot(maxima[:-1], maxima[1:], 'o')
        '''

def warming(F_hyst, amoc_hyst360q):

        ot = xr.open_dataarray("overturning_rtipSw0.nc")
        ot1 = xr.open_dataarray("overturning_rtipSw1.nc")
        ot05 = xr.open_dataarray("overturning_rtipSw05.nc")
        ot01 = xr.open_dataarray("overturning_rtipSw01.nc")
        ot005 = xr.open_dataarray("overturning_rtipSw005.nc")
        ot002 = xr.open_dataarray("overturning_rtipSw002.nc")

        t_amoc, amoc = amoc_timeseries(ot1)
        t_amoc05, amoc05 = amoc_timeseries(ot05)
        t_amoc01, amoc01 = amoc_timeseries(ot01)
        t_amoc005, amoc005 = amoc_timeseries(ot005)
        t_amoc002, amoc002 = amoc_timeseries(ot002)

        F1 = np.linspace(0.054,0.048,200)
        F1 = np.concatenate([F1,60*[0.048]])

        F05 = np.linspace(0.054,0.048,100)
        F05 = np.concatenate([F05,60*[0.048]])

        F01 = np.linspace(0.054,0.048,20)
        F01 = np.concatenate([F01,60*[0.048]])

        F005 = np.linspace(0.054,0.048,10)
        F005 = np.concatenate([F005,50*[0.048]])

        F002 = np.linspace(0.054,0.048,4)
        F002 = np.concatenate([F002,56*[0.048]])


        fig=pl.figure()
        pl.plot(t_amoc, amoc)
        pl.plot(t_amoc05, amoc05)
        pl.plot(t_amoc01, amoc01)
        pl.plot(t_amoc005, amoc005)
        pl.plot(t_amoc002, amoc002)

        fig=pl.figure()
        pl.plot(10000.*F_hyst[:27], amoc_hyst360q[:27], color='black')
        pl.plot(10000.*F_hyst[27:len(amoc_hyst360q)], amoc_hyst360q[27:], color='black')
        pl.plot(F1, amoc)
        pl.plot(F05, amoc05)
        pl.plot(F01, amoc01)
        pl.plot(F005, amoc005)
        pl.plot(F002, amoc002)

        pl.xlabel('Freshwater forcing (a.u.)')
        pl.ylabel('Max. AMOC (Sv)')

def test_eof():

        ### EOF analysis of NA stream function
        
        ot_init = xr.open_dataarray("overturning_rtipSc0i.nc")

        N = len(ot_init.time)
        fig=pl.figure()
        for i in range(16):
                pl.subplot(4,4,i+1)
                ot_init.isel(time=int(i*N/16.)).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)

        ot_init = ot_init - ot_init.mean(dim='time')

        # Create an EOF solver to do the EOF analysis. Square-root of cosine of
        # latitude weights are applied before the computation of EOFs.
        #coslat = np.cos(np.deg2rad(z_djf.coords['latitude'].values)).clip(0., 1.)
        #wgts = np.sqrt(coslat)[..., np.newaxis]
        
        solver = Eof(ot_init)#, weights=wgts

        # Retrieve the leading EOF, expressed as the covariance between the leading
        # PC time series and the input anomalies at each grid point.
        eof1 = solver.eofsAsCovariance(neofs=2)

        # leading PC time series
        pc1 = solver.pcs(npcs=2, pcscaling=1)

        fig=pl.figure()
        ax=pl.subplot(121)
        eof1[0].plot.contourf(ax=ax, levels=20, cmap=pl.cm.RdBu_r,
                             add_colorbar=False)
        ax=pl.subplot(122)
        eof1[1].plot.contourf(ax=ax, levels=20, cmap=pl.cm.RdBu_r,
                             add_colorbar=False)

        pl.figure()
        pl.subplot(211)
        pc1[:, 0].plot()
        pl.subplot(212)
        pc1[:, 1].plot()

        N = len(ot_init.time)
        fig=pl.figure()
        for i in range(16):
                pl.subplot(4,4,i+1)
                ot_init.isel(time=int(i*N/16.)).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)


def aabw_timeseries(ot):
        aabw=[]; t_aabw=[]
        for i in range(len(ot.time)):
                t_aabw.append(ot.time[i]/360.)
                aabw.append(ot.isel(time=i, depth=slice(0,32), lat=slice(0,15)).min().values/1000000.)
        return t_aabw, aabw

def amoc_timeseries(ot):
        t_amoc = []; amoc = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i].values/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max().values/1000000.)
        return t_amoc, amoc

def amoc_aabw_timeseries(ot):
        t_amoc = []; amoc = []; aabw = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i].values/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max().values/1000000.)
                aabw.append(ot.isel(time=i, depth=slice(0,32), lat=slice(0,15)).min().values/1000000.)
        return t_amoc, amoc, aabw

def amoc_aabw_timeseries_ext(ot):
        t_amoc = []; amoc = []; aabw = []; eq = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i].values/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max().values/1000000.)
                eq.append(ot.isel(time=i, depth=slice(0,27), lat=22).max().values/1000000.)
                aabw.append(ot.isel(time=i, depth=slice(25,40), lat=22).min().values/1000000.)
        return t_amoc, amoc, eq, aabw

def calc_pvalue(raw,hist,center,binwidth,data):
        hist_p = []
        center_p = []
        pvalue = 0
        for i in range(len(center)):
	        if data<np.mean(raw):
		        if center[i]<data:
			        pvalue = pvalue + hist[i]*binwidth
			        hist_p.append(hist[i])
			        center_p.append(center[i])
	        if data>np.mean(raw):
		        if center[i]>data:
			        pvalue = pvalue + hist[i]*binwidth
			        hist_p.append(hist[i])
			        center_p.append(center[i])
        return [hist_p,center_p,pvalue]  

if __name__ == '__main__':
        MAIN()

