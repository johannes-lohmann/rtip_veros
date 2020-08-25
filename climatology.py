import matplotlib.pyplot as pl
import numpy as np
import xarray as xr
import matplotlib as mpl

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

        #salt_forcing()
        #climatology()
        test()
        #forcing_all()
        pl.show()

def forcing_all():
        forc_all = xr.open_dataset("forcing_1deg_global_interp_hosing.nc")
        print forc_all

        #print forc_all.q_net

        print forc_all.isel(Time=0)['dqdt'].values
        print forc_all.isel(Time=3)['dqdt'].values

def salt_forcing():

        #N=2
        #for i in range(N):
        #        count = str(i); count = count.zfill(4)
        #        avgs = xr.open_dataset("ctrl40l.%s.averages.nc"%count)
        #        if i==0:
        #                forc = avgs['forc_salt_surface']
        #        else:
        #                forc=xr.concat((forc,avgs['forc_salt_surface']),dim='Time')

        forc = xr.open_dataarray("salt_forcing_monthly40l.nc")

        print forc

        forc0 = []; t = []
        for i in range(len(forc.Time)):
                forc0.append(forc.isel(Time=i,xt=65,yt=35))
                t.append(forc.Time[i])

        fig=pl.figure()
        pl.plot(t, forc0)
        
        

def test():

        ot = xr.open_dataarray("overturning_ctrl40l.nc")
        salt_temp = xr.open_dataset("salt_temp_ctrl40l.nc")

        ot2 = xr.open_dataarray("overturning_sflux40l.nc")
        salt_temp2 = xr.open_dataset("salt_temp_sflux40l.nc")

        ot_T0 = xr.open_dataarray("overturning_sflux40l_T0.nc")
        salt_temp_T0 = xr.open_dataset("salt_temp_sflux40l_T0.nc")

        ot_br = xr.open_dataarray("overturning_sflux40l_br.nc")
        salt_temp_br = xr.open_dataset("salt_temp_sflux40l_br.nc")

        ot_pert = xr.open_dataarray("overturning_sflux40l_pert.nc")
        salt_temp_pert = xr.open_dataset("salt_temp_sflux40l_pert.nc")

        ot_pert2 = xr.open_dataarray("overturning_sflux40l_pert2.nc")
        salt_temp_pert2 = xr.open_dataset("salt_temp_sflux40l_pert2.nc")

        ot_pert3 = xr.open_dataarray("overturning_sflux40l_pert3.nc")
        salt_temp_pert3 = xr.open_dataset("salt_temp_sflux40l_pert3.nc")

        avg = xr.open_dataset("ctrl40l.0039.averages.nc")
        #avg2 = xr.open_dataset("s40lfl.0069.averages.nc")

        print avg.xt
        print avg.yt


        t_amoc = []; amoc = []; aabw = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i]/360.)
                idx=i
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                aabw.append(ot.isel(time=i, depth=slice(0,32), lat=slice(0,15)).min()/1000000.)

        for i in range(len(ot2.time)):
                t_amoc.append(ot2.time[i]/360.)
                amoc.append(ot2.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                aabw.append(ot2.isel(time=i, depth=slice(0,32), lat=slice(0,15)).min()/1000000.)

        t_amoc_br = []; amoc_br = []; aabw_br = []
        for i in range(len(ot_br.time)):
                t_amoc_br.append(ot_br.time[i]/360.)
                amoc_br.append(ot_br.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                aabw_br.append(ot_br.isel(time=i, depth=slice(0,32), lat=slice(0,15)).min()/1000000.)

        t_amoc_T0 = []; amoc_T0 = []
        for i in range(len(ot_T0.time)):
                t_amoc_T0.append(ot_T0.time[i]/360.)
                amoc_T0.append(ot_T0.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        t_amoc_pert = []; amoc_pert = []
        for i in range(len(ot_pert.time)):
                t_amoc_pert.append(ot_pert.time[i]/360.)
                amoc_pert.append(ot_pert.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        t_amoc_pert2 = []; amoc_pert2 = []
        for i in range(len(ot_pert2.time)):
                t_amoc_pert2.append(ot_pert2.time[i]/360.)
                amoc_pert2.append(ot_pert2.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        t_amoc_pert3 = []; amoc_pert3 = []
        for i in range(len(ot_pert3.time)):
                t_amoc_pert3.append(ot_pert3.time[i]/360.)
                amoc_pert3.append(ot_pert3.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)


        fig=pl.figure()
        pl.subplot(611)
        pl.plot(t_amoc[idx:], amoc[idx:],label='Control')
        pl.ylabel('Max. AMOC (Sv.)')
        pl.xlim(4000,14000)
        pl.legend(loc='best')
        pl.subplot(612)
        pl.plot(t_amoc_pert, amoc_pert,label='Perturb T=0.2, S=0.1')
        pl.ylabel('Max. AMOC (Sv.)')
        pl.xlim(4000,14000)
        pl.legend(loc='best')
        pl.subplot(613)
        pl.plot(t_amoc_pert2, amoc_pert2,label='Perturb T=0.1, S=0.05')
        pl.ylabel('Max. AMOC (Sv.)')
        pl.xlim(4000,14000)
        pl.legend(loc='best')
        pl.subplot(614)
        pl.plot(t_amoc_pert3, amoc_pert3,label='Perturb T=0.1, S=0.05')
        pl.ylabel('Max. AMOC (Sv.)')
        pl.xlim(4000,14000)
        pl.legend(loc='best')
        pl.subplot(615)
        pl.plot(t_amoc_br, amoc_br,label='Branch < timestep')
        pl.ylabel('Max. AMOC (Sv.)')
        pl.xlim(4000,14000)
        pl.legend(loc='best')
        pl.subplot(616)
        pl.plot(t_amoc_T0, amoc_T0,label='total zero forc.')
        pl.ylabel('Max. AMOC (Sv.)')
        pl.xlim(4000,14000)
        pl.legend(loc='best')
        pl.xlabel('Time (years)')

        fig=pl.figure()
        pl.subplot(711)
        pl.plot(t_amoc, amoc)
        pl.plot(t_amoc[idx:], amoc[idx:])
        pl.plot(t_amoc_br, amoc_br)
        pl.plot(t_amoc_pert, amoc_pert)
        pl.plot(t_amoc_pert2, amoc_pert2)
        pl.plot(t_amoc_pert3, amoc_pert3)
        pl.ylabel('Max. AMOC (Sv.)')
        pl.subplot(712)
        salt_temp['salt'].plot()
        salt_temp2['salt'].plot()
        #salt_temp_br['salt'].plot()
        salt_temp_pert['salt'].plot()
        pl.subplot(713)
        salt_temp['temp'].plot()
        salt_temp2['temp'].plot()
        salt_temp_br['temp'].plot()
        salt_temp_pert['temp'].plot()
        salt_temp_pert2['temp'].plot()
        salt_temp_pert3['temp'].plot()
        pl.subplot(714)
        salt_temp['rhoN'].plot()
        salt_temp['rhoS'].plot()
        salt_temp2['rhoN'].plot()
        salt_temp2['rhoS'].plot()
        #(salt_temp['rhoN']-salt_temp['rhoS']).plot()
        pl.subplot(715)
        salt_temp['salt_tot'].plot()
        salt_temp2['salt_tot'].plot()
        salt_temp_br['salt_tot'].plot()
        salt_temp_pert['salt_tot'].plot()
        salt_temp_pert2['salt_tot'].plot()
        salt_temp_pert3['salt_tot'].plot()
        pl.subplot(716)
        salt_temp['salt_forc_tot'].plot()
        salt_temp2['salt_forc_tot'].plot()
        pl.subplot(717)
        salt_temp2['seaice'].plot()
        salt_temp_br['seaice'].plot()
        salt_temp_pert['seaice'].plot()
        salt_temp_pert2['seaice'].plot()
        salt_temp_pert3['seaice'].plot()

        fig=pl.figure()
        pl.subplot(711)
        pl.plot(t_amoc, amoc)
        pl.plot(t_amoc[idx:], amoc[idx:])
        pl.plot(t_amoc_T0, amoc_T0)
        pl.plot(t_amoc_br, amoc_br)
        pl.plot(t_amoc_pert, amoc_pert)
        pl.plot(t_amoc_pert2, amoc_pert2)
        pl.plot(t_amoc_pert3, amoc_pert3)
        pl.ylabel('Max. AMOC (Sv.)')
        pl.subplot(712)
        pl.plot(t_amoc, aabw,label='Control')
        pl.plot(t_amoc_br, aabw_br)
        pl.xlabel('Time (years)');pl.ylabel('Max. AABW (Sv.)')
        pl.subplot(713)
        salt_temp['temp'].plot()
        salt_temp2['temp'].plot()
        salt_temp_br['temp'].plot()
        salt_temp_pert['temp'].plot()
        salt_temp_pert2['temp'].plot()
        salt_temp_pert3['temp'].plot()
        pl.subplot(714)
        salt_temp['rhoN'].plot()
        salt_temp['rhoS'].plot()
        salt_temp2['rhoN'].plot()
        salt_temp2['rhoS'].plot()
        salt_temp_br['rhoN'].plot()
        salt_temp_br['rhoS'].plot()

        #(salt_temp['rhoN']-salt_temp['rhoS']).plot()
        pl.subplot(715)
        salt_temp['salt_tot'].plot()
        salt_temp2['salt_tot'].plot()
        salt_temp_br['salt_tot'].plot()
        #salt_temp_pert['salt_tot'].plot()
        #salt_temp_pert2['salt_tot'].plot()
        #salt_temp_pert3['salt_tot'].plot()
        salt_temp_T0['salt_tot'].plot()
        pl.subplot(716)
        salt_temp['salt_forc_tot'].plot()
        salt_temp2['salt_forc_tot'].plot()
        salt_temp_T0['salt_forc_tot'].plot()
        pl.subplot(717)
        salt_temp2['seaice'].plot()
        salt_temp_br['seaice'].plot()
        #salt_temp_pert['seaice'].plot()
        #salt_temp_pert2['seaice'].plot()
        #salt_temp_pert3['seaice'].plot()
        salt_temp_T0['seaice'].plot()
        pl.xlim(-300000,12000000)

        fig=pl.figure()
        pl.subplot(611)
        pl.plot(t_amoc_br, amoc_br)
        pl.ylabel('Max. AMOC (Sv.)')
        pl.subplot(612)
        pl.plot(t_amoc_br, aabw_br,label='Control')
        pl.xlabel('Time (years)');pl.ylabel('Max. AABW (Sv.)')
        pl.subplot(613)
        salt_temp_br['temp_surfN'].plot()
        pl.subplot(614)
        salt_temp_br['temp_surfS'].plot()
        pl.subplot(615)
        salt_temp_br['salt_tot'].plot()
        pl.subplot(616)
        salt_temp_br['seaice'].plot()
        pl.xlim(-300000,12000000)

        fig=pl.figure(figsize=(8,3))
        pl.plot(t_amoc[:idx], amoc[:idx], color='black')
        pl.plot(t_amoc_br, amoc_br, color='black')
        pl.ylabel('Max. AMOC (Sv)')
        pl.xlabel('Simulation time (years)')


        '''
        fig=pl.figure(figsize=(8,3))
        pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
        pl.plot(t_amoc[:idx], amoc[:idx])
        pl.ylabel('Max. AMOC (Sv.)');pl.xlabel('Time (years)')
        pl.xlim(0,4000)

        fig=pl.figure(figsize=(8,3))
        pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
        pl.plot(t_amoc[:idx], amoc[:idx])
        pl.plot(t_amoc[idx:], amoc[idx:])
        pl.xlim(0,6500)
        pl.ylabel('Max. AMOC (Sv.)');pl.xlabel('Time (years)')


        fig=pl.figure(figsize=(8,3))
        pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
        pl.plot(t_amoc[:idx], amoc[:idx])
        pl.plot(t_amoc[idx:], amoc[idx:])
        pl.plot(t_amoc_br, amoc_br)
        pl.xlim(0,6500)
        pl.ylabel('Max. AMOC (Sv.)');pl.xlabel('Time (years)')

        fig=pl.figure(figsize=(8,3))
        pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
        pl.plot(t_amoc[:idx], amoc[:idx])
        pl.plot(t_amoc[idx:], amoc[idx:])
        pl.plot(t_amoc_br, amoc_br)
        pl.plot(t_amoc_pert, amoc_pert)
        pl.plot(t_amoc_pert2, amoc_pert2)
        pl.plot(t_amoc_pert3, amoc_pert3)
        pl.xlim(0,6500)
        pl.ylabel('Max. AMOC (Sv.)');pl.xlabel('Time (years)')


        fig=pl.figure(figsize=(8,3))
        pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
        pl.plot(t_amoc[:idx], amoc[:idx])
        pl.plot(t_amoc[idx:], amoc[idx:])
        pl.plot(t_amoc_br, amoc_br)
        pl.plot(t_amoc_pert, amoc_pert)
        pl.plot(t_amoc_pert2, amoc_pert2)
        pl.plot(t_amoc_pert3, amoc_pert3)
        pl.ylabel('Max. AMOC (Sv.)');pl.xlabel('Time (years)')
        pl.xlim(0,14000)
        '''


        ###Mark regions with sea ice
        '''
        avg['temp'] = avg['temp'].fillna(10.)
        avg['temp'] = avg['temp'].where(avg['temp']>=-1.8, -15.)
        avg['temp'] = avg['temp'].where(avg['temp']!=10., np.nan)

        avg2['temp'] = avg2['temp'].fillna(10.)
        avg2['temp'] = avg2['temp'].where(avg2['temp']>=-1.8, -15.)
        avg2['temp'] = avg2['temp'].where(avg2['temp']!=10., np.nan)
        '''

        ###Count grid cells with perennial sea ice
        '''
        seaice = avg.isel(Time=-1,zt=-1)['temp'].where(avg.isel(Time=-1,zt=-1)['temp']<=-1.8, drop=True)
        seaice = seaice.fillna(10.)
        seaice = seaice.where(seaice>0., 1.)
        seaice = seaice.where(seaice!=10., np.nan)
        print seaice
        print seaice.sum(skipna=True).values
    
        fig=pl.figure()
        seaice.plot()
        '''

        '''
        zlev=-1
        fig=pl.figure()
        pl.subplot(331)
        avg.isel(zt=zlev, Time=-1)['temp'].plot()
        pl.subplot(332)
        avg.isel(zt=zlev, Time=-1)['salt'].plot()
        pl.subplot(333)
        avg.isel(zt=zlev, Time=-1)['rho'].plot()
        pl.subplot(334)
        avg2.isel(zt=zlev, Time=-1)['temp'].plot()
        pl.subplot(335)
        avg2.isel(zt=zlev, Time=-1)['salt'].plot()
        pl.subplot(336)
        avg2.isel(zt=zlev, Time=-1)['rho'].plot()
        pl.subplot(337)
        xr.plot.pcolormesh(avg2.isel(zt=zlev, Time=-1)['temp']-avg.isel(zt=zlev, Time=-1)['temp'])
        pl.subplot(338)
        xr.plot.pcolormesh(avg2.isel(zt=zlev, Time=-1)['salt']-avg.isel(zt=zlev, Time=-1)['salt'])
        pl.subplot(339)
        xr.plot.pcolormesh(avg2.isel(zt=zlev, Time=-1)['rho']-avg.isel(zt=zlev, Time=-1)['rho'])


        fig=pl.figure()
        pl.subplot(341)
        avg.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot()
        avg.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot.contour()
        pl.subplot(342)
        avg.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp'].plot()
        pl.subplot(343)
        avg.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho'].plot.contourf(levels=50,add_colorbar=True, grid=True)
        pl.subplot(344)
        ot.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        pl.subplot(345)
        avg2.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot()
        avg2.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot.contour()
        pl.subplot(346)
        avg2.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp'].plot()
        pl.subplot(347)
        avg2.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho'].plot.contourf(levels=50,add_colorbar=True, grid=True)
        pl.subplot(348)
        ot2.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        pl.subplot(349)
        xr.plot.pcolormesh(avg2.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt']-avg.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'])
        pl.subplot(3,4,10)
        xr.plot.pcolormesh(avg2.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp']-avg.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp'])
        pl.subplot(3,4,11)
        xr.plot.contourf(avg2.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho']-avg.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho'] ,levels=50,add_colorbar=True, grid=True)
        pl.subplot(3,4,12)
        xr.plot.contourf(ot2.isel(time=-1)-ot.isel(time=-1) , levels=20, add_colorbar=True)
        '''

        fig=pl.figure()
        ax1=pl.subplot(211)
        levels = np.linspace(-15., 15, 20)
        ot.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        h_contour = ax1.contour(np.asarray(ot.lat), np.asarray(ot.depth), np.asarray(ot.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        ax2=pl.subplot(212)
        ot2.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        h_contour = ax2.contour(np.asarray(ot2.lat), np.asarray(ot2.depth), np.asarray(ot2.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)


        fig=pl.figure()
        for i in range(16):
                pl.subplot(4,4,i+1)
                ot2.isel(time=i*8).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)


def climatology():

        ### monthly average salt flux of 200y run
        flux_real = xr.open_dataarray("salt_forcing_monthly40l.nc")
        ### biweekly average salt flux of 200y run
        flux_2week = xr.open_dataarray("salt_forcing_biweekly40l.nc")


        x0 = 15#65
        y0 = 15#35

        ### monthly averages data of 10y run
        avgs = xr.open_dataset("s40lfl.0000.averages.nc")

        forc1 = []; t1 = []
        for i in range(len(avgs.Time)):
                forc1.append(avgs.isel(Time=i,xt=x0,yt=y0)['forc_salt_surface'])
                t1.append(avgs.Time[i])


        salt0 = []; forc0 = []; t = []
        for i in range(len(flux_real.Time)):
                forc0.append(flux_real.isel(Time=i,xt=x0,yt=y0))
                t.append(flux_real.Time[i])

        xt = flux_real.xt
        yt = flux_real.yt
        
        t_clim = range(12)
        clim = np.zeros((12,40,90))

        forc_mean = -3.8557348258389774e-08
        
        for j in range(12):
                for i in range(len(flux_real.Time)/12):
                        clim[j,:,:] += flux_real.isel(Time=i*12+j) /(len(flux_real.Time)/12) 

        clim -= forc_mean


        forc_clim = xr.DataArray(np.asarray(clim), coords=[t_clim, yt, xt], dims=['Time', 'yt_4deg', 'xt_4deg'])

        print 'Mean salinity forcing: ', forc_clim.mean(dim='xt_4deg').mean(dim='yt_4deg').mean(dim='Time').values

        for i in range(12):
                print forc_clim.isel(Time=i).mean(dim='xt_4deg').mean(dim='yt_4deg').values

        #forc_clim = xr.DataArray(np.asarray(clim), coords=[t_clim, yt, xt], dims=['Time', 'yt', 'xt'])


        fig=pl.figure()
        forc_clim.isel(Time=0, xt_4deg=slice(40,60)).mean(dim='xt_4deg').plot()

        forc_clim.isel(Time=6, xt_4deg=slice(40,60)).mean(dim='xt_4deg').plot()

        fig=pl.figure()
        pl.subplot(221)
        forc_clim.isel(Time=10).plot()
        pl.subplot(222)
        avgs.isel(Time=10)['forc_salt_surface'].plot()
        pl.subplot(223)
        #xr.plot.pcolormesh(forc_clim.isel(Time=10)-avgs.isel(Time=10)['forc_salt_surface'])


        #clim0=[] 
        #for i in range(199):
        #        for j in range(12):
        #                clim0.append(forc_clim.isel(Time=j,xt_4deg=x0,yt_4deg=y0))
                        #clim0.append(forc_clim.isel(Time=j,xt=x0,yt=y0))

        clim0=[]
        for i in range(200):
                for j in range(12):
                        clim0.append(forc_clim.isel(Time=j,xt_4deg=x0,yt_4deg=y0))

        fig=pl.figure()
        pl.subplot(211)
        pl.plot(t, forc0)
        #pl.plot(t2, forc2, 'x-')
        #pl.plot(t3, forc3, 'x-')
        pl.subplot(212)
        pl.plot(t1, forc1,'x-')
        #pl.plot(t2, forc2,'x-')
        #pl.plot(t3, clim0,'x-')
        #pl.plot(t4, forc4)
        pl.plot(t, clim0,'x-')

        forc_all = xr.open_dataset("forcing_1deg_global_interp_hosing.nc")
        forc_all['salt_flux'] = forc_clim
        forc_all = forc_all.fillna(0.)

        forc_all.to_netcdf('forcing_1deg_global_interp_hos_flux_total0.nc')



        pl.show()

        '''
        ### make new data array with monthly average flux from biweekly data
        for i in range(len(flux_2week.Time)/2):
                if i==12:
                        data = flux_2week.isel(Time=slice(2*i-1,2*i+1)).mean(dim='Time')
                        data_newcoord = data.assign_coords(Time=flux_2week.Time[2*i-1])
                        data_expanded = data_newcoord.expand_dims('Time')
                        flux_month = data_expanded
                elif i>12:
                        data = flux_2week.isel(Time=slice(2*i-1,2*i+1)).mean(dim='Time')
                        data_newcoord = data.assign_coords(Time=flux_2week.Time[2*i-1])
                        data_expanded = data_newcoord.expand_dims('Time')
                        flux_month=xr.concat((flux_month,data_expanded),dim='Time')

        clim = np.zeros((12,40,90))
        t_clim = range(12)
        for j in range(12):
                for i in range(len(flux_month.Time)/12):
                        clim[j,:,:] += flux_month.isel(Time=i*12+j) /(len(flux_month.Time)/12)
        forc_clim = xr.DataArray(np.asarray(clim), coords=[t_clim, yt, xt], dims=['Time', 'yt_4deg', 'xt_4deg'])


        five = xr.open_dataarray("salt_forcing_1day40l_2.nc")

        forc2 = []; t2 = []
        for i in range(len(flux_2week.Time)):
                forc2.append(flux_2week.isel(Time=i,xt=x0,yt=y0))
                t2.append(flux_2week.Time[i])

        forc3 = []; t3 = []
        for i in range(len(flux_month.Time)):
                forc3.append(flux_month.isel(Time=i,xt=x0,yt=y0))
                t3.append(flux_month.Time[i])


        forc4 = []; t4 = []
        for i in range(len(five.Time)):
                forc4.append(five.isel(Time=i,xt=x0,yt=y0))
                t4.append(five.Time[i])
        '''


if __name__ == '__main__':
        MAIN()
