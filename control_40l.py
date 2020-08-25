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


def MAIN():
        comp_relax()
        #control_run()
        pl.show()

def control_run():

        ot = xr.open_dataarray("overturning_ctrl40l.nc")
        salt_temp = xr.open_dataset("salt_temp_ctrl40l.nc")        
        avgs = xr.open_dataset("ctrl40l.0039.averages.nc")

        ###Mark regions with sea ice
        avgs['temp'] = avgs['temp'].fillna(10.)
        avgs['temp'] = avgs['temp'].where(avgs['temp']>=-1.8, -15.)
        avgs['temp'] = avgs['temp'].where(avgs['temp']!=10., np.nan)

        t_amoc = []; amoc = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i]/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        fig=pl.figure()
        pl.subplot(511)
        pl.plot(t_amoc, amoc)
        pl.legend(loc='best'); pl.ylabel('Max. AMOC (Sv.)')
        pl.subplot(512)
        salt_temp['rhoN'].plot()
        salt_temp['rhoS'].plot()
        pl.subplot(513)
        salt_temp['salt_tot'].plot()
        pl.subplot(514)
        salt_temp['salt_forc_tot'].plot()
        pl.subplot(515)
        salt_temp['seaice'].plot()

        fig=pl.figure()
        ax1=pl.subplot(111)
        levels = np.linspace(-15., 15, 20)
        ot.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        h_contour = ax1.contour(np.asarray(ot.lat), np.asarray(ot.depth), np.asarray(ot.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        zlev = -1
        fig=pl.figure()
        pl.subplot(221)
        avgs.isel(Time=-1).isel(zt=zlev)['temp'].plot()
        pl.subplot(222)
        avgs.isel(Time=-1).isel(zt=zlev)['salt'].plot()
        ax1=pl.subplot(223)
        avgs.isel(Time=-1)['psi'].plot.contourf(vmin=-1.5e7,levels=20, add_colorbar=True)
        levels = np.linspace(-15., 15, 20)
        h_contour = ax1.contour(np.asarray(avgs.xu), np.asarray(avgs.yu), np.asarray(avgs.isel(Time=-1)['psi'])/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        pl.subplot(224)
        avgs.isel(Time=-1).isel(zt=zlev)['forc_salt_surface'].plot()



def comp_relax():

        ot = xr.open_dataarray("overturning_ctrl40l.nc")
        salt_temp = xr.open_dataset("salt_temp_ctrl40l.nc")

        ot360 = xr.open_dataarray("overturning_ctrl40l_s360.nc")
        salt_temp360 = xr.open_dataset("salt_temp_ctrl40l_s360.nc")

        ot360w = xr.open_dataarray("overturning_ctrl40l_s360w.nc")
        salt_temp360w = xr.open_dataset("salt_temp_ctrl40l_s360w.nc")

        avgs = xr.open_dataset("ctrl40l.0039.averages.nc")
        avg360 = xr.open_dataset("c40l360.0039.averages.nc")
        avg360w = xr.open_dataset("c40l360w.0039.averages.nc")


        ###Mark regions with sea ice
        avgs['temp'] = avgs['temp'].fillna(10.)
        avgs['temp'] = avgs['temp'].where(avgs['temp']>=-1.8, -15.)
        avgs['temp'] = avgs['temp'].where(avgs['temp']!=10., np.nan)

        avg360['temp'] = avg360['temp'].fillna(10.)
        avg360['temp'] = avg360['temp'].where(avg360['temp']>=-1.8, -15.)
        avg360['temp'] = avg360['temp'].where(avg360['temp']!=10., np.nan)


        t_amoc = []; amoc = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i]/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        t_amoc360 = []; amoc360 = []
        for i in range(len(ot360.time)):
                t_amoc360.append(ot360.time[i]/360.)
                amoc360.append(ot360.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        t_amoc360w = []; amoc360w = []
        for i in range(len(ot360w.time)):
                t_amoc360w.append(ot360w.time[i]/360.)
                amoc360w.append(ot360w.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        fig=pl.figure()
        pl.subplot(711)
        pl.plot(t_amoc, amoc, label='S-relax 30d')
        pl.plot(t_amoc360, amoc360, label='S-relax 360d')
        pl.plot(t_amoc360w, amoc360w, label='S-360d WIND')
        pl.legend(loc='best')
        pl.subplot(712)
        salt_temp['salt'].plot()
        salt_temp360['salt'].plot()
        pl.subplot(713)
        salt_temp['temp'].plot()
        salt_temp360['temp'].plot()
        pl.subplot(714)
        salt_temp['rhoN'].plot()
        salt_temp['rhoS'].plot()
        salt_temp360['rhoN'].plot()
        salt_temp360['rhoS'].plot()
        pl.subplot(715)
        salt_temp['salt_tot'].plot()
        salt_temp360['salt_tot'].plot()
        pl.subplot(716)
        salt_temp['salt_forc_tot'].plot()
        salt_temp360['salt_forc_tot'].plot()
        pl.subplot(717)
        #salt_temp['seaice'].plot()
        salt_temp360['seaice'].plot()


        zlev = -1
        fig=pl.figure()
        pl.subplot(331)
        avgs.isel(Time=-1).isel(zt=zlev)['temp'].plot()
        pl.subplot(332)
        avgs.isel(Time=-1).isel(zt=zlev)['salt'].plot()
        pl.subplot(333)
        avgs.isel(Time=-1)['psi'].plot.contourf(vmin=-1.5e7,levels=25, add_colorbar=True)

        pl.subplot(334)
        avg360.isel(Time=-1).isel(zt=zlev)['temp'].plot()
        pl.subplot(335)
        avg360.isel(Time=-1).isel(zt=zlev)['salt'].plot()
        pl.subplot(336)
        avg360.isel(Time=-1)['psi'].plot.contourf(vmin=-1.5e7,levels=25, add_colorbar=True)

        pl.subplot(337)
        xr.plot.pcolormesh(avg360.isel(Time=-1).isel(zt=zlev)['temp']-avgs.isel(Time=-1).isel(zt=zlev)['temp'])
        pl.subplot(338)
        xr.plot.pcolormesh(avg360.isel(Time=-1).isel(zt=zlev)['salt']-avgs.isel(Time=-1).isel(zt=zlev)['salt'])
        pl.subplot(339)
        xr.plot.contourf(avg360.isel(Time=-1)['psi']-avgs.isel(Time=-1)['psi'],vmin=-1.5e7,levels=25, add_colorbar=True)


        fig=pl.figure()
        pl.subplot(331)
        avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot()
        avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot.contour()
        pl.subplot(332)
        avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp'].plot()
        pl.subplot(333)
        avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho'].plot.contourf(levels=50,add_colorbar=True, grid=True)
        pl.subplot(334)
        avg360.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot()
        avg360.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot.contour()
        pl.subplot(335)
        avg360.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp'].plot()
        pl.subplot(336)
        avg360.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho'].plot.contourf(levels=50,add_colorbar=True, grid=True)
        pl.subplot(337)
        xr.plot.pcolormesh(avg360.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt']-avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'])
        pl.subplot(338)
        xr.plot.pcolormesh(avg360.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp']-avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp'])
        pl.subplot(339)
        xr.plot.contourf(avg360.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho']-avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho'] ,levels=50,add_colorbar=True, grid=True)



        fig=pl.figure()
        levels = np.linspace(-15., 15, 25)
        '''
        ax1=pl.subplot(221)
        ot.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        h_contour = ax1.contour(np.asarray(ot.lat), np.asarray(ot.depth), np.asarray(ot.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
        '''

        ax2=pl.subplot(221)
        pl.gca().set_title('No Wind')
        ot360.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=25, add_colorbar=True)
        h_contour = ax2.contour(np.asarray(ot360.lat), np.asarray(ot360.depth), np.asarray(ot360.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        ax3=pl.subplot(222)
        ot360w.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=25, add_colorbar=True)
        h_contour = ax3.contour(np.asarray(ot360w.lat), np.asarray(ot360w.depth), np.asarray(ot360w.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        pl.subplot(223)
        avg360.isel(Time=-1)['psi'].plot.contourf(levels=40, add_colorbar=True)#vmin=-1.5e7,

        pl.subplot(224)
        avg360w.isel(Time=-1)['psi'].plot.contourf(levels=40, add_colorbar=True)


        ax2=pl.subplot(211)
        ot360.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=25, add_colorbar=True)
        h_contour = ax2.contour(np.asarray(ot360.lat), np.asarray(ot360.depth), np.asarray(ot360.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
        pl.gca().set_title('No Wind')

        ax3=pl.subplot(212)
        ot360w.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=25, add_colorbar=True)
        h_contour = ax3.contour(np.asarray(ot360w.lat), np.asarray(ot360w.depth), np.asarray(ot360w.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
        pl.gca().set_title('With Wind')


        amoc = xr.open_dataarray("amoc_ctrl40l_s360w.nc")
        somoc = xr.open_dataarray("somoc_ctrl40l_s360w.nc")
        '''
        fig=pl.figure()
        ax2=pl.subplot(211)
        amoc.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=25, add_colorbar=True)
        h_contour = ax2.contour(np.asarray(amoc.lat), np.asarray(amoc.depth), np.asarray(amoc.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        ax2=pl.subplot(212)
        somoc.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=25, add_colorbar=True)
        h_contour = ax2.contour(np.asarray(somoc.lat), np.asarray(somoc.depth), np.asarray(somoc.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
        '''

        fig=pl.figure()
        ax1=pl.subplot(121)
        levels = np.arange(-15., 15, 1)
        levels0 = np.arange(-1.5e7, 1.5e7, 1000000.)
        somoc.isel(time=-1).plot.contourf(vmin=-1.1e7, levels=levels0,add_colorbar=False)
        h_contour = ax1.contour(np.asarray(somoc.lat), np.asarray(somoc.depth), np.asarray(somoc.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        ax2=pl.subplot(122)
        amoc.isel(time=-1).plot.contourf(vmin=-1.1e7, levels=levels0,add_colorbar=False)
        h_contour = ax2.contour(np.asarray(amoc.lat), np.asarray(amoc.depth), np.asarray(amoc.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)


if __name__ == '__main__':
        MAIN()
