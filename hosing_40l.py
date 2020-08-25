import matplotlib.pyplot as pl
import numpy as np
import xarray as xr
import matplotlib as mpl

pl.rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'size'   : 17})

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

mpl.rcParams['axes.xmargin'] = 0.03
mpl.rcParams['axes.ymargin'] = 0.03

mpl.rcParams['axes.unicode_minus'] = False

def MAIN():

        hosing()
        #somoc()

        pl.show()


def somoc():

        ot_ctrl_360 = xr.open_dataarray("overturning_ctrl40l_s360.nc")

        amoc = xr.open_dataarray("amoc_ctrl40l_s360.nc")
        somoc = xr.open_dataarray("somoc_ctrl40l_s360.nc")

        fig=pl.figure()
        ax1=pl.subplot(111)
        levels = np.arange(-11., 11, 1)
        levels0 = np.arange(-1.1e7, 1.1e7, 1000000.)
        ot_ctrl_360.isel(time=-1).plot.contourf(vmin=-1.1e7, levels=levels0, add_colorbar=False)#
        h_contour = ax1.contour(np.asarray(ot_ctrl_360.lat), np.asarray(ot_ctrl_360.depth), np.asarray(ot_ctrl_360.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        #ax2=pl.subplot(212)
        #ot_360.sel(time=2.736e+06).plot.contourf(vmin=-1.1e7, levels=levels0, add_colorbar=True)
        #h_contour = ax2.contour(np.asarray(ot_360.lat), np.asarray(ot_360.depth), np.asarray(ot_360.sel(time=2.736e+06))/1000000., colors='k',levels=levels)
        #h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        fig=pl.figure()
        ax1=pl.subplot(121)
        levels = np.arange(-11., 11, 1)
        levels0 = np.arange(-1.1e7, 1.1e7, 1000000.)
        somoc.isel(time=-1).plot.contourf(vmin=-1.1e7, levels=levels0,add_colorbar=False)
        h_contour = ax1.contour(np.asarray(somoc.lat), np.asarray(somoc.depth), np.asarray(somoc.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        ax2=pl.subplot(122)
        amoc.isel(time=-1).plot.contourf(vmin=-1.1e7, levels=levels0,add_colorbar=False)
        h_contour = ax2.contour(np.asarray(amoc.lat), np.asarray(amoc.depth), np.asarray(amoc.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)




def hosing():

        '''
        snap = xr.open_dataset("ctrl40l.0000.snapshot.nc")

        zlev=-1
        fig=pl.figure()
        for i in range(9):
                pl.subplot(3,3,i+1)
                snap.isel(Time=111+i*10,zt=zlev)['v'].plot()


        snap = xr.open_dataset("ctrl40l.0073.snapshot.nc")

        print snap

        zlev=-1
        fig=pl.figure()
        pl.subplot(331)
        snap.isel(Time=0,zt=zlev)['salt'].plot()
        pl.subplot(332)
        snap.isel(Time=0,zt=zlev)['temp'].plot()
        pl.subplot(333)
        snap.isel(Time=0,zt=zlev)['rho'].plot()
        pl.subplot(334)
        snap.isel(Time=0,zt=zlev)['u'].plot()
        pl.subplot(335)
        snap.isel(Time=0,zt=zlev)['v'].plot()
        pl.subplot(336)
        snap.isel(Time=0,zw=zlev)['w'].plot()
        pl.subplot(337)
        snap.isel(Time=0)['psi'].plot()

        pl.subplot(338)
        snap.isel(Time=0,zt=zlev)['Hd'].plot()

        pl.subplot(339)
        snap.isel(Time=0,zw=zlev)['Nsqr'].plot()
        '''

        #avgs = xr.open_dataset("ctrl40l.0039.averages.nc")
        avgs = xr.open_dataset("c40l360.0039.averages.nc")
        #avgs_forc = xr.open_dataset("salt40l.0029.averages.nc")
        avgs_forc = xr.open_dataset("s40l360.0035.averages.nc")

        '''
        salt_vals = np.empty((90,40,40))
        salt_vals = avgs.isel(Time=0)['salt'].values
        count=0
        for i in range(90):
                for j in range(40):
                        for k in range(40):
                                if np.isnan(salt_vals[k,j,i]):
                                        count+=1
                                        print 'garrrr'
        print count
        '''
        #y_spacing = avgs.yt.values
        #print y_spacing[1:]-y_spacing[:-1]

        print(avgs.isel(Time=0)['surface_taux'].values)

        fig=pl.figure()
        pl.subplot(221)
        avgs.isel(Time=-1)['surface_taux'].plot()
        pl.subplot(222)
        avgs.isel(Time=-1)['surface_tauy'].plot()
        pl.subplot(223)
        avgs.isel(Time=-1)['psi'].plot.contourf(levels=30, add_colorbar=True)#vmin=-1.5e7,
        
        
        ot = xr.open_dataarray("overturning_ctrl40l.nc")
        ot_forc = xr.open_dataarray("overturning_salt40l.nc")
        salt_temp = xr.open_dataset("salt_temp_ctrl40l.nc")
        salt_temp_forc = xr.open_dataset("salt_temp_salt40l.nc")

        ot_ctrl_360 = xr.open_dataarray("overturning_ctrl40l_s360.nc")
        ot_360 = xr.open_dataarray("overturning_s40l360.nc")
        salt_temp_ctrl_360 = xr.open_dataset("salt_temp_ctrl40l_s360.nc")
        salt_temp_360 = xr.open_dataset("salt_temp_s40l360.nc")

        ot_360d = xr.open_dataarray("overturning_s40l360d.nc")
        salt_temp_360d = xr.open_dataset("salt_temp_s40l360d.nc")

        ot_360u = xr.open_dataarray("overturning_s40l360u.nc")
        salt_temp_360u = xr.open_dataset("salt_temp_s40l360u.nc")

        ot_360q = xr.open_dataarray("overturning_s40l360q.nc")
        salt_temp_360q = xr.open_dataset("salt_temp_s40l360q.nc")

        #print salt_temp_360q


        t_amoc = []; amoc = []; amoc_hyst = [ot.isel(time=-1, depth=slice(0,32), lat=slice(24,40)).max()/1000000.]
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i]/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        for i in range(len(ot_forc.time)):
                t0 = ot_forc.time[i]/360.
                t_amoc.append(t0)
                amoc.append(ot_forc.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                if (t0-4000)%300==0:
                        print(t0)
                        amoc_hyst.append(ot_forc.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        t_amoc360 = []; amoc360 = [];amoc_hyst360 = [ot_360.isel(time=0, depth=slice(0,32), lat=slice(24,40)).max()/1000000.]
        for i in range(len(ot_360.time)):
                t0 = ot_360.time[i]/360.
                t_amoc360.append(t0)
                amoc360.append(ot_360.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                if (t0-4000)%300==0:
                        print(t0)
                        #amoc_hyst360.append(ot_360.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                        amoc_vals = np.asarray([ot_360.isel(time=i-j, depth=slice(0,32), lat=slice(24,40)).max()/1000000. for j in range(5)])
                        amoc_hyst360.append(np.mean(amoc_vals))

        t_amoc360d = []; amoc360d = [];amoc_hyst360d = [ot_360d.isel(time=0, depth=slice(0,32), lat=slice(24,40)).max()/1000000.]
        for i in range(len(ot_360d.time)):
                t0 = ot_360d.time[i]/360.
                t_amoc360d.append(t0)
                amoc360d.append(ot_360d.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                if (t0-5500)%300==0:
                        print(t0)
                        amoc_hyst360d.append(ot_360d.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        t_amoc360u = []; amoc360u = [];amoc_hyst360u = [ot_360u.isel(time=0, depth=slice(0,32), lat=slice(24,40)).max()/1000000.]
        for i in range(len(ot_360u.time)):
                t0 = ot_360u.time[i]/360.
                t_amoc360u.append(t0)
                amoc360u.append(ot_360u.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                if (t0-7900)%300==0:
                        print(t0)
                        amoc_hyst360u.append(ot_360u.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)

        t_amoc360q = []; amoc360q = [];amoc_hyst360q = [ot_360q.isel(time=0, depth=slice(0,32), lat=slice(24,40)).max()/1000000.]
        t_hyst = [salt_temp_360q.isel(time=0)['sst_NA']]
        for i in range(len(ot_360q.time)):
                t0 = ot_360q.time[i]/360.
                t_amoc360q.append(t0)
                amoc360q.append(ot_360q.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                if (t0-5500)%300==0:
                        print(t0)
                        #amoc_hyst360q.append(ot_360q.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
                        amoc_vals = np.asarray([ot_360q.isel(time=i-j, depth=slice(0,32), lat=slice(24,40)).max()/1000000. for j in range(5)])
                        amoc_hyst360q.append(np.mean(amoc_vals))
                        t_hyst.append(salt_temp_360q.isel(time=slice(i-10,i)).mean(dim='time')['sst_NA'])
                        #t_hyst.append(salt_temp_360q.isel(time=i)['sst_NA'])

        F_hyst = np.linspace(0.,1.7139,11)#0.,0.000006
        F_hyst = np.concatenate([F_hyst, np.linspace(1.5425,-0.3428,12)])#0.0000054,-0.0000012

        F_hyst360 = np.linspace(0.,2.0567,13)#0.,0.0000072
        F_hyst360 = np.concatenate([F_hyst360, np.linspace(1.8853,-0.3428,14)])#0.0000066,-0.0000012

        F_hyst360q = np.linspace(0.85694,2.0567,29)#0.000003,0.0000072
        F_hyst360q = np.concatenate([F_hyst360[:5],F_hyst360q, np.linspace(2.0138,0.7712,30),F_hyst360[-7:-2]])#0.00000705,0.0000027

        amoc_hyst360q = np.concatenate([amoc_hyst360[:5],np.asarray(amoc_hyst360q),amoc_hyst360[-7:-2]])

        #F_hyst360d = np.linspace(0.000003,0.0000072,15)
        #F_hyst360d = np.concatenate([F_hyst360d, np.linspace(0.0000069,0.,24)])
        #F_hyst360u = np.linspace(0.0000066,0.0,23)


        fig=pl.figure(figsize=(8,9))
        pl.subplots_adjust(left=0.1, bottom=0.09, right=0.97, top=0.96, wspace=0.0, hspace=0.32)
        ax1=pl.subplot(211)
        ax1.set_title('30-day salt relaxation time scale')
        pl.plot(F_hyst[:11], amoc_hyst[:11], 'o-', color='black')
        pl.plot(F_hyst[11:len(amoc_hyst)-2], amoc_hyst[11:-2], 'o-', color='crimson')
        ax1.spines["bottom"].set_position(("outward", 4))
        ax1.spines["left"].set_position(("outward", 4))
        pl.ylabel('Max. AMOC (Sv)')
        pl.xlabel('Freshwater forcing (PSU)')
        ax2=pl.subplot(212)
        ax2.set_title('360-day salt relaxation time scale')
        #pl.plot(F_hyst360[:13], amoc_hyst360[:13], 'x--')
        #pl.plot(F_hyst360[13:len(amoc_hyst360)-2], amoc_hyst360[13:-2], 'x--')
        #pl.plot(10000.*F_hyst360d[:len(amoc_hyst360d)], amoc_hyst360d, 's--', color='black', markerfacecolor='none')
        #pl.plot(10000.*F_hyst360u[:len(amoc_hyst360u)], amoc_hyst360u, 'x--', color='black')

        pl.plot(F_hyst360q[:34], amoc_hyst360q[:34], 'o-', color='black')
        pl.plot(F_hyst360q[34:len(amoc_hyst360q)], amoc_hyst360q[34:], 'o-', color='crimson')
        ax2.spines["bottom"].set_position(("outward", 4))
        ax2.spines["left"].set_position(("outward", 4))
        pl.xlabel('Freshwater forcing (PSU)')
        pl.ylabel('Max. AMOC (Sv)')

        fig=pl.figure()
        pl.plot(F_hyst360q[5:34], t_hyst[:29], 'x--', color='crimson')
        pl.plot(F_hyst360q[34:len(t_hyst)+5], t_hyst[29:], 'x--', color='gray')
        pl.xlabel('Freshwater forcing (PSU)')
        pl.ylabel('NA Temperature')
        
        fig=pl.figure()
        pl.subplot(511)
        pl.plot(t_amoc, amoc)
        pl.plot(t_amoc360, amoc360)
        pl.subplot(512)
        salt_temp['salt'].plot()
        salt_temp_forc['salt'].plot()
        salt_temp_360['salt'].plot()
        pl.subplot(513)
        #salt_temp['temp'].plot()
        #salt_temp_forc['temp'].plot()
        salt_temp_360q['sst_NA'].plot()
        pl.subplot(514)
        salt_temp['rhoN'].plot()
        salt_temp['rhoS'].plot()
        salt_temp_forc['rhoN'].plot()
        salt_temp_forc['rhoS'].plot()
        #(salt_temp['rhoN']-salt_temp['rhoS']).plot()
        pl.subplot(515)
        salt_temp['salt_forc'].plot()
        salt_temp_forc['salt_forc'].plot()

        #time=slice(None, 2)

        fig=pl.figure()
        salt_temp['salt'].plot()
        salt_temp_forc['salt'].plot()
        salt_temp_ctrl_360['salt'].plot()
        salt_temp_360['salt'].plot()

        zlev = -1
        fig=pl.figure()
        pl.subplot(331)
        avgs.isel(Time=-1).isel(zt=zlev)['temp'].plot()
        pl.subplot(332)
        avgs.isel(Time=-1).isel(zt=zlev)['salt'].plot()
        pl.subplot(333)
        avgs.isel(Time=-1)['psi'].plot.contourf(vmin=-1.5e7,levels=25, add_colorbar=True)

        pl.subplot(334)
        avgs_forc.isel(Time=-1).isel(zt=zlev)['temp'].plot()
        pl.subplot(335)
        avgs_forc.isel(Time=-1).isel(zt=zlev)['salt'].plot()
        pl.subplot(336)
        avgs_forc.isel(Time=-1)['psi'].plot.contourf(vmin=-1.5e7,levels=25, add_colorbar=True)

        pl.subplot(337)
        xr.plot.pcolormesh(avgs_forc.isel(Time=-1).isel(zt=zlev)['temp']-avgs.isel(Time=-1).isel(zt=zlev)['temp'])
        pl.subplot(338)
        xr.plot.pcolormesh(avgs_forc.isel(Time=-1).isel(zt=zlev)['salt']-avgs.isel(Time=-1).isel(zt=zlev)['salt'])
        pl.subplot(339)
        xr.plot.contourf(avgs_forc.isel(Time=-1)['psi']-avgs.isel(Time=-1)['psi'],vmin=-1.5e7,levels=25, add_colorbar=True)

        fig=pl.figure()
        pl.subplot(221)
        xr.plot.pcolormesh(avgs_forc.isel(Time=-1).isel(zt=zlev)['temp']-avgs.isel(Time=-1).isel(zt=zlev)['temp'])
        pl.subplot(222)
        xr.plot.pcolormesh(avgs_forc.isel(Time=-1).isel(zt=zlev)['salt']-avgs.isel(Time=-1).isel(zt=zlev)['salt'])
        pl.subplot(223)
        xr.plot.contourf(avgs_forc.isel(Time=-1)['psi']-avgs.isel(Time=-1)['psi'],vmin=-1.5e7,levels=25, add_colorbar=True)


        fig=pl.figure()
        pl.subplot(341)
        avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot()
        avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot.contour()
        pl.subplot(342)
        avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp'].plot()
        pl.subplot(343)
        avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho'].plot.contourf(levels=50,add_colorbar=True, grid=True)
        pl.subplot(344)
        ot.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        pl.subplot(345)
        avgs_forc.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot()
        avgs_forc.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'].plot.contour()
        pl.subplot(346)
        avgs_forc.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp'].plot()
        pl.subplot(347)
        avgs_forc.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho'].plot.contourf(levels=50,add_colorbar=True, grid=True)
        pl.subplot(348)
        ot_forc.isel(time=-1).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        pl.subplot(349)
        xr.plot.pcolormesh(avgs_forc.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt']-avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['salt'])
        pl.subplot(3,4,10)
        xr.plot.pcolormesh(avgs_forc.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp']-avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['temp'])
        pl.subplot(3,4,11)
        xr.plot.contourf(avgs_forc.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho']-avgs.isel(Time=-1).sel(xt=slice(284.,376.)).mean(dim='xt')['rho'] ,levels=50,add_colorbar=True, grid=True)
        pl.subplot(3,4,12)
        xr.plot.contourf(ot_forc.isel(time=-1)-ot.isel(time=-1) , levels=20, add_colorbar=True)


        fig=pl.figure()
        pl.subplot(221)
        avgs.isel(Time=-1)['forc_salt_surface'].plot()
        pl.subplot(222)
        avgs_forc.isel(Time=-1)['forc_salt_surface'].plot()
        pl.subplot(223)
        xr.plot.pcolormesh(avgs_forc.isel(Time=-1)['forc_salt_surface']-avgs.isel(Time=-1)['forc_salt_surface'])

        fig=pl.figure()
        xr.plot.pcolormesh(avgs_forc.isel(Time=-1)['forc_salt_surface']-avgs.isel(Time=-1)['forc_salt_surface'])


        fig=pl.figure()
        ax1=pl.subplot(211)
        levels = np.arange(-11., 11, 1)
        levels0 = np.arange(-1.1e7, 1.1e7, 1000000.)
        ot_ctrl_360.isel(time=-1).plot.contourf(vmin=-1.1e7, levels=levels0, add_colorbar=True)
        h_contour = ax1.contour(np.asarray(ot_ctrl_360.lat), np.asarray(ot_ctrl_360.depth), np.asarray(ot_ctrl_360.isel(time=-1))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        ax2=pl.subplot(212)
        ot_360.sel(time=2.736e+06).plot.contourf(vmin=-1.1e7, levels=levels0, add_colorbar=True)
        h_contour = ax2.contour(np.asarray(ot_360.lat), np.asarray(ot_360.depth), np.asarray(ot_360.sel(time=2.736e+06))/1000000., colors='k',levels=levels)
        h_contour.clabel(fontsize=14, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

        
        ### calibrate forcing in terms of PSU.
        '''
        forcing = xr.open_dataset("forcing_1deg_global_interp_hosing.nc")
        forcing = forcing.where(forcing.sss!=0., np.nan)
        Sref = forcing.sel(xt=slice(314.,334.),yt=slice(26.,66.)).mean(dim='xt').mean(dim='yt').mean(dim='Time')['sss']
        print Sref, salt_temp.sel(time=3517200)['salt']

        print dz/(30.*86400.)*(Sref - 0.0001*86400.*30./dz -  salt_temp.sel(time=3517200)['salt'])

        mask = np.nan_to_num(np.asarray(snap_forc['maskFrrampup']))

        pl.subplot(4,3,12)
        xr.plot.pcolormesh(dz/(30.*86400.)*((forcing.isel(Time=0)['sss']) - snap_forc.isel(Time=-1).sel(zt=-35.)['salt'] - mask*0.0001*30.*86400./dz), vmin=-0.000035)
        '''
        


        
        '''
        fig=pl.figure()
        for i in range(9):
                pl.subplot(3,3,i+1)
                #ot.isel(time=6*i).plot(vmin=-1.5e7, levels=20)
                #ot.isel(time=10*i).plot.contour(vmin=-1.5e7, levels=20, colorsdiscrete='black',add_colorbar=True)
                ot.isel(time=40*i).plot.contourf(vmin=-1.5e7, levels=20, add_colorbar=True)
        '''



if __name__ == '__main__':
        MAIN()
