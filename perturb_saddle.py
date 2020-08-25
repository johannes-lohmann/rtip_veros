import matplotlib.pyplot as pl
import numpy as np
import xarray as xr
#from eofs.xarray import Eof
#from matplotlib.collections import LineCollection
import h5py
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

        ot_c01 = xr.open_dataarray("overturning_rtipSc01.nc")
        ot_c015 = xr.open_dataarray("overturning_rtipSc015.nc")
        ot_c0055 = xr.open_dataarray("overturning_rtipSc0055.nc")
        ot_c001 = xr.open_dataarray("overturning_rtipSc001.nc")

        t_amoc_c01, amoc_c01 = amoc_timeseries(ot_c01)
        t_amoc_c015, amoc_c015 = amoc_timeseries(ot_c015)
        t_amoc_c0055, amoc_c0055 = amoc_timeseries(ot_c0055)
        t_amoc_c001, amoc_c001 = amoc_timeseries(ot_c001)


        ot_c006 = xr.open_dataarray("overturning_rtipSc006.nc")
        salt_temp_006 = xr.open_dataset("salt_temp_rtipSc006.nc")
        t_amoc_c006, amoc_c006 = amoc_timeseries(ot_c006)

        salt_temp_p0 = xr.open_dataset("salt_temp_rtipSc006_pert0.nc")
        ot_p0 = xr.open_dataarray("overturning_rtipSc006_pert0.nc")
        t_amoc_p0, amoc_p0 = amoc_timeseries(ot_p0)

        salt_temp_p1 = xr.open_dataset("salt_temp_rtipSc006_pert1.nc")
        ot_p1 = xr.open_dataarray("overturning_rtipSc006_pert1.nc")
        t_amoc_p1, amoc_p1 = amoc_timeseries(ot_p1)

        salt_temp_p2 = xr.open_dataset("salt_temp_rtipSc006_pert2.nc")
        ot_p2 = xr.open_dataarray("overturning_rtipSc006_pert2.nc")
        t_amoc_p2, amoc_p2 = amoc_timeseries(ot_p2)

        salt_temp_p3 = xr.open_dataset("salt_temp_rtipSc006_pert3.nc")
        ot_p3 = xr.open_dataarray("overturning_rtipSc006_pert3.nc")
        t_amoc_p3, amoc_p3 = amoc_timeseries(ot_p3)

        print salt_temp_006

        fig=pl.figure()
        pl.subplot(311)
        pl.plot(t_amoc_c006, amoc_c006)
        pl.plot(t_amoc_p0, amoc_p0)
        pl.plot(t_amoc_p1, amoc_p1)
        pl.plot(t_amoc_p2, amoc_p2)
        pl.plot(t_amoc_p3, amoc_p3)
        pl.subplot(312)
        salt_temp_006['salt'].plot()
        salt_temp_p0['salt'].plot()
        salt_temp_p1['salt'].plot()
        salt_temp_p2['salt'].plot()
        salt_temp_p3['salt'].plot()
        pl.subplot(313)
        salt_temp_006['sst_NA'].plot()
        salt_temp_p0['sst_NA'].plot()
        salt_temp_p1['sst_NA'].plot()
        salt_temp_p2['sst_NA'].plot()
        salt_temp_p3['sst_NA'].plot()

        #perturb_ic()
        #test_pert()



        fig=pl.figure(figsize=(10,12))
        pl.subplots_adjust(left=0.08, bottom=0.07, right=0.97, top=0.98, wspace=0.0, hspace=0.06)
        ax=pl.subplot(511)
        pl.plot(t_amoc_c015, amoc_c015, label='150y', color='black')
        pl.legend(loc='best'); pl.xlim(9000,19300);pl.ylim(2.0,7.2)
        ax.spines['bottom'].set_visible(False); ax.get_xaxis().set_visible(False)
        ax.spines["left"].set_position(("outward", 3));pl.ylabel('Max. AMOC (Sv)')
        ax=pl.subplot(512)
        pl.plot(t_amoc_c01, amoc_c01, label='100y', color='black')
        pl.legend(loc='best'); pl.xlim(9000,19300);pl.ylim(2.0,7.2)
        ax.spines['bottom'].set_visible(False); ax.get_xaxis().set_visible(False)
        ax.spines["left"].set_position(("outward", 3));pl.ylabel('Max. AMOC (Sv)')
        ax=pl.subplot(513)
        pl.plot(t_amoc_c006, amoc_c006, label='60y', color='black')
        pl.plot(t_amoc_p0, amoc_p0)
        pl.plot(t_amoc_p1, amoc_p1)
        pl.plot(t_amoc_p2, amoc_p2)
        pl.plot(t_amoc_p3, amoc_p3)
        pl.legend(loc='best'); pl.xlim(9000,19300);pl.ylim(2.0,7.2)
        ax.spines['bottom'].set_visible(False); ax.get_xaxis().set_visible(False)
        ax.spines["left"].set_position(("outward", 3));pl.ylabel('Max. AMOC (Sv)')
        ax=pl.subplot(514)
        pl.plot(t_amoc_c0055, amoc_c0055, label='55y', color='black')
        pl.legend(loc='best'); pl.xlim(9000,19300);pl.ylim(2.0,7.2)
        ax.spines['bottom'].set_visible(False); ax.get_xaxis().set_visible(False)
        ax.spines["left"].set_position(("outward", 3));pl.ylabel('Max. AMOC (Sv)')
        ax=pl.subplot(515)
        pl.plot(t_amoc_c001, amoc_c001, label='10y', color='black')
        pl.legend(loc='best'); pl.xlim(9000,19300);pl.ylim(2.0,7.2)
        ax.spines["bottom"].set_position(("outward", 3))
        #ax.spines["bottom"].set_position(("data", 0))
        #ax.spines["left"].set_position(("axes", .3))
        ax.spines["left"].set_position(("outward", 3))
        pl.xlabel('Simulation time (years)');pl.ylabel('Max. AMOC (Sv)')

        pl.show()

def test_pert():
        f = h5py.File('rtipSc006.0039.restart_pert2.h5', 'r+')
        f0 = h5py.File('rtipSc006.0039.restart.h5', 'r+')

        snap = f['snapshot']
        snap0 = f0['snapshot']

        temp = snap['temp'][:]
        temp0 = snap0['temp'][:]

        fig=pl.figure()
        pl.pcolormesh(temp[:,:,-1,0]-temp0[:,:,-1,0])
        pl.colorbar()


def perturb_ic():
        '''
        pert0 -> T+, S+
        pert1 -> T-, S+
        pert2 -> T+, S-
        pert3 -> T-, S-
        '''


        f = h5py.File('rtipSc006.0039.restart_pert3.h5', 'r+')

        print f.keys()
        snap = f['snapshot']
        print snap.keys()

        salt = snap['salt'][:]
        temp = snap['temp'][:]

        pert_salt = np.random.normal(loc=1.0, scale=1.0, size=(94,44,40,3))
        pert_temp = np.random.normal(loc=1.0, scale=1.0, size=(94,44,40,3))

        salt[...] = salt - 0.02*pert_salt
        salt[salt<=10.]=0.
        del snap['salt']
        snap.create_dataset('salt', data=salt)

        temp[...] = temp - temp/(temp+0.01)*0.05*pert_temp
        del snap['temp']
        snap.create_dataset('temp', data=temp)

        print temp[30:40,30:40,-1,0]

        fig=pl.figure()
        pl.pcolormesh(temp[:,:,-1,0])
        pl.colorbar()

        pl.show()

        f.close()


def amoc_timeseries(ot):
        t_amoc = []; amoc = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i]/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max()/1000000.)
        return t_amoc, amoc



if __name__ == '__main__':
        MAIN()
