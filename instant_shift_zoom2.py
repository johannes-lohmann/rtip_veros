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

def amoc_timeseries(ot):
        t_amoc = []; amoc = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i]/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max().values/1000000.)
        return t_amoc, amoc

ot_init = xr.open_dataarray("overturning_rtipzSpInst2.nc")
t_amoc_init, amoc_init = amoc_timeseries(ot_init)

### 10 year branch of spinup; 2-month averages
ot_init0 = xr.open_dataarray("overturning_rtipzSpInst.nc")
t_amoc_init0, amoc_init0 = amoc_timeseries(ot_init0)

### 100 year spinup; 2-year averages
ot_init1 = xr.open_dataarray("overturning_rtipSc0i.nc")
t_amoc_init1, amoc_init1 = amoc_timeseries(ot_init1)

t_amoc_init1 = np.concatenate(([9100.], t_amoc_init1[:46]))
amoc_init1 = np.concatenate(([7.05], amoc_init1[:46]))

#t_amoc_init0 = np.concatenate(([9100.], t_amoc_init0[:46]))
#amoc_init0 = np.concatenate(([7.05], amoc_init0[:46]))

'''
ics = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9','10','11','12','13','14','15', '16','17', '18', '19', '20', '21', '22', '23', '24', '25', '26','30', '31', '32', '33', '34', '35', '40', '45']

ics= ['0','1','2','3','4','5alt','6','7', '8']
ics= ['5','5alt']
'''

amoc_init1_step = [];t_amoc_init1_step = []
for i in range(len(amoc_init1)):
        amoc_init1_step.append(amoc_init1[i])
        amoc_init1_step.append(amoc_init1[i])

        t_amoc_init1_step.append(t_amoc_init1[i]-2.)
        t_amoc_init1_step.append(t_amoc_init1[i])

amoc_init0_step = [];t_amoc_init0_step = []
for i in range(len(amoc_init0)):
        amoc_init0_step.append(amoc_init0[i])
        amoc_init0_step.append(amoc_init0[i])

        t_amoc_init0_step.append(t_amoc_init0[i]-1./6.)
        t_amoc_init0_step.append(t_amoc_init0[i])

amoc_init_step = [];t_amoc_init_step = []
for i in range(len(amoc_init)):
        amoc_init_step.append(amoc_init[i])
        amoc_init_step.append(amoc_init[i])

        t_amoc_init_step.append(t_amoc_init[i]-1./52.)
        t_amoc_init_step.append(t_amoc_init[i])



fig=pl.figure(figsize=(7.5,6.5))
pl.subplots_adjust(left=0.13, bottom=0.15, right=0.97, top=0.95, wspace=0.0, hspace=0.28)
pl.subplot(311)
#pl.plot(t_amoc_init, amoc_init, 'x-', color='black')
pl.plot(t_amoc_init_step, amoc_init_step, color='green')
pl.ylabel('AMOC max. (Sv)')

pl.subplot(312)
#pl.plot(t_amoc_init0, amoc_init0, 'x-', color='black')
pl.plot(t_amoc_init0_step, amoc_init0_step, color='crimson')
pl.plot(t_amoc_init_step, amoc_init_step, color='green', alpha=0.7, linewidth=1.)
#pl.plot(t_amoc_init, amoc_init)
pl.ylabel('AMOC max. (Sv)')

pl.subplot(313)
pl.plot(t_amoc_init1_step, amoc_init1_step, color='black')
#pl.plot(t_amoc_init0, amoc_init0)
pl.plot(t_amoc_init0_step, amoc_init0_step, color='crimson', alpha=0.5, linewidth=1.)

pl.ylabel('AMOC max. (Sv)'); pl.xlabel('Branch-off time (years)')


M = 52

#fig=pl.figure()
#ax=pl.subplot(111)
#ax.set_color_cycle([pl.cm.plasma(i) for i in np.linspace(0, 1, M)])
#for ic in ics:
for i in range(37,M):
        fig=pl.figure()
        ot = xr.open_dataarray("overturning_rtipInstZ%s.nc"%i)
        #st = xr.open_dataset("salt_temp_rtipz%s.nc"%i)
        #ot = xr.open_dataarray("overturning_rtipInst%s.nc"%ic)
        t_amoc, amoc = amoc_timeseries(ot)
        pl.plot(t_amoc, amoc)

pl.show()


#tip = [1, 1, 1, 0.5, 1, 1, 1, 1, 0.5, 0, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 0, 0, 1, 0.5, 0.5, 0.5, 1, 0, 1, 1, 0.5, 1, 0.5, 0.5, 1, 0.5, 1, 0.5, 0, 0.5, 1, 0.5, 1, 0, 0.5, 0.5, 0.5, 1, 0]

#18 tip
#7 track
#21 edge...
