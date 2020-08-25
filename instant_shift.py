import matplotlib.pyplot as pl
import numpy as np
import xarray as xr


def amoc_timeseries(ot):
        t_amoc = []; amoc = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i]/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max().values/1000000.)
        return t_amoc, amoc

ics = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9','10','11','12','13','14','15', '16','17', '18', '19', '20', '21', '22', '23', '24', '25', '26','30', '31', '32', '33', '34', '35', '40', '45']
M = 46

#fig=pl.figure()
#ax=pl.subplot(111)
#ax.set_color_cycle([pl.cm.plasma(i) for i in np.linspace(0, 1, len(ics))])
#for ic in ics:
for i in range(M):
        fig=pl.figure()
        ot = xr.open_dataarray("overturning_rtipInst%s.nc"%i)
        #st = xr.open_dataset("salt_temp_rtipz%s.nc"%i)
        #ot = xr.open_dataarray("overturning_rtipInst%s.nc"%ic)
        t_amoc, amoc = amoc_timeseries(ot)
        pl.plot(t_amoc, amoc)

pl.show()

tip = [1, 1, 1, 0.5, 1, 1, 1, 1, 0.5, 0, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 0, 0, 1, 0.5, 0.5, 0.5, 1, 0, 1, 1, 0.5, 1, 0.5, 0.5, 1, 0.5, 1, 0.5, 0, 0.5, 1, 0.5, 1, 0, 0.5, 0.5, 0.5, 1, 0]

#18 tip
#7 track
#21 edge...
