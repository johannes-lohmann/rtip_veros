import matplotlib.pyplot as pl
import numpy as np
import xarray as xr

N=2
for i in range(N):
        count = str(i); count = count.zfill(4)
        avgs = xr.open_dataset("ctrl40l.%s.averages.nc"%count)
        if i==0:
                forc = avgs['forc_salt_surface']
        else:
                forc=xr.concat((forc,avgs['forc_salt_surface']),dim='Time')

forc.to_netcdf('salt_forcing_monthly40l.nc')


