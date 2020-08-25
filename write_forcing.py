import matplotlib.pyplot as pl
import numpy as np
import xarray as xr

forc = xr.open_dataset("forcing_1deg_global_interpolated.nc")

mask = forc.isel(Time=0)['sss']

mask = mask.where(mask<0., 0)

lc = mask.coords['xt']; la = mask.coords['yt']
mask.loc[dict(xt=lc[(lc>296)&(lc<362)], yt=la[(la>49)&(la<75)])] = 1
print mask

forc['freshwater_mask'] = mask

print forc

forc['freshwater_mask'].plot()

forc.to_netcdf('forcing_1deg_global_interp_hosing.nc')

pl.show()
