import matplotlib.pyplot as pl
import numpy as np
#import xarray as xr
#import h5netcdf
import h5py

f = h5py.File('ctrl40l.0039.restart_pert3.h5', 'r+')

print f.keys()

snap = f['snapshot']

print snap.keys()

salt = snap['salt'][:]
temp = snap['temp'][:]

pert_salt = np.random.normal(loc=0.0, scale=1.0, size=(94,44,40,3))
pert_temp = np.random.normal(loc=0.0, scale=1.0, size=(94,44,40,3))

salt[...] = salt + 0.05*pert_salt#0.1
salt[salt<=10.]=0.
del snap['salt']
snap.create_dataset('salt', data=salt)

temp[...] = temp + temp/(temp+0.01)*0.1*pert_temp#0.2
del snap['temp']
snap.create_dataset('temp', data=temp)

print temp[30:40,30:40,-1,0]

fig=pl.figure()
pl.pcolormesh(temp[:,:,-1,0])
pl.colorbar()

pl.show()

f.close()

