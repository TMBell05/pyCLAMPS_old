import os
from datetime import datetime, timedelta
from glob import glob

import netCDF4
import numpy as np

from pyclamps.io import concat_files

FILL_VALUE = -9999
start = datetime(2018, 7, 14)
end = datetime(2018, 7, 18)
day = start
td = timedelta(days=1)
elev = 70

while day <= end:
    filename = day.strftime('VAD_daily_%Y%m%d_{}.nc'.format(elev))
    if os.path.isfile(filename): os.remove(filename)

    print filename

    try:
        data = concat_files(glob(day.strftime('/data/lapserate/clamps/dl/vads/*%Y%m%d*.nc')))
    except IOError:
        day += td
        continue

    time = data['time']
    sort = np.argsort(time)
    time = time[sort]

    u = data['u'][sort]
    v = data['v'][sort]
    w = data['w'][sort]
    hgt = data['hgt'][0]
    r_sq = data['r_sq'][sort]
    rms = data['rms'][sort]


    # Create the netcdf
    nc = netCDF4.Dataset(filename, 'w', format="NETCDF4")

    # Create the height dimension
    nc.createDimension('height', len(hgt))
    nc.createDimension('time', len(time))

    # Add the attributes
    nc.setncattr("elev", elev)
    nc.setncattr("date", day.isoformat())

    # Create the variables
    u_var = nc.createVariable('u', 'f8', ('time', 'height',), fill_value=FILL_VALUE)
    v_var = nc.createVariable('v', 'f8', ('time', 'height',), fill_value=FILL_VALUE)
    w_var = nc.createVariable('w', 'f8', ('time', 'height',), fill_value=FILL_VALUE)
    rms_var = nc.createVariable('rms', 'f8', ('time', 'height',), fill_value=FILL_VALUE)
    r_sq_var = nc.createVariable('r_sq', 'f8', ('time', 'height',), fill_value=FILL_VALUE)

    time_var = nc.createVariable('time', 'i8', ('time',))
    time_var.setncattr('units', 'seconds since 1970-01-01 00:00:00 UTC')
    hgt_var = nc.createVariable('hgt', 'f8', ('height',), fill_value=FILL_VALUE)
    hgt_var.setncattr('units', 'Height above ground level')

    u_var[:] = np.where(np.isnan(u), FILL_VALUE, u)
    v_var[:] = np.where(np.isnan(v), FILL_VALUE, v)
    w_var[:] = np.where(np.isnan(w), FILL_VALUE, w)
    hgt_var[:] = np.where(np.isnan(hgt), FILL_VALUE, hgt)
    rms_var[:] = np.where(np.isnan(rms), FILL_VALUE, rms)
    r_sq_var[:] = np.where(np.isnan(r_sq), FILL_VALUE, r_sq)
    time_var[:] = time

    # Close the netcdf
    nc.close()

    day += td