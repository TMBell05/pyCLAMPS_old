"""
Same as lidar_to_vad.py but modified to work with dave's ppi files

Calculates VADs based on a sin fit of the data. Can specify different keynames in 'var_lookup'
to account for different qc thresholds and variable names.

TODO-Make plotting more elegant than uncommenting stuff

Author: Tyler Bell (March 2017)
"""

import argparse
import os
from datetime import datetime

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from numpy import sin, cos

from pyclamps.utils import list_to_masked_array, ray_height
from pyclamps.vad import calc_vad_3d, calc_homogeneity


# Global Values
FILL_VALUE = -9999.
VEL_LIM = (-30, 30)
HGT_LIM = (0, 1000)
PROFILES_PER_PLOT = 2


var_lookup = {'leo':
                  {
                      'vel': 'radial_wind',
                      'thresh_var': 'cnr',
                      'thresh_value': -27.,
                      'az': 'azimuth',
                      'elev': 'elevaton',
                      'range': 'range',
                  },
              'clamps1':
                  {
                      'vel': 'velocity',
                      'thresh_var': 'intensity',
                      'thresh_value': 1.003,
                      'az': 'azimuth',
                      'elev': 'elevation',
                      'range': 'range',
                  },
              'mp3':
                  {
                      'vel': 'velocity',
                      'thresh_var': 'intensity',
                      'thresh_value': 1.005,
                      'az': 'azimuth',
                      'elev': 'elevation',
                      'range': 'range',
                  },
             }


def process_file(in_file, system, height=None, sinfit_dir=None):
    """
    Processes a line of sight netcdf file and outputs the vad. If height is specified,
    it will only return values nearest that height.
    :param in_file: File to process
    :param system: System to use in lookup table
    :param height: Height to process if desired
    :param sinfit_loc: Specify if you want to output a sin fit graph. A specific height MUST be specified.
    :return:
    """
    # Open the netcdf
    nc = netCDF4.Dataset(in_file)
    dates = nc['base_time'][0]+nc['time_offset'][:]

    scans = np.unique(nc['snum'][:])
    for scan in scans:
        print "Scan {} of {}".format(scan, np.max(scans))
        scan_ind = np.where(nc['snum'][:] == scan)[0]

        date = datetime.utcfromtimestamp(dates[scan_ind][0])

        hgt = []
        u = []
        v = []
        w = []
        rmse = []
        r_sq = []

        for i, range in enumerate(nc[var_lookup[system]['range']][:]*1000.):
            # Get the required stuff for this range ring
            cnr = nc[var_lookup[system]['thresh_var']][scan_ind, i].transpose()
            vel = nc[var_lookup[system]['vel']][scan_ind, i].transpose()
            az = nc[var_lookup[system]['az']][scan_ind]
            elev = nc[var_lookup[system]['elev']][scan_ind]

            nc_name = "{prefix}_{date}_{elev}.nc"
            nc_name = nc_name.format(prefix=args.out_prefix, date=date.strftime("%Y%m%d_%H%M%S"),
                                     elev=int(np.mean(elev)))
            nc_name = os.path.join(args.out_dir, nc_name)

            if os.path.isfile(nc_name):
                continue

            # Filter out the bad values based on CNR
            az = np.where(cnr <= var_lookup[system]['thresh_value'], FILL_VALUE, az)
            vel = np.where(cnr <= var_lookup[system]['thresh_value'], FILL_VALUE, vel)

            az = list_to_masked_array(az, FILL_VALUE)
            vel = list_to_masked_array(vel, FILL_VALUE)

            # Calculate the vad and height for this range ring
            tmp_u, tmp_v, tmp_w = calc_vad_3d(az, elev, vel)

            # Calculate the RMSE
            N = float(vel.size)
            az_rad = np.deg2rad(az)
            elev_rad = np.deg2rad(elev)

            derived_vr = (sin(az_rad) * cos(elev_rad) * tmp_u) + \
                         (cos(az_rad) * cos(elev_rad) * tmp_v) + \
                         (sin(elev_rad) * tmp_w)

            tmp_E = vel - derived_vr

            # Calculate rms error
            tmp_RMSE = np.sqrt(1 / N * np.sum(tmp_E ** 2))

            tmp_r_sq = calc_homogeneity(vel, derived_vr)

            # Append to the lists for plotting
            u.append(tmp_u)
            v.append(tmp_v)
            w.append(tmp_w)
            hgt.append(ray_height(range, elev[0]))
            rmse.append(tmp_RMSE)
            r_sq.append(tmp_r_sq)

        if os.path.isfile(nc_name):
            continue
        if height is not None:
            i = (np.abs(np.asarray(hgt) - height)).argmin()
            print i
            u = u[i]
            v = v[i]
            w = w[i]
            hgt = hgt[i]
            rmse = rmse[i]

            if sinfit_dir is not None:
                filename = "sinfit_{height}m_{date}.png".format(height=height, date=date.strftime("%Y%m%d_%H%M%S"))
                filename = os.path.join(sinfit_dir, filename)
                print nc
                cnr = nc[var_lookup[system]['thresh_var']][scan_ind, i].transpose()
                vel = nc[var_lookup[system]['vel']][scan_ind, i].transpose()
                az = nc[var_lookup[system]['az']][scan_ind].transpose()
                elev = np.round(np.mean(nc[var_lookup[system]['elev']][scan_ind]), 2)

                vel = np.where(cnr <= var_lookup[system]['thresh_value'], np.nan, vel)

                az_rad = np.deg2rad(az)
                elev_rad = np.deg2rad(elev)

                print az.shape
                print vel.shape

                plt.figure()
                plt.scatter(az, vel)
                plt.plot(az, (
                u * sin(az_rad) * cos(elev_rad) + v * cos(az_rad) * cos(elev_rad) + w * sin(elev_rad)).tolist())
                plt.ylim((-20, 20))
                plt.xlim((0, 360))
                plt.title("{date} Elev={elev} Height={height}".format(date=date.strftime("%Y%m%d_%H%M%S"),
                                                                      elev=elev,
                                                                      height=hgt))
                plt.savefig(filename)
                plt.close()

        write_to_nc(nc_name, date, elev, u, v, w, hgt, rmse, r_sq)

    # Close the netcdf
    nc.close()

    # return u, v, w, hgt, rmse, r_sq, date, elev


def write_to_nc(filename, date, elev, u, v, w, hgt, rmse, r_sq):
    # Create the netcdf
    nc = netCDF4.Dataset(filename, 'w', format="NETCDF4")

    # Create the height dimension
    nc.createDimension('height', len(hgt))

    # Add the attributes
    nc.setncattr("elev", elev)
    nc.setncattr("date", date.isoformat())

    # Create the variables
    u_var = nc.createVariable('u', 'f8', ('height',), fill_value=FILL_VALUE)
    v_var = nc.createVariable('v', 'f8', ('height',), fill_value=FILL_VALUE)
    w_var = nc.createVariable('w', 'f8', ('height',), fill_value=FILL_VALUE)
    hgt_var = nc.createVariable('hgt', 'f8', ('height',), fill_value=FILL_VALUE)
    rms_var = nc.createVariable('rms', 'f8', ('height',), fill_value=FILL_VALUE)
    r_sq_var = nc.createVariable('r_sq', 'f8', ('height',), fill_value=FILL_VALUE)

    time_var = nc.createVariable('time', 'i8')
    time_var.setncattr('units', 'seconds since 1970-01-01 00:00:00 UTC')

    u_var[:] = np.where(np.isnan(u), FILL_VALUE, u)
    v_var[:] = np.where(np.isnan(v), FILL_VALUE, v)
    w_var[:] = np.where(np.isnan(w), FILL_VALUE, w)
    hgt_var[:] = np.where(np.isnan(hgt), FILL_VALUE, hgt)
    rms_var[:] = np.where(np.isnan(rmse), FILL_VALUE, rmse)
    r_sq_var[:] = np.where(np.isnan(r_sq), FILL_VALUE, r_sq)
    time_var[0] = (date - datetime(1970, 1, 1)).total_seconds()

    # Close the netcdf
    nc.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='in_files', nargs='*')
    parser.add_argument('-O', dest='out_prefix', default='vad')
    parser.add_argument('-o', dest='out_dir', default=os.getcwd())
    parser.add_argument('-s', dest='system', default=None)

    args = parser.parse_args()

    vad_ws = np.zeros(len(args.in_files))
    vad_wd = np.zeros(len(args.in_files))
    time = np.zeros(len(args.in_files), dtype=datetime)
    vad_rmse = np.zeros(len(args.in_files))

    height = 100

    for i, f in enumerate(sorted(args.in_files)):
        print(f)
        process_file(f, system=args.system)#, height=height, sinfit_dir=args.out_dir)
        # nc_name = "{prefix}_{date}_{elev}.nc"
        # nc_name = nc_name.format(prefix=args.out_prefix, date=date.strftime("%Y%m%d_%H%M%S"), elev=int(np.mean(elev)))
        # nc_name = os.path.join(args.out_dir, nc_name)

        # write_to_nc(nc_name, date, elev, u, v, w, hgt, rmse, r_sq)


        #
        # ws = np.sqrt(np.asarray(u)**2 + np.asarray(v)**2)
        # wd = np.arctan2(u, v)
        #
        # vad_ws[i] = ws
        # vad_wd[i] = np.rad2deg(wd)
        # time[i] = date
        # vad_rmse[i] = rmse



    # fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    # ind = np.logical_and(~np.isnan(vad_wd), time != 0)
    # vad_wd = np.where(vad_wd[ind] < 0, np.asarray(vad_wd[ind]) + 360., vad_wd[ind])
    # fig.suptitle("Time series %s" % time[1].strftime("%Y%m%d"))
    #
    # ax1.scatter(time[ind], vad_ws[ind], c=vad_rmse[ind])
    # ax1.set_ylim(0, 30)
    # ax1.set_xlim(time[0], time[-1])
    # ax1.set_title("Wind Speed ($m s^-1$)")
    #
    # ax2.scatter(time[ind], vad_wd[ind], c=vad_rmse[ind])
    # ax2.set_ylim(0, 360)
    # ax2.set_title("Wind Direciton")
    # plt.savefig("{date}_{height}m_timeseries.png".format(date=time[1].strftime("%Y%m%d"), height=height))
    # plt.close()






