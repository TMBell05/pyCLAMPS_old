import matplotlib
matplotlib.use('agg')

import argparse
import logging
import os

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from datetime import datetime, timedelta
from glob import glob

from pyclamps.io import concat_files, get_aeri_variables
from pyclamps.plotting import time_height
from pyclamps.calc import *

# CONSTANTS
dl_thresh = 1.01
avg_time = timedelta(minutes=5)
zmax = 1000.
fill_value = -9999

# Get the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d', required=True, action='store', dest='date', help="Date (YYYYmmdd)")
parser.add_argument('-l', required=True, action='store', dest='stare', help="Stare file dir")
parser.add_argument('-a', required=True, action='store', dest='aeri', help="AERI file dir")
parser.add_argument('-p', required=True, action='store', dest='vad', help="VAD Directory")
parser.add_argument('-o', required=True, action='store', dest='out_dir', help="Output Directory")
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true')
args = parser.parse_args()

# Set up logger
if args.verbose:
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.DEBUG)
else:
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.INFO)

# Get the date from the arguments
date = datetime.strptime(args.date, "%Y%m%d")
img_name = os.path.join(args.out_dir, date.strftime("CLAMPS_%Y%m%d.png"))

# Go ahead and get the axes to plot with
fig, (ptemp, dpt, wnd_spd, wnd_dir, vert_vel_ax, vert_std_ax) = plt.subplots(6)
fig.set_figheight(25)
fig.set_figwidth(15)

'''Stare processing'''
# Get the lidar stare data
logging.info("Getting stare data...")
stare_data = concat_files(os.path.join(args.stare, date.strftime("*%Y%m%d*.cdf")))
stare_times = np.asarray([datetime.utcfromtimestamp(d) for d in (stare_data['base_time']+stare_data['time_offset'])])

# Average the stare date
logging.info("Averaging stare data...")
stare_data['velocity'].mask = stare_data['intensity'] < dl_thresh
vert_vel_avg = []
vert_vel_std = []
time_avg = []
last = 0
for i, time in enumerate(stare_times):
    if time - stare_times[i-1] > timedelta(seconds=60):
        avg_data = stare_data['velocity'][last:i-1]

        vert_vel_avg.append(np.mean(avg_data, axis=0).filled(fill_value))
        vert_vel_std.append(np.std(avg_data, axis=0).filled(fill_value))
        time_avg.append(stare_times[(last + i)/2])

        last = i


time_avg = np.asarray(time_avg)
vert_vel_avg = np.ma.masked_where(np.asarray(vert_vel_avg) == fill_value, vert_vel_avg)
vert_vel_std = np.ma.masked_where(np.asarray(vert_vel_std) == fill_value, vert_vel_std)

# Plot the data
logging.info("Plotting stare data...")
vert_vel_ax = time_height(time_avg, stare_data['height'][2:]*1.e3, vert_vel_avg[:, 2:].transpose(), 'w',
                          ax=vert_vel_ax, datemin=date, datemax=date+timedelta(days=1), zmax=zmax, zmin=0,
                          datamin=-2, datamax=2)
vert_std_ax = time_height(time_avg, stare_data['height'][2:]*1.e3, vert_vel_std[:, 2:].transpose(), 'std',
                          ax=vert_std_ax, datemin=date, datemax=date+timedelta(days=1), zmax=zmax, zmin=0,
                          datamin=0, datamax=2)

vert_vel_ax.set_title(date.strftime("Vertical Velocity (5min avg) -- %Y%m%d"))
vert_std_ax.set_title(date.strftime("Sigma w -- %Y%m%d"))

del stare_data, vert_vel_std, vert_vel_avg, time_avg, stare_times


'''AERI plotting'''
# TODO - Make work with multiple files. Only works with one currently
logging.info("Getting AERI data")
aeri_data = get_aeri_variables(glob(os.path.join(args.aeri, date.strftime("*%Y%m%d*.cdf")))[0],
                               ['temperature', 'pressure', 'waterVapor'])

theta = potential_temp(aeri_data['temperature'], aeri_data['pressure'])

e = vapor_pressure(aeri_data['pressure'], aeri_data['waterVapor']*1e-3)
es = sat_vapor_pressure(aeri_data['temperature'])
dewpt = dewpoint(e)


logging.info("Plotting AERI Data")
ptemp = time_height(aeri_data['time'], aeri_data['height']*1e3, theta.transpose(), 'pt',
                    ax=ptemp, datemin=date, datemax=date+timedelta(days=1), zmax=zmax, zmin=0,
                    datamin=5, datamax=30)

dpt = time_height(aeri_data['time'], aeri_data['height']*1e3, (100. * e/es).transpose(), 'rh',
                  ax=dpt, datemin=date, datemax=date+timedelta(days=1), zmax=zmax, zmin=0,
                  datamin=0, datamax=100)

ptemp.set_title(date.strftime("Potential Temperature -- %Y%m%d"))
dpt.set_title(date.strftime("Relative Humidity -- %Y%m%d"))
del aeri_data

'''Vad plotting'''
logging.info("Loading VADs")
vad_data = concat_files(sorted(glob(os.path.join(args.vad, date.strftime('*%Y%m%d*_70*')))))
wspd = wind_spd(vad_data['u'], vad_data['v'])
wdir = wind_dir(vad_data['u'], vad_data['v'])
wdir = np.ma.masked_where(np.isnan(wdir), wdir)

time = np.asarray([datetime.utcfromtimestamp(d) for d in (vad_data['time'])])
# time = np.meshgrid(vad_data['hgt'][0], mdates.date2num(time))[1]

logging.info("Plotting VADs")
wnd_spd = time_height(time, vad_data['hgt'][0][2:], wspd[:, 2:].transpose(), 'ws',
                      ax=wnd_spd, datemin=date, datemax=date+timedelta(days=1), zmax=zmax, zmin=0,
                      datamin=0, datamax=20)
wnd_dir = time_height(time, vad_data['hgt'][0][2:], wdir[:, 2:].transpose(), 'wd',
                      ax=wnd_dir, datemin=date, datemax=date+timedelta(days=1), zmax=zmax, zmin=0,
                      datamin=0, datamax=360)

wnd_spd.set_title(date.strftime("Wind Speed -- %Y%m%d"))
wnd_dir.set_title(date.strftime("Wind Direction -- %Y%m%d"))


plt.tight_layout()
plt.savefig(img_name)


