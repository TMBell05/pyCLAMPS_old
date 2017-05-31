import argparse
from datetime import datetime, timedelta
import os

import matplotlib
matplotlib.use('agg')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cmocean

from pyclamps.io import read_wind_profile

# Get the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='*', action='store', dest='files')
parser.add_argument('-o', action='store', dest='out_file')
args = parser.parse_args()

# Preallocate arrays
num_files = len(args.files)
dates = np.zeros(num_files)
hgt = []
dir = []
spd = []

# Get the data
for i, f in enumerate(sorted(args.files)):
    print i, f
    # Get the date info
    date = datetime.strptime(f.split('/')[-1], 'Processed_Wind_Profile_37_%Y%m%d_%H%M%S.hpl')
    dates[i] = mdates.date2num(date)

    # Get the data from the file
    z, d, s = read_wind_profile(f)

    hgt.append(z)
    dir.append(d)
    spd.append(s)

    # Do some housekeeping to get the min day and max day
    if i == 0:
        min_date = datetime(date.year, date.month, date.day, 0)
        max_date = min_date + timedelta(days=1)

# Convert the data to arrays
hgt = np.asarray(hgt)
dir = np.asarray(dir)
spd = np.asarray(spd)

print hgt.shape

# Create the plot
hours = mdates.HourLocator()
fmt = mdates.DateFormatter('%H%M')
xlim = mdates.date2num(np.asarray([min_date, max_date]))
print min_date, max_date

_, time = np.meshgrid(hgt[0], dates)

# Create the figure
fig, (ax1, ax2) = plt.subplots(2)
fig.set_figheight(10)
fig.set_figwidth(15)

# Create the wind speed plot
c1 = ax1.pcolormesh(time, hgt, spd, vmin=0, vmax=20, cmap='jet')
ax1.set_ylim(0, 2000)
ax1.xaxis.set_major_locator(hours)
ax1.xaxis.set_major_formatter(fmt)
ax1.set_xlim(xlim)
ax1.set_title(date.strftime('Wind Speed (m/s) -- %Y%m%d'))

# Create the wind direction plot
c2 = ax2.pcolormesh(time, hgt, dir, vmin=0, vmax=360, cmap=cmocean.cm.phase)
ax2.set_ylim(0, 2000)
ax2.xaxis.set_major_locator(hours)
ax2.xaxis.set_major_formatter(fmt)
ax2.set_xlim(xlim)
ax2.set_title(date.strftime('Wind Direction (deg) -- %Y%m%d'))

# Add colorbars and save the picture
plt.colorbar(c1, ax=ax1)
plt.colorbar(c2, ax=ax2)
plt.tight_layout()
plt.savefig(args.out_file)
plt.close()


