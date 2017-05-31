import matplotlib
matplotlib.use('agg')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cmocean

from datetime import datetime, timedelta
from pyclamps.io import concat_files

start_date = datetime(2017, 5, 6)
end_date = datetime(2017, 5, 7)
date = start_date
td = timedelta(days=1)

while date <= end_date:
    vad_data = concat_files(date.strftime('/data/tbell/clamps1/dl/vads/*%Y%m%d*'))
    print vad_data

    hours = mdates.HourLocator()
    fmt = mdates.DateFormatter('%H%M')
    xlim = mdates.date2num(np.asarray([date, date+td]))

    time = np.asarray([datetime.utcfromtimestamp(d) for d in (vad_data['time'])])
    rng, time = np.meshgrid(vad_data['hgt'][0], mdates.date2num(time))

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    fig.set_figheight(15)
    fig.set_figwidth(15)
    spd = np.sqrt(vad_data['u'] ** 2 + vad_data['v'] ** 2)
    # spd = np.ma.masked_where(np.isnan(spd), spd)
    c1 = ax1.pcolormesh(time, rng, spd, vmin=0, vmax=20, cmap='jet')
    ax1.set_ylim(0, 2000)
    ax1.xaxis.set_major_locator(hours)
    ax1.xaxis.set_major_formatter(fmt)
    ax1.set_xlim(xlim)
    ax1.set_title(date.strftime('Wind Speed (m/s) -- %Y%m%d'))

    direction = np.rad2deg(np.arctan2(vad_data['v'], vad_data['u']))
    direction = -direction + 90
    direction = np.where(direction < 0, direction + 360., direction)
    direction = np.ma.masked_where(np.isnan(direction), direction)

    ax2.set_ylim(0, 2000)
    c2 = ax2.pcolormesh(time, rng, direction, vmin=0, vmax=360, cmap=cmocean.cm.phase)
    ax2.set_title(date.strftime('Wind Direction (deg) -- %Y%m%d'))

    rms = np.ma.masked_where(np.isnan(vad_data['w']), vad_data['w'])
    ax3.set_ylim(0, 2000)
    c3 = ax3.pcolormesh(time, rng, rms, vmin=-1, vmax=1, cmap=cmocean.cm.delta)
    ax2.set_title(date.strftime('Wind Direction (deg) -- %Y%m%d'))

    plt.colorbar(c1, ax=ax1)
    plt.colorbar(c2, ax=ax2)
    plt.colorbar(c3, ax=ax3)
    plt.tight_layout()
    plt.savefig(date.strftime("VAD_%Y%m%d.png"))
    plt.close()

    date += td