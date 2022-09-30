import datetime

import cdflib
import numpy as np
from matplotlib import dates
from sunpy.net import Fido
from sunpy.net import attrs as a

# # define start and end date
# start_time="2012-5-26 10:30"
# end_time="2012-5-28 15:40"
# # specify spacecraft 'ahead'/'behind'
# spacecraft = 'ahead'


def get_swaves(start_time, end_time, path=None):

    #######################
    # downloading the files
    #######################

    dataset = 'STEREO_LEVEL2_SWAVES'
    cda_dataset = a.cdaweb.Dataset(dataset)

    trange = a.Time(start_time, end_time)

    # always add 1 day to enddate because the enddate itself should be included in the data (which isn't the case)
    trange = a.Time(start_time, trange.end.to_datetime().date()+datetime.timedelta(days=1))

    result = Fido.search(trange, cda_dataset)
    downloaded_files = Fido.fetch(result, path=path)  # use Fido.fetch(result, path='/ThisIs/MyPath/to/Data/{file}')  to use a specific local folder for saving data files
    downloaded_files.sort()
    # print(downloaded_files)

    return downloaded_files


def plot_swaves(downloaded_files, spacecraft, start_time, end_time, ax, cmap='inferno'):

    ###################
    # reading the files
    ###################
    data_all = []
    time_all = []

    for i in downloaded_files:
        cdf_file = cdflib.CDF(i)

        data = cdf_file.varget("avg_intens_" + spacecraft)
        data_all.append(data)

        freq = cdf_file.varget('frequency')/1000  # in MHz

        time = cdf_file.varget('Epoch')
        time_all.append(cdflib.epochs.CDFepoch.to_datetime(time))

    # full time array for plotting
    time_arr = np.array(time_all).flatten()

    # full data array for plotting
    data_all = np.array(data_all)
    # if there are more than one 1-day file downloaded
    if data_all.shape[0] > 1:
        data_arr = np.concatenate((data_all[0], data_all[1]))
        for i in range(1, data_all.shape[0] - 1):
            data_arr = np.concatenate((data_arr, data_all[i+1]))
        # switching frequency axis
        data_arr = data_arr.T
    else:
        # switching frequency axis
        # one must choose the first entry of data_all here, because it's a list with len==1
        data_arr = data_all[0].T

    if isinstance(start_time, str):
        start = datetime.datetime.strptime(start_time, '%Y-%m-%d %H:%M')
        end = datetime.datetime.strptime(end_time, '%Y-%m-%d %H:%M')
    else:
        start, end = start_time, end_time

    ######################
    # plotting the spectra
    ######################

    colormesh = ax.pcolormesh(time_arr, freq[::-1], data_arr[::-1], vmin=0, vmax=0.5*np.max(data_arr), cmap=cmap)

    ax.set_ylabel('Frequency (MHz)')
    ax.set_xlabel('Date and time (UT)')
    ax.set_yscale('log')
    ax.set_ylim(freq[-1], freq[0])
    ax.set_yticks([0.01, 0.1, 1, 10])
    ax.set_yticklabels(['0.01', '0.1', '1', '10'])
    ax.set_xlim(start, end)

    ax.xaxis_date()
    ax.xaxis.set_major_formatter(dates.DateFormatter('%d/%m %H:%M'))
    # plt.show()

    return ax, colormesh
