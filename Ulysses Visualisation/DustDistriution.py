#interesting lon and lat for start_index = 2000
#usually lon = 250 and lat widely distributed
#velocity centred on 30

'''
#https://naif.jpl.nasa.gov/pub/naif/

# Leap Second Kernel naif0012.tls in generic_kernels/lsk
spice.furnsh("spice/generic_kernels/naif0012.tls")

# de440.bsp enthält die Bahnen der großen Planeten
# zu finden in generic_kernels/spk/planets
spice.furnsh("spice/generic_kernels/de440.bsp")

# Ulysses-Bahndaten, in Ulysses/kernels/spk
spice.furnsh("spice/ulysses_1990_2009_2050.bsp")
'''

import numpy as np
import matplotlib.pyplot as plt
from indices_dict import indices
import spiceypy as spice

spice_path = '../../spice/'
spice.furnsh(spice_path + "naif0012.tls")
one_year_et = spice.str2et('01-001T00:00')-spice.str2et('00-001T00:00')

first_data_column = 3
last_data_column = 34

min_quality_flag = 2
all_lines = False
bins = 10
PlottedQuantity = 'solar_lon'

used_cols = (i for i in range(first_data_column, last_data_column))

index = []
time = []
with open('Ulysses_Data_File_Cleaned.txt') as cleaned_ulysses_data:
    for count, line in enumerate(cleaned_ulysses_data):
        if count >= indices['first_data_line']:
            line = line.split()
            index.append(line[indices['index']])
            time.append(line[indices['date']]+'T'+line[indices['time_of_day']])

index = np.array(index)

time = spice.str2et(time)

ulysses_data = np.loadtxt('Ulysses_Data_File_Cleaned.txt', delimiter = ' ', skiprows = indices['first_data_line'], usecols = used_cols)

ulysses_data = np.concatenate((np.transpose(np.array([time])), ulysses_data), axis = 1)

if all_lines:
    start_index = 0
    end_index = len(ulysses_data)
else:
    start_et = spice.str2et('90-001T00:00')
    end_et = start_et+one_year_et
    start_index = np.argmin(np.abs(time-start_et))
    end_index = np.argmin(np.abs(time-end_et))


class CleanedDataPlotter:

    def execute(self):
        self._remove_999_values()
        self._plot_hist()

    def _remove_999_values(self):
        ulysses_data_without_999 = []
        for i in range(self.start_index, self.end_index):
            condition = self.data[i, self.rotation_angle_index] != 999 and self.data[i, self.quality_flag_index] >= self.min_quality_flag
            if condition:
                ulysses_data_without_999.append(self.data[i])
        self.data_without_999 = np.array(ulysses_data_without_999)

    def _plot_hist(self):
        plt.xlabel(self.xlabel)
        plt.ylabel('Count')
        plt.title(self.current_year)
        plt.hist(self.data_without_999[:, self.plot_index], bins = self.bins)
        plt.show()

    def __init__(self):
        self.data = None
        self.start_index = None
        self.end_index = None
        self.min_quality_flag = None
        self.rotation_angle_index = None
        self.quality_flag_index = None
        self.data_without_999 = None
        self.xlabel = None
        self.plot_index = None
        self.bins = None
        self.current_year = None


LinesPlotter = CleanedDataPlotter()
LinesPlotter.data = ulysses_data
LinesPlotter.start_index = start_index
LinesPlotter.end_index = end_index
LinesPlotter.min_quality_flag = min_quality_flag
LinesPlotter.rotation_angle_index = indices['rotation_angle_index']
LinesPlotter.quality_flag_index = indices['quality_flag_index']
LinesPlotter.xlabel = PlottedQuantity
LinesPlotter.plot_index = indices[PlottedQuantity+'_index']
LinesPlotter.bins = bins

i = 1990
while True:
    LinesPlotter.current_year = i
    LinesPlotter.execute()
    if LinesPlotter.end_index == len(ulysses_data)-1:
        break
    start_et = end_et
    end_et = end_et+one_year_et
    LinesPlotter.start_index = np.argmin(np.abs(time-start_et))
    LinesPlotter.end_index = np.argmin(np.abs(time-end_et))
    i += 1
