# %%

import numpy as np
from indices_dict import indices
import spiceypy as spice
from CleanedDataPlotter import CleanedDataPlotter

max_area = 1000

spice_path = '../../spice/'
spice.furnsh(spice_path + "naif0012.tls")
spice.furnsh(spice_path + "de440.bsp")
spice.furnsh(spice_path + "ulysses_1990_2009_2050.bsp")
one_year_et = spice.str2et('01-001T00:00')-spice.str2et('00-001T00:00')

first_data_column = 3
last_data_column = 34

min_quality_flag = 2
all_lines = True

wehry_comparison = True

used_cols = (i for i in range(first_data_column, last_data_column))

index = []
time = []
with open('Ulysses_Data_File_Cleaned.txt') as cleaned_ulysses_data:
    for count, line in enumerate(cleaned_ulysses_data):
        if count >= indices['first_data_line']:
            line = line.split()
            index.append(line[indices['index']])
            time.append(line[indices['date']]+'T'+line[indices['time_of_day']])

for i in range(len(index)):
    index[i] = int(index[i])

time = spice.str2et(time)
dist_array = []

for et in time:
    [pos, ltime] = spice.spkpos('SUN',  et,      'J2000', 'NONE', 'ULYSSES')
    dist = spice.vnorm(pos)
    dist = spice.convrt(dist, 'KM', 'AU')
    dist_array.append(dist)

ulysses_data = np.loadtxt('Ulysses_Data_File_Cleaned.txt', delimiter = ' ', skiprows = indices['first_data_line'], usecols = used_cols)

ulysses_data = np.concatenate((np.transpose(np.array([index])), ulysses_data), axis = 1)

ulysses_data = np.concatenate((np.transpose(np.array([dist_array])), ulysses_data), axis = 1)

ulysses_data = np.concatenate((np.transpose(np.array([time])), ulysses_data), axis = 1)


wehry_beta_particle_indices = [4,6,8,18,19,29,31,34,35,36,38,43,48,49,59,61,66,72,74,76,80,83,86,94,1032,1080,1145,1165,1410,1412,1421,1422,1427,1428,1429,1431,1433,1436,1438,1440,1442,1449,1450,1452,1455,1465,1984,2001,2003,2010,2012,2024,2034,2035,2048,2049,2051,2052,2053,2054,2055,2060]

wehry_beta_particle_indices = np.array(wehry_beta_particle_indices)-1

wehry_data = ulysses_data[wehry_beta_particle_indices,:]

LinesPlotter = CleanedDataPlotter()

if all_lines:
    start_index = 0
    end_index = len(ulysses_data)
else:
    start_et = spice.str2et('90-001T00:00')
    end_et = start_et+one_year_et
    start_index = np.argmin(np.abs(time-start_et))
    end_index = np.argmin(np.abs(time-end_et))

LinesPlotter.data = ulysses_data
if wehry_comparison:
    LinesPlotter.data = wehry_data
LinesPlotter.raw_data = ulysses_data
LinesPlotter.start_index = start_index
LinesPlotter.end_index = end_index
LinesPlotter.min_quality_flag = min_quality_flag
LinesPlotter.rotation_angle_index = indices['rotation_angle_index']
LinesPlotter.quality_flag_index = indices['quality_flag_index']
LinesPlotter.one_year_et = one_year_et
LinesPlotter.min_eff_area_factor = 50
LinesPlotter.eff_area_file = 'DefaultDataset.csv'
#LinesPlotter.eff_area_file = '20.dat'

i = 1990
while True:
    if not all_lines:
        LinesPlotter.current_year = i
    if LinesPlotter.end_index == len(ulysses_data)-1:
        break
    if all_lines:
        if wehry_comparison:
            LinesPlotter.execute_wehry_special()
        else:
            LinesPlotter.execute_wehry()
    else:
        try:
            LinesPlotter.execute_wehry()
        except:
            print('No beta meteoroids found in', i)
    if not all_lines:
        start_et = end_et
        end_et = end_et+one_year_et
        LinesPlotter.start_index = np.argmin(np.abs(time-start_et))
        LinesPlotter.end_index = np.argmin(np.abs(time-end_et))
        i += 1
    if all_lines:
        break
