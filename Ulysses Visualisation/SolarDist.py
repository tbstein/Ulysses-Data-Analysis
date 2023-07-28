import numpy as np
from indices_dict import indices
import spiceypy as spice
from CleanedDataPlotter import CleanedDataPlotter

spice_path = '../../spice/'
spice.furnsh(spice_path + "naif0012.tls")
spice.furnsh(spice_path + "de440.bsp")
spice.furnsh(spice_path + "ulysses_1990_2009_2050.bsp")
one_year_et = spice.str2et('01-001T00:00')-spice.str2et('00-001T00:00')

first_data_column = 3
last_data_column = 34

min_quality_flag = 0

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
    [pos, ltime] = spice.spkpos('SUN',  et,      'ECLIPJ2000', 'NONE', 'ULYSSES')
    dist = spice.vnorm(pos)
    dist = spice.convrt(dist, 'KM', 'AU')
    dist_array.append(dist)

ulysses_data = np.loadtxt('Ulysses_Data_File_Cleaned.txt', delimiter = ' ', skiprows = indices['first_data_line'], usecols = used_cols)

ulysses_data = np.concatenate((np.transpose(np.array([index])), ulysses_data), axis = 1)

ulysses_data = np.concatenate((np.transpose(np.array([dist_array])), ulysses_data), axis = 1)

ulysses_data = np.concatenate((np.transpose(np.array([time])), ulysses_data), axis = 1)


wehry_beta_particle_indices = [4,6,8,18,19,29,31,34,35,36,38,43,48,49,59,61,66,72,74,76,80,83,86,94,1032,1080,1145,1165,1410,1412,1421,1422,1427,1428,1429,1431,1433,1436,1438,1440,1442,1449,1450,1452,1455,1465,1984,2001,2003,2010,2012,2024,2034,2035,2048,2049,2051,2052,2053,2054,2055,2060]

peter = [] 
with open('Ulysses_Data_File_Cleaned.txt') as cleaned_ulysses_data:
    for count, line in enumerate(cleaned_ulysses_data):
        if count >= indices['first_data_line']:
            line0 = line.split()
            for k in wehry_beta_particle_indices:
                if int(line0[indices['index']]) == k:
                    peter.append(line0)

wehry_beta_particle_indices = np.array(wehry_beta_particle_indices)-1

wehry_data = ulysses_data[wehry_beta_particle_indices,:]

with open('wehryData.dat', 'w') as f:
        for i in peter:
            f.write(str(i) + '\n')

LinesPlotter = CleanedDataPlotter()

LinesPlotter.data = ulysses_data
    
#LinesPlotter.data = wehry_data

LinesPlotter.raw_data = ulysses_data
LinesPlotter.min_quality_flag = min_quality_flag
LinesPlotter.rotation_angle_index = indices['rotation_angle_index']
LinesPlotter.quality_flag_index = indices['quality_flag_index']
LinesPlotter.one_year_et = one_year_et
#LinesPlotter.eff_area_file = 'DefaultDataset.csv'
LinesPlotter.eff_area_file = '30.dat'

LinesPlotter.execute_wehry()
