import numpy as np
import matplotlib.pyplot as plt
from indices_dict import indices
import spiceypy as spice



class CleanedDataPlotter:

    def execute_wehry(self):
        self._remove_999_values()
        self._plot_wehry()

    def _remove_999_values(self):
        ulysses_data_without_999 = []
        for i in range(self.start_index, self.end_index):
            condition = self.data[i, self.rotation_angle_index] != 999 and self.data[i, self.quality_flag_index] >= self.min_quality_flag
            if condition:
                ulysses_data_without_999.append(self.data[i])
        self.data_without_999 = np.array(ulysses_data_without_999)

    def _plot_wehry(self):
        """
        plt.xlabel('Sector')
        plt.ylabel('Velocity')
        plt.title(self.current_year)
        plt.scatter(self.data_without_999[:, indices['sector']], self.data_without_999[:, indices['velocity_index']])
        wehry_velocity = 20
        wehry_sector = 50
        plt.plot(self.data_without_999[:, indices['sector']], np.ones(len(self.data_without_999[:, indices['sector']]))*wehry_velocity, color = 'red')
        plt.plot(np.ones(len(self.data_without_999[:, indices['sector']]))*wehry_sector, self.data_without_999[:, indices['velocity_index']], color = 'red')
        plt.show()
        
        beta_meteoroids = []
        for i in range(len(self.data)):
            condition = self.data[i, indices['sector']] <= wehry_sector and self.data[i, indices['velocity_index']] >= wehry_velocity
            if condition:
                beta_meteoroids.append(i)
        """
        
        plt.xlabel('Angle [°]')
        plt.ylabel('Velocity [km/s]')
        plt.title(self.current_year)
        
        
        velocities = self.data_without_999[:, indices['velocity_index']]
        
        plt.scatter(angles, velocities)
        wehry_velocity = 20
        wehry_angle = 50
        
        
        
        
        
        
        
        
        data, eff_area_time = self._effective_area()
        
        for i in range(10,11):
            new_data, beta_dist = self._correct_by_effective_area(beta_meteoroids, data, eff_area_time, min_eff_area=100*i/10)
            
            plt.xlabel('Distance [au]')
            plt.ylabel('Count')
            plt.title('Minimum effective area = ' + str(100*i/10) + r'cm$^2$')
            
            #plt.hist(beta_dist, weights = weights)
            plt.hist(beta_dist)
            plt.show()
            #print(sum(weights))

    def _effective_area(self):
        data = np.loadtxt('70.dat', delimiter = ',')
        plt.plot(data[:,0], data[:,1], label = '70km/s')
        plt.xlabel('Year')
        plt.ylabel('Effective Area [cm$^2$]')
        plt.legend()
        plt.show()
        eff_area_time = []
        for i in data[:,0]:        
            eff_area_time.append(spice.str2et(str(i)[:4]+'-01-001T00:00')+float(str(i)[4:])*one_year_et)
        eff_area_time = np.array(eff_area_time)
        return data, eff_area_time

    def _correct_by_effective_area(self, beta_meteoroids, data, eff_area_time, min_eff_area = 100):
        
        #weights = []     
        #for i in range(len(self.data_without_999)):
        #    if i in beta_meteoroids:
        #        closest_eff_area_idx = find_nearest_idx(eff_area_time, self.data_without_999[i,time_index])
        #        weights.append(data[closest_eff_area_idx,1]/max_area)
        #weights = np.array(weights)
        #return weights
        
        new_data = []
        beta_dist = []
        for i in range(len(self.data)):
            if i in beta_meteoroids:
                closest_eff_area_idx = find_nearest_idx(eff_area_time, self.data[i,time_index])
                if data[closest_eff_area_idx,1] > min_eff_area:
                    new_data.append(self.data[i])
                    beta_dist.append(self.dist[i])
        new_data = np.array(new_data)
        beta_dist = np.array(beta_dist)
        return new_data, beta_dist
            

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
        self.dist = None

#https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

max_area = 1000
time_index = 0

spice_path = '../../spice/'
spice.furnsh(spice_path + "naif0012.tls")
spice.furnsh(spice_path + "de440.bsp")
spice.furnsh(spice_path + "ulysses_1990_2009_2050.bsp")
one_year_et = spice.str2et('01-001T00:00')-spice.str2et('00-001T00:00')

first_data_column = 3
last_data_column = 34

min_quality_flag = 2
all_lines = True
bins = 10
PlottedQuantity = 'LAT'

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
dist_array = []

for et in time:
    [pos, ltime] = spice.spkpos('SUN',  et,      'J2000', 'NONE', 'ULYSSES')
    dist = spice.vnorm(pos)
    dist = spice.convrt(dist, 'KM', 'AU')
    dist_array.append(dist)

ulysses_data = np.loadtxt('Ulysses_Data_File_Cleaned.txt', delimiter = ' ', skiprows = indices['first_data_line'], usecols = used_cols)

ulysses_data = np.concatenate((np.transpose(np.array([time])), ulysses_data), axis = 1)

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
LinesPlotter.start_index = start_index
LinesPlotter.end_index = end_index
LinesPlotter.min_quality_flag = min_quality_flag
LinesPlotter.rotation_angle_index = indices['rotation_angle_index']
LinesPlotter.quality_flag_index = indices['quality_flag_index']
LinesPlotter.xlabel = PlottedQuantity
LinesPlotter.plot_index = indices[PlottedQuantity]
LinesPlotter.bins = bins
LinesPlotter.dist = dist_array

i = 1990
while True:
    if not all_lines:
        LinesPlotter.current_year = i
    if LinesPlotter.end_index == len(ulysses_data)-1:
        break
    LinesPlotter.execute_wehry()
    if not all_lines:
        start_et = end_et
        end_et = end_et+one_year_et
        LinesPlotter.start_index = np.argmin(np.abs(time-start_et))
        LinesPlotter.end_index = np.argmin(np.abs(time-end_et))
        i += 1
    if all_lines:
        break
