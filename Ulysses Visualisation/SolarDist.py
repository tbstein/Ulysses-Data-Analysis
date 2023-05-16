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
        plt.xlabel('Angle between detector axis and sun [°]')
        plt.ylabel('Relative Velocity [km/s]')
        plt.title(self.current_year)
        
        
        velocities = self.data_without_999[:, indices['velocity_index']]
        mass = self.data_without_999[:, indices['mass_index']]
        lon = self.data_without_999[:, indices['solar_lon_index']]
        lat = self.data_without_999[:, indices['solar_lat_index']]
        
        detector_sun_angles = []
        et = self.data_without_999[:, time_index]
        for i in range(len(self.data_without_999[:,0])):
            [stateSun, ltime] = spice.spkezr('SUN',  et[i],      'J2000', 'NONE', 'ULYSSES')
            posSun = stateSun[:3]
            posDetector = spice.latrec(np.linalg.norm(posSun), lon[i], lat[i])
            detector_sun_angles.append(spice.vsep(posDetector, posSun)*360/(2*np.pi))
            
        plt.scatter(detector_sun_angles, velocities)
        wehry_velocity = 20
        wehry_angle = 50
        plt.plot(detector_sun_angles, np.ones(len(detector_sun_angles))*wehry_velocity, color = 'red')
        plt.plot(np.ones(len(detector_sun_angles))*wehry_angle, velocities, color = 'red')
        plt.show() 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        streams1 = [[1991.727,1991.740],[1991.948,1991.953],[1991.978,1991.983],[1992.017,1992.021],[1992.048,1992.056],[1992.192,1992.196],[1992.268,1992.275],[1992.342,1992.349],[1992.417,1992.435],[1992.670,1992.685],[1992.796,1992.811]]
        streams2 = [[2002.905,2002.917],[2003.515,2003.535],[2003.643,2003.663],[2003.717,2003.742],[2003.778,2003.802],[2003.862,2003.867],[2003.919,2003.928],[2003.993,2004.007],[2004.063,2004.078],[2004.138,2004.142],[2004.207,2004.231],[2004.410,2004.442],[2004.450,2004.510],[2004.544,2004.574],[2004.582,2004.603],[2004.625,2004.654],[2004.825,2004.833],[2004.907,2004.914],[2004.989,2005.001],[2005.113,2005.135],[2005.210,2005.236],[2005.474,2005.487],[2005.565,2005.583],[2005.615,2005.642]]
        instrument = [ [2000.4905, 2000.49625], [2002.23297, 2002.2700], [2002.917692,2003.422021], [2003.5733, 2003.642642], [2004.918349, 2004.923907]]
        
        streams1 = (np.array(streams1)-2000)*one_year_et
        streams2 = (np.array(streams2)-2000)*one_year_et
        instrument = (np.array(instrument)-2000)*one_year_et
        
        streams = []
        for i in range(len(streams1)):
            streams.append(streams1[i])
        for i in range(len(streams2)):
            streams.append(streams2[i])
        
        interstellar_ecliptic_lon = 252
        interstellar_ecliptic_lat = 2.5
        tolerance = 30
        interstellar_min_vel = 14
        interstellar_min_mass = 0 #TODO convert mass in V to mass in kg or g
        
        beta_meteoroids = []
        for i, data in  enumerate(self.data_without_999):
            condition_wehry = detector_sun_angles[i] <= wehry_angle and velocities[i] >= wehry_velocity
            condition_not_interstellar_angle = lon[i] < interstellar_ecliptic_lon-tolerance or lon[i] > interstellar_ecliptic_lon+tolerance or lat[i] < interstellar_ecliptic_lat-tolerance or lat[i] > interstellar_ecliptic_lat+tolerance
            condition_not_interstellar_vel = velocities[i] > interstellar_min_vel
            condition_not_interstellar_mass = mass[i] > interstellar_min_mass
            condition_not_interstellar = condition_not_interstellar_angle and condition_not_interstellar_vel and condition_not_interstellar_mass
            condition_not_streams = True
            for interval in streams:
                if et[i] >= interval[0] and et[i] <= interval[1]:
                    condition_not_streams = False
            
            #condition_wehry = True
            #condition_not_streams = True
            #condition_not_interstellar = True
            condition = condition_wehry and condition_not_streams and condition_not_interstellar
            if condition:
                beta_meteoroids.append(True)  
            else:
                beta_meteoroids.append(False)  
        
        plt_angles = []
        plt_vel = []
        for i in range(len(beta_meteoroids)):
            if beta_meteoroids[i]:
                plt_angles.append(detector_sun_angles[i])
                plt_vel.append(velocities[i])
        
        plt.scatter(plt_angles, plt_vel)
        plt.plot(detector_sun_angles, np.ones(len(detector_sun_angles))*wehry_velocity, color = 'red')
        plt.plot(np.ones(len(detector_sun_angles))*wehry_angle, velocities, color = 'red')
        plt.xlabel('Angle between detector axis and sun [°]')
        plt.ylabel('Relative Velocity [km/s]')
        plt.show() 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        data, eff_area_time = self._effective_area()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        for i in range(1,2):
            new_data, beta_dist, eff_area_res = self._correct_by_effective_area(beta_meteoroids, data, eff_area_time, min_eff_area=min_eff_area_factor*i)
            
            plt.xlabel('Distance [au]')
            plt.ylabel('Count')
            plt.title('Minimum effective area = ' + str(min_eff_area_factor*i) + r'cm$^2$')
            
            plt.hist(beta_dist)
            plt.show()
            print(len(beta_dist))

        #Seems like 200cm2 is pretty good





























        lon = new_data[:, indices['solar_lon_index']]
        lat = new_data[:, indices['solar_lat_index']]
        detector_sun_angles1 = []
        et = new_data[:, time_index]
        for i in range(len(new_data[:,0])):
            [stateSun, ltime] = spice.spkezr('SUN',  et[i],      'J2000', 'NONE', 'ULYSSES')
            posSun = stateSun[:3]
            posDetector = spice.latrec(np.linalg.norm(posSun), lon[i], lat[i])
            detector_sun_angles1.append(spice.vsep(posDetector, posSun)*360/(2*np.pi))

        plt.scatter(detector_sun_angles1, new_data[:, indices['velocity_index']])
        plt.plot(detector_sun_angles, np.ones(len(detector_sun_angles))*wehry_velocity, color = 'red')
        plt.plot(np.ones(len(detector_sun_angles))*wehry_angle, velocities, color = 'red')
        plt.xlabel('Angle between detector axis and sun [°]')
        plt.ylabel('Relative Velocity [km/s]')
        plt.show() 





















        
        plt.scatter(new_data[:,time_index]/one_year_et+2000, new_data[:,indices['LAT']], label = 'beta', marker = 'x')
        plt.scatter(self.data[:,time_index]/one_year_et+2000, self.data[:,indices['LAT']], label = 'all', alpha = 0.5, marker = '.')
        plt.plot(new_data[:,time_index]/one_year_et+2000, np.zeros(len(new_data[:,time_index])), color = 'red')
        plt.xlabel('Time')
        plt.ylabel('Ulysses Ecliptic Latitude [°]')
        plt.legend()
        plt.show()
    















    
        #https://stackoverflow.com/questions/3843017/efficiently-detect-sign-changes-in-python
        zero_crossings = np.where(np.diff(np.signbit(self.data[:,indices['LAT']])))[0]
        zero_crossings_times = self.data[zero_crossings,time_index]
        print('Latitude = 0 at', zero_crossings_times/one_year_et+2000)
        self._effective_area_with_lat(zero_crossings_times/one_year_et+2000)
       
        count0 = 0
        count999 = 0
        countbeta = 0
        for i in range(len(self.data)):
            if self.data[i,indices['LAT']] > 0:
                count0 += 1
        for i in range(len(self.data_without_999)):
            if self.data_without_999[i,indices['LAT']] > 0:
                count999 += 1
        for i in range(len(new_data)):
            if new_data[i,indices['LAT']] > 0:
                countbeta += 1
        print('All data north fraction =', count0/len(self.data), ', no 999 north fraction =', count999/len(self.data_without_999), ', beta north fraction =', countbeta/len(new_data))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        flux_bins = []
        mean_dist_array = []
        label_array = ['Pre-Jupiter fly-by', 'South', 'North', None, None, None, None]
        color_array = ['black', 'blue', 'red', 'blue', 'red', 'blue', 'red']
        
        
        
        eff_area_array = []
        dist_array = []
        
        for i in range(len(new_data)):
            if new_data[i,time_index]<=zero_crossings_times[0]:
                eff_area_array.append(eff_area_res[i])
                dist_array.append(beta_dist[i])
        dist_array = np.array(dist_array)
        mean_dist_array.append(np.mean(dist_array))
        eff_area_array = np.array(eff_area_array)
        mean_eff_area = np.mean(eff_area_array)
        flux_bins.append(1/(mean_eff_area*(zero_crossings_times[1]-self.data[0,time_index])))
        
        for k in range(len(zero_crossings_times)-1):
            eff_area_array = []
            dist_array = []
            
            idx = find_nearest_idx(new_data[:,time_index], zero_crossings_times[k])
            
            for i in range(idx, len(new_data)):
                if new_data[i,time_index]<=zero_crossings_times[k+1]:
                    eff_area_array.append(eff_area_res[i])
                    dist_array.append(beta_dist[i])
            dist_array = np.array(dist_array)
            mean_dist_array.append(np.mean(dist_array))
            eff_area_array = np.array(eff_area_array)
            mean_eff_area = np.mean(eff_area_array)
            flux_bins.append(1/(mean_eff_area*(zero_crossings_times[k+1]-zero_crossings_times[k])))
            
            
        eff_area_array = []
        dist_array = []
        
        idx = find_nearest_idx(new_data[:,time_index], zero_crossings_times[-1])
        
        for i in range(idx, len(new_data)):
            eff_area_array.append(eff_area_res[i])
            dist_array.append(beta_dist[i])
        dist_array = np.array(dist_array)
        mean_dist_array.append(np.mean(dist_array))
        eff_area_array = np.array(eff_area_array)
        mean_eff_area = np.mean(eff_area_array)
        flux_bins.append(1/(mean_eff_area*(self.data[-1,time_index]-zero_crossings_times[-1])))
            
            
        mean_dist_array = np.array(mean_dist_array)
        flux_bins = np.array(flux_bins)
        
        for i in range(len(color_array)):
            plt.scatter(mean_dist_array[i], flux_bins[i]*10000, label = label_array[i], color = color_array[i])
        plt.xlabel('Distance [au]')
        plt.ylabel(r'Flux [1/(m$^2\cdot$s)]')
        plt.yscale('log')
        plt.legend()
        plt.show()
        
        for i in range(len(color_array)-1):
            plt.scatter(mean_dist_array[i], flux_bins[i]*10000, label = label_array[i], color = color_array[i])
        plt.xlabel('Distance [au]')
        plt.ylabel(r'Flux [1/(m$^2\cdot$s)]')
        plt.yscale('log')
        plt.legend()
        plt.show()
        

    def _effective_area(self):
        wehry = np.loadtxt('DefaultDataset.csv', delimiter = ',')
        plt.plot(wehry[:,0], wehry[:,1]*10000, label = 'Wehry', color = 'blue')
        data = np.loadtxt(eff_area_file, delimiter = ',')
        if eff_area_file =='DefaultDataset.csv':
            data[:,1] = data[:,1]*10000
        plt.plot(data[:,0], data[:,1], label = 'Own estimate', color = 'red')
        plt.xlabel('Year')
        plt.ylabel('Effective Area [cm$^2$]')
        plt.title(eff_area_file[:2]+'km/s')
        plt.legend()
        plt.show()
        eff_area_time = []
        for i in data[:,0]:        
            eff_area_time.append(spice.str2et(str(i)[:4]+'-01-001T00:00')+float(str(i)[4:])*one_year_et)
        eff_area_time = np.array(eff_area_time)
        return data, eff_area_time

    def _effective_area_with_lat(self, zero_crossing_times):
        data = np.loadtxt(eff_area_file, delimiter = ',')
        if eff_area_file =='DefaultDataset.csv':
            data[:,1] = data[:,1]*10000
        cond = np.where(data[:,0] < zero_crossing_times[0])
        cond = np.intersect1d(cond, cond)
        plt.plot(data[cond,0], data[cond,1], label = 'Pre-Jupiter fly-by', color = 'black')
        for i in range(len(zero_crossing_times)-1):
            if i == 0:
                label = 'South'
            elif i == 1:
                label = 'North'
            else:
                label = None
            if i % 2 == 0:
                color = 'blue'
            else:
                color = 'red'
            cond0 = np.where(data[:,0] > zero_crossing_times[i])
            cond1 = np.where(data[:,0] < zero_crossing_times[i+1])
            cond = np.intersect1d(cond0, cond1)
            plt.plot(data[cond,0], data[cond,1], label = label, color = color)
        cond = np.where(data[:,0] > zero_crossing_times[-1])
        cond = np.intersect1d(cond, cond)
        plt.plot(data[cond,0], data[cond,1], color = 'red')
        plt.xlabel('Year')
        plt.ylabel('Effective Area [cm$^2$]')
        plt.title(eff_area_file[:2]+'km/s')
        plt.legend()
        plt.show()

    def _correct_by_effective_area(self, beta_meteoroids, data, eff_area_time, min_eff_area = 200):
        new_data = []
        beta_dist = []
        eff_area_res = []
        for i in range(len(self.data_without_999)):
            if beta_meteoroids[i]:
                closest_eff_area_idx = find_nearest_idx(eff_area_time, self.data_without_999[i,time_index])
                if data[closest_eff_area_idx,1] > min_eff_area:
                    eff_area_res.append(data[closest_eff_area_idx,1])
                    new_data.append(self.data_without_999[i])
                    beta_dist.append(self.data_without_999[i,dist_index])
        beta_dist = np.array(beta_dist)
        new_data = np.array(new_data)
        eff_area_res = np.array(eff_area_res)
        return new_data, beta_dist, eff_area_res
            

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

#https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

min_eff_area_factor = 50
eff_area_file = 'DefaultDataset.csv'
#eff_area_file = '20.dat'

max_area = 1000
time_index = 0
dist_index = 1

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

ulysses_data = np.concatenate((np.transpose(np.array([dist_array])), ulysses_data), axis = 1)

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
