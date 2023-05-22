import numpy as np
import matplotlib.pyplot as plt
from indices_dict import indices
import spiceypy as spice
from MiscFunctions import find_nearest_idx

class CleanedDataPlotter:

    def execute_wehry(self):
        self._set_beta_parameters()
        self._remove_999_values()
        self._set_detector_sun_angles()
        self._plot_wehry()
        self._effective_area()
        self._plot_for_minimum_effective_area()
        self._plot_lat()
        self._calculate_zero_crossings_times()
        self._effective_area_with_lat()
        self._count_and_print_north_fraction()
        self._plot_flux()

    def _set_beta_parameters(self):
        self._wehry_velocity = 20
        self._wehry_angle = 50
        self._streams1 = [[1991.727,1991.740],[1991.948,1991.953],[1991.978,1991.983],[1992.017,1992.021],[1992.048,1992.056],[1992.192,1992.196],[1992.268,1992.275],[1992.342,1992.349],[1992.417,1992.435],[1992.670,1992.685],[1992.796,1992.811]]
        self._streams2 = [[2002.905,2002.917],[2003.515,2003.535],[2003.643,2003.663],[2003.717,2003.742],[2003.778,2003.802],[2003.862,2003.867],[2003.919,2003.928],[2003.993,2004.007],[2004.063,2004.078],[2004.138,2004.142],[2004.207,2004.231],[2004.410,2004.442],[2004.450,2004.510],[2004.544,2004.574],[2004.582,2004.603],[2004.625,2004.654],[2004.825,2004.833],[2004.907,2004.914],[2004.989,2005.001],[2005.113,2005.135],[2005.210,2005.236],[2005.474,2005.487],[2005.565,2005.583],[2005.615,2005.642]]
        self._instrument = [ [2000.4905, 2000.49625], [2002.23297, 2002.2700], [2002.917692,2003.422021], [2003.5733, 2003.642642], [2004.918349, 2004.923907]]
        self._interstellar_ecliptic_lon = 252
        self._interstellar_ecliptic_lat = 2.5
        self._tolerance = 30
        self._interstellar_min_vel = 14
        self._interstellar_min_mass = 2.5e-14

    def _remove_999_values(self):
        ulysses_data_without_999 = []
        for i in range(self.start_index, self.end_index):
            condition = self.data[i, self.rotation_angle_index] != 999 and self.data[i, self.quality_flag_index] >= self.min_quality_flag
            if condition:
                ulysses_data_without_999.append(self.data[i])
        self.data_without_999 = np.array(ulysses_data_without_999)

    def _set_detector_sun_angles(self):
        self._detector_sun_angles = self._calculate_angle_between_two_objects(self.data_without_999)

    def _calculate_angle_between_two_objects(self, data: list, target: str = 'SUN', observer: str = 'ULYSSES') -> list:
        lon = data[:, indices['solar_lon_index']]
        lat = data[:, indices['solar_lat_index']]
        detector_sun_angles = []
        et = data[:, self._time_index]
        for i in range(len(data[:,0])):
            [stateSun, ltime] = spice.spkezr('SUN',  et[i],      'J2000', 'NONE', 'ULYSSES')
            posSun = stateSun[:3]
            posDetector = spice.latrec(np.linalg.norm(posSun), lon[i], lat[i])
            detector_sun_angles.append(spice.vsep(posDetector, posSun)*360/(2*np.pi))
        return detector_sun_angles

    def _plot_wehry(self):       
        plt.xlabel('Angle between detector axis and sun [°]')
        plt.ylabel('Relative Velocity [km/s]')
        plt.title(self.current_year)
        
        
        self._velocities = self.data_without_999[:, indices['velocity_index']]
            
        plt.scatter(self._detector_sun_angles, self._velocities)
        plt.plot(self._detector_sun_angles, np.ones(len(self._detector_sun_angles))*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(len(self._detector_sun_angles))*self._wehry_angle, self._velocities, color = 'red')
        plt.show() 
        


        self._choose_beta_meteoroids()
        
        
         
        plt_angles = []
        plt_vel = []
        for i in range(len(self._beta_meteoroids)):
            if self._beta_meteoroids[i]:
                plt_angles.append(self._detector_sun_angles[i])
                plt_vel.append(self._velocities[i])
        
        plt.scatter(plt_angles, plt_vel)
        plt.plot(self._detector_sun_angles, np.ones(len(self._detector_sun_angles))*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(len(self._detector_sun_angles))*self._wehry_angle, self._velocities, color = 'red')
        plt.xlabel('Angle between detector axis and sun [°]')
        plt.ylabel('Relative Velocity [km/s]')
        plt.show() 

    def _choose_beta_meteoroids(self) -> list:
        
        mass = self.data_without_999[:, indices['mass_index']]
        
        streams1 = (np.array(self._streams1)-2000)*self.one_year_et
        streams2 = (np.array(self._streams2)-2000)*self.one_year_et
        #instrument = (np.array(self._instrument)-2000)*self.one_year_et
        
        streams = []
        for i in range(len(streams1)):
            streams.append(streams1[i])
        for i in range(len(streams2)):
            streams.append(streams2[i])
        
        lon = self.data_without_999[:, indices['solar_lon_index']]
        lat = self.data_without_999[:, indices['solar_lat_index']]
        et = self.data_without_999[:, self._time_index]
        
        beta_meteoroids = []
        for i, data in  enumerate(self.data_without_999):
            condition_wehry = self._detector_sun_angles[i] <= self._wehry_angle and self._velocities[i] >= self._wehry_velocity
            condition_not_interstellar_angle = lon[i] < self._interstellar_ecliptic_lon-self._tolerance or lon[i] > self._interstellar_ecliptic_lon+self._tolerance or lat[i] < self._interstellar_ecliptic_lat-self._tolerance or lat[i] > self._interstellar_ecliptic_lat+self._tolerance
            condition_not_interstellar_vel = self._velocities[i] > self._interstellar_min_vel
            condition_not_interstellar_mass = mass[i] > self._interstellar_min_mass
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
        
        self._beta_meteoroids = beta_meteoroids

    def _effective_area(self) -> (list, list):
        wehry = np.loadtxt('DefaultDataset.csv', delimiter = ',')
        plt.plot(wehry[:,0], wehry[:,1]*10000, label = 'Wehry', color = 'blue')
        eff_area_data = np.loadtxt(self.eff_area_file, delimiter = ',')
        if self.eff_area_file == 'DefaultDataset.csv':
            eff_area_data[:,1] = eff_area_data[:,1]*10000
        plt.plot(eff_area_data[:,0], eff_area_data[:,1], label = 'Own estimate', color = 'red')
        plt.xlabel('Year')
        plt.ylabel('Effective Area [cm$^2$]')
        plt.title(self.eff_area_file[:2]+'km/s')
        plt.legend()
        plt.show()
        eff_area_time = []
        for i in eff_area_data[:,0]:        
            eff_area_time.append(spice.str2et(str(i)[:4]+'-01-001T00:00')+float(str(i)[4:])*self.one_year_et)
        eff_area_time = np.array(eff_area_time)
        self._eff_area_data = eff_area_data
        self.eff_area_time = eff_area_time

    def _plot_for_minimum_effective_area(self):

        for i in range(1,2):
            self._correct_by_effective_area(min_eff_area=self.min_eff_area_factor*i)
            
            plt.xlabel('Distance [au]')
            plt.ylabel('Count')
            plt.title('Minimum effective area = ' + str(self.min_eff_area_factor*i) + r'cm$^2$')
            
            plt.hist(self._beta_dist)
            plt.show()
            print(len(self._beta_dist), 'beta meteoroids found')

        #Seems like 200cm2 is pretty good


        detector_sun_angles1 = self._calculate_angle_between_two_objects(self._new_data)

        plt.scatter(detector_sun_angles1, self._new_data[:, indices['velocity_index']])
        plt.plot(self._detector_sun_angles, np.ones(len(self._detector_sun_angles))*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(len(self._detector_sun_angles))*self._wehry_angle, self._velocities, color = 'red')
        plt.xlabel('Angle between detector axis and sun [°]')
        plt.ylabel('Relative Velocity [km/s]')
        plt.show() 

    def _correct_by_effective_area(self, min_eff_area: float = 200):
        new_data = []
        beta_dist = []
        eff_area_res = []
        for i in range(len(self.data_without_999)):
            if self._beta_meteoroids[i]:
                closest_eff_area_idx = find_nearest_idx(self.eff_area_time, self.data_without_999[i, self._time_index])
                if self._eff_area_data[closest_eff_area_idx,1] > min_eff_area:
                    eff_area_res.append(self._eff_area_data[closest_eff_area_idx,1])
                    new_data.append(self.data_without_999[i])
                    beta_dist.append(self.data_without_999[i,self._dist_index])
            
        if new_data == []:
            raise RuntimeError('No beta meteoroids found')
        
        beta_dist = np.array(beta_dist)
        new_data = np.array(new_data)
        eff_area_res = np.array(eff_area_res)
        self._new_data = new_data
        self._beta_dist = beta_dist
        self._eff_area_res = eff_area_res

    def _plot_lat(self):
        plt.scatter(self._new_data[:, self._time_index]/self.one_year_et+2000, self._new_data[:,indices['LAT']], label = 'beta', marker = 'x')
        plt.scatter(self.data[:, self._time_index]/self.one_year_et+2000, self.data[:,indices['LAT']], label = 'all', alpha = 0.5, marker = '.')
        plt.plot(self._new_data[:, self._time_index]/self.one_year_et+2000, np.zeros(len(self._new_data[:, self._time_index])), color = 'red')
        plt.xlabel('Time')
        plt.ylabel('Ulysses Ecliptic Latitude [°]')
        plt.legend()
        plt.show()

    def _calculate_zero_crossings_times(self) -> list:
        #https://stackoverflow.com/questions/3843017/efficiently-detect-sign-changes-in-python
        zero_crossings = np.where(np.diff(np.signbit(self.data[:,indices['LAT']])))[0]
        zero_crossings_times = self.data[zero_crossings, self._time_index]
        print('Latitude = 0 at', zero_crossings_times/self.one_year_et+2000)
        self._zero_crossings_times = zero_crossings_times
        
    def _effective_area_with_lat(self):
        zero_crossings_times_in_epoch = self._zero_crossings_times/self.one_year_et+2000
        data = np.loadtxt(self.eff_area_file, delimiter = ',')
        if self.eff_area_file == 'DefaultDataset.csv':
            data[:,1] = data[:,1]*10000
        cond = np.where(data[:,0] < zero_crossings_times_in_epoch[0])
        cond = np.intersect1d(cond, cond)
        plt.plot(data[cond,0], data[cond,1], label = 'Pre-Jupiter fly-by', color = 'black')
        for i in range(len(zero_crossings_times_in_epoch)-1):
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
            cond0 = np.where(data[:,0] > zero_crossings_times_in_epoch[i])
            cond1 = np.where(data[:,0] < zero_crossings_times_in_epoch[i+1])
            cond = np.intersect1d(cond0, cond1)
            plt.plot(data[cond,0], data[cond,1], label = label, color = color)
        cond = np.where(data[:,0] > zero_crossings_times_in_epoch[-1])
        cond = np.intersect1d(cond, cond)
        plt.plot(data[cond,0], data[cond,1], color = 'red')
        plt.xlabel('Year')
        plt.ylabel('Effective Area [cm$^2$]')
        plt.title(self.eff_area_file[:2]+'km/s')
        plt.legend()
        plt.show()
    
    def _count_and_print_north_fraction(self,):
        count0 = 0
        count999 = 0
        countbeta = 0
        for i in range(len(self.data)):
            if self.data[i,indices['LAT']] > 0:
                count0 += 1
        for i in range(len(self.data_without_999)):
            if self.data_without_999[i,indices['LAT']] > 0:
                count999 += 1
        for i in range(len(self._new_data)):
            if self._new_data[i,indices['LAT']] > 0:
                countbeta += 1
        print('All data north fraction =', count0/len(self.data), ', no 999 north fraction =', count999/len(self.data_without_999), ', beta north fraction =', countbeta/len(self._new_data))

    
    def _plot_flux(self):
        mean_dist_array, flux_bins = self._calculate_mean_eff_area_and_flux()
        
        
        label_array = ['Pre-Jupiter fly-by', 'South', 'North', None, None, None, None]
        color_array = ['black', 'blue', 'red', 'blue', 'red', 'blue', 'red']
        
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
        

    def _calculate_mean_eff_area_and_flux(self) -> (list, list):
        
        flux_bins = []
        mean_dist_array = []
        
        
        
        eff_area_array = []
        dist_array = []
        
        for i in range(len(self._new_data)):
            if self._new_data[i, self._time_index]<=self._zero_crossings_times[0]:
                eff_area_array.append(self._eff_area_res[i])
                dist_array.append(self._beta_dist[i])
        dist_array = np.array(dist_array)
        mean_dist_array.append(np.mean(dist_array))
        eff_area_array = np.array(eff_area_array)
        mean_eff_area = np.mean(eff_area_array)
        flux_bins.append(1/(mean_eff_area*(self._zero_crossings_times[1]-self.data[0, self._time_index])))
        
        for k in range(len(self._zero_crossings_times)-1):
            eff_area_array = []
            dist_array = []
            
            idx = find_nearest_idx(self._new_data[:, self._time_index], self._zero_crossings_times[k])
            
            for i in range(idx, len(self._new_data)):
                if self._new_data[i, self._time_index]<=self._zero_crossings_times[k+1]:
                    eff_area_array.append(self._eff_area_res[i])
                    dist_array.append(self._beta_dist[i])
            dist_array = np.array(dist_array)
            mean_dist_array.append(np.mean(dist_array))
            eff_area_array = np.array(eff_area_array)
            mean_eff_area = np.mean(eff_area_array)
            flux_bins.append(1/(mean_eff_area*(self._zero_crossings_times[k+1]-self._zero_crossings_times[k])))
            
            
        eff_area_array = []
        dist_array = []
        
        idx = find_nearest_idx(self._new_data[:, self._time_index], self._zero_crossings_times[-1])
        
        for i in range(idx, len(self._new_data)):
            eff_area_array.append(self._eff_area_res[i])
            dist_array.append(self._beta_dist[i])
        dist_array = np.array(dist_array)
        mean_dist_array.append(np.mean(dist_array))
        eff_area_array = np.array(eff_area_array)
        mean_eff_area = np.mean(eff_area_array)
        flux_bins.append(1/(mean_eff_area*(self.data[-1, self._time_index]-self._zero_crossings_times[-1])))
            
            
        mean_dist_array = np.array(mean_dist_array)
        flux_bins = np.array(flux_bins)
        
        return mean_dist_array, flux_bins
    

    def __init__(self):
        self.current_year = None
        
        self._time_index = 0
        self._dist_index = 1
        