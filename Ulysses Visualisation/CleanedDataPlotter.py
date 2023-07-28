#Things to do
#Check if Peters stream times are correct
#A future master thesis would simulate ISD fluxes, their plasma clouds, their rise times and charge yields, and therefore get simulated wall impact vs non wall impact data and compare that to the measurements
#Markus Fraenz fragen wegen symIon. Wie einfach ist den Ulysses Dust Detector zu simulieren
#Remove streams by filtering put low intensity classes
#Do some error propagation
#Do dofferent binning for flux

import numpy as np
import matplotlib.pyplot as plt
from indices_dict import indices
import spiceypy as spice
from MiscFunctions import find_nearest_idx

class CleanedDataPlotter:
    
    def execute_wehry(self):
        self._remove_999_values()
        self._set_angles_and_velocities()
        self._plot_wehry(self._detector_sun_angles, self._velocities)
        self._plot_rot(self.data_without_999)
        self._choose_beta_meteoroids()
        self._plot_wehry(self._beta_angles, self._beta_vel)
        self._effective_area()
        self._minimum_effective_area_hist()
        self._set_new_angles_and_velocities()
        #self._plot_wehry_time_separated(self._new_detector_sun_angles, self._new_velocities)
        self._plot_rot(self._new_data)
        self._compare_found_betas()
        self._plot_lat()
        self._calculate_zero_crossings_times()
        self._effective_area_with_lat()
        self._count_and_print_north_fraction()
        self._plot_effective_area_with_streams()
        self._plot_flux()

    def _remove_999_values(self):
        ulysses_data_without_999 = []
        for i in range(len(self.data[:, self.rotation_angle_index])):
            condition = self.data[i, self.rotation_angle_index] != 999 and self.data[i, self.quality_flag_index] >= self.min_quality_flag
            if condition:
                ulysses_data_without_999.append(self.data[i])
        self.data_without_999 = np.array(ulysses_data_without_999)

    def _set_angles_and_velocities(self):
        self._detector_sun_angles, self._velocities = self._calculate_angle_and_velocities_between_two_objects(self.data_without_999)
    
    def _calculate_angle_and_velocities_between_two_objects(self, data: list) -> (list, list):
        lon = data[:, indices['solar_lon_index']]
        lat = data[:, indices['solar_lat_index']]
        detector_sun_angles = []
        velDust = []
        velUlysses0 = []
        et = data[:, self._time_index]
        for i in range(len(data[:,0])):
            [stateUlysses, ltime] = spice.spkezr('ULYSSES',  et[i],      'ECLIPJ2000', 'NONE', 'SUN')
            posUlysses = stateUlysses[:3]
            velUlysses = stateUlysses[3:]
            velUlysses0.append(np.linalg.norm(velUlysses))
 
            velImpact = data[i, indices['velocity_index']]
            
            pointingDetector = spice.latrec(1, lon[i]*2*np.pi/360, lat[i]*2*np.pi/360)
            
            vel = velUlysses-velImpact*pointingDetector
            
            angle0 = spice.vsep(vel, posUlysses)*360/(2*np.pi)
            detector_sun_angles.append(angle0)
            
            """
            e_r = np.ones(3)/np.sqrt(3)
            scalar_product = np.dot(e_r, vel)
            
            velDust.append(scalar_product)
            """
            
            velDust.append(np.linalg.norm(vel))
            
        return detector_sun_angles, velDust
        

    def _plot_wehry(self, angles: list, velocities: list):       
        plt.xlabel('Angle between dust flow and sun in solar reference frame [°]')
        plt.ylabel('Radial Particle Velocity [km/s]')          
        plt.scatter(angles,velocities)
        plt.plot(self._detector_sun_angles, np.ones(len(self._detector_sun_angles))*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(len(self._detector_sun_angles))*self._wehry_angle, self._velocities, color = 'red')
        plt.show()
        
    def _plot_wehry_time_separated(self, angles: list, velocities: list):
        jupiter_flyby_et = -2.47398962e+08
        first_orbit_start_et = spice.str2et('1992-02-01')
        first_orbit_end_et = spice.str2et('1995-08-01')
        second_orbit_start_et = spice.str2et('1995-08-16')
        second_orbit_end_et = spice.str2et('2001-10-19')
        first_orbit_north_start_et = spice.str2et('1995-02-01')
        first_orbit_north_end_et = spice.str2et('1995-08-01')
        first_orbit_south_start_et = spice.str2et('1992-02-01')
        first_orbit_south_end_et = spice.str2et('1995-02-01')
        
        angles = np.array(angles)
        velocities = np.array(velocities)
        
        xlim = 120
        ylim = 100
        plt.xlabel('Angle between dust flow and sun in solar reference frame [°]')
        plt.ylabel('Radial Particle Velocity [km/s]')
        plt.title('Pre fly-by')
        pltangles = angles[np.array(np.where(self.data[:,self._time_index] < jupiter_flyby_et))[0]]
        pltvelocities = velocities[np.array(np.where(self.data[:,self._time_index] < jupiter_flyby_et))[0]]
        plt.scatter(pltangles,pltvelocities)
        plt.plot(np.linspace(0, xlim, num = ylim), np.ones(ylim)*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(xlim)*self._wehry_angle, np.linspace(0, ylim, num = xlim), color = 'red')
        plt.xlim(0,120)
        plt.ylim(0,100)
        plt.show()
        
        xlim = 200
        ylim = 100
        plt.xlabel('Angle between dust flow and sun in solar reference frame [°]')
        plt.ylabel('Radial Particle Velocity [km/s]')
        plt.title('First orbit')
        pltangles = angles[np.array(np.where((self.data[:,self._time_index] > first_orbit_start_et) & (self.data[:,self._time_index] < first_orbit_end_et)))[0]]
        pltvelocities = velocities[np.array(np.where((self.data[:,self._time_index] > first_orbit_start_et) & (self.data[:,self._time_index] < first_orbit_end_et)))[0]]
        plt.scatter(pltangles,pltvelocities)
        plt.plot(np.linspace(0, xlim, num = ylim), np.ones(ylim)*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(xlim)*self._wehry_angle, np.linspace(0, ylim, num = xlim), color = 'red')
        plt.xlim(0,200)
        plt.ylim(0,100)
        plt.show()
        
        xlim = 200
        ylim = 100
        plt.xlabel('Angle between dust flow and sun in solar reference frame [°]')
        plt.ylabel('Radial Particle Velocity [km/s]')
        plt.title('First orbit outside ecliptic')
        pltangles = angles[np.array(np.where((self.data[:,self._time_index] > first_orbit_start_et) & (self.data[:,self._time_index] < first_orbit_end_et) & (np.abs(self.data[:,indices['LAT']]) >= 18)))[0]]
        pltvelocities = velocities[np.array(np.where((self.data[:,self._time_index] > first_orbit_start_et) & (self.data[:,self._time_index] < first_orbit_end_et) & (np.abs(self.data[:,indices['LAT']]) >= 18)))[0]]
        plt.scatter(pltangles,pltvelocities)
        plt.plot(np.linspace(0, xlim, num = ylim), np.ones(ylim)*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(xlim)*self._wehry_angle, np.linspace(0, ylim, num = xlim), color = 'red')
        plt.xlim(0,200)
        plt.ylim(0,100)
        plt.show()
        
        xlim = 200
        ylim = 80
        plt.xlabel('Angle between dust flow and sun in solar reference frame [°]')
        plt.ylabel('Radial Particle Velocity [km/s]')
        plt.title('Second orbit')
        pltangles = angles[np.array(np.where((self.data[:,self._time_index] > second_orbit_start_et) & (self.data[:,self._time_index] < second_orbit_end_et)))[0]]
        pltvelocities = velocities[np.array(np.where((self.data[:,self._time_index] > second_orbit_start_et) & (self.data[:,self._time_index] < second_orbit_end_et)))[0]]
        plt.scatter(pltangles,pltvelocities)
        plt.plot(np.linspace(0, xlim, num = ylim), np.ones(ylim)*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(xlim)*self._wehry_angle, np.linspace(0, ylim, num = xlim), color = 'red')
        plt.xlim(0,200)
        plt.ylim(0,80)
        plt.show()
        
        xlim = 100
        ylim = 100
        plt.xlabel('Angle between dust flow and sun in solar reference frame [°]')
        plt.ylabel('Radial Particle Velocity [km/s]')
        plt.title('First orbit north')
        pltangles = angles[np.array(np.where((self.data[:,self._time_index] > first_orbit_north_start_et) & (self.data[:,self._time_index] < first_orbit_north_end_et)))[0]]
        pltvelocities = velocities[np.array(np.where((self.data[:,self._time_index] > first_orbit_north_start_et) & (self.data[:,self._time_index] < first_orbit_north_end_et)))[0]]
        plt.scatter(pltangles,pltvelocities)
        plt.plot(np.linspace(0, xlim, num = ylim), np.ones(ylim)*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(xlim)*self._wehry_angle, np.linspace(0, ylim, num = xlim), color = 'red')
        plt.xlim(0,100)
        plt.ylim(0,100)
        plt.show()
        
        xlim = 200
        ylim = 100
        plt.xlabel('Angle between dust flow and sun in solar reference frame [°]')
        plt.ylabel('Radial Particle Velocity [km/s]')
        plt.title('First orbit south')
        pltangles = angles[np.array(np.where((self.data[:,self._time_index] > first_orbit_south_start_et) & (self.data[:,self._time_index] < first_orbit_south_end_et)))[0]]
        pltvelocities = velocities[np.array(np.where((self.data[:,self._time_index] > first_orbit_south_start_et) & (self.data[:,self._time_index] < first_orbit_south_end_et)))[0]]
        plt.scatter(pltangles,pltvelocities)
        plt.plot(np.linspace(0, xlim, num = ylim), np.ones(ylim)*self._wehry_velocity, color = 'red')
        plt.plot(np.ones(xlim)*self._wehry_angle, np.linspace(0, ylim, num = xlim), color = 'red')
        plt.xlim(0,200)
        plt.ylim(0,100)
        plt.show()
        
        
        
    def _plot_rot(self, data):
        plt.xlabel('Time')
        plt.ylabel('Rotation angle [°]')        
        plt.scatter(data[:,self._time_index]/self.one_year_et+2000, data[:,indices['rotation_angle_index']])
        plt.show() 

    def _choose_beta_meteoroids(self) -> list:
        
        """
        altobelliTimes = [[spice.str2et('92-300T00:00'), spice.str2et('94-001T00:00')],[spice.str2et('96-150T00:00'), spice.str2et('97-250T00:00')],[spice.str2et('98-350T00:00'), spice.str2et('00-050T00:00')]]
        temp = []
        for impact in self.data_without_999:
            for i in altobelliTimes:
                if impact[self._time_index] > i[0] and impact[self._time_index] < i[1]:
                    temp.append(impact)
        self.data_without_999 = np.array(temp)
        """
        
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
        beta_angles = []
        beta_vel = []
        
        
        ISD_direction = spice.latrec(1, self._interstellar_ecliptic_lon*2*np.pi/360, self._interstellar_ecliptic_lat*2*np.pi/360)
        detector_ISD_angles = []
        velDust = []
        et = self.data_without_999[:, self._time_index]
        for i in range(len(self.data_without_999[:,0])):
            [stateUlysses, ltime] = spice.spkezr('ULYSSES',  et[i],      'ECLIPJ2000', 'NONE', 'SUN')
            velUlysses = stateUlysses[3:]
            velImpact = self.data_without_999[i, indices['velocity_index']]
            pointingDetector = spice.latrec(1, lon[i]*2*np.pi/360, lat[i]*2*np.pi/360)
            vel = velUlysses-velImpact*pointingDetector
            detector_ISD_angles.append(spice.vsep(vel, ISD_direction)*360/(2*np.pi))
            velDust.append(np.linalg.norm(vel))
        plt.scatter(et/self.one_year_et+2000, detector_ISD_angles)
        plt.plot(et/self.one_year_et+2000, self._tolerance*np.ones(len(detector_ISD_angles)), color = 'red')
        plt.xlabel('Time')
        plt.ylabel('Angle between dust velocity and ISD [°]')
        plt.show()
        
        detector_ISD_angles = np.array(detector_ISD_angles)
        velDust = np.array(velDust)
        mass = np.array(mass)
        
        condition_23 = []
        condition_01 = []
        for i in range(len(mass)):
            condition_23.append(self.data_without_999[i, indices['quality_flag_index']] > 1.5)
            condition_01.append(self.data_without_999[i, indices['quality_flag_index']] < 1.5)
            
        condition_23 = np.array(condition_23)
        conditiion_01 = np.array(condition_01)
        
        plt.hist(detector_ISD_angles[condition_23]-180, label = '2&3')
        plt.hist(detector_ISD_angles[condition_01]-180, label = '0&1', alpha = 0.5)
        plt.ylabel('Count')
        plt.xlabel('Angle between dust velocity and ISD [°]')
        plt.title('No ROT = 999')
        plt.legend()
        plt.show()
        
        condition_23 = []
        condition_01 = []
        for i in range(len(mass)):
            condition_23.append(velDust[i] > self._interstellar_min_vel and mass[i] > self._interstellar_min_mass and self.data_without_999[i, indices['quality_flag_index']] > 1.5)
            condition_01.append(velDust[i] > self._interstellar_min_vel and mass[i] > self._interstellar_min_mass and self.data_without_999[i, indices['quality_flag_index']] < 1.5)
            
        condition_23 = np.array(condition_23)
        conditiion_01 = np.array(condition_01)
        
        plt.hist(detector_ISD_angles[condition_23]-180, label = '2&3')
        plt.hist(detector_ISD_angles[condition_01]-180, label = '0&1', alpha = 0.5)
        plt.ylabel('Count')
        plt.xlabel('Angle between dust velocity and ISD [°]')
        plt.title('ISD min mass and min velocity filtering')
        plt.legend()
        plt.show()
        
        condition_23 = []
        condition_01 = []
        for i in range(len(mass)):
            condition_23.append(velDust[i] > self._interstellar_min_vel and velDust[i] < self._interstellar_max_vel and mass[i] > self._interstellar_min_mass and mass[i] < self._interstellar_max_mass and self.data_without_999[i, indices['quality_flag_index']] > 1.5)
            condition_01.append(velDust[i] > self._interstellar_min_vel and velDust[i] < self._interstellar_max_vel and mass[i] > self._interstellar_min_mass and mass[i] < self._interstellar_max_mass and self.data_without_999[i, indices['quality_flag_index']] < 1.5)
            
        condition_23 = np.array(condition_23)
        conditiion_01 = np.array(condition_01)
        
        plt.hist(detector_ISD_angles[condition_23]-180, label = '2&3')
        plt.hist(detector_ISD_angles[condition_01]-180, label = '0&1', alpha = 0.5)
        plt.ylabel('Count')
        plt.xlabel('Angle between dust velocity and ISD [°]')
        plt.title('ISD min/max mass and min/max velocity filtering')
        plt.legend()
        plt.show()
        
        condition_23 = []
        condition_01 = []
        for i in range(len(mass)):
            condition_not_streams = True
            for interval in streams:
                if et[i] >= interval[0] and et[i] <= interval[1]:
                    condition_not_streams = False
            if condition_not_streams:
                condition_23.append(velDust[i] > self._interstellar_min_vel and velDust[i] < self._interstellar_max_vel and mass[i] > self._interstellar_min_mass and mass[i] < self._interstellar_max_mass and self.data_without_999[i, indices['quality_flag_index']] > 1.5)
                condition_01.append(velDust[i] > self._interstellar_min_vel and velDust[i] < self._interstellar_max_vel and mass[i] > self._interstellar_min_mass and mass[i] < self._interstellar_max_mass and self.data_without_999[i, indices['quality_flag_index']] < 1.5)
            else:
                condition_23.append(False)
                condition_01.append(False)
            
        condition_23 = np.array(condition_23)
        conditiion_01 = np.array(condition_01)
        
        plt.hist(detector_ISD_angles[condition_23]-180, label = '2&3')
        plt.hist(detector_ISD_angles[condition_01]-180, label = '0&1', alpha = 0.5)
        plt.ylabel('Count')
        plt.xlabel('Angle between dust velocity and ISD [°]')
        plt.title('ISD min/max mass and min/max velocity and streams filtering')
        plt.legend()
        plt.show()
        
        
        plt.hist(self.data_without_999[np.where(self.data_without_999[:, indices['quality_flag_index']] > 1.5),indices['rotation_angle_index']][0], label = '2&3')
        plt.hist(self.data_without_999[np.where(self.data_without_999[:, indices['quality_flag_index']] < 1.5),indices['rotation_angle_index']][0], label = '0&1', alpha = 0.5)
        plt.ylabel('Count')
        plt.xlabel('Rotation Angle [°]')
        plt.legend()
        plt.show()
        
        condition_23 = []
        condition_01 = []
        for i in range(len(mass)):
            condition_23.append(self.data_without_999[i,self._dist_index] > 3 and self.data_without_999[i, indices['quality_flag_index']] > 1.5)
            condition_01.append(self.data_without_999[i,self._dist_index] > 3 and self.data_without_999[i, indices['quality_flag_index']] < 1.5)
        
        plt.hist(self.data_without_999[condition_23,indices['rotation_angle_index']], label = '2&3')
        plt.hist(self.data_without_999[condition_01,indices['rotation_angle_index']], label = '0&1', alpha = 0.5)
        plt.ylabel('Count')
        plt.xlabel('Rotation Angle [°]')
        plt.title('Dist > 3AU')
        plt.legend()
        plt.show()
        
        for i, data in  enumerate(self.data_without_999):
            #A minimum requirement for beta particles is a certain minimum velocity and maximum detector_sun_angle
            condition_wehry = self._detector_sun_angles[i] <= self._wehry_angle and self._velocities[i] >= self._wehry_velocity
            #Particles are classified ISD if they come from a certain direction, have a certain minimum velocity and a certain minimum mass
            condition_interstellar_angle = detector_ISD_angles[i] > 180 + self._tolerance
            condition_interstellar_vel = velDust[i] > self._interstellar_min_vel
            """
            Lennart anschreiben fuer siene IDentifikation
            """
            condition_interstellar_mass = mass[i] > self._interstellar_min_mass
            condition_interstellar_mass = True
            condition_not_interstellar = not (condition_interstellar_angle and condition_interstellar_vel and condition_interstellar_mass)
            
            #Particles are considered as jupiter steam particles if they were detected within Strub's jupiter stream times
            condition_not_streams = True
            for interval in streams:
                if et[i] >= interval[0] and et[i] <= interval[1]:
                    condition_not_streams = False
            
            #Particles are considered beta if they pass the minimum requirement and are not ISD or jupiter stream particles
            condition = condition_wehry and condition_not_streams and condition_not_interstellar
            if condition:
                beta_meteoroids.append(True)
                beta_angles.append(self._detector_sun_angles[i])
                beta_vel.append(self._velocities[i])
            else:
                beta_meteoroids.append(False)
        
        self._beta_meteoroids = beta_meteoroids
        self._beta_angles = beta_angles
        self._beta_vel = beta_vel

    def _effective_area(self) -> (list, list):
        eff_area_data = np.loadtxt(self.eff_area_file, delimiter = ',')
        if self.eff_area_file == 'DefaultDataset.csv':
            eff_area_data[:,1] = eff_area_data[:,1]*10000
        eff_area_time = []
        for i in eff_area_data[:,0]:        
            eff_area_time.append(spice.str2et(str(i)[:4]+'-01-001T00:00')+float(str(i)[4:])*self.one_year_et)
        eff_area_time = np.array(eff_area_time)
        self._eff_area_data = eff_area_data
        self.eff_area_time = eff_area_time

    def _minimum_effective_area_hist(self):
        
        GM_km3_per_s2 = 1.327e11
        au_in_km  = 149597870.7
        
        peri = []
        self._write_new_data()
        
        beta_data = np.loadtxt('silicate0.00.txt', usecols = (1,2), skiprows = 1, delimiter = ' ')
        mass = self._new_data[:,indices['mass_index']]/1000
        
        
        for i in range(len(self._new_data)):
            idx = find_nearest_idx(beta_data[:,0], mass[i])
            [state, ltime] = spice.spkezr('ULYSSES',  self._new_data[i,self._time_index],      'ECLIPJ2000', 'NONE', 'SUN')
            beta = beta_data[idx,1]
            elts = spice.oscelt(state, self._new_data[i,self._time_index], (1-beta)*GM_km3_per_s2)
            peri.append(elts[0]/au_in_km)
        
        beta_plot = []
        lat_plot = []
        for i in range(len(self._new_data)):
            idx = find_nearest_idx(beta_data[:,0], mass[i])
            beta = beta_data[idx,1]
            beta_plot.append(beta)
            lat_plot.append(self._new_data[i,indices['LAT']])
                
        plt.scatter(lat_plot, beta_plot)
        plt.xlabel('Ulysses Latitude [°]')
        plt.ylabel(r'$\beta$')
        plt.show()
        
        plt.xlabel('Perihelion distance [au]')
        plt.ylabel('Count')
        
        plt.hist(peri)
        plt.show()

    def _write_new_data(self):
        new_data = []
        beta_dist = []
        eff_area_res = []
        
        for i in range(len(self.data_without_999)):
            if self._beta_meteoroids[i]:
                closest_eff_area_idx = find_nearest_idx(self.eff_area_time, self.data_without_999[i, self._time_index])
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
        
    def _set_new_angles_and_velocities(self):
        self._new_detector_sun_angles, self._new_velocities = self._calculate_angle_and_velocities_between_two_objects(self._new_data)
        
    def _compare_found_betas(self):
        print('Found beta indices:', self._new_data[:, self._index_index])
        print(len(self._new_data[:,0]),'beta particles found')
        print('Wehry beta indices:', self._wehry_beta_particle_indices)
        print(len(self._wehry_beta_particle_indices),'beta particles found')
        my_indices = np.array(self._new_data[:, self._index_index])
        wehry_indices = np.array(self._wehry_beta_particle_indices)
        intersect_indices = np.intersect1d(my_indices, wehry_indices)
        print('Intersection beta indices:', intersect_indices)
        print(len(intersect_indices),'beta particles found')

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
        zero_crossings = np.where(np.diff(np.signbit(self.raw_data[:,indices['LAT']])))[0]
        zero_crossings_times = self.raw_data[zero_crossings, self._time_index]
        self._zero_crossings_times = zero_crossings_times
        
    def _effective_area_with_lat(self):
        zero_crossings_times_in_epoch = self._zero_crossings_times/self.one_year_et+2000
        eff_area_data = np.loadtxt(self.eff_area_file, delimiter = ',')
        if self.eff_area_file == 'DefaultDataset.csv':
            eff_area_data[:,1] = eff_area_data[:,1]*10000
            
        #Check in which hemisphere Ulysses currently is, plot that, and calculate the average effective area for that bin
        #Does the same thing three times, for pre-Jupiter, for everything but the last semi-orbit, and for the last semi-orbit
        
        #Number 1
        cond = np.where(eff_area_data[:,0] < zero_crossings_times_in_epoch[0])
        cond = np.intersect1d(cond, cond)
        avg_eff_area = [np.mean(eff_area_data[cond,1])]
        plt.plot(eff_area_data[cond,0], eff_area_data[cond,1], label = 'Pre-Jupiter fly-by', color = 'black')
        
        #Number 2 to n-1
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
            cond0 = np.where(eff_area_data[:,0] > zero_crossings_times_in_epoch[i])
            cond1 = np.where(eff_area_data[:,0] < zero_crossings_times_in_epoch[i+1])
            cond = np.intersect1d(cond0, cond1)
            avg_eff_area.append(np.mean(eff_area_data[cond,1]))
            plt.plot(eff_area_data[cond,0], eff_area_data[cond,1], label = label, color = color)
            
        #Number n
        cond = np.where(eff_area_data[:,0] > zero_crossings_times_in_epoch[-1])
        cond = np.intersect1d(cond, cond)
        avg_eff_area.append(np.mean(eff_area_data[cond,1]))
        plt.plot(eff_area_data[cond,0], eff_area_data[cond,1], color = 'red')
        plt.xlabel('Year')
        plt.ylabel('Effective Area [cm$^2$]')
        plt.title(self.eff_area_file[:2]+'km/s')
        plt.legend()
        plt.show()
        
        self._avg_eff_area = np.array(avg_eff_area)
        
    
    def _count_and_print_north_fraction(self):
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

    def _plot_effective_area_with_streams(self):
            old_eff_area_data = np.loadtxt(self.eff_area_file, delimiter = ',')
            if self.eff_area_file == 'DefaultDataset.csv':
                old_eff_area_data[:,1] = old_eff_area_data[:,1]*10000
                
            x = np.linspace(old_eff_area_data[0,0], old_eff_area_data[-1,0], num = 100000)
            eff_area_data = np.array([x, np.interp(x, old_eff_area_data[:,0], old_eff_area_data[:,1])])
            
            streams = []
            for i in range(len(self._streams1)):
                streams.append(self._streams1[i])
            for i in range(len(self._streams2)):
                streams.append(self._streams2[i])
            
            #Check when Ulysses is within a Jupiter stream
            
            one_time_legend_flag = True
            no_stream_index_list = range(len(eff_area_data[0,:]))
            for interval in streams:
                cond0 = np.where(eff_area_data[0,:] >= interval[0])
                cond1 = np.where(eff_area_data[0,:] <= interval[1])
                cond = np.intersect1d(cond0, cond1)
                if one_time_legend_flag:
                    plt.plot(eff_area_data[0,cond], eff_area_data[1,cond], label = 'Jupiter stream', color = 'black')
                    one_time_legend_flag = False
                else:
                    plt.plot(eff_area_data[0,cond], eff_area_data[1,cond], color = 'black')
                no_stream_index_list = list(set(no_stream_index_list) - set(cond))
            plt.plot(eff_area_data[0,no_stream_index_list], eff_area_data[1,no_stream_index_list], label = 'No Jupiter stream', color = 'green', alpha = 0.1)
            plt.xlabel('Year')
            plt.ylabel('Effective Area [cm$^2$]')
            plt.title(self.eff_area_file[:2]+'km/s')
            plt.legend()
            plt.show()

    
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
        
        
    #Bins over pre-fly-by, north sections an south sections and calculates fluxes and mean effective area of bin
    def _calculate_mean_eff_area_and_flux(self) -> (list, list):
        
        """
        Add number of detected particles
        """
        
        flux_bins = []
        mean_dist_array = []
        
        #Does the same thing three times, for pre-Jupiter, for everything but the last semi-orbit, and for the last semi-orbit
        
        #Number 1
        dist_array = []
        
        for i in range(len(self._new_data)):
            if self._new_data[i, self._time_index]<=self._zero_crossings_times[0]:
                dist_array.append(self._beta_dist[i])
        dist_array = np.array(dist_array)
        mean_dist_array.append(np.mean(dist_array))
        flux_bins.append(len(dist_array)/(self._avg_eff_area[0]*(self._zero_crossings_times[1]-self.raw_data[0, self._time_index])))
        
        #Number 2 to n-1
        for k in range(len(self._zero_crossings_times)-1):
            dist_array = []
            
            idx = find_nearest_idx(self._new_data[:, self._time_index], self._zero_crossings_times[k])+1
            
            for i in range(idx, len(self._new_data)):
                if self._new_data[i, self._time_index]<=self._zero_crossings_times[k+1]:
                    dist_array.append(self._beta_dist[i])
            dist_array = np.array(dist_array)
            mean_dist_array.append(np.mean(dist_array))
            flux_bins.append(len(dist_array)/(self._avg_eff_area[k+1]*(self._zero_crossings_times[k+1]-self._zero_crossings_times[k])))
            
        #Number n 
        dist_array = []
        
        idx = find_nearest_idx(self._new_data[:, self._time_index], self._zero_crossings_times[-1])
        
        for i in range(idx, len(self._new_data)):
            dist_array.append(self._beta_dist[i])
        dist_array = np.array(dist_array)
        mean_dist_array.append(np.mean(dist_array))
        flux_bins.append(len(dist_array)/(self._avg_eff_area[-1]*(self.raw_data[-1, self._time_index]-self._zero_crossings_times[-1])))
            
            
        mean_dist_array = np.array(mean_dist_array)
        flux_bins = np.array(flux_bins)
        
        
        return mean_dist_array, flux_bins
    

    def __init__(self):       
        self._time_index = 0
        self._dist_index = 1
        self._index_index = 2
        
        self._wehry_velocity = 20
        self._wehry_angle = 50
        self._streams1 = [[1991.727,1991.740],[1991.948,1991.953],[1991.978,1991.983],[1992.017,1992.021],[1992.048,1992.056],[1992.192,1992.196],[1992.268,1992.275],[1992.342,1992.349],[1992.417,1992.435],[1992.670,1992.685],[1992.796,1992.811]]
        self._streams2 = [[2002.905,2002.917],[2003.515,2003.535],[2003.643,2003.663],[2003.717,2003.742],[2003.778,2003.802],[2003.862,2003.867],[2003.919,2003.928],[2003.993,2004.007],[2004.063,2004.078],[2004.138,2004.142],[2004.207,2004.231],[2004.410,2004.442],[2004.450,2004.510],[2004.544,2004.574],[2004.582,2004.603],[2004.625,2004.654],[2004.825,2004.833],[2004.907,2004.914],[2004.989,2005.001],[2005.113,2005.135],[2005.210,2005.236],[2005.474,2005.487],[2005.565,2005.583],[2005.615,2005.642]]
        self._instrument = [ [2000.4905, 2000.49625], [2002.23297, 2002.2700], [2002.917692,2003.422021], [2003.5733, 2003.642642], [2004.918349, 2004.923907]]
        self._interstellar_ecliptic_lon = 252
        self._interstellar_ecliptic_lat = 2.5
        """
        try different tolerances
        """
        self._tolerance = 70
        self._interstellar_min_vel = 14
        self._interstellar_max_vel = 38
        self._interstellar_min_mass = 2.5e-14
        self._interstellar_max_mass = 1e-12
        self._wehry_beta_particle_indices = [4,6,8,18,19,29,31,34,35,36,38,43,48,49,59,61,66,72,74,76,80,83,86,94,1032,1080,1145,1165,1410,1412,1421,1422,1427,1428,1429,1431,1433,1436,1438,1440,1442,1449,1450,1452,1455,1465,1984,2001,2003,2010,2012,2024,2034,2035,2048,2049,2051,2052,2053,2054,2055,2060]

        