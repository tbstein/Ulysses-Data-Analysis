import numpy as np
from indices_dict import indices
import spiceypy as spice
import matplotlib.pyplot as plt

spice_path = '../../spice/'
spice.furnsh(spice_path + "naif0012.tls")
spice.furnsh(spice_path + "de440.bsp")
spice.furnsh(spice_path + "ulysses_1990_2009_2050.bsp")

first_data_column = 3
last_data_column = 34

used_cols = (i for i in range(first_data_column, last_data_column))

time = []
with open('Ulysses_Data_File_Cleaned.txt') as cleaned_ulysses_data:
    for count, line in enumerate(cleaned_ulysses_data):
        if count >= indices['first_data_line']:
            line = line.split()
            time.append(line[indices['date']]+'T'+line[indices['time_of_day']])

time = spice.str2et(time)

ulysses_data = np.loadtxt('Ulysses_Data_File_Cleaned.txt', delimiter = ' ', skiprows = indices['first_data_line'], usecols = used_cols)

ulysses_data = np.concatenate((np.transpose(np.array([time])), ulysses_data), axis = 1)


wehry_beta_particle_indices = [4,6,8,18,19,29,31,34,35,36,38,43,48,49,59,61,66,72,74,76,80,83,86,94,1032,1080,1145,1165,1410,1412,1421,1422,1427,1428,1429,1431,1433,1436,1438,1440,1442,1449,1450,1452,1455,1465,1984,2001,2003,2010,2012,2024,2034,2035,2048,2049,2051,2052,2053,2054,2055,2060]

wehry_beta_particle_indices = np.array(wehry_beta_particle_indices)-1

wehry_data = ulysses_data[wehry_beta_particle_indices,:]


data = wehry_data
time_index = 0
wehry_velocity = 20
wehry_angle = 50

def calculate_angle_and_velocities_between_two_objects(data: list) -> (list, list):
    lon = data[:, indices['solar_lon_index']]
    lat = data[:, indices['solar_lat_index']]
    detector_sun_angles = []
    velDust = []
    et = data[:, time_index]
    for i in range(len(data[:,0])):
        [stateSun, ltime] = spice.spkezr('SUN',  et[i],      'ECLIPJ2000', 'NONE', 'ULYSSES')
        posSun = stateSun[:3]
        velSun = stateSun[3:]
        velRelative = data[i, indices['velocity_index']]
        posDetector = spice.latrec(1, lon[i]*2*np.pi/360, lat[i]*2*np.pi/360)
        vel = np.linalg.norm(velSun+velRelative*posDetector)
        detector_sun_angles.append(spice.vsep(posDetector, posSun)*360/(2*np.pi))
        velDust.append(vel)
    return detector_sun_angles, velDust
    

def plot_wehry(angles: list, velocities: list):       
    plt.xlabel('Angle between detector axis and sun [°]')
    plt.ylabel('Absolute Particle Velocity [km/s]')        
    plt.scatter(angles,velocities)
    plt.plot(detector_sun_angles, np.ones(len(detector_sun_angles))*wehry_velocity, color = 'red')
    plt.plot(np.ones(len(detector_sun_angles))*wehry_angle, velocities, color = 'red')
    plt.show() 

detector_sun_angles, velocities = calculate_angle_and_velocities_between_two_objects(data)

plot_wehry(detector_sun_angles, velocities)

