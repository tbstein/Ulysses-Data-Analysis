import numpy as np
from indices_dict import indices
import spiceypy as spice
from CleanedDataPlotter import CleanedDataPlotter
import math

def find_closest_point_on_plane(antennaPointing, point_on_plane, vector_to_project):
    # Unpack the normal vector (antennaPointing)
    a, b, c = antennaPointing

    # Unpack the coordinates of a point on the plane
    x0, y0, z0 = point_on_plane

    # Unpack the vector to be projected onto the plane
    x_vec, y_vec, z_vec = vector_to_project

    # Calculate the scalar projection of the vector onto the normal vector of the plane
    scalar_projection = (a * x_vec + b * y_vec + c * z_vec) / (a**2 + b**2 + c**2)

    # Calculate the coordinates of the closest point on the plane
    closest_point_x = x_vec - scalar_projection * a + x0
    closest_point_y = y_vec - scalar_projection * b + y0
    closest_point_z = z_vec - scalar_projection * c + z0
    
    return closest_point_x, closest_point_y, closest_point_z

def rotate_point_around_axis(closest_point, angle_R, axis):
    # Convert the angle from degrees to radians
    angle_R_rad = math.radians(angle_R)

    # Normalize the axis vector to make sure it is a unit vector
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)

    # Convert the closest_point and axis to numpy arrays for easier calculations
    closest_point = np.array(closest_point)
    
    # Perform the Rodrigues' rotation formula
    rotated_vector = closest_point * math.cos(angle_R_rad) + np.cross(axis, closest_point) * math.sin(angle_R_rad) + axis * np.dot(axis, closest_point) * (1 - math.cos(angle_R_rad))

    return rotated_vector.tolist()

def rotate_vector_away_from_normal(vector, angle_degrees, normal_vector):
    # Convert the angle from degrees to radians
    angle_rad = math.radians(angle_degrees)

    # Normalize the vectors to make sure they are unit vectors
    vector = np.array(vector)
    vector = vector / np.linalg.norm(vector)

    normal_vector = np.array(normal_vector)
    normal_vector = normal_vector / np.linalg.norm(normal_vector)

    # Calculate the axis of rotation (perpendicular to the normal vector and the vector)
    axis = np.cross(normal_vector, vector)
    axis = axis / np.linalg.norm(axis)

    # Perform the Rodrigues' rotation formula with the modified angle
    rotated_vector = vector * math.cos(angle_rad) + np.cross(axis, vector) * math.sin(angle_rad) + axis * np.dot(axis, vector) * (1 - math.cos(angle_rad))

    return rotated_vector.tolist()

def calculate_detector_pointing(data: list) -> (list, list):
    rot = data[:, indices['rotation_angle_index']]
    time_index = 0
    et = data[:, time_index]
    s_lon = []
    s_lat = []
    for i in range(len(data[:,0])):
        [stateUlysses, ltime] = spice.spkezr('ULYSSES',  et[i],      'ECLIPJ2000', 'NONE', 'SUN')
        [stateEarth, ltime] = spice.spkezr('EARTH BARYCENTER',  et[i],      'ECLIPJ2000', 'NONE', 'SUN')
        posEarth = stateEarth[:3]
        posUlysses = stateUlysses[:3]
        
        posEarth = np.array(posEarth)
        posUlysses = np.array(posUlysses)
        
        antennaPointing = posEarth-posUlysses
        
        point_on_plane = posUlysses
        
        vector_to_project = (0, 0, 1)  # The vector to be projected onto the plane
        closest_point = find_closest_point_on_plane(antennaPointing, point_on_plane, vector_to_project)
        rotated_point = rotate_point_around_axis(closest_point, rot[i], antennaPointing)
        rotated_5_degrees = rotate_vector_away_from_normal(rotated_point, 5, -antennaPointing)
        
        detectorPointing = spice.reclat(rotated_5_degrees)
        s_lon.append(-detectorPointing[1])
        s_lat.append(-detectorPointing[2])
    
    return s_lon, s_lat

spice_path = '../../spice/'
spice.furnsh(spice_path + "naif0012.tls")
spice.furnsh(spice_path + "de440.bsp")
spice.furnsh(spice_path + "ulysses_1990_2009_2050.bsp")
one_year_et = spice.str2et('01-001T00:00')-spice.str2et('00-001T00:00')

first_data_column = 3
last_data_column = 34

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

ulysses_data = np.loadtxt('Ulysses_Data_File_Cleaned.txt', delimiter = ' ', skiprows = indices['first_data_line'], usecols = used_cols)

ulysses_data = np.concatenate((np.transpose(np.array([index])), ulysses_data), axis = 1)

ulysses_data = np.concatenate((np.transpose(np.array([time])), ulysses_data), axis = 1)

s_lon, s_lat = calculate_detector_pointing(ulysses_data)

cmpr_s_lon = ulysses_data[:, indices['solar_lon_index']]
cmpr_s_lat = ulysses_data[:, indices['solar_lat_index']]

print((s_lon-cmpr_s_lon)%90, (s_lat - cmpr_s_lat)%90)