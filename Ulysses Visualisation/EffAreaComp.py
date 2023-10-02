import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
from MiscFunctions import find_nearest_idx
from indices_dict import indices

#Units of spice are km and km/s

#Takes angles in degrees, not radians
def uniform(angle: float, interval: tuple = (-30,30)) -> float:
    if interval[0] <= angle and angle <= interval[1]:
        return 1/(interval[1]-interval[0])
    else:
        return 0

def avg_eff_area(eff_area: list, eff_area_angle: list, extra_angle: float) -> list:
    avg = []
    N = 360
    two_pi = np.linspace(0, 2*np.pi, num = N) #Rotation angle
    #detector_earth_angles = np.linspace(0, 95, num = 95)
    #detector_earth_angles = np.concatenate((detector_earth_angles, np.flip(detector_earth_angles)))
    detector_meteoroid_angles = np.linspace(0, extra_angle, num = extra_angle)
    detector_meteoroid_angles = np.concatenate((detector_meteoroid_angles, np.flip(detector_meteoroid_angles)))
    for i, phi in enumerate(detector_meteoroid_angles):
        A = 0
        for n, theta in enumerate(two_pi):
            dot_product = np.cos(phi/180*np.pi)+np.sin(phi/180*np.pi)*np.tan(extra_angle/180*np.pi)*np.cos(theta)
            norm = np.sqrt(1+np.tan(extra_angle/180*np.pi)**2)
            angle = np.arccos(dot_product/norm)/np.pi*180
            #index = find_nearest_idx(eff_area_angle, -2*detector_earth_angles[i]*np.cos(theta)+extra_angle+detector_earth_angles[i])
            index = find_nearest_idx(eff_area_angle, angle)
            #A += eff_area[index]/N*2*np.pi
            A += eff_area[index]/N
        avg.append(A)
    avg = np.array(avg)
    return avg

spice_path = '../../spice/'
spice.furnsh(spice_path + "naif0012.tls")
spice.furnsh(spice_path + "de440.bsp")
spice.furnsh(spice_path + "ulysses_1990_2009_2050.bsp")
one_year_et = spice.str2et('01-001T00:00')-spice.str2et('00-001T00:00')

start_et = spice.str2et('90-301T13:53')
end_et = spice.str2et('07-334T08:17')

num = 10000

index = []
time = []
lon = []
lat = []
with open('Ulysses_Data_File_Cleaned.txt') as cleaned_ulysses_data:
    for count, line in enumerate(cleaned_ulysses_data):
        if count >= indices['first_data_line']:
            line = line.split()
            time.append(line[indices['date']]+'T'+line[indices['time_of_day']])
            lon.append(float(line[indices['solar_lon_index']]))
            lat.append(float(line[indices['solar_lat_index']]))

lon = np.array(lon)
lat = np.array(lat)

time = spice.str2et(time)

delta_t = (end_et-start_et)/num

max_area = 1000
extra_angle = 95

sensitivityRAW = np.loadtxt('sensitivity_function.txt', delimiter=',')
sensitivity_angle = np.concatenate((-np.flip(sensitivityRAW[1:,0]),sensitivityRAW[:,0]))
sensitivity_value = np.concatenate((np.flip(sensitivityRAW[1:,1]),sensitivityRAW[:,1]))

sensitivity = np.array([sensitivity_angle, sensitivity_value])


eff_area_array = []

boundaries = [10,24,38]

for i in boundaries:

    uniformDist = []
    for angle in sensitivity[0]:
        uniformDist.append(uniform(angle, interval = (-i,i)))
    uniformDist = np.array(uniformDist)
    
    temp = avg_eff_area(sensitivity[1], sensitivity[0], extra_angle-10)
    
    
    convolution = np.convolve(uniformDist, temp)
    conv_plt_angle = np.linspace(-len(convolution)/2,len(convolution)/2,len(convolution))
    
    wehry = np.loadtxt('DefaultDataset.csv', delimiter = ',')
    
    velocity_dust = 30
    
    angle = []
    factor = []
    for et in time:
        [stateUlysses, ltime] = spice.spkezr('ULYSSES',  et,      'ECLIPJ2000', 'NONE', 'SUN')
        [stateEarth, ltime] = spice.spkezr('EARTH BARYCENTER',  et,      'ECLIPJ2000', 'NONE', 'SUN')
        velUlysses = stateUlysses[3:]
        posEarth = stateEarth[:3]
        posUlysses = stateUlysses[:3]
        dustVel = velocity_dust*posUlysses/(np.linalg.norm(posUlysses))
        relVel = velUlysses-dustVel
        pointingDetector = spice.latrec(1, lon[i]*2*np.pi/360, lat[i]*2*np.pi/360)
        factor.append(np.abs(np.linalg.norm(relVel)/np.linalg.norm(dustVel)))
        
        angle0 = spice.vsep(-velUlysses+dustVel, posEarth-posUlysses)
    
        angle.append(angle0)
        
        
    factor = np.array(factor)
    
    angle = np.array(angle)
    
    pltangle = angle*360/(2*np.pi)
    plttime = time/one_year_et+2000
    
    
    
    wehry_sensitivity = np.loadtxt('wehry_sensitivity.csv', delimiter = ',')
    
    eff_area = []
    for i in range(len(pltangle)):
        index = find_nearest_idx(conv_plt_angle+extra_angle, pltangle[i])
        
        eff_area.append(convolution[index])
    
    eff_area = np.array(eff_area)*factor
    eff_area_array.append(eff_area)



plt.xlabel('Time')
plt.ylabel(r'Velocity corrected effective area [cm$^2$]')
for i in range(len(boundaries)):
    plt.plot(plttime, eff_area_array[i], label = str(boundaries[i])+'°')
plt.legend()
plt.savefig('effareacompdegrees.pdf')
plt.show()

    
    
