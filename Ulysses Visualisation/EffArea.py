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

"""
def avg_eff_area(eff_area: list, eff_area_angle: list, extra_angle: float) -> list:
    avg = []
    N = 360*10
    two_pi = np.linspace(0, 2*np.pi, num = N) #Rotation angle
    detector_meteoroid_angles = np.linspace(-85, 95, num = 181)
    for i, phi in enumerate(detector_meteoroid_angles):
        A = 0
        for n, theta in enumerate(two_pi):
            #dot_product = np.cos(phi/180*np.pi)+np.sin(phi/180*np.pi)*np.tan(extra_angle/180*np.pi)*np.cos(theta)
            dot_product = np.sin((phi+85)/180*np.pi)*np.tan(85/180*np.pi)*np.cos(theta)+np.cos((phi+85)/180*np.pi)
            norm = np.sqrt(1+np.tan(85/180*np.pi)**2)
            angle = np.arccos(dot_product/norm)/np.pi*180
            #index = find_nearest_idx(eff_area_angle, -2*detector_earth_angles[i]*np.cos(theta)+extra_angle+detector_earth_angles[i])
            index = find_nearest_idx(eff_area_angle, angle)
            #A += eff_area[index]/N*2*np.pi
            A += eff_area[index]/N
        avg.append(A)
    avg = np.array(avg)
    return avg
"""


def avg_eff_area(eff_area: list, eff_area_angle: list, extra_angle: float) -> list:
    avg = []
    N = 360*10
    two_pi = np.linspace(0, 2*np.pi, num = N) #Rotation angle
    detector_meteoroid_angles = np.linspace(-85, 95, num = 181)
    for i, phi in enumerate(detector_meteoroid_angles):
        A = 0
        for n, theta in enumerate(two_pi):
            #dot_product = np.cos(phi/180*np.pi)+np.sin(phi/180*np.pi)*np.tan(extra_angle/180*np.pi)*np.cos(theta)
            dot_product = np.sin((phi+85)/180*np.pi)*np.tan(85/180*np.pi)*np.cos(theta)-np.cos((phi+85)/180*np.pi)
            norm = np.sqrt(1+np.tan(85/180*np.pi)**2)
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
#time = np.linspace(start_et, end_et, num=num)

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
#sensitivity_angle += extra_angle

#sensitivity_value = sensitivity_value/max_area*2*np.pi*np.sin(sensitivity_angle/360*(2*np.pi))

sensitivity = np.array([sensitivity_angle, sensitivity_value])


plt.ylabel(r'Area [cm$^2$]')
plt.xlabel('Angle between impactor and detector pointing [°]')
plt.title('Sensitivity function')
plt.plot(sensitivity[0,:], sensitivity[1,:])
plt.savefig('sensitivity.pdf')
plt.show()

uniformDist = []
for angle in sensitivity[0]:
    #uniformDist.append(uniform(angle-extra_angle))
    uniformDist.append(uniform(angle, interval = (-24,24)))
uniformDist = np.array(uniformDist)

plt.xlabel("Angle between impactor and detector pointing [°]")
plt.ylabel("Normalized Value")
plt.title("Uniform distribution")
plt.plot(sensitivity[0], uniformDist)
plt.savefig('uniform.pdf')
plt.show()

temp = avg_eff_area(sensitivity[1], sensitivity[0], extra_angle-10)
plt.title('Rotation averaged sensitivity function')
plt.ylabel(r'Area [cm$^2$]')
plt.xlabel('Angle between impactor and detector cone [°]')
#plt.plot(np.linspace(-85, 95, num = len(temp)), temp)
plt.plot(np.linspace(-85, 95, num = len(temp)), temp)
plt.savefig('rotationavg.pdf')
plt.show()


convolution = np.convolve(temp, uniformDist, mode='same')
conv_plt_angle = np.linspace(-len(convolution)/2,len(convolution)/2,len(convolution))
conv_plt_angle = np.linspace(-85, 95, num = len(temp))
plt.ylabel(r'Area [cm$^2$]')
plt.xlabel('Angle between impactor and detector cone [°]')
plt.title('Convolved sensitivity function')
#plt.plot(conv_plt_angle+extra_angle, convolution)
plt.plot(conv_plt_angle, convolution)
plt.savefig('convolved.pdf')
plt.show()

#convolution = sensitivity[1]
#conv_plt_angle = sensitivity[0]


wehry = np.loadtxt('DefaultDataset.csv', delimiter = ',')

velocity_dust = 30
for i in range(1):

    velocity_dust = 30+10*i
    
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
        #angle.append(spice.vsep(velUlysses+dustVel, posEarth-posUlysses))
        #angle.append(spice.vsep(-velUlysses+dustVel, pointingDetector))
        
        #angle0 = spice.vsep(-velUlysses+dustVel, pointingDetector)
        #angle1 = spice.vsep(-velUlysses+dustVel, -pointingDetector)
        
        '''
        angle0 = spice.vsep(-velUlysses+dustVel, posEarth-posUlysses)
        angle1 = spice.vsep(-velUlysses+dustVel, -posEarth+posUlysses)
        angle1 = angle0
        
        if np.abs(angle0-extra_angle) <= np.abs(angle1-extra_angle):
            angle.append(angle0)
        else:
            angle.append(angle1)
        '''
        #angle.append(95/180*np.pi-spice.vsep(-velUlysses+dustVel, posEarth-posUlysses))
        angle.append(-85/180*np.pi+spice.vsep(velUlysses-dustVel, posEarth-posUlysses))
        
        
    factor = np.array(factor)
    
    angle = np.array(angle)
    
    pltangle = angle*360/(2*np.pi)
    plttime = time/one_year_et+2000
    
    
    #plt.title('Angle between Earth-Ulysses Position Vector and Sun-Ulysses Velocity Vector + ' + str(velocity_dust) + 'km/s from Sun direction')
    plt.title('Angle between Detector cone and beta meteoroid flux')
    plt.ylabel('Angle [°]')
    plt.xlabel('Time [years]')
    plt.plot(plttime, pltangle)
    plt.savefig('angles.pdf')
    plt.show()
    
    
    
    wehry_sensitivity = np.loadtxt('wehry_sensitivity.csv', delimiter = ',')
    
    eff_area = []
    for i in range(len(pltangle)):
        index = find_nearest_idx(conv_plt_angle, pltangle[i])
        eff_area.append(convolution[index])
        #index = find_nearest_idx(np.linspace(-95, 85, num = len(temp)), pltangle[i])
        #eff_area.append(temp[index])
    
    plt.title(str(velocity_dust) + 'km/s')
    plt.xlabel('Time [years]')
    plt.ylabel(r'Effective area [cm$^2$]')
    plt.plot(plttime, eff_area, label = 'Own estimate', color = 'red', ls = 'solid')
    plt.plot(wehry[:,0], wehry[:,1]*10000, label = 'Wehry', color = 'blue', ls = 'dashed')
    plt.legend()
    plt.savefig('raweffarea.pdf')
    plt.show()
    
    eff_area = np.array(eff_area)*factor



    plt.title(str(velocity_dust) + 'km/s')
    plt.xlabel('Time [years]')
    plt.ylabel(r'Velocity corrected effective area [cm$^2$]')
    plt.plot(plttime, eff_area, label = 'Own estimate', color = 'red', ls = 'solid')
    plt.plot(wehry[:,0], wehry[:,1]*10000, label = 'Wehry', color = 'blue', ls = 'dashed')
    plt.legend()
    plt.savefig('correffarea30.pdf')
    plt.show()
    
    strub = np.loadtxt('uly_beta_sensareaearly.txt', delimiter = ',')
    
    plt.title(str(velocity_dust) + 'km/s')
    plt.xlabel('Time [years]')
    plt.ylabel(r'Velocity corrected effective area [cm$^2$]')
    plt.plot(plttime, eff_area, label = 'Own estimate', color = 'red', ls = 'solid')
    plt.plot(strub[:,0], strub[:,1], label = 'Strub', color = 'blue', ls = 'dashed')
    plt.legend()
    plt.savefig('VergleichPeterTobias.pdf')
    plt.show()
    
    plt.title(str(velocity_dust) + 'km/s')
    plt.xlabel('Time [years]')
    plt.ylabel(r'Velocity corrected effective area [cm$^2$]')
    plt.plot(plttime, eff_area, label = 'Own estimate', color = 'red', ls = 'solid')
    plt.plot(strub[:,0], strub[:,1], label = 'Strub', color = 'blue', ls = 'dashed')
    plt.plot(wehry[:,0], wehry[:,1]*10000, label = 'Wehry', color = 'green', ls = '-.')
    plt.legend()
    plt.savefig('VergleichWehry.pdf')
    plt.show()
    
    
    
    
    with open(str(velocity_dust)+'.dat', 'w') as f:
        for i in range(len(eff_area)):
            f.write(str(plttime[i]) + ', ' + str(eff_area[i]) + '\n')


    plt.title(str(velocity_dust) + 'km/s')
    plt.plot(plttime, factor)
    plt.xlabel('Time [years]')
    plt.ylabel('Velocity factor')
    plt.savefig('factor.pdf')
    plt.show()
    
    
