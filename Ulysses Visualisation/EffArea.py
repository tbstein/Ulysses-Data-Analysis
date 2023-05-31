# %%

import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
from MiscFunctions import find_nearest_idx

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
    rot_angle_equivalent = np.linspace(1, 0, num = N)
    detector_earth_angles = np.linspace(0, 95, num = 95)
    detector_earth_angles = np.concatenate((detector_earth_angles, np.flip(detector_earth_angles)))
    for i in range(len(detector_earth_angles)):
        A = 0
        for n, theta in enumerate(two_pi):
        #for n, value in enumerate(rot_angle_equivalent):
            index = find_nearest_idx(eff_area_angle, -2*detector_earth_angles[i]*np.cos(theta)+extra_angle+detector_earth_angles[i])
            #index = find_nearest_idx(eff_area_angle, -2*detector_earth_angles[i]*value+extra_angle+detector_earth_angles[i])
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
time = np.linspace(start_et, end_et, num=num)
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
plt.xlabel('Angle [°]')
plt.title('Sensitivity function')
plt.plot(sensitivity[0,:], sensitivity[1,:])
plt.show()

uniformDist = []
for angle in sensitivity[0]:
    #uniformDist.append(uniform(angle-extra_angle))
    uniformDist.append(uniform(angle))
uniformDist = np.array(uniformDist)

plt.xlabel("Angle [°]")
plt.ylabel("Normalized Value")
plt.title("Uniform distribution")
plt.plot(sensitivity[0], uniformDist)
plt.show()

temp = avg_eff_area(sensitivity[1], sensitivity[0], extra_angle)
plt.title('Rotation averaged sensitivity function')
plt.ylabel(r'Area [cm$^2$]')
plt.xlabel('Angle [°]')
plt.plot(np.linspace(-95, 95, num = len(temp)), temp)
plt.show()


convolution = np.convolve(uniformDist, temp)
conv_plt_angle = np.linspace(-len(convolution)/2,len(convolution)/2,len(convolution))
plt.ylabel(r'Area [cm$^2$]')
plt.xlabel('Angle [°]')
plt.title('Convolved sensitivity function')
#plt.plot(conv_plt_angle+extra_angle, convolution)
plt.plot(conv_plt_angle, convolution)
plt.show()

#convolution = sensitivity[1]
#conv_plt_angle = sensitivity[0]

wehry = np.loadtxt('DefaultDataset.csv', delimiter = ',')

velocity_dust = 20
for i in range(6):

    velocity_dust = 20+10*i
    
    angle = []
    factor = []
    for et in time:
        [stateSun, ltime] = spice.spkezr('SUN',  et,      'J2000', 'NONE', 'ULYSSES')
        [stateEarth, ltime] = spice.spkezr('EARTH BARYCENTER',  et,      'J2000', 'NONE', 'ULYSSES')
        velSun = stateSun[3:]
        posEarth = stateEarth[:3]
        posSun = stateSun[:3]
        dustVel = velocity_dust*posSun/(np.linalg.norm(posSun))
        relVel = velSun-dustVel
        factor.append(np.abs(np.linalg.norm(relVel)/np.linalg.norm(dustVel)))
        angle.append(spice.vsep(velSun+dustVel, posEarth))
        
    factor = np.array(factor)
    
    angle = np.array(angle)
    
    pltangle = angle*360/(2*np.pi)
    plttime = time/one_year_et+2000
    
    plt.title('Angle between Earth-Ulysses Position Vector and Sun-Ulysses Velocity Vector + ' + str(velocity_dust) + 'km/s from Sun direction')
    plt.ylabel('Angle [°]')
    plt.xlabel('Time')
    plt.plot(plttime, pltangle)
    plt.show()
    
    
    eff_area = []
    for i in range(len(pltangle)):
        #index = find_nearest_idx(conv_plt_angle+extra_angle, pltangle[i])
        index = find_nearest_idx(conv_plt_angle, extra_angle-pltangle[i])
        eff_area.append(convolution[index])
    
    
    plt.title(str(velocity_dust) + 'km/s')
    plt.xlabel('Time')
    plt.ylabel(r'Effective area [cm$^2$]')
    plt.plot(plttime, eff_area, label = 'Own estimate', color = 'red')
    plt.plot(wehry[:,0], wehry[:,1]*10000, label = 'Wehry', color = 'blue')
    plt.legend()
    plt.show()
    
    eff_area = np.array(eff_area)*factor



    plt.title(str(velocity_dust) + 'km/s')
    plt.xlabel('Time')
    plt.ylabel(r'Velocity corrected effective area [cm$^2$]')
    plt.plot(plttime, eff_area, label = 'Own estimate', color = 'red')
    plt.plot(wehry[:,0], wehry[:,1]*10000, label = 'Wehry', color = 'blue')
    plt.legend()
    plt.show()
    
    with open(str(velocity_dust)+'.dat', 'w') as f:
        for i in range(len(eff_area)):
            f.write(str(plttime[i]) + ', ' + str(eff_area[i]) + '\n')
    
    plt.title(str(velocity_dust) + 'km/s')
    plt.plot(plttime, factor)
    plt.xlabel('Time')
    plt.ylabel('Velocity factor')
    plt.show()