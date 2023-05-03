import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice

#Units of spice are km and km/s

spice_path = '../../spice/'
spice.furnsh(spice_path + "naif0012.tls")
spice.furnsh(spice_path + "de440.bsp")
spice.furnsh(spice_path + "ulysses_1990_2009_2050.bsp")
one_year_et = spice.str2et('01-001T00:00')-spice.str2et('00-001T00:00')

start_et = spice.str2et('90-301T13:53')
end_et = spice.str2et('07-334T08:17')

num = 10000
time = np.linspace(start_et, end_et, num=num)

max_area = 1000
extra_angle = 95

sensitivityRAW = np.loadtxt('sensitivity_function.txt', delimiter=',')
sensitivity_angle = np.concatenate((-np.flip(sensitivityRAW[1:,0]),sensitivityRAW[:,0]))
sensitivity_value = np.concatenate((np.flip(sensitivityRAW[1:,1]),sensitivityRAW[:,1]))
sensitivity_angle += extra_angle

#sensitivity_value = sensitivity_value/max_area*2*np.pi*np.sin(sensitivity_angle/360*(2*np.pi))

sensitivity = np.array([sensitivity_angle, sensitivity_value])

plt.ylabel(r'Area [cm$^2$]')
plt.xlabel('Angle [°]')
plt.title('Sensitivity function')
plt.plot(sensitivity[0,:], sensitivity[1,:])
plt.show()

#https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#Takes angles in degrees, not radians
def uniform(angle, interval = (-30,30)):
    if interval[0] <= angle and angle <= interval[1]:
        return 1/(interval[1]-interval[0])
    else:
        return 0

uniformDist = []
for angle in sensitivity[0]:
    uniformDist.append(uniform(angle-extra_angle))
uniformDist = np.array(uniformDist)

plt.xlabel("Angle [°]")
plt.ylabel("Normalized Value")
plt.title("Uniform distribution")
plt.plot(sensitivity[0], uniformDist)
plt.show()


convolution = np.convolve(uniformDist, sensitivity[1])
conv_plt_angle = np.linspace(-len(convolution)/2,len(convolution)/2,len(convolution))
plt.ylabel(r'Area [cm$^2$]')
plt.xlabel('Angle [°]')
plt.title('Convolved sensitivity function')
plt.plot(conv_plt_angle+extra_angle, convolution)
plt.show()


velocity_dust = 20
for i in range(6):

    velocity_dust = 20+10*i
    
    angle = []
    for et in time:
        [stateSun, ltime] = spice.spkezr('SUN',  et,      'J2000', 'NONE', 'ULYSSES')
        [stateEarth, ltime] = spice.spkezr('EARTH',  et,      'J2000', 'NONE', 'ULYSSES')
        velSun = stateSun[3:]
        posEarth = stateEarth[:3]
        posSun = stateSun[:3]
        dustVel = velocity_dust*posSun/(np.linalg.norm(posSun))
        angle.append(spice.vsep(velSun+dustVel, posEarth))
    
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
        index = find_nearest_idx(conv_plt_angle+extra_angle, pltangle[i])
        eff_area.append(convolution[index])

    plt.title(str(velocity_dust) + 'km/s')
    plt.xlabel('Time')
    plt.ylabel(r'Effective area [cm$^2$]')
    plt.plot(plttime, eff_area)
    plt.show()
    
    with open(str(velocity_dust)+'.dat', 'w') as f:
        for i in range(len(eff_area)):
            f.write(str(plttime[i]) + ', ' + str(eff_area[i]) + '\n')