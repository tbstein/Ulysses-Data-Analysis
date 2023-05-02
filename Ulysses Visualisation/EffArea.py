import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice

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

angle = []
for et in time:
    [state, ltime] = spice.spkezr('SUN',  et,      'J2000', 'NONE', 'ULYSSES')
    pos = state[:3]
    vel = state[3:]
    angle.append(spice.vsep(pos, vel))

angle = np.array(angle)

pltangle = angle*360/(2*np.pi)
plttime = time/one_year_et+2000

plt.title('Angle between Sun-Ulysses Vector and Ulysses Velocity Vector')
plt.ylabel('Angle (°)')
plt.xlabel('Time')
plt.plot(plttime, pltangle)
plt.show()

sensitivityRAW = np.loadtxt('sensitivity_function.txt', delimiter=',')
print(sensitivityRAW)
sensitivity_angle = np.concatenate((-np.flip(sensitivityRAW[1:,0]),sensitivityRAW[:,0]))
sensitivity_value = np.concatenate((np.flip(sensitivityRAW[1:,1]),sensitivityRAW[:,1]))
sensitivity_angle += 95

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

eff_area = []
for i in range(len(pltangle)):
    index = find_nearest_idx(sensitivity[0,:], pltangle[i])
    eff_area.append(sensitivity[1, index])

plt.xlabel('Time')
plt.ylabel(r'Effective area [cm$^2$]')
plt.plot(plttime, eff_area)