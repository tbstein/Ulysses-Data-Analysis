import numpy as np
import matplotlib.pyplot as plt

test0 = np.loadtxt('wehryTest.dat', delimiter = ',')
test1 = np.loadtxt('angles_to_sun_wehry_particles.txt', delimiter = ',', skiprows = 1, usecols = (i for i in range(9, 11)))

anglesStrub = test1[:,0]
anglesAlt = test1[:,1]

print(test0-anglesStrub)

plt.scatter(test0, anglesStrub-test0, label = 'test0')
#plt.scatter(anglesStrub, anglesStrub-test0, label = 'Strub')
#plt.legend()
plt.show()

a = np.where(test0-anglesStrub < -1)[:]
print(a, test0[a], anglesStrub[a], test0[a]-anglesStrub[a])