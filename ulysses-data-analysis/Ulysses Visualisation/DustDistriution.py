import numpy as np
import matplotlib.pyplot as plt

first_data_line = 6
first_data_column = 3
last_data_column = 34
rotation_angle_index = -7
used_cols = (i for i in range(first_data_column, last_data_column))

ulysses_data = np.loadtxt('Ulysses_Data_File_Cleaned.txt', delimiter = ' ', skiprows = first_data_line, usecols = used_cols)

ulysses_data_without_999 = []
for i in range(len(ulysses_data)):
    if ulysses_data[i,rotation_angle_index] != 999:
        ulysses_data_without_999.append(ulysses_data[i])

ulysses_data_without_999 = np.array(ulysses_data_without_999)

print(ulysses_data_without_999[:,rotation_angle_index])
plt.hist(ulysses_data_without_999[:,rotation_angle_index])
plt.show()