#interesting lon and lat for start_index = 2000
#usually lon = 250 and lat widely distributed
#velocity centred on 30

import numpy as np
import matplotlib.pyplot as plt
from indices_dict import indices

first_data_column = 3
last_data_column = 34

min_quality_flag = 2
all_lines = True

used_cols = (i for i in range(first_data_column, last_data_column))


index = []
date = []
time_of_day = []
with open('Ulysses_Data_File_Cleaned.txt') as cleaned_ulysses_data:
    for count, line in enumerate(cleaned_ulysses_data):
        if count >= indices['first_data_line']:
            line = line.split()
            index.append(line[indices['index']])
            date.append(line[indices['date']])
            time_of_day.append(line[indices['time_of_day']])

index = np.array(index)
year = [date[i].split('-')[0] for i in range(len(date))]
day = [date[i].split('-')[1] for i in range(len(date))]
hour = [time_of_day[i].split(':')[0] for i in range(len(time_of_day))]
minute = [time_of_day[i].split(':')[1] for i in range(len(time_of_day))]

for i in range(len(year)):
    if year[i][0] == '9':
        year[i] = '19' + year[i]
    if year[i][0] == '0':
        year[i] = '20' + year[i]

year = np.array([int(year[i]) for i in range(len(year))])
day = np.array([int(day[i]) for i in range(len(day))])
hour = np.array([int(hour[i]) for i in range(len(hour))])
minute = np.array([int(minute[i]) for i in range(len(minute))])

passed_years = np.array([year[i]-year[0] for i in range(len(year))])

def calculate_passed_time_in_minutes(passed_years: int, day: int, hour: int, minute: int) -> int:
    return minute+60*hour+60*24*day+60*24*365*passed_years

time = calculate_passed_time_in_minutes(passed_years, day, hour, minute)

ulysses_data = np.loadtxt('Ulysses_Data_File_Cleaned.txt', delimiter = ' ', skiprows = indices['first_data_line'], usecols = used_cols)

ulysses_data = np.concatenate((np.transpose(np.array([time])), ulysses_data), axis = 1)
ulysses_data = np.concatenate((np.transpose(np.array([index])), ulysses_data), axis = 1)

if all_lines:
    start_index = 0
    end_index = len(ulysses_data)
else:
    start_index = 2000
    end_index = start_index + 1000

ulysses_data_without_999 = []
for i in range(start_index, end_index):
    condition = int(ulysses_data[i,indices['rotation_angle_index']]) != 999 and float(ulysses_data[i,indices['quality_flag_index']]) >= min_quality_flag
    if condition:
        ulysses_data_without_999.append(ulysses_data[i])

ulysses_data_without_999 = np.array(ulysses_data_without_999)

print(ulysses_data_without_999[0])

plt.hist(ulysses_data_without_999[:,indices['solar_lon_index']])
plt.show()
