import re

def remove_double_spaces(data_file):
    return [re.sub(' +', ' ', line) for line in data_file]

def remove_initial_spaces(data):
    for count, value in enumerate(data):
        if value[0] == ' ':
            data[count] = data[count][1:]

with open('dpf90-07.out') as ulysses_data_file:
    ulysses_data_file_cleaned = remove_double_spaces(ulysses_data_file)

remove_initial_spaces(ulysses_data_file_cleaned)

with open('Ulysses_Data_File_Cleaned.txt', 'w') as output_file:
    for line in ulysses_data_file_cleaned:
        output_file.write(line)