import csv
import numpy as np

with open('IC_BH_velocity.csv', 'r') as csvfile:
	csvreader = csv.reader(csvfile)
	header = next(csvreader)
	v_values = next(csvreader)

v_values = [np.float32(j) for j in v_values]

print(header)
print(v_values)
print(v_values[header.index('Mean velocity [km/s]')])
