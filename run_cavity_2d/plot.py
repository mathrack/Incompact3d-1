#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Read channel.dat
tmp = np.loadtxt("cavity_bulk.dat", dtype=str, delimiter="\n", skiprows=1)

# Number of columns (<u>, <v>, <T>, <u'²>, <v'²>, <T'²>)
ncol = len(tmp[0].split())

# Create numpy array
array=np.zeros((tmp.shape[0], ncol))
for i in range(tmp.shape[0]):
    array[i,:] = np.array([float(item) for item in tmp[i].split()])

# Plot
for i in range(ncol):
    if i<ncol/2:
        plt.plot(array[:,i])
        plt.title("Column " + str(i))
    else:
        plt.plot(np.log10(array[:,i]))
        plt.title("Column " + str(i) + ", log scale")
    plt.show()
