#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Subroutine to extract / plot the histogram of a variable
def extract_and_plot(f, name, prec):
    # Matplotlib figure and axis
    if (plot):
        fig, ax = plt.subplots()
    #
    # Size of the histogram
    tmp = f.read(4)
    n = np.fromfile(f, dtype=np.int32, count=1)
    tmp = f.read(4)
    #
    # First bin
    tmp = f.read(4)
    xmin = np.fromfile(f, dtype=prec, count=1)
    xmax = np.fromfile(f, dtype=prec, count=1)
    xnumr = np.fromfile(f, dtype=np.int32, count=1)
    np.fromfile(f, dtype=prec, count=1)
    tmp = f.read(4)
    if (plot):
        xmin = xmin[0]
        xmax = xmax[0]
        xnumr = max(xnumr[0], 0.5)
        if xnumr>=1:
            tmp = ax.add_patch(Rectangle((xmin, 0.), xmax-xmin, float(xnumr)))
        if xnumr<0:
            print("Error in first bin")
        tmp = ax.set_xlim(left=xmin)
    #
    # Read the size bins
    for ibin in range(n[0]-2):
        tmp = f.read(4)
        xmin = np.fromfile(f, dtype=prec, count=1)
        xmax = np.fromfile(f, dtype=prec, count=1)
        xnum = np.fromfile(f, dtype=np.int32, count=1)
        tmp = f.read(4)
        if (plot):
            xmin = xmin[0]
            xmax = xmax[0]
            xnum = xnum[0]
            xnumr = max(xnum, xnumr)
            if xnum>0:
                tmp = ax.add_patch(Rectangle((xmin, 0.), xmax-xmin, float(xnum)))
                tmp = ax.set_ylim(top=xnumr)
            if xnum<0:
                print("Error in bin "+str(ibin))
    #
    # Last bin
    tmp = f.read(4)
    xmin = np.fromfile(f, dtype=prec, count=1)
    xmax = np.fromfile(f, dtype=prec, count=1)
    xnum = np.fromfile(f, dtype=np.int32, count=1)
    np.fromfile(f, dtype=prec, count=1)
    tmp = f.read(4)
    if (plot):
        xmin = xmin[0]
        xmax = xmax[0]
        xnum = xnum[0]
        xnumr = max(xnum, xnumr)
        if xnum>0:
            tmp = ax.add_patch(Rectangle((xmin, 0.), xmax-xmin, float(xnum)))
            tmp = ax.set_ylim(top=xnumr)
        if xnum<0:
            print("Error in last bin")
        tmp = ax.set_xlim(right=xmax)
    #
    # Display
    if (plot):
        tmp = ax.set_ylim(bottom=1.)
        tmp = ax.set_xlabel("Value")
        tmp = ax.set_ylabel("Number of samples")
        tmp = ax.set_yscale("log")
        tmp = ax.set_title("Histogram for "+name)
        tmp = fig.show()

# Flag to activate the plot
plot = True

# Name of the file to process
fname = "out/dp_histogram_100_100.bin"

# Open the file
fdp = open(fname, "rb")

# Name of the file to process
fname = "out/sp_histogram_100_100.bin"

# Open the file
fsp = open(fname, "rb")

extract_and_plot(fsp, "u (SP)", np.float32)
extract_and_plot(fdp, "u (DP)", np.float64)
extract_and_plot(fsp, "v (SP)", np.float32)
extract_and_plot(fdp, "v (DP)", np.float64)
extract_and_plot(fsp, "t (SP)", np.float32)
extract_and_plot(fdp, "t (DP)", np.float64)

input("Input for x-derivative")

extract_and_plot(fsp, "dudx (SP)", np.float32)
extract_and_plot(fdp, "dudx (DP)", np.float64)
extract_and_plot(fsp, "dvdx (SP)", np.float32)
extract_and_plot(fdp, "dvdx (DP)", np.float64)
extract_and_plot(fsp, "dtdx (SP)", np.float32)
extract_and_plot(fdp, "dtdx (DP)", np.float64)

input("Input for y-derivative")

extract_and_plot(fsp, "dudy (SP)", np.float32)
extract_and_plot(fdp, "dudy (DP)", np.float64)
extract_and_plot(fsp, "dvdy (SP)", np.float32)
extract_and_plot(fdp, "dvdy (DP)", np.float64)
extract_and_plot(fsp, "dtdy (SP)", np.float32)
extract_and_plot(fdp, "dtdy (DP)", np.float64)

# Close the file
fsp.close()
fdp.close()
