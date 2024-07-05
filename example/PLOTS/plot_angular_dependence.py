###=========================================================================
### Plots sfg spectra for a series of selected (theta,psi) angles
###=========================================================================

import os
import sys

import numpy as np
from matplotlib import pyplot as plt

fct = 1e6  # scaling factor

###=========================================================================
### Define angles to plot

theta_list = [10, 30, 60]
psi_list = [0, 30, 60, 90]

ndim1 = len(theta_list)
ndim2 = len(psi_list)

print("ndim = ", ndim1, ndim2)

###=========================================================================
### Read *.sfg files in path folder

path = "../ppp/azimuthal_avg/"

sfg_files = []
for temp in sorted(os.listdir(path)):
    if temp.endswith(".sfg"):
        sfg_files.append(path + temp)
# print('files = ',sfg_files)

nfiles = len(sfg_files)
print("nfiles = ", nfiles)

###=========================================================================
### Determine angle list

angle_list = []

for i in range(nfiles):
    theta = float(sfg_files[i].split("_")[-4])
    psi = float(sfg_files[i].split("_")[-3])
    angle_list.append([theta, psi])

###=========================================================================
### Define function to read data from files


def get_sfg_data(theta, psi):

    ### get index of file for specific value of (theta,psi)
    idx = -1
    for i in range(nfiles):
        if angle_list[i][0] == theta and angle_list[i][1] == psi:
            idx = i
            break
    if idx == -1:
        sys.exit(
            "[get_sfg_data] Error. Angles ({},{}) not found in sfg files.".format(
                theta, psi
            )
        )

    ### read sfg data: x= freq, y= sfg intensity
    temp = np.loadtxt(sfg_files[idx])
    x = temp[:, 0]
    y = temp[:, 1]

    return x, y


###=========================================================================
### Plot data

### define some options
plt.rcParams["axes.linewidth"] = 2.0

labelsz = 12
ticksz = 10

figsize = [8, 8]

### create plot
fig, ax = plt.subplots(ndim2, ndim1, figsize=figsize, sharex=True, sharey=True)
fig.subplots_adjust(wspace=0.0, hspace=0.0)

### plot data
for i in range(ndim2):
    psi = psi_list[i]
    for j in range(ndim1):
        theta = theta_list[j]

        ### load data
        x, y = get_sfg_data(theta, psi)

        ax[i, j].plot(x, y * fct, "k", lw=2, label="{}{}".format(theta, psi))

### set axis limits
ax[0, 0].set_xlim([1950, 2149])

### set axis ticks
ax[0, 0].tick_params(axis="both", which="major", labelsize=ticksz)

### set x/y labels
for i in range(ndim1):
    ax[-1, i].set_xlabel(r"Wavenumber $\mathregular{(cm^{-1})}$", size=labelsz)

for i in range(ndim2):
    ax[i, 0].set_ylabel(r"$\mathregular{|\chi^{(2)}|^2}$ (arb. u.)", size=labelsz)

### set additional labels
for i in range(ndim1):
    theta = theta_list[i]
    ax[0, i].set_xlabel(
        r"$\mathregular{\theta = }$" + "{}$\degree$".format(theta), size=14.0
    )
    ax[0, i].xaxis.set_label_position("top")

for i in range(ndim2):
    psi = psi_list[i]
    ax[i, -1].set_ylabel(
        r"$\mathregular{\psi = }$" + "{}$\degree$".format(psi),
        size=14.0,
        rotation=270,
        labelpad=20,
    )
    ax[i, -1].yaxis.set_label_position("right")

### show plots
plt.show()

### save plots
fig.savefig("angular_dependence.pdf", format="pdf", dpi=1200)

###=========================================================================
### End program
print("done")
