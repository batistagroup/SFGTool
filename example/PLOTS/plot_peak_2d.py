###=========================================================================
### Plots sfg peaks as a function of (theta,psi) angles
###=========================================================================

import os
import sys

import numpy as np
from matplotlib import pyplot as plt

fct = 24.5  # scaling factor

###=========================================================================
### Read *.line files in path folder

path = "../ppp/azimuthal_avg/"

sfg_files = []
for temp in sorted(os.listdir(path)):
    if temp.endswith(".line"):
        sfg_files.append(path + temp)
# print('files = ',sfg_files)

nfiles = len(sfg_files)
print("nfiles = ", nfiles)

###=========================================================================
### Determine angle list

theta_list = []
psi_list = []

angle_list = []

for i in range(nfiles):
    theta = float(sfg_files[i].split("_")[-4])
    psi = float(sfg_files[i].split("_")[-3])

    angle_list.append([theta, psi])

    if theta not in theta_list:
        theta_list.append(theta)

    if psi not in psi_list:
        psi_list.append(psi)

theta_list = np.sort(theta_list)
psi_list = np.sort(psi_list)

ndim1 = len(theta_list)
ndim2 = len(psi_list)

###=========================================================================
### Read data from files


def get_sfg_data(theta, psi):

    ### get index of file for specific value of (theta,psi)
    idx = -1
    for i in range(nfiles):
        if angle_list[i][0] == theta and angle_list[i][1] == psi:
            idx = i
            break
    if idx == -1:
        sys.exit("[get_sfg_data] Error. Angles not found in sfg files.")

    ### read sfg data: x= freq, y= ampl
    temp = np.loadtxt(sfg_files[idx])
    x = temp[:, 1]
    y = temp[:, 2]

    return x, y


peak1 = np.zeros([ndim1, ndim2])
peak2 = np.zeros([ndim1, ndim2])
peak3 = np.zeros([ndim1, ndim2])

for i in range(ndim1):
    for j in range(ndim2):

        theta = theta_list[i]
        psi = psi_list[j]

        x, y = get_sfg_data(theta, psi)

        peak1[i, j] = y[0] * fct
        peak2[i, j] = y[1] * fct
        peak3[i, j] = y[2] * fct

print("max peak = ", np.max(peak1), np.max(peak2), np.max(peak3))

###=========================================================================
### Plot data (combined)

### define some options

plt.rcParams["axes.linewidth"] = 2.0

lw = 1.5

labelsz = 14.0
titlesz = 14.0
ticksz = 12

cmap = "bwr"
cmap = plt.cm.get_cmap(cmap).copy()
extend = "neither"
nlevels = 50

### create plot
figsize = [14, 6]
fig, ax = plt.subplots(1, 3, figsize=figsize, sharex=True, sharey=True)
fig.subplots_adjust(wspace=0.05)

### define levels
lmax = 1.0
lmin = -1.0
ldel = (lmax - lmin) / nlevels
levels = np.arange(lmin, lmax + ldel, ldel)

### plot data
CS0 = ax[0].contourf(
    theta_list, psi_list, np.transpose(peak1), cmap=cmap, levels=levels, extend=extend
)
CS1 = ax[1].contourf(
    theta_list, psi_list, np.transpose(peak2), cmap=cmap, levels=levels, extend=extend
)
CS2 = ax[2].contourf(
    theta_list, psi_list, np.transpose(peak3), cmap=cmap, levels=levels, extend=extend
)

### add lines
CL0 = ax[0].contour(
    theta_list,
    psi_list,
    np.transpose(peak1),
    colors="k",
    linewidths=lw,
    levels=levels[::2],
)
CL1 = ax[1].contour(
    theta_list,
    psi_list,
    np.transpose(peak2),
    colors="k",
    linewidths=lw,
    levels=levels[::2],
)
CL2 = ax[2].contour(
    theta_list,
    psi_list,
    np.transpose(peak3),
    colors="k",
    linewidths=lw,
    levels=levels[::2],
)

### set colorbar
cax = fig.add_axes([0.91, 0.11, 0.015, 0.77])
CB = fig.colorbar(CS0, cax=cax)
CB.add_lines(CL0)
CB.set_ticks(levels[::5])
CB.ax.tick_params(labelsize=ticksz)

### set axis ticks
for i in range(3):
    ax[i].tick_params(axis="both", which="major", labelsize=ticksz)

### set labels
for i in range(3):
    ax[i].set_xlabel(r"$\mathregular{\theta \ (\degree)}$", size=labelsz)
ax[0].set_ylabel(r"$\mathregular{\psi \ (\degree)}$", size=labelsz)

### show plots
# plt.show()

###=========================================================================
### Plot data (individual)


def plot_2d(x, y, z):

    ### create plot
    figsize = [8, 8]
    fig, ax = plt.subplots(figsize=figsize)

    ### define levels
    lmax = np.max(z)
    lmin = np.min(z)
    lmax = 1.0
    lmin = -1.0
    ldel = (lmax - lmin) / nlevels
    levels = np.arange(lmin, lmax + ldel, ldel)

    ### plot data
    CS = ax.contourf(x, y, np.transpose(z), cmap=cmap, levels=levels, extend=extend)
    CL = ax.contour(
        x, y, np.transpose(z), colors="k", linewidths=lw, levels=levels[::2]
    )

    ### set colorbar
    CB = fig.colorbar(CS)
    CB.add_lines(CL)
    CB.set_ticks(levels[::5])
    CB.ax.tick_params(labelsize=ticksz)

    ### set axis ticks
    ax.tick_params(axis="both", which="major", labelsize=ticksz)

    ### set labels
    ax.set_xlabel(r"$\mathregular{\theta \ (\degree)}$", size=labelsz)
    ax.set_ylabel(r"$\mathregular{\psi \ (\degree)}$", size=labelsz)

    return fig


fig1 = plot_2d(theta_list, psi_list, peak1)
fig2 = plot_2d(theta_list, psi_list, peak2)
fig3 = plot_2d(theta_list, psi_list, peak3)

### show plots
plt.show()

###=========================================================================
### Save plots

fig.savefig("peaks_2d.pdf", format="pdf", dpi=1200)

fig1.savefig("peak1_2d.pdf", format="pdf", dpi=1200)
fig2.savefig("peak2_2d.pdf", format="pdf", dpi=1200)
fig3.savefig("peak3_2d.pdf", format="pdf", dpi=1200)

###=========================================================================
### End program
print("done")
