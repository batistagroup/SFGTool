###=========================================================================
### Plots polar (phi) dependence of sfg spectra for selected (theta,psi) angles
###=========================================================================

import os
import sys

import numpy as np
from matplotlib import pyplot as plt

fct = 40  # normalization factor
itype = 1  # 0=amplitude ; 1=|amplitude| ; 2=amplitude**2

###=========================================================================
### Define angles to plot

theta_psi_list = [[0, 0], [30, 30], [30, 0], [60, 0]]

ndim1 = len(theta_psi_list)

print("theta/psi angles: ", theta_psi_list)

###=========================================================================
### Read *.line files in path folder

path = "../ppp/azimuthal_anisotropy/"

sfg_files = []
for temp in sorted(os.listdir(path)):
    if temp.endswith(".line"):
        sfg_files.append(path + temp)
print("files = ", sfg_files)

nfiles = len(sfg_files)
print("nfiles = ", nfiles)

###=========================================================================
### Determine angle list

phi_list = []
angle_list = []
for i in range(nfiles):
    theta = int(sfg_files[i].split("_")[-4])
    psi = int(sfg_files[i].split("_")[-3])
    phi = int(sfg_files[i].split("_")[-2])

    if phi not in phi_list:
        phi_list.append(phi)

    angle_list.append([theta, psi, phi])

phi_list = np.sort(phi_list)
print("phi_list = ", phi_list)

ndim2 = len(phi_list)

###=================================================
### Read data from files


def get_sfg_data(theta, psi, phi):

    ### get index of file for specific value of (theta,psi,phi)
    idx = -1
    for i in range(nfiles):
        if (
            angle_list[i][0] == theta
            and angle_list[i][1] == psi
            and angle_list[i][2] == phi
        ):
            idx = i
            break
    if idx == -1:
        sys.exit(
            "[get_sfg_data] Error. Angles ({},{},{}) not found in sfg files.".format(
                theta, psi, phi
            )
        )

    ### read sfg data: x= freq, y= sfg amplitude
    temp = np.loadtxt(sfg_files[idx])
    x = temp[:, 1]
    y = temp[:, 2]

    return x, y


peak1 = np.zeros([ndim1, ndim2])
peak2 = np.zeros([ndim1, ndim2])
peak3 = np.zeros([ndim1, ndim2])

for i in range(ndim1):

    theta = theta_psi_list[i][0]
    psi = theta_psi_list[i][1]
    for j in range(ndim2):
        phi = phi_list[j]
        x, y = get_sfg_data(theta, psi, phi)
        if itype == 0:
            peak1[i, j] = y[0] * fct
            peak2[i, j] = y[1] * fct
            peak3[i, j] = y[2] * fct
        elif itype == 1:
            peak1[i, j] = np.abs(y[0]) * fct
            peak2[i, j] = np.abs(y[1]) * fct
            peak3[i, j] = np.abs(y[2]) * fct
        elif itype == 2:
            peak1[i, j] = y[0] ** 2 * fct
            peak2[i, j] = y[1] ** 2 * fct
            peak3[i, j] = y[2] ** 2 * fct

###=================================================
### Plot data

### define some options
plt.rcParams["axes.linewidth"] = 2.0

titlesz = 16
labelsz = 10
ticksz = 10

figsize = [8, 8]

### loop over (theta,psi) angles
for i in range(ndim1):

    ### create plot
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=figsize)

    ### plot data
    ax.plot(np.radians(phi_list), peak1[i, :], label="a`(2)")
    ax.plot(np.radians(phi_list), peak2[i, :], label="a``")
    ax.plot(np.radians(phi_list), peak3[i, :], label="a`(1)")

    ### customize lines and ticks
    ax.set_rmax(1.1)
    if itype != 0:
        ax.set_rmin(0.0)
    else:
        ax.set_rmin(-1.1)
    if itype != 0:
        ax.set_rticks(np.arange(0.0, 1.1, 0.1))
    else:
        ax.set_rticks(np.arange(-1.0, 1.1, 0.2))
    # 	ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line

    ### set title/labels/legend
    theta = str(theta_psi_list[i][0]) + r"$\mathregular{\degree}$"
    psi = str(theta_psi_list[i][1]) + r"$\mathregular{\degree}$"
    title = r"$\mathregular{(\theta,\psi) =}$" + "({},{})".format(theta, psi)
    ax.set_title(title, size=titlesz)
    ax.legend(fontsize=12, frameon=False)

    ### save figure
    theta = theta_psi_list[i][0]
    psi = theta_psi_list[i][1]
    filename = "azimuthal_dependende_{}_{}.pdf".format(theta, psi)
    fig.savefig(filename, format="pdf", dpi=1200)

### show plots
plt.show()
