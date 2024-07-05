###=================================================
### Plots sfg spectra for different polarizations
###=================================================

import numpy as np
from matplotlib import pyplot as plt

fct = 1e6  # normalization factor

###=================================================
### Define files

theta = 30
psi = 30

sfg_files = []

sfg_files.append("../ppp/azimuthal_avg/Rebpy_{}_{}_avg_ppp.sfg".format(theta, psi))
sfg_files.append("../ssp/Rebpy_{}_{}_avg_ssp.sfg".format(theta, psi))
sfg_files.append("../sps/Rebpy_{}_{}_avg_sps.sfg".format(theta, psi))

label = []
label.append("ppp")
label.append("ssp")
label.append("sps")

nfiles = len(sfg_files)
print("nfiles = ", nfiles)

###=================================================
### Plot data

### define some options
plt.rcParams["axes.linewidth"] = 2.0

titlesz = 16
legendsz = 14
labelsz = 12
ticksz = 12

figsize = [8, 8]

### create plot
fig, ax = plt.subplots(nfiles, figsize=figsize, sharex=True)

### plot data
for i in range(nfiles):

    ### read sfg data: x= freq, y= sfg intensity
    temp = np.loadtxt(sfg_files[i])
    x = temp[:, 0]
    y = temp[:, 1]

    ### plot
    ax[i].plot(x, y * fct, c="k", lw=2, label=label[i])

### set axis limits
for i in range(nfiles):
    ax[i].set_xlim([1950, 2150])

### set axis ticks
for i in range(nfiles):
    ax[i].tick_params(axis="both", which="major", labelsize=ticksz)


### set title
theta = str(theta) + r"$\mathregular{\degree}$"
psi = str(psi) + r"$\mathregular{\degree}$"
title = r"$\mathregular{(\theta,\psi) =}$" + "({},{})".format(theta, psi)
ax[0].set_title(title, size=titlesz)

### set labels
ax[-1].set_xlabel(r"Wavenumber $\mathregular{(cm^{-1})}$", size=labelsz)

for i in range(nfiles):
    ax[i].set_ylabel(r"$\mathregular{|\chi^{(2)}|^2}$ (arb. u.)", size=labelsz)

for i in range(nfiles):
    ax[i].legend(loc=1, fontsize=legendsz, numpoints=10, frameon=False)

### show plots
plt.show()

### save plots
fig.savefig("polarization_dependence.pdf", format="pdf", dpi=1200)
