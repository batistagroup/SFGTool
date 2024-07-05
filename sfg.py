#!/usr/bin/python
import sys

import numpy as np

Ha2cm = 2.1947463136320e5  # Hartree to cm-1 (Eh/hc , CODATA 2018)
Bohr2Ang = 0.529177210903  # Bohr to Angstrom (CODATA 2018)
me2amu = (
    9.1093837015e-31 / 1.66053906660e-27
)  # electron mass to atomic mass unit (me/mu , CODATA 2018)


###====================================================================
### Main program
###====================================================================
def main():

    global sfgtype, sfgpol
    global theta, psi, phi
    global thetatype, psitype, phitype
    global rotconv

    global iverbose

    global visangle, visfreq, irangle, irfreq
    global n1, n2, n3

    global freq
    global dipole
    global polar
    global modes
    global chi, chi_amplitude
    global dip, dip_amplitude

    global thetamin, thetamax, thetastep
    global psimin, psimax, psistep
    global phimin, phimax, phistep

    global scalef, width
    global freqmin, freqmax, freqstep
    global freqaxis
    global Lijk
    global freqfile, dipfile, polarfile, outfile, expfile
    global nr_phase, nr_ampl

    ###=================================================
    ### Usage
    usage = "\n Usage: python {} sfg.inp \n".format(sys.argv[0])
    if len(sys.argv) != 2:
        sys.exit(usage)
    ###=================================================

    ###=================================================
    ### Initialize variables with default values
    set_defaults()
    ###=================================================

    ###=================================================
    ### Read input file
    inputfile = sys.argv[1]
    read_input(inputfile)
    ###=================================================

    ###=================================================
    ### Define Fresnel factors
    Lijk = Fresnel_factors()
    ###=================================================

    ###=================================================
    ### Read frequency, dipole and polarizability

    ### frequency
    if freqfile == "null":
        sys.exit(
            '\n Error. File for reading freq. information ("freqfile") must be set. \n'
        )
    freq = read_freq(freqfile)

    ### dipole
    if dipfile == "null":
        sys.exit(
            '\n Error. File for reading dipole information ("dipfile") must be set. \n'
        )
    dipole = read_dipole(dipfile)

    ### polar
    if sfgtype != "pmirras":
        if polarfile == "null":
            sys.exit(
                '\n Error. File for reading polar information ("polarfile") must be set. \n'
            )
        polar = read_polar(polarfile, sfgtype)

    ###=================================================
    ### Select only data in the freq. range

    ftemp = np.copy(freq)
    dtemp = np.copy(dipole)
    if sfgtype != "pmirras":
        ptemp = np.copy(polar)

    freq = []
    dipole = []
    polar = []
    modes = []
    for i in range(len(ftemp)):
        if ftemp[i] > freqmin and ftemp[i] < freqmax:
            modes.append(i + 1)
            freq.append(ftemp[i])
            dipole.append(dtemp[i])
            if sfgtype != "pmirras":
                polar.append(ptemp[i])

    ###=================================================
    ### Print data to file

    if iverbose != 0:
        if outfile == "null":
            sys.exit(
                '\n Error. File for outputing data information ("outfile") must be set. \n'
            )
        print_data(outfile)

    ###=================================================
    ### Define spectral grid

    freqaxis = define_axis()

    ###=================================================
    ### Compute SFG spectra

    print("Computing SFG spectra...")

    ### determine initial angle (if not fixed in input)
    if thetatype == "scan":
        theta = thetamin
    if psitype == "scan" or psitype == "avg":
        psi = psimin
    if phitype == "scan" or phitype == "avg":
        phi = phimin

    if sfgtype == "pmirras":
        phitype = "fix"
        phi = 0.0

    while theta < thetamax:
        while psi < psimax:
            while phi < phimax:
                ### define label
                label1, label2, label3 = get_label()
                print("for angles: ({},{},{})".format(label1, label2, label3))
                ### create chi tensor in laboratory frame for each normal mode
                if (
                    psitype == "avg" and phitype == "avg"
                ):  # twist+azimuthal average (psi+phi angles)
                    print("performing twist+azimuthal average (psi+phi avg)")
                    if sfgtype == "pmirras":
                        dip = dip_psiphi_avg(theta)
                    else:
                        chi = chi_psiphi_avg(theta)
                elif phitype == "avg":  # azimuthal average (phi angle)
                    print("performing azimuthal average (phi avg)")
                    if sfgtype == "pmirras":
                        dip = dip_phi_avg(theta, psi)
                    else:
                        chi = chi_phi_avg(theta, psi)
                else:  # fix angles
                    print("fix angles (no avg)")
                    rotmat = Euler_matrix(theta, psi, phi, rotconv)
                    if sfgtype == "pmirras":
                        dip = dip_tensor(rotmat)
                    else:
                        chi = chi_tensor(rotmat)
                ### compute spectra
                if sfgtype == "pmirras":
                    dip_amplitude = dip_amplitudes()
                else:
                    chi_amplitude = sfg_amplitudes()
                ### save data to file
                if sfgtype == "pmirras":
                    save_dip_data()
                else:
                    save_chi_data()
                ### go to next azimuthal (phi) angle if scanning
                if phitype == "scan":
                    phi += phistep
                else:
                    break

            ### reset azimuthal (phi) angle if scanning other angles
            if (thetatype == "scan" or psitype == "scan") and phitype == "scan":
                phi = phimin

            ### go to next twist (psi) angle if scanning
            if psitype == "scan":
                psi += psistep
            else:
                break

        ### reset twist (psi) angle if scanning other angles
        if thetatype == "scan" and psitype == "scan":
            psi = psimin

        ### go to next tilt (theta) angle if scanning
        if thetatype == "scan":
            theta += thetastep
        else:
            break

    print("Computing SFG spectra...done")

    print("\n THANKS FOR USING THE PROGRAM. HAVE A GOOD DAY! \n")
    return


###====================================================================
### End Main program
###====================================================================

###====================================================================
### Auxiliary functions
###====================================================================


### ---------------------------------------------------------------------
def set_defaults():
    """Initialize variables and set default values"""

    global sfgtype, sfgpol
    global theta, psi, phi
    global thetatype, psitype, phitype
    global scalef, width
    global thetamin, thetamax, thetastep
    global psimin, psimax, psistep
    global phimin, phimax, phistep
    global freqmin, freqmax, freqstep
    global nr_phase, nr_ampl
    global visangle, visfreq, irangle, irfreq
    global n1, n2, n3
    global iverbose
    global freqfile, dipfile, polarfile, outfile, expfile
    global rotconv

    ### Options controlling verbosity of output
    iverbose = 1

    ### Options controlling type of SFG
    sfgtype = "sfg"
    sfgpol = "ppp"

    ### Options controlling molecular orientation
    theta = 0.0
    psi = 0.0
    phi = 0.0

    thetatype = "fix"  # for internal use only
    psitype = "fix"  # for internal use only
    phitype = "fix"  # for internal use only

    ### Options controlling Euler angles
    rotconv = "zyz"
    thetamin = np.radians(0.0)
    thetamax = np.radians(180.0)
    thetastep = np.radians(5.0)
    psimin = np.radians(0.0)
    psimax = np.radians(360.0)
    psistep = np.radians(5.0)
    phimin = np.radians(0.0)
    phimax = np.radians(360.0)
    phistep = np.radians(5.0)

    ### Options controlling input/output files
    freqfile = "null"
    dipfile = "null"
    polarfile = "null"
    outfile = "null"
    expfile = "null"

    ### Options controlling spectral output
    scalef = 1.0
    width = 1.0
    freqmin = 0.0
    freqmax = 4000.0
    freqstep = 0.1

    ### Options controlling non-resonant background
    nr_ampl = 0.0
    nr_phase = 1.0

    ### Options controlling Fresnel factors (3-layer model)
    visangle = 45.0
    visfreq = 12500.0  # 800 nm
    irangle = 45.0
    irfreq = 3000.0
    n1 = np.zeros(3)
    n1[0] = 1.000  # air refractive index
    n1[1] = 1.000  # air refractive index
    n1[2] = 1.000  # air refractive index
    n2 = np.zeros(3)
    n2[0] = 1.333  # water refractive index
    n2[1] = 1.333  # water refractive index
    n2[2] = 1.333  # water refractive index
    n3 = np.zeros(3)
    n3[0] = 1.180  # from Ref. PRB 59, 12632 (1999)
    n3[1] = 1.180  # from Ref. PRB 59, 12632 (1999)
    n3[2] = 1.180  # from Ref. PRB 59, 12632 (1999)

    return


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def read_input(inputfile):
    """Read inputfile"""

    global sfgtype, sfgpol
    global theta, psi, phi
    global thetatype, psitype, phitype
    global rotconv
    global scalef, width
    global thetamin, thetamax, thetastep
    global psimin, psimax, psistep
    global phimin, phimax, phistep
    global freqmin, freqmax, freqstep
    global iverbose

    global visangle, visfreq, irangle, irfreq
    global n1, n2, n3

    global freqfile, dipfile, polarfile, outfile, expfile
    global nr_phase, nr_ampl

    print("Reading input file...")

    ifile = open(inputfile, "r")
    for line in ifile:

        if line.startswith("#"):  # skip commented lines
            continue
        if line.startswith("\n"):  # skip blank lines
            continue
        # 		line = line.lower() 								# convert to lower case

        keyword = line.split("=")[0].strip()  # get keyword
        value = line.split("=")[1].split("#")[0].strip()  # get value without comment

        ### Options controlling verbosity of output
        if keyword == "iverbose":
            iverbose = int(value)

        ### Options controlling type of SFG
        if keyword == "sfgtype":
            sfgtype = value
        if keyword == "sfgpol":
            sfgpol = value

        ### Options controlling molecular orientation
        if keyword == "theta":
            try:
                theta = np.radians(float(value))
            except ValueError:
                if value == "scan":
                    theta = value
                    thetatype = value
                else:
                    sys.exit(
                        '\n [read_input] Error. theta should be a number or "scan" \n'
                    )
        if keyword == "psi":
            try:
                psi = np.radians(float(value))
            except ValueError:
                if value == "scan" or value == "avg":
                    psi = value
                    psitype = value
                else:
                    sys.exit(
                        '\n [read_input] Error. psi should be a number or "scan" or "avg" \n'
                    )
        if keyword == "phi":
            try:
                phi = np.radians(float(value))
            except ValueError:
                if value == "scan" or value == "avg":
                    phi = value
                    phitype = value
                else:
                    sys.exit(
                        '\n [read_input] Error. phi should be a number or "scan" or "avg" \n'
                    )

        ### Options controlling Euler angles
        if keyword == "thetamin":
            thetamin = np.radians(float(value))
        if keyword == "thetamax":
            thetamax = np.radians(float(value))
        if keyword == "thetastep":
            thetastep = np.radians(float(value))
        if keyword == "psimin":
            psimin = np.radians(float(value))
        if keyword == "psimax":
            psimax = np.radians(float(value))
        if keyword == "psistep":
            psistep = np.radians(float(value))
        if keyword == "phimin":
            phimin = np.radians(float(value))
        if keyword == "phimax":
            phimax = np.radians(float(value))
        if keyword == "phistep":
            phistep = np.radians(float(value))
        if keyword == "rotconv":
            rotconv = value

        ### Options controlling input/output files
        if keyword == "freqfile":
            freqfile = value
        if keyword == "dipfile":
            dipfile = value
        if keyword == "polarfile":
            polarfile = value
        if keyword == "outfile":
            outfile = value
        if keyword == "expfile":
            expfile = value

        ### Options controlling spectral output
        if keyword == "scalef":
            scalef = float(value)
        if keyword == "width":
            width = float(value)
        if keyword == "freqmin":
            freqmin = float(value)
        if keyword == "freqmax":
            freqmax = float(value)
        if keyword == "freqstep":
            freqstep = float(value)

        ### Options controlling non-resonant background
        if keyword == "nr_ampl":
            nr_ampl = float(value)
        if keyword == "nr_phase":
            nr_phase = np.radians(float(value))

        ### Options controlling Fresnel factors (3-layer model)
        if keyword == "visangle":
            visangle = np.radians(float(value))
        if keyword == "irangle":
            irangle = np.radians(float(value))
        if keyword == "visbfreq":
            visfreq = float(value)
        if keyword == "irfreq":
            irfreq = float(value)
        if keyword == "n1":
            values = value.split()
            if len(values) == 1:  # freq.-independent refractive index
                n1[0] = float(values[0])
                n1[1] = float(values[0])
                n1[2] = float(values[0])
            elif len(values) == 3:  # freq.-dependent refractive index
                n1[0] = float(values[0])
                n1[1] = float(values[1])
                n1[2] = float(values[2])
            else:
                sys.exit(
                    "\n Error. Refractive index n1 should be 1 (freq.-independent) or 3 (freq.-dependent) numbers. \n"
                )
        if keyword == "n2":
            values = value.split()
            if len(values) == 1:  # freq.-independent refractive index
                n2[0] = float(values[0])
                n2[1] = float(values[0])
                n2[2] = float(values[0])
            elif len(values) == 3:  # freq.-dependent refractive index
                n2[0] = float(values[0])
                n2[1] = float(values[1])
                n2[2] = float(values[2])
            else:
                sys.exit(
                    "\n Error. Refractive index n2 should be 1 (freq.-independent) or 3 (freq.-dependent) numbers. \n"
                )
        if keyword == "n3":
            values = value.split()
            if len(values) == 1:  # freq.-independent refractive index
                n3[0] = float(values[0])
                n3[1] = float(values[0])
                n3[2] = float(values[0])
            elif len(values) == 3:  # freq.-dependent refractive index
                n3[0] = float(values[0])
                n3[1] = float(values[1])
                n3[2] = float(values[2])
            else:
                sys.exit(
                    "\n Error. Refractive index n3 should be 1 (freq.-independent) or 3 (freq.-dependent) numbers. \n"
                )

    print("Reading input file...done")
    return


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def Fresnel_factors():
    """
    Compute Fresnel factors using three-layer model.
    Note1: Frequency order is SFG,Vis,IR
    Note2: Implementation correspond to co-propagation in reflection geometry.
    Ref.: International Reviews in Physical Chemistry, 24:2, 191-256 (2005).
    """
    # FIXME TODO:include other options

    print("Computing Fresnel factors using three-layer model...")

    ### sfg freq
    sfgfreq = visfreq + irfreq

    ### incident/reflection angles (beta_i angles in Ref.)
    beta = np.zeros(3)
    beta[0] = np.arcsin(
        (visfreq * np.sin(visangle) + irfreq * np.sin(irangle)) / sfgfreq
    )  # phase-matching condition [Eq. 4]
    beta[1] = visangle
    beta[2] = irangle

    ### Fresnel factors [Eq. 5]
    Lxx = np.zeros(3)
    Lyy = np.zeros(3)
    Lzz = np.zeros(3)
    for i in range(3):
        gamma_i = np.arcsin(n1[i] * np.sin(beta[i]) / n2[i])  # refractive angle
        Lxx[i] = (
            2.0
            * n1[i]
            * np.cos(gamma_i)
            / (n1[i] * np.cos(gamma_i) + n2[i] * np.cos(beta[i]))
        )
        Lyy[i] = (
            2.0
            * n1[i]
            * np.cos(beta[i])
            / (n1[i] * np.cos(beta[i]) + n2[i] * np.cos(gamma_i))
        )
        Lzz[i] = (
            2.0
            * n2[i]
            * np.cos(beta[i])
            / (n1[i] * np.cos(gamma_i) + n2[i] * np.cos(beta[i]))
            * (n1[i] / n3[i]) ** 2
        )

    ### include incident/reflection angles in Fresnel factors
    for i in range(3):
        Lxx[i] *= np.cos(beta[i])
        Lzz[i] *= np.sin(beta[i])

    ### L_ijk
    Lijk = np.zeros([3, 3, 3])

    Lijk[0, 0, 2] = Lxx[0] * Lxx[1] * Lzz[2]  # Lxxz
    Lijk[0, 2, 0] = Lxx[0] * Lzz[1] * Lxx[2]  # Lxzx
    Lijk[0, 1, 2] = Lxx[0] * Lyy[1] * Lzz[2]  # Lxyz
    Lijk[0, 2, 1] = Lxx[0] * Lzz[1] * Lyy[2]  # Lxzy

    Lijk[1, 0, 2] = Lyy[0] * Lxx[1] * Lzz[2]  # Lyxz
    Lijk[1, 2, 0] = Lyy[0] * Lzz[1] * Lxx[2]  # Lyzx
    Lijk[1, 1, 2] = Lyy[0] * Lyy[1] * Lzz[2]  # Lyyz
    Lijk[1, 2, 1] = Lyy[0] * Lzz[1] * Lyy[2]  # Lyzy

    Lijk[2, 0, 0] = Lzz[0] * Lxx[1] * Lxx[2]  # Lzxx
    Lijk[2, 1, 0] = Lzz[0] * Lyy[1] * Lxx[2]  # Lzyx
    Lijk[2, 0, 1] = Lzz[0] * Lxx[1] * Lyy[2]  # Lzxy
    Lijk[2, 1, 1] = Lzz[0] * Lyy[1] * Lyy[2]  # Lzyy
    Lijk[2, 2, 2] = Lzz[0] * Lzz[1] * Lzz[2]  # Lzzz

    print("Computing Fresnel factors...done")
    return Lijk


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def read_freq(filename):
    """
    Read frequencies from file
    Units: cm^{-1}
    """

    try:
        lines = open(filename, "r").readlines()
    except IOError:
        sys.exit('\n [read_freq] Error. Could not open file "{}". \n'.format(filename))

    print('Reading frequencies from file "{}"...'.format(filename))

    freq = read_freq_logfile(lines)

    print('Reading frequencies from file "{}"...done'.format(filename))

    return freq


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def read_dipole(filename):
    """
    Read dipole from file
    Units: Debye
    """
    # FIXME TODO: include additional ways of reading dipole info

    try:
        lines = open(filename, "r").readlines()
    except IOError:
        sys.exit(
            '\n [read_dipole] Error. Could not open file "{}". \n'.format(filename)
        )

    print('Reading dipoles from file "{}"...'.format(filename))

    dipole = read_dipole_logfile(lines)

    print('Reading dipoles from file "{}"...done'.format(filename))

    return dipole


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def read_polar(filename, sfgtype):
    """
    Read polarizabilities from file
    Units: \AA^{3}
    """
    # FIXME TODO: include additional ways of reading polar info

    try:
        lines = open(filename, "r").readlines()
    except IOError:
        sys.exit('\n [read_polar] Error. Could not open file "{}". \n'.format(filename))

    print('Reading polarizabilities from file "{}"...'.format(filename))

    if sfgtype == "sfg":
        polar = read_polar_logfile_sfg(lines)
    elif sfgtype == "drsfg":
        polar = read_polar_logfile_drsfg(lines)
    else:
        sys.exit('\n [read_polar] Error. sfgtype "{}" not defined. \n'.format(sfgtype))

    print('Reading polarizabilities from file "{}"...done'.format(filename))

    return polar


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def define_axis():
    """Define spectral axis"""

    if expfile != "null":  # read axis from experimental spectra
        ifile = open(expfile, "r")
        axis = []
        for line in ifile:
            temp = line.split()
            axis.append(float(temp[0]))
        ifile.close()
    else:  # define using (freqmin,freqmax)
        axis = []
        temp = freqmin
        while temp < freqmax:
            axis.append(temp)
            temp += freqstep
    return axis


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def init_chi():
    """Initialize complex chi tensor"""
    return np.zeros([len(modes), 3, 3, 3], dtype=complex)


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def Euler_matrix(theta, psi, phi, rotconv):
    """
    Compute Euler matrix for given angles
    """

    R = np.zeros([3, 3])

    if rotconv == "zyz":  # ZYZ convention

        R[0, 0] = np.cos(phi) * np.cos(theta) * np.cos(psi) - np.sin(phi) * np.sin(psi)
        R[0, 1] = -np.cos(phi) * np.cos(theta) * np.sin(psi) - np.sin(phi) * np.cos(psi)
        R[0, 2] = np.cos(phi) * np.sin(theta)

        R[1, 0] = np.sin(phi) * np.cos(theta) * np.cos(psi) + np.cos(phi) * np.sin(psi)
        R[1, 1] = -np.sin(phi) * np.cos(theta) * np.sin(psi) + np.cos(phi) * np.cos(psi)
        R[1, 2] = np.sin(phi) * np.sin(theta)

        R[2, 0] = -np.sin(theta) * np.cos(psi)
        R[2, 1] = np.sin(theta) * np.sin(psi)
        R[2, 2] = np.cos(theta)

    elif rotconv == "zxz":  # ZXZ convention

        R[0, 0] = -np.sin(phi) * np.cos(theta) * np.sin(psi) + np.cos(phi) * np.cos(psi)
        R[1, 0] = np.cos(phi) * np.cos(theta) * np.sin(psi) + np.sin(phi) * np.cos(psi)
        R[2, 0] = np.sin(theta) * np.sin(psi)

        R[0, 1] = -np.sin(phi) * np.cos(theta) * np.cos(psi) - np.cos(phi) * np.sin(psi)
        R[1, 1] = np.cos(phi) * np.cos(theta) * np.cos(psi) - np.sin(phi) * np.sin(psi)
        R[2, 1] = np.sin(theta) * np.cos(psi)

        R[0, 2] = np.sin(phi) * np.sin(theta)
        R[1, 2] = -np.cos(phi) * np.sin(theta)
        R[2, 2] = np.cos(theta)

    else:  # error
        sys.exit(
            '\n [Euler_matrix] Error. rotconv "{}" not defined. \n'.format(rotconv)
        )

    return R


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def chi_tensor(rotmat):
    """Construct chi tensor for particular rotmat using the hyperpolarizability
            obtained from the dipole and polarizability tensors
    Units: Debye*\AA^3
    """

    chi = init_chi()

    for i in range(len(modes)):

        ### rotate dipole
        dtemp = np.matmul(rotmat, dipole[i])

        ### rotate polar
        ptemp = np.matmul(rotmat, np.matmul(polar[i], np.transpose(rotmat)))

        ### compute hyperpolarizability
        chi[i] = np.tensordot(ptemp, dtemp, axes=0)  # tensor product

    return chi


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def sfg_amplitudes():
    """Compute SFG amplitudes for each normal modes"""
    # TODO: include mix polarizations

    amplitude = []

    for i in range(len(modes)):
        if sfgpol == "ppp":
            temp = (
                Lijk[2, 0, 0] * chi[i][2, 0, 0]
                + Lijk[2, 2, 2] * chi[i][2, 2, 2]
                - Lijk[0, 0, 2] * chi[i][0, 0, 2]
                - Lijk[0, 2, 0] * chi[i][0, 2, 0]
            )
        elif sfgpol == "ssp":
            temp = Lijk[1, 1, 2] * chi[i][1, 1, 2]
        elif sfgpol == "sps":
            temp = Lijk[1, 2, 1] * chi[i][1, 2, 1]
        elif sfgpol == "psp":
            temp = Lijk[2, 1, 0] * chi[i][2, 1, 0] - Lijk[0, 1, 2] * chi[i][0, 1, 2]
        elif sfgpol == "pss":
            temp = Lijk[2, 1, 1] * chi[i][2, 1, 1]
        elif sfgpol == "zzz":
            temp = Lijk[2, 2, 2] * chi[i][2, 2, 2]
        else:
            sys.exit(
                '\n [sfg_amplitudes] Error. sfgpol "{}" not implemented. \n'.format(
                    sfgpol
                )
            )

        amplitude.append(temp)

    return amplitude


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def get_label():
    """Define label for saving data"""

    if thetatype == "fix" or thetatype == "scan":
        label1 = round(np.degrees(theta))
    else:
        sys.exit(
            '\n [get_label] Error. Label not defined for thetatype "{}". \n'.format(
                thetatype
            )
        )

    if psitype == "fix" or psitype == "scan":
        label2 = round(np.degrees(psi))
    elif psitype == "avg":
        label2 = psitype
    else:
        sys.exit(
            '\n [get_label] Error. Label not defined for psitype "{}". \n'.format(
                psitype
            )
        )

    if phitype == "fix" or phitype == "scan":
        label3 = round(np.degrees(phi))
    elif phitype == "avg":
        label3 = phitype
    else:
        sys.exit(
            '\n [get_label] Error. Label not defined for phitype "{}". \n'.format(
                phitype
            )
        )

    return label1, label2, label3


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def save_chi_data():
    """Save spectra to file"""

    ### define label
    label1, label2, label3 = get_label()

    ### line spectra
    filename = outfile + "_{}_{}_{}_{}.line".format(label1, label2, label3, sfgpol)
    ofile = open(filename, "w")
    ofile.write(
        "# {:4} {:10} {:18} {:18} \n".format("mode", "freq", "A.real", "A.imag")
    )
    for i in range(len(modes)):
        ofile.write(
            "{:4} {:8.2f} {:12.6f} {:12.6f} \n".format(
                modes[i], freq[i], chi_amplitude[i].real, chi_amplitude[i].imag
            )
        )
    ofile.close()

    ### spectra
    filename = outfile + "_{}_{}_{}_{}.sfg".format(label1, label2, label3, sfgpol)
    ofile = open(filename, "w")
    ofile.write("# {:10} {:14} {:10} {} \n".format("freq", "|X|^2", "X.imag", "X.real"))
    nr_fct = nr_ampl * np.exp(nr_phase * 1j)  # non-resonant part
    for i in range(len(freqaxis)):  # loop over w_ir
        temp = nr_fct
        for j in range(len(modes)):  # loop over q
            temp += chi_amplitude[j] / (freqaxis[i] - float(freq[j]) + width * (1j))
        ofile.write(
            "{:8.2f}   {:1.4e}   {:1.4e}   {:1.4e} \n".format(
                freqaxis[i], abs(temp) ** 2, temp.imag, temp.real
            )
        )
    ofile.close()

    return


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def chi_phi_avg(theta, psi):
    """Perform azimuthal average (phi angle) of chi"""

    chi = init_chi()

    phicount = 0
    phi = phimin
    while phi < phimax:
        ### define Euler matrix
        rotmat = Euler_matrix(theta, psi, phi, rotconv)
        ### create chi tensor for each normal mode
        chi += chi_tensor(rotmat)
        ### go to next angle
        phi += phistep
        phicount += 1
    ### normalize average
    chi /= phicount

    return chi


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def chi_psiphi_avg(theta):
    """Perform twist and azimuthal average (psi+phi angle) of chi"""

    chi = init_chi()

    psicount = 0
    psi = psimin
    while psi < psimax:  # average twist angle (psi)
        chi += chi_phi_avg(theta, psi)
        ### go to next twist angle (psi)
        psi += psistep
        psicount += 1
    ### normalize average
    chi /= psicount

    return chi


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def print_data(filename):

    try:
        ofile = open(outfile + ".dat", "w")
    except IOError:
        sys.exit(
            '\n [print_data] Error. Could not open file "{}".\n '.format(
                filename + ".dat"
            )
        )

    ### compute hyperpolarizability for reference frame
    if sfgtype != "pmirras":
        rotmat = Euler_matrix(0, 0, 0, rotconv)
        beta = chi_tensor(rotmat)

    ### program version
    ofile.write("-----------------------\n")
    ofile.write("-    SFG script       -\n")
    ofile.write("-   Version 0.1       -\n")
    ofile.write("-----------------------\n")

    ### print molecular info
    for i in range(len(modes)):
        ### label
        ofile.write("\n Mode {}".format(modes[i]))
        ### freq
        ofile.write("\n Frequency {:.2f} cm-1".format(freq[i]))
        ### dip
        ofile.write("\n Dipole: ")
        ofile.write("\n {} {} {}".format(dipole[i][0], dipole[i][1], dipole[i][2]))
        ### polar
        if sfgtype != "pmirras":
            ofile.write("\n Polarizability: ")
            for j in range(3):
                ofile.write(
                    "\n {} {} {}".format(polar[i][j][0], polar[i][j][1], polar[i][j][2])
                )
        ### hyperpolar
        if iverbose >= 3 and sfgtype != "pmirras":
            ofile.write("\n HyperPolarizability: ")
            for j in range(3):
                for k in range(3):
                    ofile.write(
                        "\n {} {} {}".format(
                            beta[i][j][k][0], beta[i][j][k][1], beta[i][j][k][2]
                        )
                    )
        ofile.write("\n")
    ofile.write("\n")

    ### print Fresnel factors info
    if iverbose >= 2 and sfgtype != "pmirras":
        ofile.write("----------------------")
        ofile.write("\n SFG Fresnel Factors\n")
        ofile.write("----------------------")
        for i in range(3):
            for j in range(3):
                ofile.write(
                    "\n {:.4f} {:.4f} {:.4f}".format(
                        Lijk[i][j][0], Lijk[i][j][1], Lijk[i][j][2]
                    )
                )
            ofile.write("\n")
        ofile.write("\n")
        ofile.write("----------------------")

    ### close file
    ofile.close()

    return


### ---------------------------------------------------------------------
def read_freq_logfile(lines):
    """
    Read harmonic frequencies from Gaussian logfile
    Units: Gaussian prints freq. in cm^{-1}
    Gaussian keyword: freq
    """

    label = "Frequencies"

    freq = []

    ierr = 1
    for i in range(len(lines)):
        if label in lines[i]:
            ierr = 0
            temp = lines[i].split()
            if len(temp) == 5:  # non-linear molecule
                freq.append(float(temp[2]) * scalef)
                freq.append(float(temp[3]) * scalef)
                freq.append(float(temp[4]) * scalef)
            elif len(temp) == 3:  # linear molecule
                freq.append(float(temp[2]) * scalef)

    if ierr == 1:
        sys.exit(
            '\n [read_freq] Error. Label "{}" not found in logfile. \n'.format(label)
        )

    return freq


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def read_dipole_logfile(lines):
    """
    Read dipole from Gaussian logfile
    Gaussian keyword: freq(PrintDerivatives)
    Note: Gaussian actually computes the dipole derivatives wrt mass-weigthed normal modes
          in (km/mol)^1/2 units [Source: Gaussian support private email].
              Therefore we convert from dipole derivative to dipole using the harmonic formula
                                            mu = sqrt(hbar/2.*omega) * dmu/dq
              and convert to Debye units.
    """

    label = "Dipole derivative"

    dipole = []

    ierr = 1
    for i in range(len(lines)):
        if label in lines[i]:
            ierr = 0
            temp = lines[i].split()
            temp1 = float(
                temp[5].replace("D", "E")
            )  # replace exponential 'D' notation with 'E'
            temp2 = float(temp[6].replace("D", "E"))
            temp3 = float(temp[7].replace("D", "E"))
            temp = [temp1, temp2, temp3]
            dipole.append(temp)

    if ierr == 1:
        sys.exit(
            '\n [read_dipole] Error. Label "{}" not found in logfile. \n'.format(label)
        )

    ### convert dipole derivatives to dipole
    for i in range(len(dipole)):
        dipole[i] /= np.sqrt(2.0 * freq[i])

    ### convert units
    fct = 1.0 / 42.2561  # 42.2561 km/mol = Debye^2 / \AA^2 / amu
    for i in range(len(dipole)):
        dipole[i] *= np.sqrt(fct * Ha2cm * me2amu) * Bohr2Ang

    return dipole


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def read_polar_logfile_sfg(lines):
    """
    Read polarizability from Gaussian logfile
    Gaussian keyword: freq(Polar,PrintDerivatives)
    Note: Gaussian actually computes the polar. derivatives wrt mass-weigthed normal modes
          in \AA^2/amu^1/2 units [Source: Gaussian support private email].
              Therefore we convert from polar. derivative to polar. using the harmonic formula
                                            alpha = sqrt(hbar/2.*omega) * dalpha/dq
              and convert to \AA^3 units.
    """

    label = "Polarizability derivatives"
    offset = 2
    polar = []
    ierr = 1
    for i in range(len(lines)):
        if label in lines[i]:
            ierr = 0
            ptemp = []
            for j in range(3):
                temp = lines[i + offset + j].split()
                temp1 = float(
                    temp[1].replace("D", "E")
                )  # replace exponential 'D' notation with 'E'
                temp2 = float(temp[2].replace("D", "E"))
                temp3 = float(temp[3].replace("D", "E"))
                temp = [temp1, temp2, temp3]
                ptemp.append(temp)
            polar.append(ptemp)
    if ierr == 1:
        sys.exit(
            '\n [read_polar_logfile_sfg] Error. Label "{}" not found in logfile. \n'.format(
                label
            )
        )

    ### convert polar derivatives to polar
    for i in range(len(polar)):
        polar[i] /= np.sqrt(2.0 * freq[i])

    ### convert units
    for i in range(len(polar)):
        polar[i] *= np.sqrt(Ha2cm * me2amu) * Bohr2Ang

    return polar


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def read_polar_logfile_drsfg(lines):
    """
    Read resonant polarizability from Gaussian logfile
    Gaussian keyword: freq(readFCHT)
    Units: Gaussian prints polar. in \AA^3.
    """

    label = "ELECTRIC DIPOLE-ELECTRIC DIPOLE"
    offset = 3
    polar = []
    ierr = 1
    for i in range(len(lines)):
        if label in lines[i]:
            ierr = 0
            ptemp = np.zeros([3, 3], dtype=complex)
            a = 0
            b = 0
            for j in range(9):
                temp = lines[i + offset + j].split()
                tempr = float(temp[1].replace("D", "E"))
                tempi = float(temp[2].replace("D", "E"))
                ptemp[a, b] = tempr + 1j * tempi
                b += 1
                if b == 3:
                    a += 1
                    b = 0
            polar.append(ptemp)
    if ierr == 1:
        sys.exit(
            '\n [read_polar_logfile_drsfg] Error. Label "{}" not found in logfile. \n'.format(
                label
            )
        )

    return polar


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def dip_amplitudes():
    """Compute PM-IRRAS amplitudes for each normal modes"""

    amplitude = []

    for i in range(len(modes)):
        temp = dip[i][2] ** 2

        amplitude.append(temp)

    return amplitude


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def dip_tensor(rotmat):
    """Construct pm-irras dip tensor for particular rotmat using dipole tensors
    Units: Debye
    """

    dip = np.zeros([len(modes), 3], dtype=float)

    for i in range(len(modes)):
        ### rotate dipole
        dip[i] = np.matmul(rotmat, dipole[i])

    return dip


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def dip_phi_avg(theta, psi):
    """Perform azimuthal average (phi angle) of dip"""

    dip = np.zeros([len(modes), 3], dtype=float)

    phicount = 0
    phi = phimin
    while phi < phimax:
        ### define Euler matrix
        rotmat = Euler_matrix(theta, psi, phi, rotconv)
        ### create pmirras dip tensor for each normal mode
        dip += dip_tensor(rotmat)
        ### go to next angle
        phi += phistep
        phicount += 1
    ### normalize average
    dip /= phicount

    return dip


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def dip_psiphi_avg(theta):
    """Perform twist and azimuthal average (psi+phi angle) of dip"""

    dip = np.zeros([len(modes), 3], dtype=float)

    psicount = 0
    psi = psimin
    while psi < psimax:  # average twist angle (psi)
        dip += dip_phi_avg(theta, psi)
        ### go to next twist angle (psi)
        psi += psistep
        psicount += 1
    ### normalize average
    dip /= psicount

    return dip


### ---------------------------------------------------------------------


### ---------------------------------------------------------------------
def save_dip_data():
    """Save spectra to file"""

    ### define label
    label1, label2, label3 = get_label()

    ### line spectra
    filename = outfile + "_{}_{}_pmirras.line".format(label1, label2)
    ofile = open(filename, "w")
    ofile.write("# {:4} {:10} {:16} \n".format("mode", "freq", "I"))
    for i in range(len(modes)):
        ofile.write(
            "{:4} {:10.4f} {:16.10f} \n".format(modes[i], freq[i], dip_amplitude[i])
        )
    ofile.close()

    ### spectra
    filename = outfile + "_{}_{}_pmirras.pmirras".format(label1, label2)
    ofile = open(filename, "w")
    ofile.write("# {:10} {:14} \n".format("freq", "I"))
    for i in range(len(freqaxis)):  # loop over w_ir
        temp = 0.0
        for j in range(len(modes)):  # loop over q
            temp += dip_amplitude[j] * real_Lorentzian(freq[j], freqaxis[i], width)
        ofile.write("{:8.2f}   {:1.4e}   \n".format(freqaxis[i], temp))
    ofile.close()

    return


### ---------------------------------------------------------------------


def real_Lorentzian(x, x0, gamma):
    return (1.0 / np.pi) * (0.5 * gamma) / ((x - x0) ** 2 + (0.5 * gamma) ** 2)


###====================================================================
### __name__
###====================================================================

if __name__ == "__main__":
    main()
