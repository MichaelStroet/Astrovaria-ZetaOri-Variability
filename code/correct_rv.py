
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})

import matplotlib.cm as cm
from matplotlib.colors import Normalize

from directories import getPaths

# ======================================================================================
# Spectra
# ======================================================================================

def getObservationSpectra(dirs):

    observed_lines = {}
    for line in os.listdir(dirs["zetaOri_obs"]):
        line_dir = os.path.join(dirs["zetaOri_obs"], line)
        observed_lines[line] = {}
        for epoch in os.listdir(line_dir):
            name = epoch[:-4]#.replace("BJD", "")
            spectrum = np.genfromtxt(os.path.join(line_dir, epoch))
            observed_lines[line][name] = spectrum

    return observed_lines

def reverseShiftedSpectrum(spectrum):
    wave, flux = [spectrum[:,0], spectrum[:,1]]
    shift = wave * RV / C_KM # Doppler shift: Dlam / lam = v / c
    new_wave = wave - shift
    new_spectrum = np.dstack((new_wave, flux))[0]
    return new_spectrum

# ======================================================================================
# Saving files
# ======================================================================================

def getHeader(file):
    with open(file, "r") as f:
        header = f"Corrected for radial velocity of {RV} km/s\n"
        line = "# \n"
        while line.startswith("#"):
            header += line.replace("# ", "")
            line = f.readline()
    return header[:-1]

def saveNewSpectrum(spectrum, line, epoch):
    if not os.path.exists(os.path.join(dirs["zetaOri_rv"], line)):
        os.mkdir(os.path.join(dirs["zetaOri_rv"], line))

    old_file = os.path.join(dirs["zetaOri_obs"], line, f"{epoch}.dat")
    header = getHeader(old_file)

    new_file = os.path.join(dirs["zetaOri_rv"], line, f"{epoch}.dat")
    np.savetxt(new_file, spectrum, header = header)

def correctRV(obs_spectra):
    for line in obs_spectra:
        for epoch in obs_spectra[line]:
            spectrum = obs_spectra[line][epoch]
            new_spectrum = reverseShiftedSpectrum(obs_spectra[line][epoch])
            saveNewSpectrum(new_spectrum, line, epoch)

            # plt.figure()
            # plt.plot(spectrum[:,0], spectrum[:,1], color = "black", label = "observed")
            # plt.plot(new_spectrum[:,0], new_spectrum[:,1], color = "red", label = "corrected")
            #
            # plt.title(f"{line}: {epoch}")
            # plt.show()

# ======================================================================================

# Speed of light
C_KM = 299792.458 # km/s

# Zeta Ori radial velocity
RV = 20.55 # km/s

if __name__ == "__main__":

    dirs = getPaths()
    obs_spectra = getObservationSpectra(dirs)
    correctRV(obs_spectra)

    for line in obs_spectra:

        fig, axes = plt.subplots(4, 4, num = line, figsize = (12,8), sharex=True, sharey=True, constrained_layout=True)

        fig.suptitle(f"RV corrections for {line}")
        fig.supxlabel(r"Wavelength [$\AA$]")
        fig.supylabel(r"Flux [normed]")

        mdl_spectrum = np.genfromtxt(os.path.join(dirs["results"], "single_standard_thin", f"{line}.dat"))
        mdl_min = mdl_spectrum[:,0][np.argmin(mdl_spectrum[:,1])]

        for epoch, ax in zip(obs_spectra[line], fig.axes):
            old_spectrum = obs_spectra[line][epoch]
            old_min = old_spectrum[:,0][np.argmin(old_spectrum[:,1])]
            ax.plot(old_spectrum[:,0], old_spectrum[:,1], linewidth = 1, color = "black")
            ax.axvline(old_min, linewidth = 0.5, color = "black")

            new_spectrum = np.genfromtxt(os.path.join(dirs["zetaOri_rv"], line, f"{epoch}.dat"))
            new_min = new_spectrum[:,0][np.argmin(new_spectrum[:,1])]
            ax.plot(new_spectrum[:,0], new_spectrum[:,1], linewidth = 1, color = "blue")
            ax.axvline(new_min, linewidth = 0.5, color = "blue")

            ax.plot(mdl_spectrum[:,0], mdl_spectrum[:,1], linewidth = 1, alpha=0.8, color = "red")
            ax.axvline(mdl_min, linewidth = 0.5, color = "red")

            ax.set_title(epoch)

            window = 3
            if line == "HALPHA":
                ax.set_xlim(left = mdl_min - 10, right = mdl_min + 10)
            else:
                ax.set_xlim(left = mdl_min - window, right = mdl_min + window)

        plt.savefig(os.path.join(dirs["main"], "zetaOri", f"{line}_RV.png"))

    plt.show()
