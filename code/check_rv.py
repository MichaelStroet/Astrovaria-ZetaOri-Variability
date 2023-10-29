
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})

import matplotlib.cm as cm
from matplotlib.colors import Normalize

from directories import getPaths

# ======================================================================================
# Auxilary functions
# ======================================================================================

def getObservationSpectra(dirs):

    observed_lines = {}
    for line in os.listdir(dirs["zetaOri_rv"]):
        line_dir = os.path.join(dirs["zetaOri_rv"], line)
        observed_lines[line] = {}
        for epoch in os.listdir(line_dir):
            name = epoch[:-4]#.replace("BJD", "")
            spectrum = np.genfromtxt(os.path.join(line_dir, epoch))
            observed_lines[line][name] = spectrum

    return observed_lines

def getModelSpectrum(model_dir, observed_lines):
    """
    Extract the spectrum for each line of each model
    """

    spectrum = {}
    for line in observed_lines:
        spectrum[line] = np.genfromtxt(os.path.join(model_dir, f"{line}.dat"))

    return spectrum

def getLinesData(file):
    """
    From the list of lines file, get the spectrum data of the lines to fit
    using a list of unique indeces.
    Returns the lines from the file and the corresponding spectrum data.
    """
    # Get the lines from the file
    contents = np.genfromtxt(file, dtype="str")

    # Store the contents as a dictionary
    lines = {}
    for line in contents:
        lines[line[0]] = {}                               # Line species
        lines[line[0]]["rest"] = line[1].astype("float")  # Rest wavelengths (Angstrom)
        lines[line[0]]["range"] = line[2].astype("float") # Range around rest (Angstrom)
        lines[line[0]]["norm"] = line[3].astype("float")  # Line amplitude
        lines[line[0]]["sigma"] = line[4].astype("float") # Line standard deviation (Angstrom)

    return lines

# ======================================================================================

def gaussian(x, center, norm, sigma, velocity = 0, continuum = 0):
    """
    Model of a gaussian line with doppler shift given by radial velocity in km/s
    Note: Normalisations are inverted. Positive norms become absorption lines.
    """
    wave_shift = velocity * center / C_KM # Doppler shift: Dlam / lam = v / c
    return -norm * np.exp(np.power((x - center - wave_shift) / sigma, 2) / -2) + continuum


def fitSingleLine(data, rest, guesses):
    """
    Fit a single line with a given model function and guess parameters.
    Returns the fit parameters and RV error
    """
    from scipy.optimize import curve_fit
    fit, cov = curve_fit(lambda x, *parameters : gaussian(x, rest, *parameters, continuum = 1),
                         xdata = data[:,0],
                         ydata = data[:,1],
                         p0 = guesses)
    rv_error = np.sqrt(np.diag(cov))[-1]

    return fit, rv_error

# ======================================================================================

def plotComparisonFits(names, fits, errs):

    plt.figure(figsize = (5,4))

    # Plot the seperate line fits\
    plt.errorbar(fits, range(len(fits)), xerr=errs, color = "blue", fmt="o", capsize = 3)

    plt.title(f"Gaussian radial velocity fits for Zeta Ori")
    plt.xlabel("Radial velocity [km/s]")
    plt.yticks(ticks = range(len(fits)), labels = names)

    legend = plt.legend(prop = {"size": 14})
    plt.grid()
    plt.tight_layout()

    # image_path = os.path.join(dir, f"fit_comparison_{target}.png")
    # plt.savefig(image_path, bbox_inches='tight')

# Speed of light
C_KM = 299792.458 # km/s

if __name__ == "__main__":

    dirs = getPaths()
    model_name = "single_standard_thin"
    dirs["output"] = os.path.join(dirs["results"], model_name)

    obs_spectra = getObservationSpectra(dirs)
    model_spectrum = getModelSpectrum(dirs["output"], obs_spectra)

    lines_txt = getLinesData("lines.txt")
    line_names = []
    fits = []
    errs = []
    print(lines_txt)


    for line in obs_spectra:

        plt.figure(line)

        # Plot observed spectra

        for epoch in obs_spectra[line]:
            spectrum = obs_spectra[line][epoch]
            if epoch == "AVERAGE":
                plt.plot(spectrum[:,0], spectrum[:,1], color = "black", label = "Observed", zorder = 2)
            else:
                plt.plot(spectrum[:,0], spectrum[:,1], color = "silver", zorder = 1)

        # Plot guess model
        rest = lines_txt[line]["rest"]
        norm = lines_txt[line]["norm"]
        sigma = lines_txt[line]["sigma"]
        line_range = lines_txt[line]["range"]
        vel = 0

        average = obs_spectra[line]["AVERAGE"]
        line_data = average[np.where((average[:,0] > rest - line_range) & (average[:,0] < rest + line_range))]

        fit, err = fitSingleLine(line_data, rest, [norm, sigma, vel])
        line_names.append(line)
        fits.append(fit[-1])
        errs.append(err)

        model_x = line_data[:,0]
        model_y = gaussian(model_x, rest, *fit) + 1

        plt.plot(model_x, model_y, color = "red", label = f"RV: {fit[-1]:.1f} km/s", zorder = 3)

        # Plot fastwind model
        # plt.plot(model_spectrum[line][:,0], model_spectrum[line][:,1], color = "red", zorder = 3)

        plt.title(line)
        plt.legend()
        plt.grid()

    print(line_names)
    for i, name in enumerate(line_names):
        print(name)
        print(f"{fits[i]:.1f} +- {errs[i]:.1f} km/s\n")

    plotComparisonFits(line_names, fits, errs)


    plt.show()
