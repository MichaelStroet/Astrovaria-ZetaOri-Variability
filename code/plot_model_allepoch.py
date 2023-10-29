
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

import matplotlib.cm as cm
from matplotlib.colors import Normalize

from directories import getPaths

# =================================================================================================
# Auxilary functions
# =================================================================================================

def getObservationSpectra(dirs):

    left_peak = "BJD2458406.5453007.dat"
    average = "AVERAGE.dat"
    right_peak = "BJD2458404.599795.dat"

    observed_lines = {}
    for line in os.listdir(dirs["zetaOri_rv"]):
        observed_lines[line] = []
        epochs_dir = os.path.join(dirs["zetaOri_rv"], line)
        for epoch in os.listdir(epochs_dir):

            if not epoch.startswith("AVERAGE"):
                observed_lines[line].append(np.genfromtxt(os.path.join(epochs_dir, epoch)))

    return observed_lines

def getModelSpectra(model_dir, observed_lines):
    """
    Extract the spectrum for each line of each model
    """

    spectrum = {}
    for line in observed_lines:
        spectrum[line] = np.genfromtxt(os.path.join(model_dir, f"{line}.dat"))

    return spectrum

# =================================================================================================
# Plotting functions
# =================================================================================================

def plotSingleLine(model_dir, line, observed_lines, model_spectrum):
    """

    """

    obs_l = observed_lines[line]["left"]
    obs_m = observed_lines[line]["mean"]
    obs_r = observed_lines[line]["right"]

    plt.figure(figsize = (16,9))

    # Continuum line
    plt.axhline(1, color = "red", linestyle = "dashed")

    # Observed spectra
    obs_size = 2
    plt.plot(obs_l[:,0], obs_l[:,1], linewidth = obs_size, color = "silver", label = "left")
    plt.plot(obs_r[:,0], obs_r[:,1], linewidth = obs_size, color = "silver", label = "right")
    plt.plot(obs_m[:,0], obs_m[:,1], linewidth = obs_size, color = "black", label = "mean")

    plt.plot(model_spectrum[line][:,0], model_spectrum[line][:,1], linewidth = 2)

    plt.title(line)
    plt.xlabel(r"Wavelength [$\AA$]")
    plt.ylabel("Flux [arbitrary]")

    plt.xlim(left = min(obs_m[:,0]), right = max(obs_m[:,0]))

    plt.legend()
    plt.tight_layout()

    plt.savefig(os.path.join(model_dir, line + ".png"))
    plt.close()

def plotAllLines(model_dir, model_name, observed_lines, model_spectrum):

    fig, axes = plt.subplots(2, 5, figsize = (18,8), constrained_layout=True)

    for line, ax in zip(observed_lines, fig.axes):
        print(line)
        # Plot the observed spectra
        obs_size = 2
        if line != "HALPHA":

            spectra = np.array(observed_lines[line])

            for spectrum in spectra:
                # print(f"{line} {len(spectrum)}\n{spectrum}\n\n")
                ax.plot(spectrum[:,0], spectrum[:,1], linewidth = obs_size, color = "silver")

            avg_flux = np.mean(spectra[:, :, 1], axis = 0)
            ax.plot(spectrum[:,0], avg_flux, linewidth = obs_size, color = "black")

            ax.set_xlim(left = max([min(spectrum[:,0]), min(model_spectrum[line][:,0])]),
                        right = min([max(spectrum[:,0]), max(model_spectrum[line][:,0])]))

        else:

            annelotte_Ha = []
            for epoch in os.listdir(dirs["zetaOri_Annelotte"]):
                if epoch.startswith("zetOri_Ha"):
                    spectrum = np.genfromtxt(os.path.join(dirs["zetaOri_Annelotte"], epoch))
                    HALPHA = 6562.790
                    C_KM = 299792.458
                    spectrum[:,0] = HALPHA * (1 + spectrum[:,0] / C_KM)
                    ax.plot(spectrum[:,0], spectrum[:,1], linewidth = obs_size, color = "silver")
                    annelotte_Ha.append(spectrum)

            annelotte_Ha = np.array(annelotte_Ha)
            avg_flux = np.mean(annelotte_Ha[:, :, 1], axis = 0)
            ax.plot(spectrum[:,0], avg_flux, linewidth = obs_size, color = "black")

            ax.set_xlim(left = 6540, right = 6590)

        # Plot the models
        model_size = 3
        ax.plot(model_spectrum[line][:,0], model_spectrum[line][:,1], linewidth = model_size, color = "red")

        titles = {"HALPHA" : r"H$\alpha$",
                  "HBETA" : r"H$\beta$",
                  "HGAMMA" : r"H$\gamma$",
                  "HDELTA" : r"H$\delta$",
                  "HEI4026" : r"HeI $\lambda$4026",
                  "HEI4471" : r"HeI $\lambda$4471",
                  "HEI4713" : r"HeI $\lambda$4713",
                  "HEII4200" : r"HeII $\lambda$4200",
                  "HEII4686" : r"HeII $\lambda$4686",
                  "HEII5411" : r"HeII $\lambda$5411"}
        ax.set_title(titles[line])
        ax.axhline(1, color = "black", linestyle = "dashed")

    fig.supxlabel(r"Wavelength [$\AA$]")
    fig.supylabel(r"Normalised Flux")

    plt.savefig(os.path.join(model_dir,"StandardModel.pdf"))
    # plt.savefig(os.path.join(model_dir,"AllLines.pdf"))

# =================================================================================================
# Main plot funciton
# =================================================================================================

def plotModel(dirs, model_name):

    # Get the ovserved spectra
    observed_lines = getObservationSpectra(dirs)

    # Get the model spectra
    model_spectra = getModelSpectra(dirs["model"], observed_lines)

    # Make a plot for each line individually
    # for line in observed_lines:
    #     plotSingleLine(dirs["model"], line, observed_lines, model_spectra)

    # Make a big plot for all lines simultaneously
    plotAllLines(dirs["model"], model_name, observed_lines, model_spectra)

# =================================================================================================

if __name__ == "__main__":

    dirs = getPaths()
    # model_name = "single_M-6,47_B1,10_g3,30_T31_Y0,08_Cf50,00_Cs0,10"
    # model_dir = os.path.join(dirs["results"], "04 - bunch of random models", "variations_lots", model_name)
    model_name = "standard_model_01feb"
    model_dir = os.path.join(dirs["results"], "02 - finding new standard", model_name)
    dirs["model"] = model_dir

    plotModel(dirs, model_name)

    plt.show()
