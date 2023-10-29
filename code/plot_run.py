
import os
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as cm
from matplotlib.colors import Normalize

from directories import getPaths
from plot_epochs import plotEpochsLine

plt.rcParams.update({'font.size': 16})

# ======================================================================================
# Plotting functions
# =====================================================================================

def getAnnelotteSpectra(dirs):
    annelotte_Ha = []
    for epoch in os.listdir(dirs["zetaOri_Annelotte"]):
        if epoch.startswith("zetOri_Ha"):
            spectrum = np.genfromtxt(os.path.join(dirs["zetaOri_Annelotte"], epoch))
            HALPHA = 6562.790
            C_KM = 299792.458
            spectrum[:,0] = HALPHA * (1 + spectrum[:,0] / C_KM)
            # ax.plot(spectrum[:,0], spectrum[:,1], linewidth = obs_size, color = "silver")
            annelotte_Ha.append(spectrum)
    annelotte_Ha = np.array(annelotte_Ha)
    avg_flux = np.mean(annelotte_Ha[:, :, 1], axis = 0)

    return annelotte_Ha, avg_flux

def getRunSpectra(run_dir, line):
    spectra = {}
    values = []
    files = os.listdir(run_dir)
    for file in files:
        if os.path.isdir(os.path.join(run_dir, file)):
            value = file.split("_")[-1].replace(",", ".")
            values.append(float(value))

            spectrum = np.genfromtxt(os.path.join(run_dir, file, f"{line}.dat"))
            spectra[value] = spectrum

    return spectra, values

# ======================================================================================

if __name__ == "__main__":

    dirs = getPaths()

    observed_lines, observed_avg = getAnnelotteSpectra(dirs)

    Mdot_dir = os.path.join(dirs["results"], "03 - 2nd parameter study", "Run_M_-6,70_-6,22")
    Fcl_dir = os.path.join(dirs["results"], "03 - 2nd parameter study", "Run_Cf_1,10_200")

    Mdot_run, Mdots = getRunSpectra(Mdot_dir, "HALPHA")
    Fcl_run, Fcls = getRunSpectra(Fcl_dir, "HALPHA")
    print(Mdot_run)

    fig, (axM, axF) = plt.subplots(1, 2, figsize = (16,5), sharex=True, sharey=True, constrained_layout=True)

    # fig.supxlabel(r"Wavelength [$\AA$]")
    fig.supylabel(r"Normalised Flux")

    for ax in [axM, axF]:
        for spectrum in observed_lines:
            ax.plot(spectrum[:,0], spectrum[:,1], linewidth = 2.5, color = "silver")
        ax.plot(spectrum[:,0], observed_avg, linewidth = 2.5, color = "black")
        ax.axhline(1, color = "black", linestyle = "dashed")
        ax.set_xlim(left = 6546, right = 6584)
        ax.set_xlabel(r"Wavelength [$\AA$]", fontsize = 18)

    # Define a color map for the Mdot plot
    print(min(Mdots))
    norm = Normalize(min(Mdots), max(Mdots))
    scalarmap = cm.ScalarMappable(norm=norm, cmap="viridis")
    scalarmap._A = []

    for value in Mdot_run:
        spectrum = Mdot_run[value]
        color = scalarmap.to_rgba(float(value))
        axM.plot(spectrum[:,0], spectrum[:,1], linewidth = 1, color = color)

    plt.colorbar(scalarmap, ax = axM, ticks = Mdots, label = r"$\dot{M} [M_\odot/yr]$")

    # Define a color map for the Fcl plot
    norm = Normalize(min(Fcls), max(Fcls))
    scalarmap = cm.ScalarMappable(norm=norm, cmap="viridis")
    scalarmap._A = []

    for value in Fcl_run:
        spectrum = Fcl_run[value]
        color = scalarmap.to_rgba(float(value))
        axF.plot(spectrum[:,0], spectrum[:,1], linewidth = 1, color = color)

    plt.colorbar(scalarmap, ax = axF, ticks = Fcls, label = r"$f_{\rm cl}$")

    plt.savefig(os.path.join(dirs["results"], "plots", "MdotFcl_variations.pdf"))

    plt.show()
