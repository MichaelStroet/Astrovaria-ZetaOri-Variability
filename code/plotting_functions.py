
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})

import matplotlib.cm as cm
from matplotlib.colors import Normalize

from directories import getPaths

# ======================================================================================
# Data reading
# ======================================================================================

def getObservationSpectra(observed_dir, lines = []):
    """
    Get the observed epochs for all or a specific set of lines.
    Returns dictionary of dictionaries for each line each epoch.
    """
    if len(lines) < 1:
        lines = os.listdir(observed_dir)

    observed_lines = {}
    for line in lines:
        line_dir = os.path.join(observed_dir, line)
        observed_lines[line] = {}
        for epoch in os.listdir(line_dir):
            name = epoch[:-4].replace("BJD", "")
            spectrum = np.genfromtxt(os.path.join(line_dir, epoch))
            observed_lines[line][name] = spectrum

    return observed_lines

def getAnnelotteSpectra(dir):
    """
    Gets Annelotte's corrected HALPHA spectra and converts them from km/s to Angstrom.
    """
    files = [file for file in os.listdir(dir) if file.endswith(".txt")]

    annelotte_spectra = {}
    for file in files:
        HALPHA = 6562.790
        C_KM = 299792.458
        spectrum = np.genfromtxt(os.path.join(dir, file))
        spectrum[:,0] = HALPHA * (1 + spectrum[:,0] / C_KM)

        bjd = float(file[:-4].split("_")[-1])
        annelotte_spectra[f"{bjd:.3f}"] = spectrum

    return annelotte_spectra

def getModelSpectrum(model_dir, observed_lines):
    """
    Extract the spectrum of each line for the given model
    """

    spectrum = {}
    for line in observed_lines:
        spectrum[line] = np.genfromtxt(os.path.join(model_dir, f"{line}.dat"))

    return spectrum

# ======================================================================================
# Plotting runs
# ======================================================================================

def plotRunAllLines():
    pass

def plotRunOneLine():
    pass

def plotRun():
    plotRunAllLines()
    plotRunOneLine()

# ======================================================================================
# Single Models
# ======================================================================================

def plotModelAllLines(dirs, model_name):

    # Get the ovserved spectra
    observed_lines = getObservationSpectra(dirs["zetaOri_rv"], os.listdir(dirs["zetaOri_rv"]))

    # Set up figure for all lines plot
    fig, axes = plt.subplots(2, 5, num = "All lines", figsize = (18,8), constrained_layout=True)
    fig.suptitle(f"Lines of {model_name}".replace("_", " ").replace(",", "."))

    for line, ax in zip(observed_lines, fig.axes):
        # Plot the observed spectra
        observed_line = observed_lines[line]
        obs_average = observed_line["AVERAGE"]
        obs_size = 1
        for epoch in observed_line:
            if epoch != "AVERAGE":
                spectrum = observed_line[epoch]
                ax.plot(spectrum[:,0], spectrum[:,1], linewidth = obs_size, color = "silver")
        ax.plot(obs_average[:,0], obs_average[:,1], linewidth = obs_size, color = "black", label = "AVERAGE")

        # Plot the models
        model_size = 2
        spectrum = np.genfromtxt(os.path.join(dirs["model"], f"{line}.dat"))
        ax.plot(spectrum[:,0], spectrum[:,1], linewidth = model_size, color="red")

        ax.set_title(line)
        ax.axhline(1, color = "black", linestyle = "dashed")
        ax.set_xlim(left = min(obs_average[:,0]), right = max(obs_average[:,0]))

    plt.savefig(os.path.join(dirs["output"], f"{model_name}.png"))
    plt.close()

def plotModelOneLine():
    pass

def plotModel():
    plotModelAllLines()
    plotModelOneLine()

# ======================================================================================
# Observations
# ======================================================================================

def plotObservedLine(line):

    observed_lines = getObservationSpectra(dirs["zetaOri_rv"], os.listdir(dirs["zetaOri_rv"]))
    observed_line = observed_lines[line]
    average = observed_line["AVERAGE"]
    size = 1

    plt.figure(figsize = (16,8))
    for epoch in observed_line:
        if epoch != "AVERAGE":
            spectrum = observed_line[epoch]
            plt.plot(spectrum[:,0], spectrum[:,1], linewidth = size)
    plt.plot(average[:,0], average[:,1], linewidth = size*2, color = "black", label = "AVERAGE")

    plt.title(line)
    plt.axhline(1, color = "black", linestyle = "dashed")
    plt.xlim(left = min(average[:,0]), right = max(average[:,0]))

    # plt.show()

# ======================================================================================

if __name__ == "__main__":

    dirs = getPaths()

    plotObservedLine("HALPHA")

    # dirs["run"] = os.path.join(dirs["results"], "Mdot_Cf_invariant_T31")
    #
    # for model_name in [file for file in os.listdir(dirs["run"]) if os.path.isdir(os.path.join(dirs["run"], file))]:
    #     dirs["model"] = os.path.join(dirs["run"], model_name)
    #     dirs["output"] = dirs["run"]
    #
    #     plotModelAllLines(dirs, model_name)

    plt.show()
