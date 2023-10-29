
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

def getObservationSpectra(dirs, lines):

    observed_lines = {}
    for line in lines:
        line_dir = os.path.join(dirs["zetaOri_rv"], line)
        observed_lines[line] = {}
        for epoch in os.listdir(line_dir):
            name = epoch[:-4].replace("BJD", "")
            spectrum = np.genfromtxt(os.path.join(line_dir, epoch))
            observed_lines[line][name] = spectrum

    return observed_lines

def getModelSpectrum(dirs, observed_lines):
    """
    Extract the spectrum for each line of each model
    """

    spectrum = {}
    for line in observed_lines:
        spectrum[line] = np.genfromtxt(os.path.join(dirs["output"], f"{line}.dat"))

    return spectrum

def calculateRollingAverage(data, average_radius):
    """
    Calculate the rolling average of a timeseries from the points in a given radius.
    """

    if average_radius < 1:
        return data

    # Calculate the convolution kernal size and values
    kernel_size = 2 * average_radius + 1
    kernel = np.ones(kernel_size) / kernel_size

    # Points falling inside the radius from the edges are removed
    time = data[average_radius:-average_radius, 0]

    # Convolve the data with the kernel to get the rolling average
    value = np.convolve(data[:,1], kernel, mode = "valid")
    return np.dstack((time, value))[0]

# ======================================================================================
# Plotting functions
# ======================================================================================

def plotEpochsLine(model_dir, model_name, model_spectrum, epoch_spectra, line):

    fig, axes = plt.subplots(4, 4, num = line, figsize = (12,8), constrained_layout=True, sharex=True, sharey=True)
    fig.suptitle(f"Epochs of {line} with {model_name}".replace("_", " ").replace(",", "."))

    for epoch, ax in zip(epoch_spectra, fig.axes):
        spectrum = epoch_spectra[epoch]

        # Plot the observed spectra
        obs_size = 2
        spectrum = calculateRollingAverage(spectrum, 0)
        ax.plot(spectrum[:,0], spectrum[:,1], linewidth = obs_size, color = "black")

        # Plot the models
        model_size = 2
        ax.plot(model_spectrum[:,0], model_spectrum[:,1], linewidth = model_size, color = "red")

        ax.set_title(epoch)
        ax.axhline(1, color = "black", linestyle = "dashed")

        ax.set_xlim(left = min(spectrum[:,0]), right = max(spectrum[:,0]))
        y_pad = 0.20 * (max(model_spectrum[:,1]) - min(model_spectrum[:,1]))
        ax.set_ylim(bottom = min(model_spectrum[:,1]) - y_pad, top = max(model_spectrum[:,1]) + y_pad)

    plt.savefig(os.path.join(model_dir, f"Epochs_{line}.png"))
    plt.close()

def plotResidualEpochs(model_dir, model_name, model_spectrum, epoch_spectra, line, rest_wave):

    fig, axes = plt.subplots(4, 4, num = line, figsize = (12,8), constrained_layout=True, sharex=True, sharey=True)
    fig.suptitle(f"Residual epochs of {line} with {model_name}".replace("_", " ").replace(",", "."))

    for epoch, ax in zip(epoch_spectra, fig.axes):
        spectrum = epoch_spectra[epoch]
        interpolated_model = np.interp(spectrum[:,0], model_spectrum[:,0], model_spectrum[:,1])
        residuals = interpolated_model / spectrum[:,1]

        size = 1
        ax.plot(spectrum[:,0], residuals, linewidth = size, color = "black")
        ax.set_title(epoch)
        ax.axhline(1, color = "red", linestyle = "dashed")
        ax.axvline(rest_wave, linewidth = 0.5, color = "red")

        ax.set_xlim(left = min(spectrum[:,0]), right = max(spectrum[:,0]))

        if line == "HALPHA":
            ax.set_ylim(bottom = 0.9, top = 1.1)

    plt.savefig(os.path.join(model_dir, f"Residuals_{line}.png"))
    plt.close()

# ======================================================================================
# Main plot funciton
# ======================================================================================

def plotEpochs(dirs, model_name, lines = None):

    line_dict = {
    "HEI4026":       4026.1914,
    "HDELTA":        4101.734,
    "HEII4200":      4199.8390,
    "HGAMMA":        4340.472,
    "HEI4471":       4471.4802,
    "HEII4686":      4685.698,
    "HEI4713":       4713.1457,
    "HBETA":         4861.350,
    "HEII5411":      5411.516,
    "HALPHA":        6562.790
    }

    if type(lines) == type(None):
        lines = os.listdir(dirs["zetaOri_rv"])

    # Get the ovserved spectra
    observed_lines = getObservationSpectra(dirs, lines)

    # Get the model spectrum
    model_spectrum = getModelSpectrum(dirs, observed_lines)

    # Make a big plot
    for line in observed_lines:
        plotEpochsLine(dirs["output"], model_name, model_spectrum[line], observed_lines[line], line)
        # plotResidualEpochs(dirs["output"], model_name, model_spectrum[line], observed_lines[line], line, line_dict[line])

def plotAllLinesSingleModel(dirs, model_name):

    # Get the ovserved spectra
    observed_lines = getObservationSpectra(dirs, os.listdir(dirs["zetaOri_rv"]))

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
        spectrum = np.genfromtxt(os.path.join(dirs["output"], f"{line}.dat"))
        ax.plot(spectrum[:,0], spectrum[:,1], linewidth = model_size, color="red")

        ax.set_title(line)
        ax.axhline(1, color = "black", linestyle = "dashed")
        ax.set_xlim(left = min(obs_average[:,0]), right = max(obs_average[:,0]))

    plt.savefig(os.path.join(dirs["output"], "AllLines.png"))
    plt.close()

def plotAllLinesRun(dirs, model_name):

    # Get the ovserved spectra
    observed_lines = getObservationSpectra(dirs, os.listdir(dirs["zetaOri_rv"]))

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
        spectrum = np.genfromtxt(os.path.join(dirs["output"], f"{line}.dat"))
        ax.plot(spectrum[:,0], spectrum[:,1], linewidth = model_size, color="red")

        ax.set_title(line)
        ax.axhline(1, color = "black", linestyle = "dashed")
        ax.set_xlim(left = min(obs_average[:,0]), right = max(obs_average[:,0]))

    plt.savefig(os.path.join(os.path.dirname(dirs["output"]), f"{model_name}.png"))
    plt.close()

# ======================================================================================

if __name__ == "__main__":

    dirs = getPaths()

    # variations_dir = os.path.join(dirs["results"], "Mdot-Cf")
    # for model_name in [file for file in os.listdir(variations_dir) if os.path.isdir(os.path.join(variations_dir, file))]:
    #     print(model_name)
    #     dirs["output"] = os.path.join(variations_dir, model_name)
    #     plotAllLinesRun(dirs, model_name)
    #     # plotEpochs(dirs, model_name)

    model_name = "standard_model"
    dirs["output"] = os.path.join(dirs["results"], model_name)
    plotEpochs(dirs, model_name)

    plt.show()
