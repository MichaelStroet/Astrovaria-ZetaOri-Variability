import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

from directories import getPaths
from plotting_functions import getObservationSpectra, getModelSpectrum


def calculateEquivalentWidth(spectrum, line_l = 6540, line_r = 6590, cont_l = 6538, cont_r = 6542):
    """
    Calculates the equivalent width of a line within a given wavelength range.
    Calculates the error on the EW using Vollmann&Eversberg 2006
        => err = sqrt(1 + avg_continuum_flux / avg_flux) * (wavelength_range - EW) / SNR
    Returns the equivalent width and error.
    """

    # Calcuate continuum values between the given bounds
    continuum = spectrum[np.where((spectrum[:,0] >= cont_l) & (spectrum[:,0] <= cont_r))]
    cont_avg = np.mean(continuum[:,1])
    cont_std = np.std(continuum[:,1])

    # Get the spectrum between the given bounds
    line = spectrum[np.where((spectrum[:,0] >= line_l) & (spectrum[:,0] <= line_r))]
    flux_avg = np.mean(line[:,1])

    # Calculate the equivalent width and its error
    equivalent_width = integrate.simpson(y = (1 - line[:,1]/cont_avg), x = line[:,0])
    error = np.sqrt(1 + (cont_avg/flux_avg)) * ((line_r - line_l) - equivalent_width) / (1/cont_std)

    return equivalent_width, error

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

def cutSpectrum(spectrum, line):
    """
    Correct a spectrum for telluric lines by lazily cutting out their regions.
    """

    with open("EW_wavelengths.txt", "r") as file:
        contents = file.read().split("\n")

    for row in contents:
        if row.startswith(line):
            ranges = row[len(line)+1:].split()

    if type(ranges) == type("none"):
        return spectrum

    indeces = []
    for range in ranges:
        left, right = range.split(":")
        indeces += np.where((spectrum[:,0] >= float(left)) & (spectrum[:,0] <= float(right)))[0].tolist()

    new_spectrum = spectrum[indeces]
    return new_spectrum

# ======================================================================================

if __name__ == "__main__":

    dirs = getPaths()

    observed_lines = getObservationSpectra(dirs["zetaOri_rv"], os.listdir(dirs["zetaOri_rv"]))

    model_spectrum = getModelSpectrum(os.path.join(dirs["results"], "standard_model"), observed_lines)
    annelotte_spectra = getAnnelotteSpectra(dirs["zetaOri_Annelotte"])

    output_dir = os.path.join(dirs["EW"], "corrected_spectra_plots")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    line = "HALPHA"
    bound_left = 6540
    bound_right = 6590

    original_widths = []
    cut_widths = []
    annelotte_widths = []

    original_errs = []
    cut_errs = []
    annelotte_errs = []

    dates = []
    epoch_0 = 0

    # EW of the model
    model_width = calculateEquivalentWidth(model_spectrum[line], line_l=bound_left, line_r=bound_right)

    for epoch in observed_lines[line]:

        if epoch != "AVERAGE":

            # Add the epoch date to the list
            date = float(epoch)
            dates.append(date)

            if epoch_0 == 0:
                epoch_0 = date

            # Calculate EW of the original spectrum with telluric lines
            original_spectrum = observed_lines[line][epoch]
            original_EW, original_err = calculateEquivalentWidth(original_spectrum, line_l=bound_left, line_r=bound_right)
            original_widths.append(original_EW)
            original_errs.append(original_err)

            # Calculate EW of the 'cut-corrected' spectrum
            cut_spectrum = cutSpectrum(original_spectrum, line)
            cut_EW, cut_err = calculateEquivalentWidth(cut_spectrum, line_l=bound_left, line_r=bound_right)
            cut_widths.append(cut_EW)
            cut_errs.append(cut_err)

            try:
                annelotte_spectrum = annelotte_spectra[f"{date:.3f}"]
            except:
                annelotte_spectrum = np.ones(original_spectrum.shape)
                annelotte_spectrum[:,0] = original_spectrum[:,0]

            plt.figure(f"{date - epoch_0:.2f}")
            plt.plot(original_spectrum[:,0], original_spectrum[:,1], color = "silver", label = "original")
            plt.plot(cut_spectrum[:,0], cut_spectrum[:,1], color = "tomato", label = "lazy cut")
            plt.plot(annelotte_spectrum[:,0], annelotte_spectrum[:,1], color = "royalblue", label = "Annelotte")

            plt.title(f"{epoch} {line}")
            plt.xlabel(r"Wavelength [$\AA$]")
            plt.ylabel(r"Flux")

            plt.legend()
            plt.grid()

            plt.savefig(os.path.join(output_dir, f"Correction_{epoch}_{line}.png"))

    annelotte_dates = []
    for epoch in annelotte_spectra:

        # Add the epoch date to the list
        date = float(epoch)
        annelotte_dates.append(date)

        # Calculate EW of Annelotte's corrected spectrum
        annelotte_spectrum = annelotte_spectra[epoch]
        annelotte_spectrum = annelotte_spectrum[np.where((annelotte_spectrum[:,0] > original_spectrum[0,0]) & (annelotte_spectrum[:,0] < original_spectrum[-1,0]))]
        annelotte_EW, annelotte_err = calculateEquivalentWidth(annelotte_spectrum, line_l=bound_left, line_r=bound_right)
        annelotte_widths.append(annelotte_EW)
        annelotte_errs.append(annelotte_err)


    dates = np.array(dates) - epoch_0
    annelotte_dates = np.array(annelotte_dates) - epoch_0

    plt.figure(figsize=(14,8))
    size = 5
    fmt = "s"
    plt.errorbar(dates, original_widths, yerr = original_errs, capsize = size, fmt = fmt, linewidth = 1, color = "black", label = "original")
    plt.errorbar(dates, cut_widths, yerr = cut_errs, capsize = size, fmt = fmt, linewidth = 1, color = "orangered", label = "the lazy way")
    plt.errorbar(annelotte_dates, annelotte_widths, yerr = annelotte_errs, capsize = size, fmt = fmt, linewidth = 1, color = "mediumblue", label = "corrected")

    plt.axhline(0, linestyle = "dashed", color = "black")

    plt.title(f"{line} EW over time")
    plt.xlabel(f"time since {epoch_0:.3f} [days]")
    plt.ylabel(r"EW [$\AA$]")

    plt.legend()
    plt.grid()

    plt.savefig(os.path.join(dirs["EW"], "EW_over_time.png"))
    plt.show()
