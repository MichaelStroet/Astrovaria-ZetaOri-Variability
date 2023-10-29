import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

from directories import getPaths
from plotting_functions import getObservationSpectra, getModelSpectrum
from equivalent_width import calculateEquivalentWidth

dirs = getPaths()
original_dir = os.path.join(dirs["zetaOri_Annelotte"], "originals")
annelotte_files = [file for file in os.listdir(dirs["zetaOri_Annelotte"]) if file.endswith(".txt")]

for file in annelotte_files:

    # Get the observation date
    bjd = file[:-4].split("_")[-1]

    # Check if there is a matching original spectrum
    for original_file in os.listdir(original_dir):
        if original_file.startswith(f"BJD{bjd[:-1]}"):

            # Get the original spectrum
            original_spectrum = np.genfromtxt(os.path.join(original_dir, original_file))

            # Get Annelotte's corrected spectrum
            spectrum = np.genfromtxt(os.path.join(dirs["zetaOri_Annelotte"], file))
            HALPHA = 6562.790
            C_KM = 299792.458
            spectrum[:,0] = HALPHA * (1 + spectrum[:,0] / C_KM)
            spectrum = spectrum[np.where((spectrum[:,0] >= original_spectrum[0,0]) & (spectrum[:,0] <= original_spectrum[-1,0]))]

            # Interpolate the corrected spectrum to match the original spectrum x-values
            interpolated_flux = np.interp(x = original_spectrum[:,0], xp = spectrum[:,0], fp = spectrum[:,1])
            interpolated_spectrum = np.dstack((original_spectrum[:,0], interpolated_flux))[0]

            # Calculate the difference as a new spectrum
            telluric_spectrum = np.dstack((original_spectrum[:,0], 1 + original_spectrum[:,1] - interpolated_spectrum[:,1]))[0]

            # Calculate the equivalent widths
            original_EW, original_err = calculateEquivalentWidth(original_spectrum)
            interpolated_EW, interpolated_err = calculateEquivalentWidth(interpolated_spectrum)
            telluric_EW, telluric_err = calculateEquivalentWidth(telluric_spectrum)

            plt.figure(figsize = (12,6))

            size = 2
            original_label = f"EW={original_EW:.2f} +- {original_err:.2f} A"
            plt.plot(original_spectrum[:,0], original_spectrum[:,1], linewidth = size, color = "dimgray", label=original_label)

            telluric_label = f"EW={telluric_EW:.2f} +- {telluric_err:.2f} A"
            plt.plot(telluric_spectrum[:,0], telluric_spectrum[:,1], linewidth = size, color = "forestgreen", label=telluric_label)

            interpolated_label = f"EW={interpolated_EW:.2f} +- {interpolated_err:.2f} A"
            plt.plot(interpolated_spectrum[:,0], interpolated_spectrum[:,1], linewidth = size, color = "mediumblue", label=interpolated_label)

            plt.axhline(1, linestyle = "solid", color = "black", label = "Continuum")
            plt.title(f"HALPHA {bjd}")
            plt.xlabel(r"wavelength [$\AA$]")
            plt.ylabel("Normalised flux")

            plt.ylim(bottom = 0.7, top = 1.1)
            plt.legend()
            plt.grid()

            plt.savefig(os.path.join(dirs["EW"], "telluric_spectra", f"telluric_HALPHA_{bjd}.png"))
