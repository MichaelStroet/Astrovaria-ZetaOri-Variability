import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

HALPHA = 6562.790
C_KM = 299792.458
files = os.listdir(os.getcwd())

original_dir = "originals"
plot_dir = "plots"
all_spectra = []
for file in files:
    if file.endswith(".txt"):# and file.startswith("zetaOri_Ha"):

        spectrum = np.genfromtxt(file)
        spectrum[:,0] = HALPHA * (1 + spectrum[:,0] / C_KM)
        all_spectra.append(spectrum)
        bjd = file[:-4].split("_")[-1]

        # plt.figure(figsize = (12,6))
        #
        # for original_file in os.listdir(original_dir):
        #     if original_file.startswith(f"BJD{bjd[:-1]}"):
        #         original_spectrum = np.genfromtxt(os.path.join(original_dir, original_file))
        #         plt.plot(original_spectrum[:,0], original_spectrum[:,1], color = "red", alpha = 0.7, label = "Original")
        #
        # plt.plot(spectrum[:,0], spectrum[:,1], color = "black", label = "Annelotte")
        # plt.axhline(1, linestyle = "dashed", color = "black", label = "Continuum")
        # plt.title(f"HALPHA {bjd}")
        # plt.xlabel(r"wavelength [$\AA$]")
        # plt.ylabel("Normalised flux")
        #
        # plt.ylim(bottom = 0.65, top = 1.15)
        # plt.legend()
        # plt.grid()
        #
        # plt.tight_layout()
        # plt.savefig(os.path.join(plot_dir, file[:-4] + ".pdf"))
        #
        # plt.close()


plt.figure(figsize = (12,6))
all_spectra = np.array(all_spectra)
average_flux = np.mean(all_spectra[:,:,1], axis = 0)
average_spectrum = np.dstack((all_spectra[0][:,0], average_flux))[0]
np.savetxt("AVERAGE_HALPHA_ANNELOTTE.dat", average_spectrum)

for spectrum in all_spectra:
    plt.plot(spectrum[:,0], spectrum[:,1], linewidth = 3)


plt.axvline(HALPHA, linestyle = "dashed", linewidth = 3, color = "black", label = r"H$\alpha\ \lambda$" + f"{HALPHA:.3f}")
# plt.plot(spectrum[:,0], average_flux, linewidth = 3, color = "black", label = "Average")
plt.axhline(1, linestyle = "solid", color = "black")

# plt.title(f"HALPHA")
plt.xlabel(r"Wavelength [$\AA$]")
plt.ylabel("Normalised flux")
plt.xlim(left = 6540, right = 6590)

plt.legend()
plt.grid()

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "HALPHA_spectra.pdf"), dpi=300)


plt.show()
