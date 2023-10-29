
import os
import numpy as np
import matplotlib.pyplot as plt


import matplotlib.cm as cm
from matplotlib.colors import Normalize, TwoSlopeNorm, LogNorm

from directories import getPaths
from plotting_functions import getObservationSpectra, getAnnelotteSpectra
from equivalent_width import calculateEquivalentWidth

plt.rcParams.update({'font.size': 18})
# plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
# plt.rcParams.update({'figure.autolayout': True})

# def calculateModelEW(dirs, models, lines):
#
#     model_EW = {}
#     for model in models:
#         model_dir = os.path.join(dirs["grid"], model)
#         model_EW[model] = {}
#         for line in lines:
#             spectrum = np.genfromtxt(os.path.join(model_dir, f"{line}.dat"))
#             model_EW[model][line], EW_error = calculateEquivalentWidth(spectrum)
#
#     return model_EW

def getGridModels(grid_dir, line):
    """
    Gets the line spectra of all models in the grid.
    Returns spectra in a dictionary and the unique Mdot and CLF values in arrays.
    """

    # Name of the grid model dictionaries
    models = [file for file in os.listdir(grid_dir) if file.startswith("grid")]

    spectra = {}
    Mdots = []
    CLFs = []

    # Loop over each model and get the spectrum and Mdot/CLF values
    for i, model in enumerate(models):
        spectrum = np.genfromtxt(os.path.join(grid_dir, model, f"{line}.dat"))

        Mdot = float(model.split("_")[1][1:].replace(",", "."))
        if not Mdot in Mdots:
            Mdots.append(Mdot)

        CLF = float(model.split("_")[2][2:].replace(",", "."))
        if not CLF in CLFs:
            CLFs.append(CLF)

        spectra[f"{Mdot}_{CLF}"] = spectrum

    return spectra, np.sort(Mdots), np.sort(CLFs)

def calcChiSquared(observed, expected):
    """
    Calculates The Chi^2 goodness-of-fit statistic
    observed: The model spectrum
    expected: The real spectrum
    """
    return np.sum((observed - expected)**2 / expected)

def createChiSquaredHeatmap(observed_spectrum, model_spectra, Mdots, CLFs):

    # Create empty heatmap
    heatmap = np.zeros((len(Mdots), len(CLFs)))

    # Loop over each model
    best_chi2 = np.inf
    for i, Mdot in enumerate(Mdots):
        for j, CLF in enumerate(CLFs):
            # Get and interpolate the model spectrum
            model = model_spectra[f"{Mdot}_{CLF}"]
            interpolated_model = np.interp(observed_spectrum[:,0], model[:,0], model[:,1])
            new_model = np.dstack((observed_spectrum[:,0], interpolated_model))[0]

            # Calculate the chi squared statistic
            chi2 = calcChiSquared(new_model, observed_spectrum)

            # Add the chi^2 to the heatmap
            heatmap[i,j] = chi2

            if chi2 < best_chi2:
                best_chi2 = chi2
                best_Mdot = Mdot
                best_CLF = CLF

    return heatmap, [best_chi2, best_Mdot, best_CLF]

def plotHeatmap(dirs, heatmap, Mdots, factors, title = ""):

    norm = Normalize(vmin = np.min(heatmap), vmax = 5)

    plt.figure()
    plt.imshow(heatmap, norm = norm, cmap = "binary")
    plt.colorbar()
    # plt.colorbar(label = r"$\Delta$$W_{\rm eq}$ ($W_{model}$ - $W_{obs}$) [$\AA$]")

    plt.title(title)
    plt.yticks(range(len(Mdots)), Mdots)
    plt.ylabel(r"$\log_{10}$($\dot{M}$) [$\log_{10}$(M$\odot$/yr)]")

    plt.xticks(range(len(factors)), factors)
    plt.xlabel(r"Clumping factor")

    # plt.savefig(os.path.join(dirs["grid"], f"EW_{line}.pdf"))



if __name__ == "__main__":



    dirs = getPaths()
    dirs["grid"] = os.path.join(dirs["results"], "05 - grid_M-5,5_-6,5_Cf1,1_100")
    output_dir = os.path.join(dirs["results"], "fitting_grid")

    # Get the observed rv-corrected Halpha spectra for all epochs
    # observed_spectra = getObservationSpectra(dirs["zetaOri_rv"])5
    observed_spectra = getAnnelotteSpectra(dirs["zetaOri_Annelotte"])
    stellar_abs = np.genfromtxt(os.path.join(dirs["results"], "04 - bunch of random models", "single_M-9_T31_Cf1","HALPHA.dat"))

    # Get the grid model spectra and Mdot and CLF values
    line = "HALPHA"
    bound_left = 6540
    bound_right = 6590
    model_spectra, Mdots, CLFs = getGridModels(dirs["grid"], line)

    all_heatmaps = {}
    best_fits = {}
    best_fits_list = []
    observed_EWs = []
    observed_errs = []
    absorption_spectra = {}
    for epoch in observed_spectra:
        observed_spectrum = observed_spectra[epoch]

        heatmap, best_fit = createChiSquaredHeatmap(observed_spectrum, model_spectra, Mdots, CLFs)
        all_heatmaps[epoch] = heatmap
        best_fits[epoch] = best_fit
        best_fits_list.append(best_fit)

        # Calculate the epoch's equivalent width
        obs_EW, obs_err = calculateEquivalentWidth(observed_spectrum, line_l=bound_left, line_r=bound_right)
        observed_EWs.append(obs_EW)
        observed_errs.append(obs_err)

        absorption_spectrum = np.interp(observed_spectrum[:,0], stellar_abs[:,0], stellar_abs[:,1])
        absorption_spectra[epoch] = absorption_spectrum

    # --------------------------------------------------------------------------
    # Chi^2 heatmaps
    # --------------------------------------------------------------------------
    fig, axes = plt.subplots(3, 8, figsize = (18,7), sharex=True, sharey=True, layout = "constrained")
    # fig.suptitle(r"$\chi^{2}$ of grid models")


    norm = Normalize(vmin = np.min(heatmap), vmax = 1)
    for ax, epoch in zip(fig.axes, observed_spectra):
        heatmap = all_heatmaps[epoch]
        img = ax.imshow(heatmap, norm = norm, cmap = "magma_r")
        ax.set_title(epoch, fontsize = 13)
        ax.set_xticks(range(len(CLFs))[::2],CLFs[::2], fontsize = 11, rotation = 45)
        ax.set_yticks(range(len(Mdots))[::2], Mdots[::2], fontsize = 11)

    for ax in fig.axes[len(all_heatmaps):]:
        ax.axis("off")

    fig.supxlabel(r"$f_{\rm cl}$")
    fig.supylabel(r"log($\dot{M}$) [M$\odot$/yr]")
    fig.colorbar(img, ax=fig.axes, pad = 0.02, label = r"$\chi^2$")

    plt.savefig(os.path.join(output_dir, "chi2_heatmaps.pdf"), dpi=300)
    plt.close()

    # --------------------------------------------------------------------------
    # Best fit spectra
    # --------------------------------------------------------------------------
    fig, axes = plt.subplots(3, 8, figsize = (16,6), sharex=True, sharey=True, constrained_layout=True)
    # fig.suptitle(f"'Best' fits of Halpha")
    fig.supxlabel(r"Wavelength [$\AA$]")
    fig.supylabel(r"Normalised Flux")

    epochs = []
    best_Mdots = []
    best_CLFs = []
    for ax, epoch in zip(fig.axes, observed_spectra):
        chi2, Mdot, CLF = best_fits[epoch]
        epochs.append(float(epoch))
        best_Mdots.append(Mdot)
        best_CLFs.append(CLF)

        observed_spectrum = observed_spectra[epoch]
        model_spectrum = model_spectra[f"{Mdot}_{CLF}"]

        ax.plot(observed_spectrum[:,0], observed_spectrum[:,1], color = "black")
        ax.plot(model_spectrum[:,0], model_spectrum[:,1], color = "red")
        ax.axhline(1, linestyle = "dashed", color = "blue")
        ax.set_title(f"{epoch}", fontsize = 16)
        ax.set_xlim(left = bound_left, right = bound_right)
        ax.grid()

    for ax in fig.axes[len(all_heatmaps):]:
        ax.axis("off")

    plt.savefig(os.path.join(output_dir, "best_fits.pdf"))
    plt.close()

    # --------------------------------------------------------------------------
    # spectra with absorption subtracted
    # --------------------------------------------------------------------------
    fig, axes = plt.subplots(4, 6, figsize = (12,6), sharex=True, sharey=True, constrained_layout=True)
    # fig.suptitle(f"'Best' fits of Halpha")
    fig.supxlabel(r"Wavelength [$\AA$]")
    fig.supylabel(r"Normalised Flux")

    abs_spectra = {}
    for ax, epoch in zip(fig.axes, observed_spectra):

        observed_spectrum = observed_spectra[epoch]
        ax.plot(observed_spectrum[:,0], observed_spectrum[:,1], color = "gray")

        absorption_spectrum = absorption_spectra[epoch]
        ax.plot(observed_spectrum[:,0], observed_spectrum[:,1] - absorption_spectrum + 1, color = "black")
        # ax.plot(observed_spectrum[:,0], absorption_spectrum, color = "red")

        # chi2, Mdot, CLF = best_fits[epoch]
        # model_spectrum = model_spectra[f"{Mdot}_{CLF}"]
        # absorption_spectrum = np.interp(model_spectrum[:,0], stellar_abs[:,0], stellar_abs[:,1])
        # ax.plot(model_spectrum[:,0], model_spectrum[:,1] - absorption_spectrum + 1, color = "red")

        ax.axhline(1, linestyle = "dashed", color = "blue")
        ax.axvline(6562.790, linestyle = "dashed", linewidth = 1, color = "black", label = r"H$\alpha\ \lambda$" + f"{6562.790:.3f}")
        ax.set_title(epoch, fontsize = 16)
        ax.set_xlim(left = bound_left, right = bound_right)
        ax.grid()

    for ax in fig.axes[len(all_heatmaps):]:
        ax.axis("off")

    plt.savefig(os.path.join(output_dir, "spectra_absorption_subtracted.pdf"))
    plt.close()

    # --------------------------------------------------------------------------
    # spectra with absorption subtracted WITH MODEL
    # --------------------------------------------------------------------------
    # fig, axes = plt.subplots(4,6, figsize = (12,6), sharex=True, sharey=True, constrained_layout=True)
    fig, axes = plt.subplots(3, 8, figsize = (16,6), sharex=True, sharey=True, constrained_layout=True)
    # fig.suptitle(f"'Best' fits of Halpha")
    fig.supxlabel(r"Wavelength [$\AA$]")
    fig.supylabel(r"Normalised Flux")


    for ax, epoch in zip(fig.axes, observed_spectra):

        observed_spectrum = observed_spectra[epoch]
        ax.plot(observed_spectrum[:,0], observed_spectrum[:,1], color = "gray")

        absorption_spectrum = absorption_spectra[epoch]
        ax.plot(observed_spectrum[:,0], observed_spectrum[:,1] - absorption_spectrum + 1, color = "black")
        # ax.plot(observed_spectrum[:,0], absorption_spectrum, color = "red")

        chi2, Mdot, CLF = best_fits[epoch]
        model_spectrum = model_spectra[f"{Mdot}_{CLF}"]
        absorption_spectrum = np.interp(model_spectrum[:,0], stellar_abs[:,0], stellar_abs[:,1])

        ax.plot(model_spectrum[:,0], model_spectrum[:,1] - absorption_spectrum + 1, color = "red")

        ax.axhline(1, linestyle = "dashed", color = "blue")
        ax.axvline(6562.790, linestyle = "dashed", linewidth = 1, color = "black", label = r"H$\alpha\ \lambda$" + f"{6562.790:.3f}")

        absorption_spectrum = np.interp(observed_spectrum[:,0], stellar_abs[:,0], stellar_abs[:,1])
        corrected_fluxes = observed_spectrum[:,1] - absorption_spectrum + 1
        corrected_spectrum = np.dstack((observed_spectrum[:,0], corrected_fluxes))[0]
        obs_EW, obs_err = calculateEquivalentWidth(corrected_spectrum, line_l=bound_left, line_r=bound_right)

        absorption_spectrum = np.interp(model_spectrum[:,0], stellar_abs[:,0], stellar_abs[:,1])
        corrected_fluxes = model_spectrum[:,1] - absorption_spectrum + 1
        corrected_spectrum = np.dstack((model_spectrum[:,0], corrected_fluxes))[0]
        mod_EW, mod_err = calculateEquivalentWidth(corrected_spectrum, line_l=bound_left, line_r=bound_right)

        # ax.set_title(f"{epoch}\nWobs: {obs_EW:.1f}\nWmod: {mod_EW:.1f}", fontsize = 12)
        ax.set_title(f"{epoch}", fontsize = 16)
        ax.set_xlim(left = bound_left, right = bound_right)
        ax.grid()

    for ax in fig.axes[len(all_heatmaps):]:
        ax.axis("off")

    plt.savefig(os.path.join(output_dir, "spectra_absorption_subtracted_with_model.pdf"))
    # plt.close()

    # --------------------------------------------------------------------------
    # Best fit values over time
    # --------------------------------------------------------------------------
    epoch_0 = min(epochs)
    epochs_0 = np.array(epochs) - epoch_0
    best_Mdots = np.array(best_Mdots)
    best_CLFs = np.array(best_CLFs)
    Mdot_sqrtCLF = 10**best_Mdots * np.sqrt(best_CLFs)

    fig, (axM, axF, axMF) = plt.subplots(3, 1, figsize = (12,6), sharex=True, constrained_layout=True)

    axM.plot(epochs_0, best_Mdots, linewidth=2, marker="o", markersize=8, color = "black")
    axM.set_ylim(bottom=-6.55, top=-5.95)
    # axM.set_title("Best fit values over time")
    axM.set_ylabel(r"log$(\dot{M})$ [M$\odot$/yr]",fontsize=16)
    axM.grid()

    axF.plot(epochs_0, best_CLFs, linewidth=2, marker="o", markersize=8, color = "black")
    axF.set_ylim(bottom=0, top=70)
    axF.set_ylabel(r"$f_{\rm cl}$")
    axF.grid()

    axMF.plot(epochs_0, np.log10(Mdot_sqrtCLF), linewidth=2, marker="o", markersize=8, color = "black")
    # axMF.set_ylim(bottom=1.95e-6, top=2.65e-6)
    axMF.set_xlabel(f"Time since first epoch [days]")
    axMF.set_ylabel(r"$\dot{M}$ * $\sqrt{f_{\rm cl}}$")
    axMF.grid()

    plt.savefig(os.path.join(output_dir, "fits_over_time.pdf"))
    plt.close()

    # --------------------------------------------------------------------------
    # Equivalent widths best fits over time
    # --------------------------------------------------------------------------

    plt.figure(figsize = (15,4.5), constrained_layout=True)

    model_EWs = []
    model_errs = []
    for epoch in observed_spectra:
        chi2, Mdot, CLF = best_fits[epoch]

        model_spectrum = model_spectra[f"{Mdot}_{CLF}"]
        model_EW, model_err = calculateEquivalentWidth(model_spectrum, line_l=bound_left, line_r=bound_right)

        model_EWs.append(model_EW)
        model_errs.append(model_err)

    size = 5
    fmt = "s"
    plt.errorbar(epochs_0, observed_EWs, yerr = observed_errs, capsize = size, fmt = fmt, linewidth = 1, color = "black", label = "Observed")
    # plt.errorbar(epochs_0, model_EWs, yerr = model_errs, capsize = size, fmt = fmt, linewidth = 1, color = "red", label = "Best fit")

    plt.axhline(0, linestyle = "dashed", color = "black")

    plt.xlabel(f"Time since first epoch [days]")
    plt.ylabel(r"$W_{\rm eq}$ [$\AA$]")
    plt.grid()

    plt.savefig(os.path.join(output_dir, "EW_over_time.pdf"))
    # plt.close()
    # --------------------------------------------------------------------------
    # Equivalent widths ABSORPTION CORRECTED best fits over time
    # --------------------------------------------------------------------------


    emisison_EWs = []
    emission_errs = []
    model_EWs = []
    model_errs = []
    for epoch in observed_spectra:
        chi2, Mdot, CLF = best_fits[epoch]

        observed_spectrum = observed_spectra[epoch]
        absorption_spectrum = np.interp(observed_spectrum[:,0], stellar_abs[:,0], stellar_abs[:,1])
        corrected_fluxes = observed_spectrum[:,1] - absorption_spectrum + 1
        corrected_spectrum = np.dstack((observed_spectrum[:,0], corrected_fluxes))[0]
        emission_EW, emission_err = calculateEquivalentWidth(corrected_spectrum, line_l=bound_left, line_r=bound_right)
        emisison_EWs.append(-emission_EW)
        emission_errs.append(emission_err)

    MdotCLFs = 10**best_Mdots * best_CLFs**0.5
    model_Mdots = 10**np.arange(-6.55, -5.95, 0.01)
    model_CLFs = np.linspace(1, 65, 1000)
    model_MdotCLFs = np.linspace(2.0e-6, 2.7e-6, 1000)

    Nactive = 5
    Mdots_quiet = 10**best_Mdots[:-Nactive]
    Mdots_active = 10**best_Mdots[-Nactive:]
    CLFs_quiet = best_CLFs[:-Nactive]
    CLFs_active = best_CLFs[-Nactive:]
    MdotCLFs_quiet = MdotCLFs[:-Nactive]
    MdotCLFs_active = MdotCLFs[-Nactive:]

    EW_quiet = emisison_EWs[:-Nactive]
    EW_active = emisison_EWs[-Nactive:]
    EW_err_quiet = emission_errs[:-Nactive]
    EW_err_active = emission_errs[-Nactive:]

    from scipy.optimize import curve_fit

    def powerlaw(x, alpha, m):
        return m*x**alpha



    plt.figure(figsize = (15,6), constrained_layout=True)

    q_fit, q_err = curve_fit(powerlaw, Mdots_quiet, EW_quiet, sigma=EW_err_quiet)#, p0 = [1, 1e7])
    q_model = powerlaw(model_Mdots, *q_fit)
    print("\nquiet Mdots")
    print(f"α: {q_fit[0]:.3f} ± {np.sqrt(np.diag(q_err))[0]:.3f}")
    print(f"m: {q_fit[1]:.3f} ± {np.sqrt(np.diag(q_err))[1]:.3f}")

    plt.errorbar(Mdots_quiet, EW_quiet, yerr=EW_err_quiet, color="black", linestyle="None", fmt="s", capsize=5, label="Quiet")
    plt.plot(model_Mdots, q_model, color="black", linestyle="solid", label = f"α: {q_fit[0]:.3f} ± {np.sqrt(np.diag(q_err))[0]:.3f}")

    a_fit, a_err = curve_fit(powerlaw, Mdots_active, EW_active, sigma=EW_err_active)#, p0 = [1, 1e7])
    a_model = powerlaw(model_Mdots, *a_fit)
    print("\nactive Mdots")
    print(f"α: {a_fit[0]:.3f} ± {np.sqrt(np.diag(a_err))[0]:.3f}")
    print(f"m: {a_fit[1]:.3f} ± {np.sqrt(np.diag(a_err))[1]:.3f}")

    plt.errorbar(Mdots_active, EW_active, yerr=EW_err_active ,color="red", linestyle="None", fmt="s", capsize=5, label="Active")
    plt.plot(model_Mdots, a_model, color="red", linestyle="solid", label = f"α: {a_fit[0]:.3f} ± {np.sqrt(np.diag(a_err))[0]:.3f}")

    fit, err = curve_fit(powerlaw, 10**best_Mdots, emisison_EWs, sigma=emission_errs)#, p0 = [1, 1e7])
    model = powerlaw(model_Mdots, *fit)
    print("\nMdots_all")
    print(f"α: {fit[0]:.3f} ± {np.sqrt(np.diag(err))[0]:.3f}")
    print(f"m: {fit[1]:.3f} ± {np.sqrt(np.diag(err))[1]:.3f}")

    plt.plot(model_Mdots, model, color="blue", linestyle="solid", label = f"α: {fit[0]:.3f} ± {np.sqrt(np.diag(err))[0]:.3f}")

    plt.xlabel(r"$\dot{M}$ [M$\odot$/yr]")
    plt.ylabel(r"$W_{\rm eq,emission}$ [$\AA$]")
    plt.grid()
    plt.legend()

    plt.savefig(os.path.join(output_dir, "Alpha_fitting_Mdot.pdf"))
    plt.close()


    plt.figure(figsize = (15,6), constrained_layout=True)

    q_fit, q_err = curve_fit(powerlaw, CLFs_quiet, EW_quiet, sigma=EW_err_quiet)#, p0 = [1, 1e7])
    q_model = powerlaw(model_CLFs, *q_fit)
    print("\nquiet CLFs")
    print(f"α: {q_fit[0]:.3f} ± {np.sqrt(np.diag(q_err))[0]:.3f}")
    print(f"m: {q_fit[1]:.3f} ± {np.sqrt(np.diag(q_err))[1]:.3f}")

    plt.errorbar(CLFs_quiet, EW_quiet, yerr=EW_err_quiet, color="black", linestyle="None", fmt="s", capsize=5, label="Quiet")
    plt.plot(model_CLFs, q_model, color="black", linestyle="solid", label = f"α: {q_fit[0]:.3f} ± {np.sqrt(np.diag(q_err))[0]:.3f}")

    a_fit, a_err = curve_fit(powerlaw, CLFs_active, EW_active, sigma=EW_err_active)#, p0 = [1, 1e7])
    a_model = powerlaw(model_CLFs, *a_fit)
    print("\nactive CLFs")
    print(f"α: {a_fit[0]:.3f} ± {np.sqrt(np.diag(a_err))[0]:.3f}")
    print(f"m: {a_fit[1]:.3f} ± {np.sqrt(np.diag(a_err))[1]:.3f}")

    plt.errorbar(CLFs_active, EW_active, yerr=EW_err_active ,color="red", linestyle="None", fmt="s", capsize=5, label="Active")
    plt.plot(model_CLFs, a_model, color="red", linestyle="solid", label = f"α: {a_fit[0]:.3f} ± {np.sqrt(np.diag(a_err))[0]:.3f}")

    fit, err = curve_fit(powerlaw, best_CLFs, emisison_EWs, sigma=emission_errs)#, p0 = [1, 1e7])
    model = powerlaw(model_CLFs, *fit)
    print("\nCLFs_all")
    print(f"α: {fit[0]:.3f} ± {np.sqrt(np.diag(err))[0]:.3f}")
    print(f"m: {fit[1]:.3f} ± {np.sqrt(np.diag(err))[1]:.3f}")
    plt.plot(model_CLFs, model, color="blue", linestyle="solid", label = f"α: {fit[0]:.3f} ± {np.sqrt(np.diag(err))[0]:.3f}")

    plt.xlabel(r"$f_{\rm cl}$")
    plt.ylabel(r"$W_{\rm eq,emission}$ [$\AA$]")
    plt.grid()
    plt.legend()

    plt.savefig(os.path.join(output_dir, "Alpha_fitting_Fcl.pdf"))
    plt.close()


    plt.figure(figsize = (15,5), constrained_layout=True)

    model_q = np.linspace(2e-6, 2.1e-6, 100)
    q_fit, q_err = curve_fit(powerlaw, MdotCLFs_quiet, EW_quiet, sigma=EW_err_quiet, p0 = [8, 1e46], bounds=((5,1e45),(10,1e80)))
    q_model = powerlaw(model_q, *q_fit)
    print("\nquiet MdotCLFs")
    print(f"α: {q_fit[0]:.3f} ± {np.sqrt(np.diag(q_err))[0]:.3f}")
    print(f"m: {q_fit[1]:.3f} ± {np.sqrt(np.diag(q_err))[1]:.3f}")

    plt.errorbar(MdotCLFs_quiet, EW_quiet, yerr=EW_err_quiet, color="black", linestyle="None", fmt="s", capsize=5, label="Quiet")
    # plt.plot(model_q, q_model, color="black", linestyle="solid", label = f"α: {q_fit[0]:.3f} ± {np.sqrt(np.diag(q_err))[0]:.3f}")

    model_a = np.linspace(2.5e-6, 2.6e-6, 100)
    a_fit, a_err = curve_fit(powerlaw, MdotCLFs_active, EW_active, sigma=EW_err_active, p0 = [8, 1e45], bounds=((5,1e45),(10,1e80)))
    a_model = powerlaw(model_a, *a_fit)
    print("\nactive MdotCLFs")
    print(f"α: {a_fit[0]:.3f} ± {np.sqrt(np.diag(a_err))[0]:.3f}")
    print(f"m: {a_fit[1]:.3f} ± {np.sqrt(np.diag(a_err))[1]:.3f}")

    plt.errorbar(MdotCLFs_active, EW_active, yerr=EW_err_active ,color="red", linestyle="None", fmt="s", capsize=5, label="Active")
    # plt.plot(model_a, a_model, color="red", linestyle="solid", label = f"α: {a_fit[0]:.3f} ± {np.sqrt(np.diag(a_err))[0]:.3f}")

    fit, err = curve_fit(powerlaw, MdotCLFs, emisison_EWs, sigma=emission_errs, p0 = [1, 1e7])
    model = powerlaw(model_MdotCLFs, *fit)
    print("\nMdotCLFs_all")
    print(f"α: {fit[0]:.3f} ± {np.sqrt(np.diag(err))[0]:.3f}")
    print(f"m: {fit[1]:.3f} ± {np.sqrt(np.diag(err))[1]:.3f}")
    plt.plot(model_MdotCLFs, model, color="blue", linestyle="solid", label = f"α: {fit[0]:.3f} ± {np.sqrt(np.diag(err))[0]:.3f}")

    plt.xlabel(r"$\dot{M}\sqrt{f_{\rm cl}}$")
    plt.ylabel(r"$W_{\rm eq,emission}$ [$\AA$]")
    plt.grid()
    plt.legend(loc = "lower right")
    plt.savefig(os.path.join(output_dir, "Alpha_fitting_MdotFcl.pdf"))
    plt.close()

    # --------------------------------------------------------------------------
    # Jump zoom in
    # --------------------------------------------------------------------------

    fig, axes = plt.subplots(ncols=3, nrows=2, sharey = "row", figsize = (12,5), gridspec_kw={"height_ratios": [3,5]}, constrained_layout=True)
    gs = axes[0,0].get_gridspec()

    # remove the underlying axes
    for ax in axes[0,:]:
        ax.remove()

    # Add big subplot
    axbig = fig.add_subplot(gs[0,:])
    axbig.plot(epochs_0, Mdot_sqrtCLF, marker="s", color = "black")
    axbig.set_xlabel(f"Time since first epoch [days]")
    axbig.set_ylabel(r"$\dot{M}$ * $\sqrt{f_{\rm cl}}$")

    axbig.set_xlim(left = 4.7)
    axbig.grid()

    epoch_ids = [-7,-5,-2]
    for id in epoch_ids:
        axbig.plot(epochs_0[id], 10**best_Mdots[id] * np.sqrt(best_CLFs)[id], marker="s", color = "red", ms = 15, markeredgewidth = 3, fillstyle = "none")
        axbig.set_ylim(bottom=1.95e-6, top=2.65e-6)
    axes[1,0].set_ylabel("Normalised Flux")
    fig.supxlabel(r"Wavelength [$\AA$]")
    for ax, id in zip(axes[1,:], epoch_ids):
        epoch = f"{epochs[id]}"
        chi2, Mdot, CLF = best_fits[epoch]

        observed_spectrum = observed_spectra[epoch]
        model_spectrum = model_spectra[f"{Mdot}_{CLF}"]

        ax.plot(observed_spectrum[:,0], observed_spectrum[:,1], color = "black")
        ax.plot(model_spectrum[:,0], model_spectrum[:,1], color = "red")
        ax.axhline(1, linestyle = "dashed", color = "blue")
        ax.set_xlim(left = bound_left, right = bound_right)
        ax.grid()

    plt.savefig(os.path.join(output_dir, "fits_over_time_zoomed.pdf"))
    plt.close()

    plt.show()
