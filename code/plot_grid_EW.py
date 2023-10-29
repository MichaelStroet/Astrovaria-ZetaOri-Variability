
import os
import numpy as np
import matplotlib.pyplot as plt

# plt.rcParams['font.size'] = 18

# plt.rcParams.update({'figure.autolayout': True})

import matplotlib.cm as cm
from matplotlib.colors import Normalize, TwoSlopeNorm, LogNorm

from directories import getPaths
from plotting_functions import getObservationSpectra
from equivalent_width import calculateEquivalentWidth

plt.rcParams.update({'font.size': 18})
# def calculateObservedEW(observed_lines):
#     """
#     Calculates the equivalent width + error of the observed spectra
#     """
#     observed_EW = {}
#     for line in observed_lines:
#         observed_EW[line] = {}
#         for epoch in observed_lines[line]:
#             observed_EW[line][epoch] = {}
#
#             EW, error = calculateEquivalentWidth(observed_lines[line][epoch])
#             observed_EW[line][epoch]["value"] = EW
#             observed_EW[line][epoch]["error"] = error
#
#     return observed_EW

def calculateModelEW(dirs, models, lines):

    model_EW = {}
    for model in models:
        model_dir = os.path.join(dirs["grid"], model)
        model_EW[model] = {}
        for line in lines:
            spectrum = np.genfromtxt(os.path.join(model_dir, f"{line}.dat"))
            model_EW[model][line], EW_error = calculateEquivalentWidth(spectrum)

    return model_EW

def getGridModels(grid_dir):

    models = [file for file in os.listdir(grid_dir) if file.startswith("grid")]
    values = np.zeros((len(models), 2))
    for i, model in enumerate(models):
        for split in model.split("_"):
            if split.startswith("M"):
                values[i,0] = float(split[1:].replace(",", "."))
            elif split.startswith("Cf"):
                values[i,1] = float(split[2:].replace(",", "."))

    return models, values

def createHeatmapObsModel(obs_eq_width, model_EW, Mdots, factors, line):
    heatmap = np.zeros((len(Mdots), len(factors)))
    for i, Mdot in enumerate(Mdots):
        for j, CLF in enumerate(factors):
            model = f"grid_M{Mdot:.1f}_Cf{CLF:.2f}".replace(".", ",")
            model_eq_width = model_EW[model][line]
            heatmap[i,j] = model_eq_width - obs_eq_width

    return heatmap

def createHeatmapModel(model_EW, Mdots, factors, line):
    heatmap = np.zeros((len(Mdots), len(factors)))
    for i, Mdot in enumerate(Mdots):
        for j, CLF in enumerate(factors):
            model = f"grid_M{Mdot:.1f}_Cf{CLF:.2f}".replace(".", ",")
            heatmap[i,j] = model_EW[model][line]

    return heatmap

def plotHeatmap(dirs, heatmap, Mdots, factors, title = "", cbar_label = "", two_norm = False):

            if two_norm:
                # min_distance = np.min(np.abs([np.min(heatmap), np.max(heatmap)]))
                # norm = TwoSlopeNorm(vcenter = 0, vmin = -min_distance, vmax = +min_distance)
                norm = TwoSlopeNorm(vcenter = -0.75, vmin = -2.5, vmax = 1.0)
            else:
                norm = Normalize(vmin = -np.max(heatmap), vmax = np.max(heatmap))
            plt.figure(figsize = (9,7), constrained_layout=True)
            plt.imshow(heatmap, norm = norm, cmap = "seismic")
            plt.colorbar(label = cbar_label)
            # plt.colorbar(label = r"$\Delta$$W_{\rm eq}$ (W$_{model}$ - W$_{obs}$) [$\AA$]")

            # plt.title(title)
            plt.yticks(range(len(Mdots)), Mdots, fontsize = 16)
            plt.ylabel(r"$\log$($\dot{M}$) [M$\odot$/yr]", fontsize = 22)

            plt.xticks(range(len(factors)), factors, fontsize = 16, rotation = 25)
            plt.xlabel(r"$f_{\rm cl}$", fontsize = 22)

            plt.savefig(os.path.join(dirs["grid"], f"Grid_EWs_{line}.pdf"), dpi=300)

def plot1DParameter(dirs, line, parameter, values, opposite_values, scalarmap, EW_min, EW_max, heatmap):

    fig, axes = plt.subplots(2, 1, figsize=(13,6), sharex = True, gridspec_kw={"height_ratios": [7,3]}, constrained_layout=True)
    ax, ax_zoom = axes
    for i, par in enumerate(values):
        if parameter == "CLF":
            EWs = heatmap[i,:]
        elif parameter == "Mdot":
            EWs = heatmap[:,i]
        color = scalarmap.to_rgba(par)
        ax.plot(opposite_values, EWs, linewidth = 1, marker = "s", ms = 5, color=color, label = f"{par}")
        ax_zoom.plot(opposite_values, EWs, linewidth = 3, marker = "s", ms = 10, color=color)

    ax.axhline(EW_max, linestyle = "dashed", color = "black")#, label = "observed")
    ax.axhline(EW_min, linestyle = "dashed", color = "black")
    ax_zoom.axhline(EW_max, linestyle = "dashed", color = "black")
    ax_zoom.axhline(EW_min, linestyle = "dashed", color = "black")
    ax_zoom.axvline(-6.22, linestyle = "solid", linewidth = 3, color = "black")
    ax_zoom.axvline(-6.06, linestyle = "solid", linewidth = 3, color = "black")

    # ax.set_title(f"{line} equivalent widths for varying {parameter}")
    if parameter == "CLF":
        ax_zoom.set_xlabel(r"log($\dot{M}$) [M$\odot$/yr]")
    elif parameter == "Mdot":
        ax_zoom.set_xlabel(r"Clumping factor")
    ax.set_ylabel(r"$W_{\rm eq}$ [$\AA$]")
    ax_zoom.set_ylabel(r"$W_{\rm eq}$ [$\AA$]")

    EW_range = EW_max - EW_min
    # ax_zoom.set_ylim(bottom = EW_min - 0.25*EW_range, top = EW_max + 0.25*EW_range)
    ax_zoom.set_ylim(bottom = EW_min - 0.5, top = EW_max + 0.5)
    ax.legend(fontsize = 12)
    ax.grid()
    ax_zoom.grid()

    # fig.colorbar(scalarmap, ax=axes, ticks = values, label = parameter)
    plt.savefig(os.path.join(dirs["grid"], f"Constant_{parameter}_{line}.pdf"))


if __name__ == "__main__":

    dirs = getPaths()
    dirs["grid"] = os.path.join(dirs["results"], "05 - grid_M-5,5_-6,5_Cf1,1_100")
    lines = ["HALPHA"]

    # # Get the observed rv-corrected spectra
    # observed_lines = getObservationSpectra(dirs["zetaOri_rv"], lines=lines)
    #
    # # Calculate the observed EW
    # observed_EW = calculateObservedEW(observed_lines)

    # Get the grid model files and Mdot-CLF values
    model_names, model_vals = getGridModels(dirs["grid"])
    Mdots = np.unique(model_vals[:,0])
    factors = np.unique(model_vals[:,1])

    # Calculate the model EW
    model_EW = calculateModelEW(dirs, model_names, lines)

    for line in lines:
        heatmap = createHeatmapModel(model_EW, Mdots, factors, line)
        plotHeatmap(dirs, heatmap, Mdots, factors, title=f"Equivalent widths of {line} models", cbar_label=r"$W_{\rm eq}$ [$\AA$]", two_norm = True)

        # Define the observed EW limits
        EW_max = 0.0
        EW_min = -1.5

        # Define a color map for the CLF plot
        norm_CLF = LogNorm(vmin = min(factors), vmax = max(factors))
        scalarmap_CLF = cm.ScalarMappable(norm=norm_CLF, cmap="plasma")
        # scalarmap_CLF._A = []

        # Define a color map for the Mdot plot
        norm_Mdot = Normalize(vmin = min(Mdots), vmax = max(Mdots))
        scalarmap_Mdot = cm.ScalarMappable(norm=norm_Mdot, cmap="plasma")
        scalarmap_Mdot._A = []

        # Plot the 1D parameter graphs
        plot1DParameter(dirs, line, "Mdot", Mdots, factors, scalarmap_Mdot, EW_min, EW_max, heatmap)
        plot1DParameter(dirs, line, "CLF", factors, Mdots, scalarmap_CLF, EW_min, EW_max, heatmap)

        for i, Mdot in enumerate(Mdots):
            for j, CLF in enumerate(factors):
                model = f"grid_M{Mdot:.1f}_Cf{CLF:.2f}".replace(".", ",")
                heatmap[i,j] = model_EW[model][line]

        # Determine masslosses between EW maxima
        CLF_table = {}
        CLF_EW = {}
        for CLF in factors:
            CLF_table[f"{CLF}"] = []
            CLF_EW[f"{CLF}"] = []
            for Mdot in Mdots:
                model = f"grid_M{Mdot:.1f}_Cf{CLF:.2f}".replace(".", ",")
                EW = model_EW[model][line]
                CLF_EW[f"{CLF}"].append(EW)
                if EW >= EW_min and EW <= EW_max:
                    CLF_table[f"{CLF}"].append(Mdot)

        interp_mdot = np.arange(-6.7, -5.4, 0.00001)
        for CLF in CLF_table:
            EWs = CLF_EW[CLF]
            interp_EW = np.interp(x = interp_mdot, xp = Mdots, fp = EWs)
            print(f"{CLF}: {CLF_table[CLF]}")
            interp_inside = interp_mdot[np.where((interp_EW >= EW_min) & (interp_EW <= EW_max))]
            print(f"Interpolated: {interp_inside[0]:.3f} <-> {interp_inside[-1]:.3f} = {np.abs(interp_inside[0] - interp_inside[-1]):.3f}")

    plt.show()
