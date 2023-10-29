"""
Creates models for a grid of two parameters (Massloss and clumping factor)
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt

from broaden import broaden
from directories import getPaths
from plot_epochs import plotEpochs, plotAllLinesSingleModel
from running_model import runModel, saveLines

def runModelAndSaveLines(dirs, catalogue, **parameters):
    """
    Runs a model with the given parameters,
    applies rotational broadening to the lines
    and saves them to file.
    """

    # Check if the catalogue directory exists, if not, create it.
    dirs["catalogue"] = os.path.join(dirs["fw"], catalogue)
    if not os.path.exists(dirs["catalogue"]):
        os.mkdir(dirs["catalogue"])

    # Check if the output lines directory exists, if not, create it.
    dirs["output"] = os.path.join(dirs["results"], catalogue)
    if not os.path.exists(dirs["output"]):
        os.mkdir(dirs["output"])

    # Run a new model
    print(f"\n------------------------------\nRunning {catalogue}\n------------------------------\n")
    runModel(dirs, catalogue, **parameters)

    # Save the lines to file
    print("\n------------------------------\nSaving lines to results\n------------------------------\n")
    saveLines(dirs)

    # Create plots of the run
    print("\n------------------------------\nPlotting results\n------------------------------\n")
    plotEpochs(dirs, catalogue, lines = ["HALPHA"])
    plotAllLinesSingleModel(dirs, catalogue)


if __name__ ==  "__main__":

    # Get the paths to the directories
    dirs = getPaths()

    par_labels = {"Mdot":     "M",
                  "Beta":     "B",
                  "logg":     "g",
                  "Teff":     "T",
                  "Rs":       "R",
                  "CLF":      "Cf",
                  "VCLSTART": "Cs",
                  "VCLMAX":   "Cm",
                  "Vinf":     "Vi",
                  "Y":        "Y"}

    # Grid parameters
    masslosses = np.arange(5.9, 6.6, 0.1)
    factors = np.logspace(np.log10(1.1), np.log10(100), 11)

    # Calculate a model for each grid point
    for Mdot in masslosses:
        for CLF in factors:
            input_parameters = {}
            input_parameters["Mdot"] = -Mdot
            input_parameters["CLF"] = CLF

            catalogue = f"grid_M-{Mdot:.1f}_Cf{CLF:.2f}"

            # Run the model and save the output
            runModelAndSaveLines(dirs, catalogue, **input_parameters)
