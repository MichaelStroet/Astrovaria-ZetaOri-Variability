"""
Creates models for a run of a single parameter
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt

from broaden import broaden
from directories import getPaths
from plot_1param_run import plotRun
from running_model import runModel, saveLines

def runModelAndSaveLines(dirs, catalogue, run_name, **parameters):
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
    dirs["output"] = os.path.join(dirs["run"], catalogue)
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
    plotRun(dirs, run_name)

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

    # Define the run parameters
    par = "Teff"
    min = 28000
    max = 35000
    scale = "lin"
    num = 12
    
    # Set the name of the run and create an output directory
    run_name = f"Run_{par_labels[par]}_{min:.2f}_{max:.2f}".replace('.', ',')
    dirs["run"] = os.path.join(dirs["results"], run_name)
    if not os.path.exists(dirs["run"]):
        os.mkdir(dirs["run"])

    # Get the parameter values
    if scale == "lin":
        values = np.linspace(min, max, num)
    elif scale == "log":
        values = np.log10(np.logspace(min, max, num))
    else:
        print(f"Error: Unknown scale {scale}")
        exit(1)

    print(f"{par}:\n{values}\n")

    # Run the model for each parameter value
    for value in values:
        catalogue = f"{par_labels[par]}_{value:.2f}".replace(".", ",")
        print(catalogue)
        # Define the INDAT parameters
        input_parameters = {par: value}

        input_parameters["VCLSTART"] = 0.05
        if par == "VCLSTART":
            input_parameters["VCLMAX"] = 0.8

        # Run the model and save the output
        runModelAndSaveLines(dirs, catalogue, run_name, **input_parameters)
