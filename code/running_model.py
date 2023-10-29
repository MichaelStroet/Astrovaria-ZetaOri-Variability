
import os, sys
import numpy as np
import matplotlib.pyplot as plt

from broaden import broaden
from directories import getPaths
from plot_1param_run import plotRun

def writeIndat(dirs,catalogue,Teff=30000,logg=3.3,Rs=22,Rmax=120,Tmin=0.6,Mdot=-6.47,Vmin=0.1,Vinf=1850,Beta=1.1,Y=0.08,Vturb=10,Z=1,CLF=50,VCLSTART=0.1,VCLMAX=0.2,CLUMPING="THIN"):
    """
    Writes a new INDAT.DAT file for the given parameters.
    """
    NO_HOPF = "T"

    # Write INDAT.DAT
    indat_file = os.path.join(dirs["fw"], "INDAT.DAT")
    with open(indat_file, "w") as file:
        string = f"'{catalogue}'\n"
        string += f"T T 0 100\n"
        string += f"0.0\n"
        string += f"{Teff} {logg} {Rs}\n"
        string += f"{Rmax} {Tmin}\n"
        string += f"{10**Mdot} {Vmin} {Vinf} {Beta} 0.1\n"
        string += f"{Y} 2\n"
        string += f"F {NO_HOPF} F T T\n"
        string += f"{Vturb} {Z} T T\n"
        string += f"T F 1 2\n"
        if CLUMPING == "THIN":
            string += f"{CLF} {VCLSTART} {VCLMAX}\n"
            string += f"T"
        else:
            string += f"THICK\n"
            string += f"{CLF} {VCLSTART} {VCLMAX}\n"
            string += f"0.1 0.05 0.1\n"
            string += f"0.5 0.05 0.1\n"
            string += f"1.0 0.05 0.1\n"
            string += f"XRAYS 0.17765126\n"
            string += f"0.75 25.0 1.5 460.0 -1000\n"
        file.write(string)


def runModel(dirs, catalogue, **parameters):
    """
    Runs the fastwind model and calculates the resulting line profiles.
    """

    # Write a new INDAT.DAT for the given parameters
    writeIndat(dirs, catalogue, **parameters)

    # Determine options for pformalsol
    iescat = 0
    if "Vturb" in parameters:
        Vturb = parameters["Vturb"]
    else:
        Vturb = 10

    # Change to the fastwind directory and run the model
    print("\n------------------------------\nRunning pnlte_A10HHe\n------------------------------\n")
    os.chdir(dirs["fw"])
    os.system("./pnlte_A10HHe.eo")

    # Calculate the lines of the model and return to the code directory
    print("\n------------------------------\nRunning pformalsol_A10HHe\n------------------------------\n")
    os.system(f"yes '{catalogue}\n{Vturb}\n{iescat}' | ./pformalsol_A10HHe.eo")
    os.chdir(dirs["code"])

def broadenLines(dirs, Vsini, Vmacro):
    """
    Extract the spectra from the OUT. files and applies rotational broadening.
    Returns a dictionary of lines with broadened spectra.
    """

    # Get the line output files
    line_files = [file for file in os.listdir(dirs["catalogue"]) if file.startswith("OUT.")]

    # Get the spectra of the lines
    lines = {}
    for file in line_files:

        # Extract the spectrum from the file
        content = np.genfromtxt(os.path.join(dirs["catalogue"], file), max_rows = 161)
        wave, flux = [content[:,2], content[:,4]]

        # Apply rotational broadening
        spec_res = wave[1] / (wave[1] - wave[0])
        new_wave, new_flux = broaden(wave, flux, Vsini, spec_res, Vmacro)

        # Add the spectrum to the dictionary
        line = file[4:-6]
        lines[line] = np.dstack((new_wave, new_flux))[0]

    return lines

def saveLines(dirs, Vsini=110, Vmacro=0):
    """
    Saves the broadened lines to file
    """
    # Rotationally broaden the lines
    lines = broadenLines(dirs, Vsini, Vmacro)

    # Save the lines to file
    for line in lines:
        path = os.path.join(dirs["output"], f"{line}.dat")
        np.savetxt(path, lines[line])
