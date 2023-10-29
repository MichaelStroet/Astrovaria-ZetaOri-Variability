import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

from directories import getPaths
from plotting_functions import getObservationSpectra, getModelSpectrum, plotObservedLine

dirs = getPaths()
grid_dir = os.path.join(dirs["results"], "grid_M-5,5_-6,5_Cf1,1_100")
models = os.listdir(grid_dir)

for model in models:
    split = model.split("_")
    Mdot = float(split[1][1:].replace(",", "."))
    CLF = float(split[2][2:].replace(",", "."))
    new_name = f"grid_M{Mdot:.1f}_Cf{CLF:.2f}".replace(".", ",")
    print(model)
    print(new_name)
    print()

    if model != new_name:
        os.rename(os.path.join(grid_dir, model), os.path.join(grid_dir, new_name))
