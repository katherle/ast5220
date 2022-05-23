import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy import constants as const
from astropy.table import QTable  # to use tables with units
from astropy.visualization import quantity_support # plots with units
from matplotlib_inline.backend_inline import set_matplotlib_formats
from scipy.optimize import fsolve

#making plots look nicer:
quantity_support()
set_matplotlib_formats('svg')
from matplotlib import cm
from cycler import cycler
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(7, 7/1.25)) # Larger figure sizes
plt.rc('font', size=14)

tmp = np.asarray(pd.read_csv('perturbations_k0.1.txt', sep = " ", header = None))
columns = ["x", "delta_cdm", "v_cdm", "delta_b", "v_b", "Theta0", "Theta1", "Phi", "Psi",
            "Pi", "ST_0", "ST_5", "ST_50", "ST_500"]
perturb_1 = QTable(tmp, names = columns)

#plot the source function to test
fig, ax = plt.subplots()
custom_cycler = cycler("color", ['xkcd:turquoise', 'xkcd:light purple', 'xkcd:dark blue'])
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], perturb_1["ST_50"], label = "k = 0.1/Mpc")
ax.set_xlim([-7.5, 0])
ax.set_xlabel("x")
ax.set_title("Source function")
ax.grid()
ax.legend()
plt.tight_layout()
plt.show()
