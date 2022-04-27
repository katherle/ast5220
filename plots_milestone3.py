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

tmp = np.asarray(pd.read_csv('perturbations_k0.01.txt', sep = " ", header = None))
columns = ["x", "delta_cdm", "v_cdm", "delta_b", "v_b", "Theta0", "Theta1", "Theta2"]
#"Phi", "Psi", "Pi", "ST1", "ST5", "ST50", "ST500"]
perturb_1e2 = QTable(tmp, names = columns)
#unit_names = [" ", " ", " ", " ", " ", " ", " ", " ", " "]
#for key, unit in zip(recombination.keys(), unit_names):
#    recombination[key].unit = unit
#    recombination[key] = recombination[key].si

#for convenience
x_1e2 = perturb_1e2["x"]
delta_cdm_1e2 = perturb_1e2["delta_cdm"]
v_cdm_1e2 = perturb_1e2["v_cdm"]
delta_b_1e2 = perturb_1e2["delta_b"]
v_b_1e2 = perturb_1e2["v_b"]
theta0_1e2 = perturb_1e2["Theta0"]
theta1_1e2 = perturb_1e2["Theta1"]
theta2_1e2 = perturb_1e2["Theta2"]

#plot delta_cdm and v_cdm
fig, ax = plt.subplots()
custom_cycler = cycler("color", ['xkcd:turquoise', 'xkcd:lavender', 'xkcd:dark blue'])
ax.set_prop_cycle(custom_cycler)
ax.plot(x_1e2, delta_cdm_1e2, label = "k = 0.01")
ax.plot(x_1e2, delta_b_1e2, ls = "--", color = 'xkcd:turquoise')
ax.set_xlim([-18, 0])
#ax.set_yscale("log")
ax.set_title(r"$\delta_{\rm CDM}$, $\delta_b$")
ax.grid()
ax.legend()
#plt.show()

#plot theta0
fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(x_1e2, theta0_1e2, label = "k = 0.01")
ax.set_xlim([-18, 0])
ax.set_title(r"$\Theta_0$")
ax.grid()
ax.legend()
#plt.show()

for i in np.array([1]):
    print(i)
