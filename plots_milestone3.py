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
columns = ["x", "delta_cdm", "v_cdm", "delta_b", "v_b", "Theta0", "Theta1", "Phi"]
perturb_1 = QTable(tmp, names = columns)

tmp = np.asarray(pd.read_csv('perturbations_k0.01.txt', sep = " ", header = None))
perturb_01 = QTable(tmp, names = columns)

tmp = np.asarray(pd.read_csv('perturbations_k0.001.txt', sep = " ", header = None))
perturb_001 = QTable(tmp, names = columns)


#plot delta_cdm and v_cdm
fig, ax = plt.subplots()
custom_cycler = cycler("color", ['xkcd:turquoise', 'xkcd:light purple', 'xkcd:dark blue'])
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], np.abs(perturb_1["delta_cdm"]), label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], np.abs(perturb_01["delta_cdm"]), label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], np.abs(perturb_001["delta_cdm"]), label = "k = 0.001/Mpc")
ax.plot(perturb_1["x"], np.abs(perturb_1["delta_b"]), ls = "--")
ax.plot(perturb_01["x"], np.abs(perturb_01["delta_b"]), ls = "--")
ax.plot(perturb_001["x"], np.abs(perturb_001["delta_b"]), ls = "--")
ax.set_xlim([-18, 0])
#ax.set_ylim([10**(-18), 10**5])
ax.set_yscale("log")
ax.set_title(r"$\delta_{\rm CDM}$, $\delta_b$")
ax.grid()
ax.legend()
plt.show()

fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], np.abs(perturb_1["v_cdm"]), label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], np.abs(perturb_01["v_cdm"]), label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], np.abs(perturb_001["v_cdm"]), label = "k = 0.001/Mpc")
ax.plot(perturb_1["x"], np.abs(perturb_1["v_b"]), ls = "--")
ax.plot(perturb_01["x"], np.abs(perturb_01["v_b"]), ls = "--")
ax.plot(perturb_001["x"], np.abs(perturb_001["v_b"]), ls = "--")
ax.set_xlim([-18, 0])
ax.set_yscale("log")
ax.set_title(r"$v_{\rm CDM}$, $v_b$")
ax.grid()
ax.legend()
plt.show()

#plot theta0
fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], perturb_1["Theta0"], label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], perturb_01["Theta0"], label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], perturb_001["Theta0"], label = "k = 0.001/Mpc")
ax.set_xlim([-18, 0])
ax.set_title(r"$\Theta_0$")
ax.grid()
ax.legend()
plt.show()

#plot theta1
fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], perturb_1["Theta1"], label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], perturb_01["Theta1"], label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], perturb_001["Theta1"], label = "k = 0.001/Mpc")
ax.set_xlim([-18, 0])
ax.set_title(r"$\Theta_1$")
ax.grid()
ax.legend()
plt.show()

#plot phi
fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], perturb_1["Phi"], label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], perturb_01["Phi"], label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], perturb_001["Phi"], label = "k = 0.001.Mpc")
ax.set_xlim([-18, 0])
ax.set_title(r"$\Phi$")
ax.grid()
ax.legend()
plt.show()
