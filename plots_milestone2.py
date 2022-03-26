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

tmp = np.asarray(pd.read_csv('recombination.txt', sep = " ", header = None))
columns = ["x", "Xe", "ne", "tau", "dtaudx", "ddtauddx", "g_tilde", "dgdx_tilde", "ddgddx_tilde"]
recombination = QTable(tmp, names = columns)
unit_names = [" ", " ", "m-3", " ", " ", " ", " ", " ", " "]
for key, unit in zip(recombination.keys(), unit_names):
    recombination[key].unit = unit
    recombination[key] = recombination[key].si

#for convenience
x = recombination["x"]
Xe = recombination["Xe"]
ne = recombination["ne"]
tau = recombination["tau"]
dtau = recombination["dtaudx"]
ddtau = recombination["ddtauddx"]
g = recombination["g_tilde"]
dg = recombination["dgdx_tilde"]
ddg = recombination["ddgddx_tilde"]

#plot tau and its first two derivatives
fig, ax = plt.subplots()
custom_cycler = cycler("color", cm.viridis(np.linspace(0, 1, 3)))
ax.set_prop_cycle(custom_cycler)
ax.plot(x, tau, label = r"$\tau(x)$")
ax.plot(x, -dtau, label = r"$-\tau'(x)$")
ax.plot(x, ddtau, label = r"$\tau''(x)$")
ax.axvline(-6.98684, ls = '--', lw = 1, color = "limegreen", label = "decoupling")
ax.axvline(-7.16378, ls = '--', lw = 1, color = "darkblue", label = "recombination")
ax.set_yscale("log")
ax.set_xlim([-12, 0])
ax.set_ylim([1e-8, 1e8])
ax.set_xlabel("x")
ax.set_title("Optical depth")
plt.legend()
plt.grid()
plt.savefig("tau.pdf")
plt.show()

#plot g_tilde and its first two derivatives
fig, ax = plt.subplots(2, 1, figsize = (7, 10))
ax[0].set_prop_cycle(custom_cycler)
ax[0].set_xlim([-12, 0])
ax[0].plot(x, g, label = r"$\tilde{g}(x)$")
ax[0].plot(x, dg, label = r"$\tilde{g}'(x)$")
ax[0].plot(x, ddg, label = r"$\tilde{g}''(x)$")
ax[0].set_title("Unscaled")

ax[1].set_prop_cycle(custom_cycler)
ax[1].set_xlim([-8, -6])
ax[1].plot(x, g/np.std(g), label = r"$\tilde{g}(x)$")
ax[1].plot(x, dg/np.std(dg), label = r"$\tilde{g}'(x)$")
ax[1].plot(x, ddg/np.std(ddg), label = r"$\tilde{g}''(x)$")
ax[1].axvline(-6.98684, ls = '--', lw = 1, color = "limegreen", label = "decoupling")
ax[1].axvline(-7.16378, ls = '--', lw = 1, color = "darkblue", label = "recombination")
ax[1].set_title("Scaled")

for axi in ax:
    axi.set_xlabel("x")
    axi.set_ylabel("Visibility function")
    axi.legend(loc=3, fontsize=12)
    axi.grid()
plt.tight_layout()
plt.savefig("g_tilde.pdf")
plt.show()

#plot Xe
fig, ax = plt.subplots()
ax.plot(x, Xe, color = "grey")
ax.axvline(-6.98684, ls = '--', lw = 1, color = "limegreen", label = "decoupling")
ax.axvline(-7.16378, ls = '--', lw = 1, color = "darkblue", label = "recombination")
ax.axhline(0.00019854, ls = '--', lw = 1, color = "k", label = "freezeout")
ax.set_xlabel("x")
ax.set_title("Xe(x)")
ax.set_xlim([-8, -6])
ax.set_yscale("log")
plt.legend()
plt.grid()
plt.savefig("Xe.pdf")
plt.show()
