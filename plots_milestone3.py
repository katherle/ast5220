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

tmp = np.asarray(pd.read_csv('perturbations_k0.01.txt', sep = " ", header = None))
perturb_01 = QTable(tmp, names = columns)

tmp = np.asarray(pd.read_csv('perturbations_k0.001.txt', sep = " ", header = None))
perturb_001 = QTable(tmp, names = columns)

H0 = 67*u.km/u.s/u.Mpc
TCMB = 2.7255*u.K

omegaM0 = 0.05+0.267
omegaR0 = 2*np.pi**2/30*(const.k_B*TCMB)**4/(const.hbar**3 * const.c**5)*8*np.pi*const.G/(3*H0**2)
omegaLambda0 = 1 - omegaM0 - omegaR0

xrm = np.log(omegaR0/omegaM0)
func = lambda x: omegaLambda0*np.exp(x) - 1/2*omegaM0*np.exp(-2*x) - omegaR0*np.exp(-3*x)
x_acc = fsolve(func, 0)

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
ax.axvline(xrm, ls = "--", lw = 1, color = "tab:blue", label = r"$x_{\rm rm}$")
ax.axvline(x_acc, ls = "--", lw = 1, color = "tab:green", label = r"$x_{\rm acc}$")
ax.axvline(-6.99, ls = '--', lw = 1, color = "gray", label = "decoupling")
ax.set_xlim([-15, 0])
#ax.set_ylim([10**(-18), 10**5])
ax.set_yscale("log")
ax.set_xlabel("x")
ax.set_ylabel(r"$\delta_{\rm CDM}$, $\delta_b$")
ax.set_title("Matter density perturbations")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("delta.pdf")
plt.show()

fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], np.abs(perturb_1["v_cdm"]), label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], np.abs(perturb_01["v_cdm"]), label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], np.abs(perturb_001["v_cdm"]), label = "k = 0.001/Mpc")
ax.plot(perturb_1["x"], np.abs(perturb_1["v_b"]), ls = "--")
ax.plot(perturb_01["x"], np.abs(perturb_01["v_b"]), ls = "--")
ax.plot(perturb_001["x"], np.abs(perturb_001["v_b"]), ls = "--")
ax.axvline(xrm, ls = "--", lw = 1, color = "tab:blue", label = r"$x_{\rm rm}$")
ax.axvline(x_acc, ls = "--", lw = 1, color = "tab:green", label = r"$x_{\rm acc}$")
ax.axvline(-6.98684, ls = '--', lw = 1, color = "gray", label = "decoupling")
ax.set_xlim([-15, 0])
ax.set_yscale("log")
ax.set_xlabel("x")
ax.set_ylabel(r"$v_{\rm CDM}$, $v_b$")
ax.set_title("Matter velocity perturbations")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("v.pdf")
plt.show()

#plot theta0
fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], perturb_1["Theta0"], label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], perturb_01["Theta0"], label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], perturb_001["Theta0"], label = "k = 0.001/Mpc")
ax.set_xlim([-15, 0])
ax.set_xlabel("x")
ax.set_ylabel(r"$\Theta_0$")
ax.set_title("First photon multipole moment")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("Theta0.pdf")
plt.show()

#plot theta1
fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], perturb_1["Theta1"], label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], perturb_01["Theta1"], label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], perturb_001["Theta1"], label = "k = 0.001/Mpc")
ax.set_xlim([-15, 0])
ax.set_xlabel("x")
ax.set_ylabel(r"$\Theta_1$")
ax.set_title("Second photon multipole moment")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("Theta1.pdf")
plt.show()

#plot phi
fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], perturb_1["Phi"], label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], perturb_01["Phi"], label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], perturb_001["Phi"], label = "k = 0.001/Mpc")
"""
ax.plot(perturb_1["x"], perturb_1["Psi"], ls = "--")
ax.plot(perturb_01["x"], perturb_01["Psi"], ls = "--")
ax.plot(perturb_001["x"], perturb_001["Psi"], ls = "--")
"""
ax.axvline(xrm, ls = "--", lw = 1, color = "tab:blue", label = r"$x_{\rm rm}$")
ax.axvline(x_acc, ls = "--", lw = 1, color = "tab:green", label = r"$x_{\rm acc}$")
ax.set_xlim([-15, 0])
ax.set_xlabel("x")
ax.set_ylabel(r"$\Phi$")
ax.set_title("Gravitational potentials")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("Phi.pdf")
plt.show()

#plot psi
"""
fig, ax = plt.subplots()
ax.set_prop_cycle(custom_cycler)
ax.plot(perturb_1["x"], perturb_1["Psi"], label = "k = 0.1/Mpc")
ax.plot(perturb_01["x"], perturb_01["Psi"], label = "k = 0.01/Mpc")
ax.plot(perturb_001["x"], perturb_001["Psi"], label = "k = 0.001/Mpc")
ax.axvline(xrm, ls = "--", lw = 1, color = "tab:blue", label = r"$x_{\rm rm}$")
ax.axvline(x_acc, ls = "--", lw = 1, color = "tab:green", label = r"$x_{\rm acc}$")
ax.set_xlim([-18, 0])
ax.set_title(r"$\Psi$")
ax.grid()
ax.legend()
plt.savefig("Psi.pdf")
plt.show()
"""
