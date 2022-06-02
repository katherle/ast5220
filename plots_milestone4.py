import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy import constants as const
from astropy.table import QTable  # to use tables with units
from astropy.visualization import quantity_support # plots with units
from matplotlib_inline.backend_inline import set_matplotlib_formats
import scipy as sp

#making plots look nicer:
quantity_support()
set_matplotlib_formats('svg')
from matplotlib import cm
from cycler import cycler
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(7, 7/1.25)) # Larger figure sizes
plt.rc('font', size=14)


#plot the source function to test
tmp = np.asarray(pd.read_csv('perturbations_k0.1.txt', sep = " ", header = None))
columns = ["x", "delta_cdm", "v_cdm", "delta_b", "v_b", "Theta0", "Theta1", "Phi", "Psi",
            "Pi", "ST_0", "ST_5", "ST_100", "ST_500"]
ST = QTable(tmp, names = columns)

fig, ax = plt.subplots(2, 1, figsize = (5, 8))
#custom_cycler = cycler("color", ['xkcd:turquoise', 'xkcd:light purple', 'xkcd:dark blue'])
#ax.set_prop_cycle(custom_cycler)
for axi in ax:
    axi.plot(ST["x"], ST["ST_100"]/1e-3)
    axi.set_xlabel("x")
    axi.set_title("Source function")
    axi.grid()

ax[0].set_xlim([-7.5, 0])
ax[1].set_xlim([-1.5, 0])
plt.tight_layout()
#lt.show()

#plot the bessel function
tmp = np.asarray(pd.read_csv("bessel_test.txt", sep = " ", header = None))
columns = ["x", "j_50"]
Bessel = QTable(tmp, names = columns)

fig, ax = plt.subplots()
ax.plot(Bessel["x"], Bessel["j_50"])
#ax.set_ylim([-0.15, 0.2])
plt.tight_layout()
plt.grid()
#plt.show()

#plot the CMB power spectrum
tmp     = np.asarray(pd.read_csv('cells.txt', sep = " ", header = None))
columns = ["ell", "C_ell"]
Cell    = QTable(tmp, names = columns)
unit_names = [" ", "uK2"]
for key, unit in zip(Cell.keys(), unit_names):
    Cell[key].unit = unit

#for convenience
l      = Cell["ell"]
C_ell  = Cell["C_ell"] #already normalized

fig, ax = plt.subplots()
#ax.plot(l, abs(C_ell)*l/100, label = r"$C_\ell\left(\frac{\ell}{100}\right)$")
ax.plot(l, C_ell, label = r"$C_\ell$")
#ax.set_ylim([-1000, 8000])
ax.set_ylim([1e-3, 1e7])
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel(r"$\ell$")
ax.set_ylabel(r"$\frac{\ell(\ell+1)C_\ell}{2\pi}$ ($\mu$K$^2$)")
ax.set_title("CMB power spectrum")
ax.legend()
ax.grid()
plt.tight_layout()
#plt.show()

# generate cmb map
import healpy as hp
lmax = 2000
Cl = C_ell*2*np.pi/(l*(l+1))/1e12 #undo normalization
Cl = np.concatenate((np.array([0, 0]), Cl), axis = None) #hp expects l to start at 0

np.random.seed(583)
alm = hp.synalm(Cl, lmax = lmax, new = True)
high_nside = 1024
cmb_map = hp.alm2map(alm, nside = high_nside, lmax = lmax)
hp.mollview(cmb_map, min=-300*1e-6, max=300*1e-6, unit = "K", title = "CMB Temperature", bgcolor = "k")
#plt.show()

#matter power spectrum
tmp       = np.asarray(pd.read_csv('matter.txt', sep = " ", header = None))
columns   = ["k", "P(k)", "theta6", "theta30", "theta100", "theta200", "theta500", "theta1000"]
Matter    = QTable(tmp, names = columns)
unit_names = ["Mpc-1", "Mpc3", " ", " ", " ", " ", " ", " "]
for key, unit in zip(Matter.keys(), unit_names):
    Matter[key].unit = unit

#for convenience
k     = Matter["k"]*3.08567758e22
P_k   = Matter["P(k)"]/(3.08567758e22**3)
theta_6  = Matter["theta6"]
theta_30 = Matter["theta30"]
theta_100 = Matter["theta100"]
theta_200 = Matter["theta200"]
theta_500 = Matter["theta500"]
theta_1000 = Matter["theta1000"]

h    = 0.67
H0   = h*10**5*u.m/u.s/u.Mpc
k_eq = 0.05

fig, ax = plt.subplots()
ax.plot(k/h, np.abs(P_k)*(h**3))
ax.axvline(k_eq, ls = "--", color = "k")
ax.set_yscale("log")
ax.set_xscale("log")
ax.grid()
plt.tight_layout()
#plt.show()

eta0 = (46.5729*u.Gyr*const.c).to("Mpc")

fig, ax = plt.subplots()
ax.plot(k*eta0, theta_6, lw = 0.5, label = "ell=6")
ax.plot(k*eta0, theta_30, lw = 0.5, label = "ell=30")
ax.plot(k*eta0, theta_100, lw = 0.5, label = "ell=100")
#ax.plot(k*eta0, theta_200, lw = 0.5, label = "ell=200")
#ax.plot(k*eta0, theta_500, lw = 0.5, label = "ell=500")
#ax.plot(k*eta0, theta_1000, lw = 0.5, label = "ell=1000")
#ax.set_ylim([-0.015, 0.015])
ax.set_xlim([0, 1000])
ax.legend()
plt.grid()
plt.tight_layout()
#plt.show()

fig, ax = plt.subplots()
ax.plot(k*eta0, theta_6**2, lw = 0.5, label = "ell=6")
ax.plot(k*eta0, theta_30**2, lw = 0.5, label = "ell=30")
ax.plot(k*eta0, theta_100**2, lw = 0.5, label = "ell=100")
#ax.plot(k*eta0, theta_200, lw = 0.5, label = "ell=200")
#ax.plot(k*eta0, theta_500, lw = 0.5, label = "ell=500")
#ax.plot(k*eta0, theta_1000, lw = 0.5, label = "ell=1000")
#ax.set_ylim([-0.015, 0.015])
ax.set_xlim([0, 1000])
ax.legend()
plt.grid()
plt.tight_layout()
#plt.show()

fig, ax = plt.subplots()
#ax.plot(const.c*k/H0, theta_6, lw = 0.5)
ax.plot(k*eta0, sp.special.spherical_jn(6, (k*eta0).value), lw = 0.5, label = r"$j_{\ell}(k\eta_0)$")
ax.plot(k*eta0, theta_6, lw = 0.5, label = "ell=6")
#ax.set_ylim([-0.015, 0.015])
ax.set_xlim([0, 1000])
ax.legend()
plt.grid()
plt.tight_layout()
plt.show()
