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
    axi.plot(ST["x"], ST["ST_100"]/1e-3, lw = 1, color = "grey")
    axi.set_xlabel("x")
    axi.set_title("Source function")
    axi.grid()

ax[0].set_xlim([-7.5, 0])
ax[1].set_xlim([-1.5, 0])
plt.tight_layout()
plt.savefig("ST.pdf")
plt.show()

"""
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
"""

#plot the CMB power spectrum
tmp          = np.asarray(pd.read_csv('cells.txt', sep = " ", header = None))
tmp_sw       = np.asarray(pd.read_csv('cells_sw.txt', sep = " ", header = None))
tmp_isw      = np.asarray(pd.read_csv('cells_isw.txt', sep = " ", header = None))
tmp_doppler  = np.asarray(pd.read_csv('cells_doppler.txt', sep = " ", header = None))
tmp_quad     = np.asarray(pd.read_csv('cells_quad.txt', sep = " ", header = None))
columns      = ["ell", "C_ell"]
Cell         = QTable(tmp, names = columns)
Cell_sw      = QTable(tmp_sw, names = columns)
Cell_isw     = QTable(tmp_isw, names = columns)
Cell_doppler = QTable(tmp_doppler, names = columns)
Cell_quad    = QTable(tmp_quad, names = columns)
unit_names = [" ", "uK2"]
for key, unit in zip(Cell.keys(), unit_names):
    Cell[key].unit         = unit
    Cell_sw[key].unit      = unit
    Cell_isw[key].unit     = unit
    Cell_doppler[key].unit = unit
    Cell_quad[key].unit    = unit

#for convenience
l             = Cell["ell"]
C_ell         = Cell["C_ell"] #already normalized
C_ell_sw      = Cell_sw["C_ell"]
C_ell_isw     = Cell_isw["C_ell"]
C_ell_doppler = Cell_doppler["C_ell"]
C_ell_quad    = Cell_quad["C_ell"]

#Planck data
colnames_low = ["ell", "D_ell", "err_up", "err_down"]
planck_low = pd.read_csv("https://cmb.wintherscoming.no/data/planck_cell_low.txt", sep='\s+', skiprows = [0], names = colnames_low)
colnames_high = ["ell", "D_ell", "err_down", "err_up", "best_fit"]
planck_high = pd.read_csv("https://cmb.wintherscoming.no/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt", sep = '\s+', skiprows = [0], names = colnames_high)

fig, ax = plt.subplots()
ax.plot(l, C_ell, color = "grey", label = r"$C_\ell$")
#ax.plot(l**1.018, C_ell*np.exp(-0.05*(l/200)**1.5), color = "b", label = r"$C_\ell$ (corrected)")
ax.errorbar(planck_low["ell"], planck_low["D_ell"], yerr = [planck_low["err_down"], planck_low["err_up"]], fmt = ".", color = "k", label = "Planck data")
ax.errorbar(planck_high["ell"], planck_high["D_ell"], yerr = [planck_high["err_down"], planck_high["err_up"]], fmt = ".", color = "k")
ax.set_ylim([-1000, 7000])
ax.set_xscale("log")
ax.set_xlabel(r"$\ell$")
ax.set_ylabel(r"$\frac{\ell(\ell+1)C_\ell}{2\pi}$ ($\mu$K$^2$)")
ax.set_title("CMB power spectrum")
ax.legend()
ax.grid()
plt.tight_layout()
plt.savefig("CMB_spectrum.pdf")
plt.show()

fig, ax = plt.subplots()
custom_cycler = cycler("color", ['xkcd:light purple', 'xkcd:dark blue', 'xkcd:turquoise', 'xkcd:green'])
ax.set_prop_cycle(custom_cycler)
ax.plot(l, C_ell, color = "grey", label = "Total")
ax.plot(l, C_ell_sw, lw = 1, label = "SW term")
ax.plot(l, C_ell_isw, lw = 1, label = "ISW term")
ax.plot(l, C_ell_doppler, lw = 1, label = "Doppler term")
ax.plot(l, C_ell_quad, lw = 1, label = "Quadrupole term")
ax.set_ylim([10**(-2), 10**7])
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_title("Comparison of source function terms")
ax.set_xlabel(r"$\ell$")
ax.set_ylabel(r"$\frac{\ell(\ell+1)C_\ell}{2\pi}$ ($\mu$K$^2$)")
ax.legend()
plt.grid()
plt.tight_layout()
plt.savefig("CMB_terms.pdf")
plt.show()

# generate cmb map
import healpy as hp
lmax = 2000
Cl = C_ell*2*np.pi/(l*(l+1))/1e12 #undo normalization
Cl = np.concatenate((np.array([0, 0]), Cl), axis = None) #hp expects l to start at 0

np.random.seed(583)
alm = hp.synalm(Cl, lmax = lmax, new = True)
high_nside = 1024
cmb_map = hp.alm2map(alm, nside = high_nside, lmax = lmax)
cmap = cm.RdYlBu
hp.mollview(cmb_map, min=-300*1e-6, max=300*1e-6, cmap = cmap, unit = "K", title = "CMB Temperature")
plt.savefig("CMB_map.png")
plt.show()

#matter power spectrum
tmp       = np.asarray(pd.read_csv('matter.txt', sep = " ", header = None))
columns   = ["k", "P(k)", "theta6", "theta30", "theta100", "theta200", "theta500", "theta1000"]
Matter    = QTable(tmp, names = columns)
unit_names = ["Mpc-1", "Mpc3", " ", " ", " ", " ", " ", " "]
for key, unit in zip(Matter.keys(), unit_names):
    Matter[key].unit = unit

#for convenience
Mpc = 3.08567758e22
k     = Matter["k"]*Mpc
P_k   = Matter["P(k)"]/(Mpc**3)
theta_6  = Matter["theta6"]
theta_30 = Matter["theta30"]
theta_100 = Matter["theta100"]
theta_200 = Matter["theta200"]
theta_500 = Matter["theta500"]
theta_1000 = Matter["theta1000"]

h    = 0.67
H0   = h*10**5*u.m/u.s/u.Mpc

TCMB = 2.7255*u.K
Neff = 0.0

omegaM0 = 0.05+0.267
omegarad0 = 2*np.pi**2/30*(const.k_B*TCMB)**4/(const.hbar**3 * const.c**5)*8*np.pi*const.G/(3*H0**2)
omegaR0 = omegarad0*(1 + Neff*7/8*(4/11)**(4/3))
omegak0 = 0
omegaLambda0 = 1 - omegaM0 - omegak0 - omegaR0

#matter radiation equality
a_eq = omegaR0/omegaM0
H_eq = H0 * np.sqrt(omegaM0/a_eq**3 + omegaR0/a_eq**4 + omegak0/a_eq**2 + omegaLambda0)
k_eq = a_eq*H_eq/const.c


#data
colnames_matter = ["k/h", "P(k)/h**3", "Error"]
data_matter = pd.read_csv("https://cmb.wintherscoming.no/data/reid_DR7.txt", sep = "\s+", skiprows = [0], names = colnames_matter)

fig, ax = plt.subplots()
ax.plot(k, np.abs(P_k)*(h**3), color = "grey", label = "Theory")
ax.errorbar(data_matter["k/h"], data_matter["P(k)/h**3"], yerr = data_matter["Error"], fmt = ".", color = "k", label = "SDSS data")
ax.axvline(k_eq, ls = "--", lw = 1, color = "k", label = r"k$_{\rm eq}$")
ax.set_xlabel(r"k [h Mpc$^{-1}$]")
ax.set_ylabel(r"P(k) [Mpc$^3$/h$^3$]")
ax.set_title("Matter power spectrum")
ax.set_yscale("log")
ax.set_xscale("log")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("matter_spectrum.pdf")
plt.show()

eta0 = (46.5729*u.Gyr*const.c).to("Mpc")

fig, ax = plt.subplots(2, 1, figsize = (6, 8))
ax[0].set_prop_cycle(custom_cycler)
ax[0].plot(const.c*k/H0, theta_6, lw = 0.5, label = r"$\ell$=6")
#ax.plot(k, theta_30, lw = 0.5, label = "ell=30")
ax[0].plot(const.c*k/H0, theta_100, lw = 0.5, label = r"$\ell$=100")
ax[0].plot(const.c*k/H0, theta_200, lw = 0.5, label = r"$\ell$=200")
ax[0].plot(const.c*k/H0, theta_500, lw = 0.5, label = r"$\ell$=500")
#ax.plot(const.c*k/H0, theta_1000, lw = 0.5, label = r"$\ell$=1000")
ax[0].set_ylim([-0.015, 0.015])
ax[0].set_xlim([0, 450])
ax[0].set_xlabel(r"ck/H$_0$")
ax[0].set_ylabel(r"$\Theta_\ell$")
ax[0].set_title("Transfer function")
ax[0].legend()
ax[0].grid()

ax[1].set_prop_cycle(custom_cycler)
ax[1].plot(const.c*k/H0, theta_6**2*H0/k/const.c, lw = 0.5, label = r"$\ell$=6")
#ax.plot(k, theta_30**2/k, lw = 0.5, label = "ell=30")
ax[1].plot(const.c*k/H0, theta_100**2*H0/k/const.c, lw = 0.5, label = r"$\ell$=100")
ax[1].plot(const.c*k/H0, theta_200**2*H0/k/const.c, lw = 0.5, label = r"$\ell$=200")
ax[1].plot(const.c*k/H0, theta_500**2*H0/k/const.c, lw = 0.5, label = r"$\ell$=500")
#ax.plot(const.c*k/H0, theta_1000**2*H0/k/const.c, lw = 0.5, label = "ell=1000")
ax[1].set_ylim([0, 2.5e-6])
ax[1].set_xlim([0, 200])
ax[1].set_xlabel(r"ck/H$_0$")
ax[1].set_ylabel(r"$\Theta_\ell^2 H_0/ck$")
ax[1].set_title("Spectrum integrand")
ax[1].grid()
ax[1].legend()

plt.tight_layout()
plt.savefig("transfer.pdf")
plt.show()

"""
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
"""
