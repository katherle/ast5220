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

# import txt file
tmp = np.asarray(pd.read_csv('cosmology.txt', sep = " ", header = None))
columns = ["x", "eta", "t", "H", "Hp", "dHpdx", "ddHpddx", "OmegaB", "OmegaCDM", "OmegaLambda", "OmegaR", "OmegaNu", "d_lumin"]
background = QTable(tmp, names = columns)
unit_names = [" ", "m", "s", "s-1", "s-1", "s-1", "s-1", " ", " ", " ", " ", " ", "m"]
for key, unit in zip(background.keys(), unit_names):
    background[key].unit = unit
    background[key] = background[key].si

#put things in sensible units
background["eta"] = background["eta"].to("Mpc")
background["t"] = background["t"].to("Gyr")
background["H"] = background["H"].to("km/(s Mpc)")
background["Hp"] = background["Hp"].to("km/(s Mpc)")
background["dHpdx"] = background["dHpdx"].to("km/(s Mpc)")
background["ddHpddx"] = background["ddHpddx"].to("km/(s Mpc)")
background["d_lumin"] = background["d_lumin"].to("Gpc")

x = background["x"] #for convenience

#find today's values:
H0 = 67*u.km/u.s/u.Mpc
TCMB = 2.7255*u.K
Neff = 0

omegaM0 = 0.05+0.267
omegarad0 = 2*np.pi**2/30*(const.k_B*TCMB)**4/(const.hbar**3 * const.c**5)*8*np.pi*const.G/(3*H0**2)
omegaR0 = omegarad0*(1 + Neff*7/8*(4/11)**(4/3))
omegak0 = 0
omegaLambda0 = 1 - omegaM0 - omegak0 - omegaR0

#equalities
x_rm = np.log(omegaR0/omegaM0)
x_mde = np.log(omegaM0/omegaLambda0)/3

#start of acceleration: dda/dt2 = 0
func = lambda x: omegaLambda0*np.exp(x) - 1/2*omegaM0*np.exp(-2*x) - omegaR0*np.exp(-3*x)
x_acc = fsolve(func, 0)

print(x_rm, x_mde, x_acc)

#plot H and associates
fig, ax = plt.subplots(2, 2, figsize = (10, 10/1.25))
ax[0, 0].plot(x, background["H"]/100, color = "grey")
ax[0, 0].set_yscale("log")
ax[0, 0].set_ylim([10**(-1), 10**8])
ax[0, 0].set_title(r"$\frac{H(x)}{100}$")

ax[0, 1].plot(x, background["Hp"]/100, color = "grey")
ax[0, 1].set_yscale("log")
ax[0, 1].set_ylim([10**(-1), 10**3])
ax[0, 1].set_title(r"$\frac{\mathcal{H}(x)}{100}$")

ax[1, 0].plot(x, background["dHpdx"]/background["Hp"], color = "grey")
ax[1, 0].set_title(r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$")

#ax[1, 1].plot(x, (background["dHpdx"]/background["Hp"])**2, "k--", lw = 1.2, label = r"$(\mathcal{H}'/\mathcal{H})^2$")
ax[1, 1].plot(x, background["ddHpddx"]/background["Hp"], color = "grey")
ax[1, 1].set_title(r"$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$")

for row in ax:
    for axis in row:
        axis.grid()
        axis.set_xlim([-12, 5])
        axis.set_xlabel("x")
        axis.axvline(x_rm, ls = '--', lw = 1, color = "tab:blue", label = r"x$_{rm}$")
        axis.axvline(x_mde, ls = '--', lw = 1, color = "tab:orange", label = r"x$_{m\Lambda}$")
        axis.axvline(x_acc, ls = '--', lw = 1, color = "tab:green", label = r"x$_{\rm acc}$")
        axis.legend()
plt.tight_layout()
plt.savefig("H_etc.pdf")
#plt.show()

#plot eta and t
fig, ax = plt.subplots(1, 2, figsize = (10, 5))
ax[0].plot(x, (background["eta"]/const.c).to("Gyr"), color = "grey")
ax[0].set_ylabel(r"$\eta$(x)/c [Gyr]")
ax[0].set_title("Conformal time")
ax[0].set_ylim(bottom=1e-2)

ax[1].plot(x, background["t"], color = "grey")
ax[1].set_ylabel(r"t(x) [Gyr]")
ax[1].set_title("Proper time")
ax[1].set_ylim(bottom=1e-8)

for axis in ax:
    axis.grid()
    axis.set_xlim([-12, 5])
    #axis.set_ylim([-1, 50])
    axis.set_yscale("log")
    axis.set_xlabel("x")
    axis.axvline(x_rm, ls = '--', lw = 1, color = "tab:blue", label = r"x$_{rm}$")
    axis.axvline(x_mde, ls = '--', lw = 1, color = "tab:orange", label = r"x$_{m\Lambda}$")
    axis.axvline(x_acc, ls = '--', lw = 1, color = "tab:green", label = r"x$_{\rm acc}$")
    axis.legend()
plt.tight_layout()
plt.savefig("times.pdf")
#plt.show()

#plot eta*Hp/c
plt.figure()
plt.plot(x, (background["eta"]*background["Hp"]/const.c).si, color = "grey")
plt.title(r"$\frac{\eta(x)\mathcal{H}(x)}{c}$")
plt.xlabel("x")
plt.yscale("log")
plt.xlim([-12, 5])
#plt.ylim([0.75, 3])
plt.grid()
plt.axvline(x_rm, ls = '--', lw = 1, color = "tab:blue", label = r"x$_{rm}$")
plt.axvline(x_mde, ls = '--', lw = 1, color = "tab:orange", label = r"x$_{m\Lambda}$")
plt.axvline(x_acc, ls = '--', lw = 1, color = "tab:green", label = r"x$_{\rm acc}$")
plt.legend()
plt.savefig("etaHp_c.pdf")
#plt.show()

#plot omegas
omegaR = background["OmegaR"] + background["OmegaNu"]
omegaM = background["OmegaB"] + background["OmegaCDM"]
plt.figure()
plt.plot(x, omegaR, label = r"$\Omega_R$")
plt.plot(x, omegaM, label = r"$\Omega_M$")
plt.plot(x, background["OmegaLambda"], label = r"$\Omega_\Lambda$")
plt.axvline(x_rm, ls = '--', lw = 1, color = "tab:blue", label = r"x$_{rm}$")
plt.axvline(x_mde, ls = '--', lw = 1, color = "tab:orange", label = r"x$_{m\Lambda}$")
plt.axvline(x_acc, ls = '--', lw = 1, color = "tab:green", label = r"x$_{\rm acc}$")
plt.legend()
plt.xlim([-20, 5])
plt.ylabel(r"$\Omega_i$")
plt.xlabel("x")
plt.title("Relative densities")
plt.grid()
plt.savefig("omegas.pdf")
#plt.show()

#get supernova observations
colnames = ["z", "d_L (Gpc)", "Error (Gpc)"]
supernova = pd.read_csv("https://cmb.wintherscoming.no/data/supernovadata.txt", sep='\s+', skiprows = [0], names = colnames)

#plot luminosity distance v redshift
z = 1/np.exp(x) - 1
plt.figure()
plt.plot(z, background["d_lumin"], color = "grey", label = "Luminosity distance")
plt.errorbar(supernova["z"], supernova["d_L (Gpc)"], yerr = supernova["Error (Gpc)"], fmt = ".", color = "k", label = "Supernova data")
plt.yscale("log")
plt.xscale("log")
plt.xlim(right = 2)
plt.ylim(top = 10**2)
#plt.axvline(1/np.exp(x_rm)-1, ls = '--', lw = 1, color = "tab:blue", label = r"z$_{rm}$")
plt.axvline(1/np.exp(x_mde)-1, ls = '--', lw = 1, color = "tab:orange", label = r"z$_{m\Lambda}$")
plt.axvline(1/np.exp(x_acc)-1, ls = '--', lw = 1, color = "tab:green", label = r"z$_{\rm acc}$")
plt.legend()
plt.title("Luminosity distance")
plt.xlabel("z")
plt.ylabel(r"d$_L$ [Gpc]")
plt.grid()
plt.savefig("lumin.pdf")
plt.show()
