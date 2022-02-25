import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy import constants as const
from astropy.table import QTable  # to use tables with units
from astropy.visualization import quantity_support # plots with units
from matplotlib_inline.backend_inline import set_matplotlib_formats

#making plots look nicer:
quantity_support()
set_matplotlib_formats('svg')
from matplotlib import cm
from cycler import cycler
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(7, 7/1.25)) # Larger figure sizes
plt.rc('font', size=10)

# import txt file
tmp = np.asarray(pd.read_csv('cosmology.txt', sep = " ", header = None))
columns = ["x", "eta", "t", "H", "Hp", "dHpdx", "ddHpddx", "OmegaB", "OmegaCDM", "OmegaLambda", "OmegaR", "OmegaNu", "d_lumin"]
background = QTable(tmp, names = columns)
unit_names = [" ", "m", "s", "s-1", "s-1", "s-1", "s-1", " ", " ", " ", " ", " ", " ", " ", " "]
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

#plot H and associates
fig, ax = plt.subplots(2, 2, figsize = (10, 10/1.25))
ax[0, 0].plot(background["x"], background["H"]/100)
ax[0, 0].set_yscale("log")
ax[0, 0].set_ylim([10**(-1), 10**8])
ax[0, 0].set_ylabel("")
ax[0, 0].set_title(r"H(x) [100 km s$^{-1}$ Mpc$^{-1}$]")

ax[0, 1].plot(background["x"], background["Hp"]/100)
ax[0, 1].set_yscale("log")
ax[0, 1].set_ylim([10**(-1), 10**3])
ax[0, 1].set_ylabel("")
ax[0, 1].set_title(r"$\mathcal{H}$(x) [100 km s$^{-1}$ Mpc$^{-1}$]")

ax[1, 0].plot(background["x"], background["dHpdx"]/background["Hp"])
ax[1, 0].set_title(r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$")

ax[1, 1].plot(background["x"], background["ddHpddx"]/background["Hp"])
ax[1, 1].set_title(r"$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$")

for row in ax:
    for axis in row:
        axis.grid()
        axis.set_xlim([-12, 0])
plt.savefig("H_etc.pdf")
plt.show()

#plot eta and t
fig, ax = plt.subplots(1, 2, figsize = (10, 5))
ax[0].plot(background["x"], background["eta"])
ax[0].set_ylabel(r"$\eta$(x) [Mpc]")
ax[0].set_ylim([1, 10**5])
ax[0].set_title("Conformal time")

ax[1].plot(background["x"], background["t"])
ax[1].set_ylabel(r"t(x) [Gyr]")
ax[1].set_ylim([10**(-8), 10**2])
ax[1].set_title("Proper time")

for axis in ax:
    axis.grid()
    axis.set_xlim([-12, 0])
    axis.set_yscale("log")
plt.savefig("times.pdf")
plt.show()

#plot eta*Hp/c
plt.figure()
plt.plot(background["x"], (background["eta"]*background["Hp"]/const.c).si)
plt.title(r"$\frac{\eta(x)\mathcal{H}(x)}{c}$")
plt.xlim([-15, 0])
plt.ylim([0.75, 3])
plt.grid()
plt.savefig("etaHp_c.pdf")
plt.show()

#plot omegas
omegaM = background["OmegaB"] + background["OmegaCDM"]
omegaR = background["OmegaR"] + background["OmegaNu"]
plt.figure()
plt.plot(background["x"], omegaR, label = r"$\Omega_R$")
plt.plot(background["x"], omegaM, label = r"$\Omega_M$")
plt.plot(background["x"], background["OmegaLambda"], label = r"$\Omega_\lambda$")
plt.title(r"$\Omega_i$")
plt.legend()
plt.xlim([-20, 5])
plt.grid()
plt.savefig("omegas.pdf")
plt.show()
