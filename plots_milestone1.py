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
plt.rc('figure', figsize=(10, 10/1.25)) # Larger figure sizes
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
#print(background["H"])
#print(background["dHpdx"]/background["Hp"])

#plot H and associates
fig, ax = plt.subplots(2, 2)
ax[0, 0].plot(background["x"], background["H"]*100)
ax[0, 0].set_yscale("log")
ax[0, 0].set_title(r"H(x) [100 km s$^{-1}$ Mpc$^{-1}$]")

ax[0, 1].plot(background["x"], background["Hp"]*100)
ax[0, 1].set_yscale("log")
ax[0, 1].set_title(r"$\mathcal{H}$(x) [100 km s$^{-1}$ Mpc$^{-1}$]")

ax[1, 0].plot(background["x"], background["dHpdx"]/background["Hp"])
ax[1, 0].set_title(r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$ []")

ax[1, 1].plot(background["x"], background["ddHpddx"]/background["Hp"])
ax[1, 1].set_title(r"$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$ []")

plt.savefig("H_etc.pdf")
plt.show()
