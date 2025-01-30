"""Compares probe-based in-plane statistics with volume-averaged ones"""

import matplotlib.pyplot as plt
from os.path import join
import numpy as np
import h5py

plt.style.use("tableau-colorblind10")
plt.rc("text", usetex=True)
plt.rcParams.update({"font.size": 8})
plt.rcParams.update(
    {"axes.grid": True, "grid.color": "#DDDDDD", "grid.linestyle": "dashed"}
)


# Path with the simulation cases
DATA = "."

# Case directories
case_names = ["M1", "M2", "M3"]

# Path for saving figures
SAVE_PATH = "figures"

# Viscosity
nu = 5e-5

# %%


def read_data():
    times = ["750", "750", "750"]
    cases = {}

    for i, name in enumerate(case_names):
        cases[name] = {}
        time = times[i]

        cases[name]["y"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "UMean_X.xy",
            )
        )[:, 0]
        cases[name]["u"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "UMean_X.xy",
            )
        )[:, 1]

        cases[name]["uu"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "UPrime2Mean_XX.xy",
            )
        )[:, 1]
        cases[name]["vv"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "UPrime2Mean_YY.xy",
            )
        )[:, 1]
        cases[name]["ww"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "UPrime2Mean_ZZ.xy",
            )
        )[:, 1]
        cases[name]["uv"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "UPrime2Mean_XY.xy",
            )
        )[:, 1]
        cases[name]["tau"] = (
            0.5
            * nu
            / cases[name]["y"][1]
            * (cases[name]["u"][1] + cases[name]["u"][-2])
        )

        cases[name]["utau"] = np.sqrt(cases[name]["tau"])
        cases[name]["delta_nu"] = nu / cases[name]["utau"]
        cases[name]["y+"] = cases[name]["y"] / cases[name]["delta_nu"]

    return cases


cases = read_data()

# %% The plane data

data_m2 = np.genfromtxt(join(DATA, "M2", "plane_stats.txt"))

umean = data_m2[:, 0]
uvar = data_m2[:, 1]
vvar = data_m2[:, 2]
wvar = data_m2[:, 3]


ny = umean.size
# %%

plt.figure(figsize=(5, 2.5))
plt.subplot(121)
for i, name in enumerate(["M2"]):
    plt.semilogx(
        cases[name]["y+"], cases[name]["uu"] / 0.05**2, lw=1, color="C" + str(i)
    )

plt.plot(cases["M2"]["y+"][1 : ny + 1], uvar / 0.05**2, "C1", alpha=0.5)


plt.xlim(1, 1000)
plt.ylim(-1, 9)
plt.ylabel(r"$\overline {u'u'}^+$")
plt.xlabel(r"$y^+$")


plt.subplot(122)
for i, name in enumerate(["M2"]):
    plt.semilogx(cases[name]["y+"], cases[name]["vv"] / 0.05**2, ":", lw=1, color="C0")
    plt.plot(cases[name]["y+"], cases[name]["ww"] / 0.05**2, "-.", lw=1, color="C0")

plt.plot(cases["M2"]["y+"][1 : ny + 1], vvar / 0.05**2, ":C1")
plt.plot(cases["M2"]["y+"][1 : ny + 1], wvar / 0.05**2, "-.C1")
plt.xlim(1, 1000)
plt.ylim(0, 2)
plt.ylabel(r"$\overline {v'v'}^+$, $\overline {w'w'}^+$")
plt.xlabel(r"$y^+$")
plt.tight_layout()
plt.savefig(
    join(SAVE_PATH, "plane_vs_volume.pdf"), bbox_inches="tight", pad_inches=0.02
)
