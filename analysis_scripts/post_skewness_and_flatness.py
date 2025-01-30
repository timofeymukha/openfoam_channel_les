"""Plots skewness and flatness. DNS data available upon request."""

import matplotlib.pyplot as plt
from os.path import join
import numpy as np
import h5py
from scipy.stats import skew

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

FIGSIZE = (5, 2.5)

# Viscosity
nu = 5e-5


# Available upon request
# simson_ho = np.genfromtxt(join(DATA, "dns", "simson_high_order.txt"))


# %%


def read_data():
    times = ["750", "750", "750"]
    cases = {}

    for i, name in enumerate(case_names):
        print(name)
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

        cases[name]["tau"] = (
            0.5
            * nu
            / cases[name]["y"][1]
            * (cases[name]["u"][1] + cases[name]["u"][-2])
        )

        cases[name]["utau"] = np.sqrt(cases[name]["tau"])
        cases[name]["delta_nu"] = nu / cases[name]["utau"]

        cases[name]["y"] = cases[name]["y"][1 : cases[name]["y"].size // 2]

        cases[name]["y+"] = cases[name]["y"] / cases[name]["delta_nu"]

        ho = np.genfromtxt(join(DATA, name, "plane_stats.txt"))

        cases[name]["s_u"] = ho[:, 4]
        cases[name]["s_v"] = ho[:, 5]
        cases[name]["s_w"] = ho[:, 6]

        cases[name]["f_u"] = ho[:, 7]
        cases[name]["f_v"] = ho[:, 8]
        cases[name]["f_w"] = ho[:, 9]
    return cases


cases = read_data()

# %% skewness

fig = plt.figure(figsize=FIGSIZE)

ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

for i, name in enumerate(case_names):
    ax1.semilogx(cases[name]["y+"], cases[name]["s_u"], lw=1, color="C" + str(i + 4))

    ax2.semilogx(cases[name]["y+"], cases[name]["s_v"], lw=1, color="C" + str(i + 4))

    ax3.semilogx(
        cases[name]["y+"],
        cases[name]["s_w"],
        lw=1,
        label=name,
        color="C" + str(i + 4),
    )

# ax1.plot(simson_ho_sd[:, 1], simson_ho[:, 2], 'k', lw=1,)
# ax2.plot(simson_ho_sd[:, 1], simson_ho[:, 3], 'k', lw=1,)
# ax3.plot(simson_ho_sd[:, 1], simson_ho[:, 4], 'k', lw=1,label="DNS")
ax3.legend()

ax1.set_xlim(1, 1000)
ax2.set_xlim(1, 1000)
ax3.set_xlim(1, 1000)
ax1.set_title(r"$S_u$")
ax2.set_title(r"$S_v$")
ax3.set_title(r"$S_w$")
ax1.set_xlabel(r"$y^+$")
ax2.set_xlabel(r"$y^+$")
ax3.set_xlabel(r"$y^+$")
plt.tight_layout()

plt.savefig(join(SAVE_PATH, "skewness.pdf"), bbox_inches="tight", pad_inches=0.02)

# %% flatness

fig = plt.figure(figsize=FIGSIZE)

ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

for i, name in enumerate(case_names):
    ax1.semilogx(cases[name]["y+"], cases[name]["f_u"], lw=1, color="C" + str(i + 4))

    ax2.semilogx(cases[name]["y+"], cases[name]["f_v"], lw=1, color="C" + str(i + 4))

    ax3.semilogx(
        cases[name]["y+"],
        cases[name]["f_w"],
        lw=1,
        label=name,
        color="C" + str(i + 4),
    )


# ax1.plot(simson_ho_sd[:, 1], simson_ho[:, 5], 'k', lw=1,)
# ax2.plot(simson_ho_sd[:, 1], simson_ho[:, 6], 'k', lw=1,)
# ax3.plot(simson_ho_sd[:, 1], simson_ho[:, 7], 'k', label="DNS", lw=1,)

ax3.legend()

ax1.set_xlim(0.1, 1000)
ax2.set_xlim(0.1, 1000)
ax3.set_xlim(0.1, 1000)
ax1.set_title(r"$F_u$")
ax2.set_title(r"$F_v$")
ax3.set_title(r"$F_w$")
ax1.set_xlabel(r"$y^+$")
ax2.set_xlabel(r"$y^+$")
ax3.set_xlabel(r"$y^+$")
plt.tight_layout()

plt.savefig(join(SAVE_PATH, "flatness.pdf"), bbox_inches="tight", pad_inches=0.02)
