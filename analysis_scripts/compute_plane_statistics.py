"""Computes statistics of probe data, includng skewness and flatness."""

import matplotlib.pyplot as plt
from os.path import join
import numpy as np
import h5py
from scipy.stats import skew


# Path with the simulation cases
DATA = ".."

# Case directories
case_names = ["M1", "M2", "M3"]

# Path for saving figures
SAVE_PATH = "figures"

FIGSIZE = (5, 2.5)

# Viscosity
nu = 5e-5

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
        cases[name]["y+"] = cases[name]["y"] / cases[name]["delta_nu"]
    return cases


cases = read_data()


# %%

for mesh in ["M1", "M2", "M3"]:
    print(mesh)
    data = h5py.File(join(DATA, mesh, f"probes_u.hdf5"), "r")

    t = data["times"][:].flatten()
    locs = data["locations"][:]

    y = np.unique(locs[:, 1])
    z = np.unique(locs[:, 2])
    ny = y.size
    nz = z.size
    nt = t.size

    stats = {}
    for compi, comp in enumerate(["u", "v", "w"]):

        stats[comp] = {}

        u_data = data["data"][:, :, compi].reshape(nt, ny, nz)

        u = np.mean(u_data, axis=0)
        u = np.mean(u, axis=1)

        if compi == 0:
            u_data -= u[None, :, None]

        var_u = np.var(u_data, axis=0)
        var_u = np.mean(var_u, axis=1)

        s_u = skew(u_data, axis=0)
        s_u = np.mean(s_u, axis=1)

        f_u = np.mean(u_data**4, axis=0)
        f_u = np.mean(f_u, axis=1)
        f_u /= var_u**2

        stats[comp]["mean"] = u
        stats[comp]["variance"] = var_u
        stats[comp]["skewness"] = s_u
        stats[comp]["flatness"] = f_u

    np.savetxt(
        join(DATA, mesh, "plane_stats.txt"),
        np.column_stack(
            (
                stats["u"]["mean"],
                stats["u"]["variance"],
                stats["v"]["variance"],
                stats["w"]["variance"],
                stats["u"]["skewness"],
                stats["v"]["skewness"],
                stats["w"]["skewness"],
                stats["u"]["flatness"],
                stats["v"]["flatness"],
                stats["w"]["flatness"],
            )
        ),
    )
