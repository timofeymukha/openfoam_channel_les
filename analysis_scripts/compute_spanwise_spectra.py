"""Shows how to compute the 1D energy spectrum in z"""

from os.path import join
import numpy as np
import h5py
from scipy.signal import periodogram

# %%
DATA = "."

SAVE_PATH = "figures"

# Case directories
case_names = ["M1", "M2", "M3", "M2_sigma"]
# Viscosity
nu = 5e-5


def read_data():
    times = ["750", "750", "750", "750"]
    cases = {}

    for i, name in enumerate(case_names):
        cases[name] = {}
        time = times[i]

        cases[name]["y"] = np.genfromtxt(
            join(DATA, name, "postProcessing", "collapsedFields", time, "UMean_X.xy")
        )[:, 0]
        cases[name]["u"] = np.genfromtxt(
            join(DATA, name, "postProcessing", "collapsedFields", time, "UMean_X.xy")
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

# %%

# Select the case
mesh = "M2_sigma"

for compi, comp in enumerate(["u", "v", "w"]):
    uhat = []
    for f in [join(DATA, mesh, "probes_u.hdf5")]:
        data = h5py.File(f, "r")

        t = data["times"][:]
        print(t[0], t[-1])
        locs = data["locations"][:]

        y = np.unique(locs[:, 1])
        z = np.unique(locs[:, 2])
        ny = y.size
        nz = z.size
        nt = t.size
        lz = 6

        # Make sure we are divisable by 10 with no remainder
        nt = nt - (nt % 10)

        chunk_size = nt // 10  # Size of each chunk

        # Iterate over the dataset in strides of 1/10th the size
        for i in range(0, nt, chunk_size):
            print(i)
            # Define the slice to load a chunk of the dataset
            uprime = data["data"][i : i + chunk_size, :, compi].reshape(
                chunk_size, ny, nz
            )

            if compi == 0:
                uprime -= cases[mesh]["u"][np.newaxis, :ny, np.newaxis]

            f, ux = periodogram(
                uprime, fs=nz / lz, scaling="density", axis=2, detrend=False
            )

            uhat.append(np.mean(ux, axis=0))

    # Average across the chunks
    result = np.zeros(uhat[0].shape)
    for i in range(len(uhat)):
        result += uhat[i]

    result = result / len(uhat)
    np.savetxt(join(DATA, mesh, f"z_spectrum_{comp}.txt"), result)
