"""
    Demonstrates a few ways to compute two-point correlations in space.
"""

import matplotlib.pyplot as plt
from os.path import join
import numpy as np
import h5py
from tqdm import tqdm

plt.style.use("tableau-colorblind10")
plt.rc("text", usetex=True)
plt.rcParams.update({"font.size": 8})
plt.rcParams.update(
    {"axes.grid": True, "grid.color": "#DDDDDD", "grid.linestyle": "dashed"}
)


# Path with the simulation cases
DATA = "."

# Case directories
case_names = [
    "M1",
    "M2",
    "M3",
]

# Path for saving figures
SAVE_PATH = "figures"

# Viscosity
nu = 5e-5

# Nominal friction velocity
utau = 0.05

# %%
# Case directories
case_names = ["M1", "M2", "M3"]


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
    return cases


cases = read_data()

# %% Using FFT, one time step at a time

mesh = "M1"

for compi, comp in enumerate(["u", "v", "w"]):

    data = h5py.File(join(DATA, mesh, f"probes_u.hdf5"), "r")

    t = data["times"][:].flatten()
    print("Time range", t[0], t[-1])
    locs = data["locations"][:]

    y = np.unique(locs[:, 1])
    yplus = y * utau / nu
    z = np.unique(locs[:, 2])
    ny = y.size
    nz = z.size
    nt = t.size

    tpcorr = np.zeros((ny, int(nz / 2)))
    d = nz * np.ones(2 * nz - 1)

    for ti in tqdm(range(nt)):
        # Reshape
        uprime = data["data"][ti, :, compi].reshape(ny, nz)

        # Demean
        if compi == 0:
            uprime -= cases[mesh]["u"][1 : ny + 1, None]

        # Take fft in z
        Frf = np.fft.fft(uprime, n=nz, axis=1)

        # Compute power spectrum and inverse fft to get the autocovariance
        acov = np.fft.ifft(Frf * np.conjugate(Frf), axis=1)[:nz] / nz
        acov = acov.real
        tpcorr += acov[:, : nz // 2]

    tpcorr = tpcorr / tpcorr[:, 0][:, None]

    np.savetxt(join(DATA, mesh, f"z_tpcorr_{comp}_fft.txt"), tpcorr)

# %% Using FFT loading the whole data in memory.
# Faster, but needs a lot of RAM (Terabytes).

mesh = "M1"

for compi, comp in enumerate(["u", "v", "w"]):

    data = h5py.File(join(DATA, mesh, f"probes_u.hdf5"), "r")

    t = data["times"][:].flatten()
    print("Time range", t[0], t[-1])
    locs = data["locations"][:]

    y = np.unique(locs[:, 1])
    z = np.unique(locs[:, 2])
    ny = y.size
    nz = z.size
    nt = t.size

    uprime = data["data"][:, :, compi].reshape(nt, ny, nz)

    if compi == 0:
        uprime -= cases[mesh]["u"][None, 1 : ny + 1, None]

    Frf = np.fft.fft(uprime, n=nz, axis=2)

    tpcorr = np.fft.ifft(Frf * np.conjugate(Frf), axis=2)[..., :nz]
    tpcorr = tpcorr.real
    tpcorr = tpcorr[:, :, : nz // 2]
    tpcorr = np.sum(tpcorr, axis=0)
    tpcorr /= tpcorr[:, 0][:, None]

    np.savetxt(join(DATA, mesh, f"z_tpcorr_{comp}_fullfft.txt"), tpcorr)

# %% Does not use FFT, instead rolls the array in z and averages among the
# realizations with different reference points.
# Supports dividing the data into (big) chunks along time when processing,
# but still loads the entire component velocity array in memory.

mesh = "M1"

# Number of chunks to dvidide the data in. Ideally set to one to put everything
# in memory
n_chunks = 1

for compi, comp in enumerate(["u"]):  # , "v", "w"]):

    # Data from different chunks

    # Will hold the computed correlation for each chunk in a list
    tpcorr_list = []

    for f in [join(DATA, mesh, "probes_u.hdf5")]:
        data = h5py.File(f, "r")

        t = data["times"][:]
        print("Time range", t[0], t[-1])
        locs = data["locations"][:]

        y = np.unique(locs[:, 1])
        z = np.unique(locs[:, 2])
        ny = y.size
        nz = z.size
        nt = t.size
        lz = 6

        # Cut off some data to make sure we are perfrectly divisable by n_chunks
        nt = nt - (nt % n_chunks)
        print("Number of timesteps", nt)

        chunk_size = nt // n_chunks  # Size of each chunk

        uprime = data["data"][:, :, compi].reshape(nt, ny, nz)

        if compi == 0:
            # demean
            uprime -= cases[mesh]["u"][None, 1 : ny + 1, None]

        # Iterate over the dataset in chunk-size strides
        for chnki in range(0, nt, chunk_size):

            # Define the slice to load a chunk of the dataset

            rolled = uprime[chnki : chnki + chunk_size, :, :]
            tpcorr = np.zeros((ny, nz))

            for i in tqdm(range(nz - 1)):

                # Broadcast the multiplication of the 0th column across all
                # and take the mean in time.
                tpcorr += np.mean(rolled[:, :, 0][:, :, None] * rolled[:, :, :], axis=0)

                # Roll to get a new reference point
                rolled = np.roll(rolled, 1, axis=2)

            tpcorr /= nz - 1

            tpcorr_list.append(tpcorr)

    # Average across the chunks
    result = np.zeros(tpcorr_list[0].shape)
    for i in range(len(tpcorr_list)):
        result += tpcorr_list[i]

    result = result / len(tpcorr_list)

    result[:, 1 : int(nz / 2)] += result[:, -1 : -int(nz / 2) : -1]
    result[:, 1 : int(nz / 2)] /= 2
    result = result[:, : int(nz / 2)]

    # Normalize to get correlation and not covariance
    for i in range(result.shape[0]):
        result[i, :] /= result[i, 0]

    np.savetxt(join(DATA, mesh, f"z_tpcorr_{comp}_roll.txt"), result)
