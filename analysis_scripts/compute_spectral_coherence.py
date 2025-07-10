""" Shows how to compute the linear spectral coherence.

"""

import matplotlib.pyplot as plt
from os.path import join
import numpy as np
import h5py

# %%
DATA = ".."

nu = 5e-5

data = h5py.File(join(DATA, "M1", "probes_u.hdf5"), "r")

t = data["times"][:]
locs = data["locations"][:]

y = np.unique(locs[:, 1])
z = np.unique(locs[:, 2])

ny = y.size
nz = z.size
nt = t.size

utau = 5.00256e-02
yplus = y * utau / nu

lz = 6
dz = lz / nz


# Load data and reshape to time x y x z
data = data["data"][:, :, 1].reshape(nt, ny, nz)

# demean u
umean = np.mean(data, axis=0)
umean = np.mean(umean, axis=1)
data = data - umean[np.newaxis, :, np.newaxis]

# %% FFT in z
uhat = np.fft.fft(data, axis=2)

# The frequency and wave lengths. We use rrftfreq to not include the negative
# ones. This will be used for plotting.
freq = np.fft.rfftfreq(nz, dz)

# We can assert that the largest lengthscale is = lz and smallest = 2*dz
# as per Nyquist
lambda_z = 1 / freq[1:]

lambda_zplus = lambda_z * utau / nu

# %%

# Reference y+ value
yp_ref = 2

# Index of ref y+
ref_index = np.argmin(np.abs(yplus - yp_ref))

# Reference uhat
ref_uhat = uhat[:, ref_index, :]

# Nominator for the coherence. Broadcasts multiplication by ref_uhat* to all
# y values and takes the mean in time
nominator = np.mean(uhat * np.conjugate(ref_uhat)[:, np.newaxis, :], axis=0)

nominator = np.abs(nominator) ** 2

denominator = (
    np.mean(np.abs(uhat) ** 2, axis=0)
    * np.mean(np.abs(ref_uhat) ** 2, axis=0)[np.newaxis, :]
)

coherence = nominator / denominator

# The coherence at this point is symmetric around
coherence = coherence[:, : int(nz / 2)] + coherence[:, int(nz / 2) :][:, ::-1]
coherence *= 0.5
