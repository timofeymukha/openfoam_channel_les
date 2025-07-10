""" Provides an example of how to use ofreaders to convert the ascii probe file
to hdf5.

"""

import numpy as np
from os.path import join
import ofreaders
import h5py

# %%

DATA = ".."


# %%

u = ofreaders.Probes.read(join(DATA, "M1", "postProcessing", "probes", "250", "U"))
u.save_to_hdf5(join(DATA, "M1", "probes_u.hdf5"))

p = ofreaders.Probes.read(join(DATA, "M1", "postProcessing", "probes", "250", "p"))
p.save_to_hdf5(join(DATA, "M1", "probes_p.hdf5"))
