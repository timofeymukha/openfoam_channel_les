import numpy as np
import matplotlib.pyplot as plt
from os.path import join

# Set paths
CASE = "E:/OpenFOAM/timofey-v2306/run/channel_flow/openfoam_channel_les"
PATH = join(CASE, "postProcessing/surfaces/0/inlet/faceCentres")


with open(PATH) as points_file:
    points = [line.rstrip(")\n") for line in points_file]

points = [line.lstrip("(") for line in points]
points = points[3:-1]
points = np.genfromtxt(points)[:, 1:]

# Remove upper half
points = points[np.where(points[:, 0] < 1)[0], :]

# Get y and z values for the mesh
y = np.unique(points[:, 0])
z = np.unique(points[:, 1])

Y, Z = np.meshgrid(y, z, indexing="ij")


points = np.column_stack((Y.flatten(), Z.flatten()))


with open(join(CASE, "system/probes"), "w", newline="\n") as f:
    f.write("probes" + "\n")
    f.write("{\n")
    f.write("    type probes;\n")
    f.write("    name probes;\n")
    f.write('    libs ("sampling");\n')
    f.write("    executionInterval   5;\n")
    f.write("    executionControl   timeStep;\n")
    f.write("    writeInterval   5;\n")
    f.write("    writeControl   timeStep;\n")
    f.write("    probeLocations\n")
    f.write("    (\n")
    for i in range(points.shape[0]):
        f.write(
            "    ("
            + str(1e-3)
            + " "
            + str(points[i, 0])
            + " "
            + str(points[i, 1])
            + ")\n"
        )
    f.write("    );\n")
    f.write("    fields (U p);\n")
    f.write("}\n")
