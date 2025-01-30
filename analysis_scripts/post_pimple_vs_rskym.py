import matplotlib.pyplot as plt
from os.path import join
import numpy as np
import ofreaders

plt.style.use("tableau-colorblind10")
plt.rc("text", usetex=True)
plt.rcParams.update({"font.size": 8})
plt.rcParams.update(
    {"axes.grid": True, "grid.color": "#DDDDDD", "grid.linestyle": "dashed"}
)

FIGSIZE = (5, 2.5)

# Path with the simulation cases
DATA = "."

# Case directories
case_names = ["M2", "M2_rk"]

labels = ["pimpleFoam", "RKSymFoam"]

# Path for saving figures
SAVE_PATH = "figures"

# Viscosity
nu = 5e-5
# %% Load DNS data
dns_mean = np.genfromtxt(
    join(DATA, "dns", "LM_Channel_1000_mean_prof.dat"), comments="%"
)
dns_fluct = np.genfromtxt(
    join(DATA, "dns", "LM_Channel_1000_vel_fluc_prof.dat"), comments="%"
)
dns_omega = np.genfromtxt(
    join(DATA, "dns", "LM_Channel_1000_vor_pres_fluc_prof.dat"), comments="%"
)
dns_budget = np.genfromtxt(
    join(DATA, "dns", "LM_Channel_1000_RSTE_k_prof.dat"), comments="%"
)

dns_utau = 5.00256e-02
# %%


def read_data():
    times = ["750", "750"]
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

        cases[name]["prod"] = np.genfromtxt(
            join(DATA, name, "graphs", time, "prod.xy")
        )[:, 1]
        cases[name]["transport"] = np.genfromtxt(
            join(DATA, name, "graphs", time, "trans.xy")
        )[:, 1]
        cases[name]["dissipation"] = np.genfromtxt(
            join(DATA, name, "graphs", time, "diss.xy")
        )[:, 1]
        cases[name]["diffusion"] = np.genfromtxt(
            join(DATA, name, "graphs", time, "vDiff.xy")
        )[:, 1]
        cases[name]["p_diffusion"] = np.genfromtxt(
            join(DATA, name, "graphs", time, "pDiff.xy")
        )[:, 1]

        cases[name]["eps_num"] = (
            cases[name]["prod"][:]
            + cases[name]["dissipation"]
            + cases[name]["diffusion"]
            + cases[name]["transport"]
            + cases[name]["p_diffusion"]
        )
    return cases


cases = read_data()

# %% utau errors
for i, name in enumerate(case_names):
    error = (cases[name]["utau"] - dns_utau) / dns_utau * 100
    print(f"{name}: {error}")


for i, name in enumerate(case_names):
    error = (cases[name]["utau"] ** 2 - dns_utau**2) / dns_utau**2 * 100
    print(f"{name}: {error}")


# %%
plt.figure(figsize=FIGSIZE)
plt.subplot(121)
for i, name in enumerate(case_names):
    plt.semilogx(
        cases[name]["y+"],
        cases[name]["uu"] / cases[name]["utau"] ** 2,
        lw=1,
        color="C" + str(i + 4),
        label=labels[i],
    )

    plt.plot(
        cases[name]["y+"],
        cases[name]["uv"] / cases[name]["utau"] ** 2,
        "--",
        lw=1,
        color="C" + str(i + 4),
    )

plt.plot(dns_mean[:, 1], dns_fluct[:, 2], "k", lw=1, label="DNS")
plt.plot(dns_mean[:, 1], dns_fluct[:, 5], "--k", lw=1)
plt.xlim(1, 1000)
plt.ylim(-1, 12)
plt.ylabel(r"$\overline {u'u'}$, $\overline {u'v'}$")
plt.xlabel(r"$y^+$")
plt.legend()

plt.subplot(122)
for i, name in enumerate(case_names):
    plt.semilogx(
        cases[name]["y+"],
        cases[name]["vv"] / cases[name]["utau"] ** 2,
        "--",
        lw=1,
        color="C" + str(i + 4),
    )
    plt.plot(
        cases[name]["y+"],
        cases[name]["ww"] / cases[name]["utau"] ** 2,
        "-.",
        lw=1,
        color="C" + str(i + 4),
    )


plt.plot(dns_mean[:, 1], dns_fluct[:, 3], "--k", lw=1)
plt.plot(dns_mean[:, 1], dns_fluct[:, 4], "-.k", lw=1)
plt.xlim(1, 1000)
plt.ylim(0, 2.5)
plt.ylabel(r"$\overline {v'v'}$, $\overline {w'w'}$")
plt.xlabel(r"$y^+$")

plt.tight_layout()
plt.savefig(join(SAVE_PATH, "re_rk.pdf"), bbox_inches="tight", pad_inches=0.02)
