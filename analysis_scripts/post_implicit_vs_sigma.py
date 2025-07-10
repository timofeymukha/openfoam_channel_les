import matplotlib.pyplot as plt
from os.path import join
import numpy as np
import ofreaders
import h5py

plt.style.use("tableau-colorblind10")
plt.rc("text", usetex=True)
plt.rcParams.update({"font.size": 8})
plt.rcParams.update(
    {"axes.grid": True, "grid.color": "#DDDDDD", "grid.linestyle": "dashed"}
)

FIGSIZE = (5, 2.5)

# Path with the simulation cases
DATA = ".."

# Case directories
case_names = ["M1", "M1_sigma", "M2", "M2_sigma", "M1_rk", "M1_rk_sigma"]

colors = ["maroon", "salmon", "mediumblue", "deepskyblue", "darkgreen", "palegreen"]

labels = [
    "M1, Implicit",
    "M1, Sigma",
    "M2, Implicit",
    "M2, Sigma",
    "M1, RKSymFoam, Implicit",
    "M1, RKSymFoam, Sigma",
]

# Path for saving figures
SAVE_PATH = "../figures"

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

dns_spectra = h5py.File(join(DATA, "dns", "LM_Channel_1000_1d_energy_spectra.h5"), "r")

dns_lz = dns_spectra["Lz"][0]
dns_lx = dns_spectra["Lx"][0]
dns_yplus = dns_spectra["Y_plus"][:]
dns_kz = dns_spectra["kz"][:]
dns_z_spectrum_u = dns_spectra["Euu_kz"][:]
dns_z_spectrum_v = dns_spectra["Evv_kz"][:]
dns_z_spectrum_w = dns_spectra["Eww_kz"][:]

dns_utau = 5.00256e-02

dns_utau = 5.00256e-02
# %%


def read_data():
    times = ["750", "1500", "750", "750", "750", "750"]
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

    cases["M2_sigma"]["nut"] = np.genfromtxt(
        join(DATA, "M2_sigma", "postProcessing", "collapsedFields", "750", "nutMean.xy")
    )[:, 1]
    cases["M1_sigma"]["nut"] = np.genfromtxt(
        join(
            DATA, "M1_sigma", "postProcessing", "collapsedFields", "1500", "nutMean.xy"
        )
    )[:, 1]
    cases["M1_rk_sigma"]["nut"] = np.genfromtxt(
        join(
            DATA,
            "M1_rk_sigma",
            "postProcessing",
            "collapsedFields",
            "750",
            "nutMean.xy",
        )
    )[:, 1]
    return cases


cases = read_data()

# %% utau errors


def print_utau_errors():
    print("u_tau errors")
    for i, name in enumerate(case_names):
        error = (cases[name]["utau"] - dns_utau) / dns_utau * 100
        print(f"{name}: {error}")

    print("tau_w errors")
    for i, name in enumerate(case_names):
        error = (cases[name]["utau"] ** 2 - dns_utau**2) / dns_utau**2 * 100
        print(f"{name}: {error}")


print_utau_errors()


# %% Reynolds stresses
def re_stresses():
    global labels
    fig = plt.figure(figsize=(FIGSIZE))
    # First subplot
    ax1 = plt.subplot(121)
    for i, name in enumerate(case_names):
        plt.semilogx(
            cases[name]["y+"],
            cases[name]["uu"] / cases[name]["utau"] ** 2,
            lw=1,
            color=colors[i],
            label=labels[i],
        )
        plt.plot(
            cases[name]["y+"],
            cases[name]["uv"] / cases[name]["utau"] ** 2,
            "--",
            lw=1,
            color=colors[i],
        )

    plt.plot(dns_mean[:, 1], dns_fluct[:, 2], "k", lw=1, label="DNS")
    plt.plot(dns_mean[:, 1], dns_fluct[:, 5], "--k", lw=1)

    plt.xlim(1, 1000)
    plt.ylim(-1, 13)
    plt.ylabel(r"$\overline {u'u'}^+$, $\overline {u'v'}^+$")
    plt.xlabel(r"$y^+$")

    # Second subplot
    ax2 = plt.subplot(122)
    for i, name in enumerate(case_names):
        plt.semilogx(
            cases[name]["y+"],
            cases[name]["vv"] / cases[name]["utau"] ** 2,
            "--",
            lw=1,
            color=colors[i],
        )
        plt.plot(
            cases[name]["y+"],
            cases[name]["ww"] / cases[name]["utau"] ** 2,
            "-.",
            lw=1,
            color=colors[i],
        )

    plt.plot(dns_mean[:, 1], dns_fluct[:, 3], "--k", lw=1)
    plt.plot(dns_mean[:, 1], dns_fluct[:, 4], "-.k", lw=1)

    plt.xlim(1, 1000)
    plt.ylim(0, 2.5)
    plt.ylabel(r"$\overline {v'v'}^+$, $\overline {w'w'}^+$")
    plt.xlabel(r"$y^+$")
    plt.tight_layout()

    # Combined legend outside of the subplots
    handles, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.1))

    plt.subplots_adjust(bottom=0.5)  # Make space for the legend

    plt.savefig(join(SAVE_PATH, "re_sigma.pdf"), bbox_inches="tight", pad_inches=0.02)


re_stresses()
# %% Numerical dissipation


def dissipation():
    fig = plt.figure(figsize=(5, 3))
    labels_ = ["M1, Implicit", "M2, Implicit", "M1, RKSymFoam"]
    for i, name in enumerate(["M1", "M2", "M1_rk"]):

        nu_num = -nu * cases[name]["eps_num"] / cases[name]["dissipation"]

        plt.subplot(111)
        plt.semilogx(
            cases[name]["y+"][1:-1],
            np.abs(nu_num / nu),
            lw=1,
            color=f"C{i+4}",
            linestyle="-",
            label=labels_[i] + r", $|\nu_\mathrm{num}| / \nu$",
        )

    plt.plot(
        cases["M1_sigma"]["y+"],
        cases["M1_sigma"]["nut"] / nu,
        "--",
        color="C4",
        lw=1,
        label=r"M1, Sigma, $\nu_\mathrm{sgs} / \nu$",
    )
    plt.semilogy(
        cases["M2_sigma"]["y+"],
        cases["M2_sigma"]["nut"] / nu,
        "--",
        color=f"C5",
        lw=1,
        label=r"M2 Sigma, $\nu_\mathrm{sgs} / \nu$",
    )
    plt.plot(
        cases["M1_rk_sigma"]["y+"],
        cases["M1_rk_sigma"]["nut"] / nu,
        "--",
        color="C6",
        lw=1,
        label=r"M1, RKSymFoam, Sigma, $\nu_\mathrm{sgs} / \nu$",
    )
    plt.xlim(1, 1000)
    plt.ylim(1e-5, 2)
    plt.xlabel(r"$y^+$")
    plt.ylabel(r"Viscosity ratio")
    plt.tight_layout()

    handles, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.1))

    plt.subplots_adjust(bottom=0.45)


dissipation()
