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
DATA = "."

# Case directories
case_names = ["M2", "M2_sigma"]

labels = ["Implicit LES", "Sigma LES"]

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

    cases["M2_sigma"]["nut"] = np.genfromtxt(
        join(DATA, name, "postProcessing", "collapsedFields", "750", "nutMean.xy")
    )[:, 1]
    return cases


cases = read_data()

# %% utau errors

print("u_tau errors")
for i, name in enumerate(case_names):
    error = (cases[name]["utau"] - dns_utau) / dns_utau * 100
    print(f"{name}: {error}")

print("tau_w errors")
for i, name in enumerate(case_names):
    error = (cases[name]["utau"] ** 2 - dns_utau**2) / dns_utau**2 * 100
    print(f"{name}: {error}")


# %% Reynolds stresses
def re_stresses():

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
    plt.savefig(join(SAVE_PATH, "re_sigma.pdf"), bbox_inches="tight", pad_inches=0.02)


re_stresses()
# %% Numerical dissipation


def dissipation():
    plt.figure(figsize=FIGSIZE)

    for i, name in enumerate(["M2"]):

        plt.subplot(121)
        plt.semilogx(
            cases[name]["y+"][1:-1],
            cases[name]["eps_num"] / cases[name]["utau"] ** 4 * nu,
            lw=1,
            color=f"C{i+4}",
            linestyle="-",
            label=labels[i],
        )

    plt.xlim(1, 1000)
    plt.legend()
    plt.xlabel(r"$y^+$")
    plt.ylabel(r"$k$ budget residual in inner scaling, $\epsilon^+_\mathrm{num}$")

    for i, name in enumerate(["M2"]):

        nu_num = -nu * cases[name]["eps_num"] / cases[name]["dissipation"]

        plt.subplot(122)
        plt.semilogx(
            cases[name]["y+"][1:-1],
            nu_num / nu,
            lw=1,
            color=f"C{i+4}",
            linestyle="-",
            label=labels[i] + r", $\nu_\mathrm{num} / \nu$",
        )

    plt.plot(
        cases["M2_sigma"]["y+"],
        cases["M2_sigma"]["nut"] / nu,
        color=f"C{i+5}",
        lw=1,
        label=r"Sigma, $\nu_\mathrm{sgs} / \nu$",
    )

    plt.xlim(0.2, 1000)
    plt.ylim(-0.01, 0.15)
    plt.xlabel(r"$y^+$")
    plt.ylabel(r"Viscosity ratio")
    plt.legend()
    plt.tight_layout()

    plt.savefig(join(SAVE_PATH, "nut_sigma.pdf"), bbox_inches="tight", pad_inches=0.02)


dissipation()


# %%
def tke_budget():
    plt.figure(figsize=FIGSIZE)

    for i, case in enumerate(case_names):

        plt.semilogx(
            cases[case]["y+"][1:-1],
            cases[case]["prod"] / cases[case]["utau"] ** 4 * nu,
            color=f"C{i+4}",
            lw=1,
            label=labels[i],
        )
        plt.semilogx(
            cases[case]["y+"][1:-1],
            -cases[case]["dissipation"] / cases[case]["utau"] ** 4 * nu,
            "--",
            color=f"C{i+4}",
            lw=1,
        )
        plt.semilogx(
            cases[case]["y+"][1:-1],
            cases[case]["diffusion"] / cases[case]["utau"] ** 4 * nu,
            "-.",
            color=f"C{i+4}",
            lw=1,
        )
        plt.semilogx(
            cases[case]["y+"][1:-1],
            cases[case]["transport"] / cases[case]["utau"] ** 4 * nu,
            lw=1,
            color=f"C{i+4}",
            linestyle=(5, (10, 3)),
        )
        plt.semilogx(
            cases[case]["y+"][1:-1],
            cases[case]["p_diffusion"] / cases[case]["utau"] ** 4 * nu,
            lw=2,
            color=f"C{i+4}",
            linestyle=":",
        )

    plt.semilogx(dns_mean[:, 1], dns_budget[:, 2], "k", lw=1)
    plt.plot(dns_mean[:, 1], dns_budget[:, -2], "--k", lw=1)
    plt.plot(dns_mean[:, 1], dns_budget[:, 4], "-.k", lw=1)
    plt.plot(dns_mean[:, 1], dns_budget[:, 3], "k", lw=1, linestyle=(5, (10, 3)))

    plt.plot(dns_mean[:, 1], dns_budget[:, 6], ":k", lw=2)
    plt.xlim(1, 1000)
    plt.ylim(-0.15, 0.25)
    plt.xlabel(r"$y^+$")
    plt.ylabel(r"TKE budget terms")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        join(SAVE_PATH, "budget_sigma.pdf"),
        bbox_inches="tight",
        pad_inches=0.02,
    )


tke_budget()
