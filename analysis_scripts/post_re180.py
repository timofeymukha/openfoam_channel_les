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

DATA = "../re180_simulations"

case_names = ["M1_re180", "M1_re180_sigma", "M1_re180_smag"]
colors = ["C4", "C5", "C6"]

SAVE_PATH = "../figures"

nu = 3.5e-4

labels = ["Implicit", "Sigma", "vD Smagorinsky"]

# %% Load DNS data

dns = {}
dns["mean"] = np.genfromtxt(
    join(DATA, "..", "dns", "LM_Channel_0180_mean_prof.dat"), comments="%"
)
dns["fluct"] = np.genfromtxt(
    join(DATA, "..", "dns", "LM_Channel_0180_vel_fluc_prof.dat"), comments="%"
)
dns["omega"] = np.genfromtxt(
    join(DATA, "..", "dns", "LM_Channel_0180_vor_pres_fluc_prof.dat"), comments="%"
)
dns["budget"] = np.genfromtxt(
    join(DATA, "..", "dns", "LM_Channel_0180_RSTE_k_prof.dat"), comments="%"
)

dns["utau"] = 6.37309e-02

# %%


def read_data():
    times = ["3100", "1450", "1850"]
    cases = {}

    for i, name in enumerate(case_names):
        print(name)

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

        cases[name]["oxox"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "VorticityPrime2Mean_XX.xy",
            )
        )[:, 1]
        cases[name]["oyoy"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "VorticityPrime2Mean_YY.xy",
            )
        )[:, 1]
        cases[name]["ozoz"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "VorticityPrime2Mean_ZZ.xy",
            )
        )[:, 1]
        cases[name]["oxoy"] = np.genfromtxt(
            join(
                DATA,
                name,
                "postProcessing",
                "collapsedFields",
                time,
                "VorticityPrime2Mean_XY.xy",
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

        if "_sigma" in name or "vd" in name:
            cases[name]["nut"] = np.genfromtxt(
                join(DATA, name, "postProcessing", "collapsedFields", time, "nut.xy")
            )[:, 1]

    return cases


cases = read_data()


# %% utau errors


def print_u_tau_errors():
    print("u_tau errors")
    for i, name in enumerate(case_names):
        error = (cases[name]["utau"] - dns["utau"]) / dns["utau"] * 100
        print(f"{name}: {error}")

    print("tau_w errors")
    for i, name in enumerate(case_names):
        error = (cases[name]["utau"] ** 2 - dns["utau"] ** 2) / dns["utau"] ** 2 * 100
        print(f"{name}: {error}")


print_u_tau_errors()

# %%


def plot_re_stresses():
    global labels
    fig = plt.figure(figsize=FIGSIZE)
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

    plt.plot(dns["mean"][:, 1], dns["fluct"][:, 2], "k", lw=1, label="DNS")
    plt.plot(dns["mean"][:, 1], dns["fluct"][:, 5], "--k", lw=1)
    plt.xlim(1, 180)
    plt.ylim(-1, 8)
    plt.ylabel(r"$\overline {u'u'}^+$, $\overline {u'v'}^+$")
    plt.xlabel(r"$y^+$")

    plt.subplot(122)
    for i, name in enumerate(case_names):
        plt.semilogx(
            cases[name]["y+"],
            cases[name]["vv"] / cases[name]["utau"] ** 2,
            ":",
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

    plt.plot(dns["mean"][:, 1], dns["fluct"][:, 3], ":k", lw=1)
    plt.plot(dns["mean"][:, 1], dns["fluct"][:, 4], "-.k", lw=1)
    plt.xlim(1, 180)
    # plt.ylim(0, 2.5)
    plt.ylabel(r"$\overline {v'v'}^+$, $\overline {w'w'}^+$")
    plt.xlabel(r"$y^+$")

    plt.tight_layout()

    handles, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.1))

    plt.subplots_adjust(bottom=0.35)  # Make space for the legend

    plt.savefig(join(SAVE_PATH, "re180_re.pdf"), bbox_inches="tight", pad_inches=0.02)


plot_re_stresses()


# %%
def plot_dissipation():
    plt.figure(figsize=FIGSIZE)

    for i, name in enumerate(["M1_re180"]):
        plt.subplot(121)
        plt.plot(
            cases[name]["y+"][1:-1],
            cases[name]["eps_num"] / cases[name]["utau"] ** 4 * nu,
            lw=1,
            color=f"C{i+4}",
            linestyle="-",
            label=labels[i],
        )

    plt.plot(dns["mean"][:, 1], dns["budget"][:, -1], "k", label="DNS", lw=1)
    plt.legend()
    plt.xlabel(r"$y^+$")
    plt.xlim(1, 180)
    plt.ylabel(r"$\epsilon^+_\mathrm{num}$")

    for i, name in enumerate(["M1_re180"]):
        nu_num = -nu * cases[name]["eps_num"] / cases[name]["dissipation"]

        plt.subplot(122)
        plt.plot(
            cases[name]["y+"][1:-1],
            nu_num / nu,
            lw=1,
            color=f"C{i+4}",
            linestyle="-",
            label=labels[i] + r", $\nu_\mathrm{num} / \nu$",
        )

    plt.plot(
        cases["M1_re180_smag"]["y+"][1:-1],
        cases["M1_re180_smag"]["nut"][1:-1] / nu,
        lw=1,
        color=f"C6",
        linestyle="-",
        label=r"Smagorinsky, $\nu_\mathrm{sgs} / \nu$",
    )

    plt.plot(
        cases["M1_re180_sigma"]["y+"][1:-1],
        cases["M1_re180_sigma"]["nut"][1:-1] / nu,
        lw=1,
        color=f"C5",
        linestyle="-",
        label=r"Sigma, $\nu_\mathrm{sgs} / \nu$",
    )
    plt.ylabel("Viscosity ratio")
    plt.xlim(1, 180)
    plt.legend()
    plt.xlabel(r"$y^+$")
    plt.tight_layout()

    plt.savefig(
        join(SAVE_PATH, "eps_num_re180.pdf"), bbox_inches="tight", pad_inches=0.02
    )


plot_dissipation()
