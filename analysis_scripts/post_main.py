import matplotlib.pyplot as plt
from os.path import join
import numpy as np
import ofreaders
import h5py
from scipy.interpolate import interp1d

plt.style.use("tableau-colorblind10")
plt.rc("text", usetex=True)
plt.rcParams.update({"font.size": 8})
plt.rcParams.update(
    {"axes.grid": True, "grid.color": "#DDDDDD", "grid.linestyle": "dashed"}
)


# Path with the simulation cases
DATA = "."

# Case directories
case_names = ["M1", "M2", "M3"]

# Path for saving figures
SAVE_PATH = "figures"

FIGSIZE = (5, 2.5)

# Viscosity
nu = 5e-5
# %% Load DNS data
dns_mean = np.genfromtxt(
    join(DATA, "dns", "LM_Channel_1000_mean_prof.dat"), comments="%"
)
dns_yplus = dns_mean[:, 1]
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
dns_nz = dns_spectra["nz"][0]
dns_lx = dns_spectra["Lx"][0]
dns_yplus = dns_spectra["Y_plus"][:]
dns_kz = dns_spectra["kz"][:]
dns_z_spectrum_u = dns_spectra["Euu_kz"][:]
dns_z_spectrum_v = dns_spectra["Evv_kz"][:]
dns_z_spectrum_w = dns_spectra["Eww_kz"][:]

dns_utau = 5.00256e-02


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

        cases[name]["yp_evolution"] = ofreaders.collect_genfromtxt_over_time(
            join(DATA, name, "postProcessing", "averageYPlus"),
            "surfaceFieldValue.dat",
        )

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

        cases[name]["z_spectrum_u"] = np.genfromtxt(
            join(DATA, name, "z_spectrum_u.txt")
        )
        cases[name]["z_spectrum_v"] = np.genfromtxt(
            join(DATA, name, "z_spectrum_v.txt")
        )
        cases[name]["z_spectrum_w"] = np.genfromtxt(
            join(DATA, name, "z_spectrum_w.txt")
        )

        cases[name]["eps_num"] = (
            cases[name]["prod"][:]
            + cases[name]["dissipation"]
            + cases[name]["diffusion"]
            + cases[name]["transport"]
            + cases[name]["p_diffusion"]
        )

        cases[name]["z_tpcorr_u"] = np.genfromtxt(join(DATA, name, "z_tpcorr_u.txt"))
        cases[name]["z_tpcorr_v"] = np.genfromtxt(join(DATA, name, "z_tpcorr_v.txt"))
        cases[name]["z_tpcorr_w"] = np.genfromtxt(join(DATA, name, "z_tpcorr_w.txt"))

        cases[name]["coherence"] = np.genfromtxt(join(DATA, name, "coherence.txt"))

    return cases


cases = read_data()

# %% Grid resolution in the wall-normal direction


def delta_y():
    intrp = interp1d(dns_mean[:, 0], dns_kolmogorov)

    fig = plt.figure(figsize=(2.5, 2.5))

    ax1 = fig.add_subplot(111)
    ax1.semilogx(
        cases["M1"]["y+"][2:],
        cases["M1"]["y+"][2:] - cases["M1"]["y+"][1:-1],
        "-ok",
        ms=2,
    )
    ax1.set_xlim(0, 1000)
    ax1.set_xlabel(r"$y^+$")
    ax1.set_ylabel(r"$\Delta y^+$")
    ax1.set_ylim(0, 40)
    plt.tight_layout()

    plt.savefig(join(SAVE_PATH, "delta_y.pdf"), bbox_inches="tight", pad_inches=0.02)


delta_y()
# %% utau errors


def print_utau_errors():
    print("Error in u_tau")
    for i, name in enumerate(case_names):
        error = (cases[name]["utau"] - dns_utau) / dns_utau * 100
        print(f"{name}: {error}")

    print("Error in tau_w")
    for i, name in enumerate(case_names):
        error = (cases[name]["utau"] ** 2 - dns_utau**2) / dns_utau**2 * 100
        print(f"{name}: {error}")


print_utau_errors()

# %%


def utau_convergence():

    for i, name in enumerate(["M3"]):
        t_cut = cases[name]["yp_evolution"][:, 0]
        ind = np.argmin(np.abs(t_cut - 250))
        t_cut = t_cut[ind:]
        utau_signal = cases[name]["yp_evolution"][ind:, 1] * nu / cases[name]["y"][1]
        cases[name]["yp_mean"] = np.zeros(t_cut.size)

        for j in range(cases[name]["yp_mean"].size):
            cases[name]["yp_mean"][j] = np.mean(utau_signal[:j]) / dns_utau

    plt.figure(figsize=FIGSIZE)
    for i, name in enumerate(["M3"]):
        t = cases[name]["yp_evolution"][:, 0]
        ind = np.argmin(np.abs(t - 250))
        t_cut = t[ind:]
        utau_signal = (
            cases[name]["yp_evolution"][:, 1] * nu / cases[name]["y"][1] / dns_utau
        )

        plt.plot(
            t[100:] * dns_utau,
            utau_signal[100:],
            ":",
            lw=1,
            color="C" + str(i + 4),
            label="LES",
        )
        plt.plot(
            t_cut * dns_utau,
            cases[name]["yp_mean"],
            lw=1,
            label="LES, running average",
            color="C" + str(i + 4),
        )

    plt.hlines([1], xmin=0, xmax=40, colors="Black", lw=1)
    plt.vlines(
        [12.5],
        ymin=0.9,
        ymax=1.2,
        label="Averging start",
        colors="Black",
        lw=1,
        linestyles="dashed",
    )
    plt.legend()
    plt.xlim(0, 37)
    plt.ylim(0.94, 1.02)
    plt.xlabel(r"$t u^{dns}_\tau / \delta$")
    plt.tight_layout()
    plt.savefig(
        join(SAVE_PATH, "utau_convergence.pdf"),
        bbox_inches="tight",
        pad_inches=0.02,
    )


utau_convergence()

# %% Mean velocity


def mean_velocity():
    plt.figure(figsize=FIGSIZE)
    plt.subplot(121)
    for i, name in enumerate(case_names):
        plt.plot(
            cases[name]["y"],
            cases[name]["u"],
            lw=1,
            color="C" + str(i + 4),
            label=name,
        )

    plt.plot(dns_mean[:, 0], dns_mean[:, 2] * dns_utau, "k", label="DNS", lw=1)
    plt.xlim(0, 1)
    plt.ylim(0, 1.25)
    plt.ylabel(r"$\bar u/U_b$")
    plt.xlabel(r"$y/\delta$")
    plt.legend()

    plt.subplot(122)
    for i, name in enumerate(case_names):
        plt.semilogx(
            cases[name]["y+"],
            cases[name]["u"] / cases[name]["utau"],
            lw=1,
            color="C" + str(i + 4),
        )

    plt.plot(dns_mean[:, 1], dns_mean[:, 2], "k", label="DNS", lw=1)
    plt.xlim(1, 1000)
    plt.ylim(0, 25)
    plt.xlabel(r"$y^+$")
    plt.ylabel(r"$\bar u^+$")
    plt.tight_layout()

    plt.savefig(join(SAVE_PATH, "u.pdf"), bbox_inches="tight", pad_inches=0.02)


mean_velocity()

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
            label=name,
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

    plt.plot(dns_mean[:, 1], dns_fluct[:, 3], ":k", lw=1)
    plt.plot(dns_mean[:, 1], dns_fluct[:, 4], "-.k", lw=1)
    plt.xlim(1, 1000)
    plt.ylim(0, 2.5)
    plt.ylabel(r"$\overline {v'v'}$, $\overline {w'w'}$")
    plt.xlabel(r"$y^+$")

    plt.tight_layout()
    plt.savefig(join(SAVE_PATH, "re.pdf"), bbox_inches="tight", pad_inches=0.02)


re_stresses()

# %% Vorticity fluctuations


def vorticity_fluctuations():
    plt.figure(figsize=FIGSIZE)
    plt.subplot(121)
    for i, name in enumerate(case_names):
        utau = cases[name]["utau"]
        delta_nu = cases[name]["delta_nu"]
        plt.semilogx(
            cases[name]["y+"],
            cases[name]["oxox"] * (delta_nu / utau) ** 2,
            lw=1,
            color="C" + str(i + 4),
            label=name,
        )

        plt.plot(
            cases[name]["y+"],
            cases[name]["oxoy"] * (delta_nu / utau) ** 2,
            "--",
            lw=1,
            color="C" + str(i + 4),
        )

    plt.plot(dns_mean[:, 1], dns_omega[:, 2], "k", lw=1, label="DNS")
    plt.plot(dns_mean[:, 1], dns_omega[:, 5], "--k", lw=1)
    plt.xlim(1, 1000)
    plt.ylim(-0.01, 0.05)
    plt.ylabel(r"$\overline {\omega_x'\omega_x'}$, $\overline {\omega_x'\omega_y'}$")
    plt.xlabel(r"$y^+$")
    plt.legend()

    plt.subplot(122)
    for i, name in enumerate(case_names):
        utau = cases[name]["utau"]
        delta_nu = cases[name]["delta_nu"]
        plt.semilogx(
            cases[name]["y+"],
            cases[name]["oyoy"] * (delta_nu / utau) ** 2,
            ":",
            lw=1,
            color="C" + str(i + 4),
            label=name,
        )
        plt.semilogx(
            cases[name]["y+"],
            cases[name]["ozoz"] * (delta_nu / utau) ** 2,
            "-.",
            lw=1,
            color="C" + str(i + 4),
        )

    plt.plot(dns_mean[:, 1], dns_omega[:, 3], ":k", lw=1, label="DNS")
    plt.plot(dns_mean[:, 1], dns_omega[:, 4], "-.k", lw=1)
    plt.xlim(1, 1000)
    plt.ylim(0, 0.2)
    plt.ylabel(r"$\overline {\omega_y'\omega_y'}$, $\overline {\omega_z'\omega_z'}$")
    plt.xlabel(r"$y^+$")
    plt.tight_layout()
    plt.savefig(join(SAVE_PATH, "vorticity.pdf"), bbox_inches="tight", pad_inches=0.02)


vorticity_fluctuations()

# %% k budget


def k_budget():

    plt.figure(figsize=FIGSIZE)

    plt.plot(
        cases["M3"]["y+"][1:-1],
        cases["M3"]["prod"] / cases["M3"]["utau"] ** 4 * nu,
        color="C5",
        lw=1,
        label="Production",
    )
    plt.plot(
        cases["M3"]["y+"][1:-1],
        cases["M3"]["dissipation"] / cases["M3"]["utau"] ** 4 * nu,
        "--",
        color="C5",
        lw=1,
        label="Dissipation",
    )
    plt.plot(
        cases["M3"]["y+"][1:-1],
        cases["M3"]["diffusion"] / cases["M3"]["utau"] ** 4 * nu,
        "-.",
        color="C5",
        lw=1,
        label="Viscous transport",
    )
    plt.plot(
        cases["M3"]["y+"][1:-1],
        cases["M3"]["transport"] / cases["M3"]["utau"] ** 4 * nu,
        lw=1,
        color="C5",
        linestyle=":",
        label="Turbulent transport",
    )
    plt.plot(
        cases["M3"]["y+"][1:-1],
        cases["M3"]["p_diffusion"] / cases["M3"]["utau"] ** 4 * nu,
        lw=1,
        color="C5",
        dashes=[8, 4, 2, 4, 2, 4],
        label="Pressure transport",
    )

    plt.plot(dns_mean[:, 1], dns_budget[:, 2], "k", lw=1)
    plt.plot(dns_mean[:, 1], -dns_budget[:, -2], "--k", lw=1)
    plt.plot(dns_mean[:, 1], dns_budget[:, 4], "-.k", lw=1)
    plt.plot(dns_mean[:, 1], dns_budget[:, 3], ":k", lw=1)

    plt.plot(dns_mean[:, 1], dns_budget[:, 6], "k", lw=1, dashes=[8, 4, 2, 4, 2, 4])
    plt.xlim(1, 50)
    plt.ylim(-0.25, 0.25)
    plt.xlabel(r"$y^+$")
    plt.ylabel(r"$k$ budget terms")
    plt.legend()
    plt.tight_layout()
    plt.savefig(join(SAVE_PATH, "budget.pdf"), bbox_inches="tight", pad_inches=0.02)


k_budget()


# %% Log -log spectra


def log_log_spectra():
    nz = 1000
    dz = 6 / nz
    freq = np.fft.rfftfreq(nz, dz)
    kz = freq * 2 * np.pi

    plt.figure(figsize=(6, 2.5))

    plt.subplot(131)
    for yi, yp in enumerate([12, 200, 850]):
        ind = np.argmin(np.abs(cases["M3"]["y+"] - yp))

        val = cases["M3"]["y+"][ind]

        dns_ind = np.argmin(np.abs(val - dns_yplus))

        # Additional 2pi factor so we divide by delta kz, and not delta f.
        plt.loglog(
            kz,
            cases["M3"]["z_spectrum_u"][ind, :]
            / (2 * np.pi)
            / cases["M3"]["utau"] ** 2,
            label=r"$y^+ =$" + " " + str(yp),
            lw=0.5,
            color=f"C{yi + 4}",
        )

        plt.semilogx(dns_kz, dns_z_spectrum_u[:, dns_ind] / dns_utau**2, "k", lw=0.5)

        plt.semilogx(dns_kz, 0.1 * dns_kz ** (-5 / 3), ":k", lw=1)

        plt.xlabel(r"$k_z$")
        plt.ylabel(r"$E^+_{uu}$")
        plt.ylim(1e-12, 5e-0)

    plt.legend()

    plt.subplot(132)
    for yi, yp in enumerate([12, 200, 850]):
        ind = np.argmin(np.abs(cases["M3"]["y+"] - yp))

        val = cases["M3"]["y+"][ind]

        dns_ind = np.argmin(np.abs(val - dns_yplus))

        plt.loglog(
            kz,
            cases["M3"]["z_spectrum_v"][ind, :]
            / (2 * np.pi)
            / cases["M3"]["utau"] ** 2,
            label=r"$y^+ =$" + " " + str(yp),
            lw=0.5,
            color=f"C{yi + 4}",
        )

        plt.semilogx(dns_kz, dns_z_spectrum_v[:, dns_ind] / dns_utau**2, "k", lw=0.5)

        plt.xlabel(r"$k_z$")
        plt.ylabel(r"$E^+_{vv}$")
        plt.ylim(1e-12, 1e-0)

    plt.subplot(133)
    for yi, yp in enumerate([12, 200, 850]):
        ind = np.argmin(np.abs(cases["M3"]["y+"] - yp))

        val = cases["M3"]["y+"][ind]

        dns_ind = np.argmin(np.abs(val - dns_yplus))

        plt.loglog(
            kz,
            cases["M3"]["z_spectrum_w"][ind, :]
            / (2 * np.pi)
            / cases["M3"]["utau"] ** 2,
            label=r"$y^+ =$" + " " + str(yp),
            lw=0.5,
            color=f"C{yi + 4}",
        )

        plt.semilogx(dns_kz, dns_z_spectrum_w[:, dns_ind] / dns_utau**2, "k", lw=0.5)

        plt.xlabel(r"$k_z$")
        plt.ylabel(r"$E^+_{ww}$")
        plt.ylim(1e-12, 1e-0)

    plt.tight_layout()
    plt.savefig(
        join(SAVE_PATH, "spectra_loglog.pdf"),
        bbox_inches="tight",
        pad_inches=0.02,
    )


log_log_spectra()

# %% Premultiplied spectra 2D and coherence


def premultiplied_spectra_and_coherence():
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    nz = 1000
    dz = 6 / nz
    freq = np.fft.rfftfreq(nz, dz)

    fig, ax = plt.subplots(ncols=2, figsize=FIGSIZE)
    uxt = cases["M3"]["z_spectrum_u"].T
    lzplus = 1 / freq[1:] * dns_utau / nu

    # since we premultiply, here the 2pi factor cancels out.
    p = ax[0].pcolormesh(
        cases["M3"]["y+"][:116],
        lzplus,
        freq[1:, np.newaxis] * uxt[1:, :] / dns_utau**2,
        edgecolor="face",
        vmin=0,
        vmax=4,
    )

    ax[0].contour(
        cases["M3"]["y+"][:116],
        lzplus,
        cases["M3"]["coherence"].T,
        colors="White",
        linewidths=0.5,
        linestyles="dashed",
        levels=[0.2, 0.4, 0.6, 0.95],
    )

    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlim(1, dns_utau / nu)
    ax[0].set_ylabel(r"$\lambda_z^+$")
    ax[0].set_xlabel(r"$y^+$")
    ax[0].set_ylim(25, 5e3)
    ax[0].grid()
    ax[0].set_title("LES")

    dns_fz = dns_kz / (2 * np.pi)
    dns_lambdaz = 1 / dns_fz[1:]

    p = ax[1].pcolormesh(
        dns_yplus,
        dns_lambdaz * dns_utau / nu,
        dns_kz[1:, np.newaxis] * dns_z_spectrum_u[1:, :] / dns_utau**2,
        edgecolor="face",
        vmin=0,
        vmax=4,
    )

    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlim(1, dns_utau / nu)
    ax[1].set_ylabel(r"$\lambda_z^+$")
    ax[1].set_xlabel(r"$y^+$")
    ax[1].set_ylim(25, 5e3)
    ax[1].set_title("DNS")
    ax[1].grid()
    plt.tight_layout()
    plt.savefig(
        join(SAVE_PATH, "spectra_2d.png"),
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.02,
    )


premultiplied_spectra_and_coherence()
# %% Spatial two-point correaltions in z


def two_point_correlations():
    styles = ["-", "--", "-"]
    nz = [250, 500, 1000]

    dns_dz = dns_lz / dns_nz
    dns_z = np.arange(0, dns_nz) * dns_dz

    fig = plt.figure(figsize=(6, 2.5))
    ax1 = fig.add_subplot(131)
    ax1.set_xlim(0, 3)
    ax2 = fig.add_subplot(132)
    ax2.set_xlim(0, 3)
    ax3 = fig.add_subplot(133)
    ax3.set_xlim(0, 3)

    for i, name in enumerate(["M1", "M2", "M3"]):

        if i < 2:
            continue

        dz = 6 / nz[i]
        z = np.linspace(0, 3, int(nz[i] / 2))

        for yi, yp in enumerate([15, 200, 850]):

            dns_acf = np.fft.irfft(dns_z_spectrum_u.T, axis=1, n=dns_nz)
            dns_acf = dns_acf / dns_acf[:, 0][:, None]

            ind = np.argmin(np.abs(cases[name]["y+"] - yp))

            val = cases[name]["y+"][ind]

            dns_ind = np.argmin(np.abs(val - dns_yplus))

            ax1.plot(
                z,
                cases[name]["z_tpcorr_u"][ind, :],
                ls=styles[i],
                label=r"$y^+ =$" + " " + str(yp),
                lw=1,
                color=f"C{yi + 4}",
            )

            ax1.set_xlabel(r"$z / \delta$")
            ax1.set_ylabel(r"$R_{uu}$")

        if i == 0:
            ax1.legend()

        for yi, yp in enumerate([15, 200, 850]):

            dns_acf = np.fft.irfft(dns_z_spectrum_v.T, axis=1, n=dns_nz).real
            dns_acf = dns_acf / dns_acf[:, 0][:, None]

            ind = np.argmin(np.abs(cases[name]["y+"] - yp))

            val = cases[name]["y+"][ind]

            dns_ind = np.argmin(np.abs(val - dns_yplus))

            ax2.plot(
                z,
                cases[name]["z_tpcorr_v"][ind, :],
                ls=styles[i],
                label=r"$y^+ =$" + " " + str(yp),
                lw=1,
                color=f"C{yi + 4}",
            )

            ax2.set_xlabel(r"$z / \delta$")
            ax2.set_ylabel(r"$R_{vv}$")
            ax2.legend()

        for yi, yp in enumerate([15, 200, 850]):

            dns_acf = np.fft.irfft(dns_z_spectrum_w.T, axis=1, n=dns_nz).real
            dns_acf = dns_acf / dns_acf[:, 0][:, None]

            ind = np.argmin(np.abs(cases[name]["y+"] - yp))

            val = cases[name]["y+"][ind]
            print(cases[name]["y"][ind])

            dns_ind = np.argmin(np.abs(val - dns_yplus))

            ax3.plot(
                z,
                cases[name]["z_tpcorr_w"][ind, :],
                ls=styles[i],
                label=r"$y^+ =$" + " " + str(yp),
                lw=1,
                color=f"C{yi + 4}",
            )

            ax3.set_xlabel(r"$z / \delta$")
            ax3.set_ylabel(r"$R_{ww}$")

    plt.tight_layout()

    plt.savefig(
        join(SAVE_PATH, "two_point_corr.pdf"),
        bbox_inches="tight",
        pad_inches=0.02,
    )


# %% Re stress invariant map (Lumely tirangle)


def compute_invariants(tensor):
    # Ensure the tensor is a NumPy array
    tensor = np.array(tensor)

    # Compute trace of T and T^2
    trace_T = np.trace(tensor)
    trace_T2 = np.trace(np.matmul(tensor, tensor))

    # Compute the second invariant
    I2 = 0.5 * (trace_T**2 - trace_T2)

    # Compute the third invariant (determinant of T)
    I3 = np.linalg.det(tensor)

    return I2, I3


def lumely_triangle():
    second = np.zeros(cases["M3"]["y"].size // 2)
    third = np.zeros(cases["M3"]["y"].size // 2)

    second_dns = np.zeros(dns_yplus.size)
    third_dns = np.zeros(dns_yplus.size)

    for i in range(cases["M3"]["y"].size // 2):
        if i == 0:
            continue
        re_tensor = [
            [cases["M3"]["uu"][i], cases["M3"]["uv"][i], 0],
            [cases["M3"]["uv"][i], cases["M3"]["vv"][i], 0],
            [0, 0, cases["M3"]["ww"][i]],
        ]

        k = 0.5 * (cases["M3"]["uu"][i] + cases["M3"]["vv"][i] + cases["M3"]["ww"][i])

        re_tensor = np.array(re_tensor)

        # normalize
        re_tensor = re_tensor / (2 * k)

        # Extract deviatoric part
        re_tensor -= 1 / 3 * np.eye(re_tensor.shape[0]) * np.trace(re_tensor)
        second[i], third[i] = compute_invariants(re_tensor)

    for i in range(dns_fluct.shape[0]):
        re_tensor = [
            [dns_fluct[i, 2], dns_fluct[i, 5], 0],
            [dns_fluct[i, 5], dns_fluct[i, 3], 0],
            [0, 0, dns_fluct[i, 4]],
        ]
        re_tensor = np.array(re_tensor)

        k = 0.5 * (dns_fluct[i, 2] + dns_fluct[i, 3] + dns_fluct[i, 4])

        # normalize
        re_tensor = re_tensor / (2 * k)

        # Extract deviatoric part
        re_tensor -= 1 / 3 * np.eye(re_tensor.shape[0]) * np.trace(re_tensor)

        second_dns[i], third_dns[i] = compute_invariants(re_tensor)

    eta = (-1 / 3 * second[1:]) ** 0.5
    xi = (0.5 * third[1:]) ** (1 / 3)

    eta_dns = (-1 / 3 * second_dns[1:]) ** 0.5
    xi_dns = (0.5 * third_dns[1:]) ** (1 / 3)

    xi_tri = np.linspace(-1 / 6, 1 / 3)
    third_tri = np.linspace(-1 / 36, 2 / 9)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    ax.scatter(
        xi_dns[::5],
        eta_dns[::5],
        alpha=1,
        color="r",
        facecolors="none",
        s=25,
        label="DNS",
        marker="s",
    )

    p = ax.scatter(
        xi,
        eta,
        c=cases["M3"]["y"][1 : second.size],
        cmap="viridis",
        s=5,
        label="M3",
    )

    ax.plot(xi_tri, xi_tri, "--k")
    ax.plot(xi_tri, -xi_tri, "--k")
    ax.plot(xi_tri, np.sqrt(1 / 27 + 2 * xi_tri**3), "--k")

    ax.annotate(
        "one-component",  # Text to display
        xy=(1 / 3, 1 / 3),  # Point to annotate
        xytext=(0.3, 0.12),  # Position of the text
        arrowprops=dict(
            facecolor="black",  # Color of the arrow
            shrink=0.05,  # Shrinking factor
            width=0.5,  # Width of the arrow
            headwidth=4,  # Width of the arrow head
        ),
        horizontalalignment="right",  # Horizontal alignment of the text
        verticalalignment="top",
        fontsize=10,
    )

    ax.annotate(
        "isotropic",
        xy=(0, 0),
        xytext=(0.2, 0.05),
        arrowprops=dict(facecolor="black", shrink=0.05, width=0.5, headwidth=4),
        horizontalalignment="right",
        verticalalignment="top",
        fontsize=10,
    )

    ax.text(
        1 / 8 + 0.015,
        1 / 8 - 0.015,
        "axisymmetric",
        rotation=35,
        rotation_mode="anchor",
        fontsize=10,
    )

    ax.text(
        0,
        0.2,
        "two-component",
        fontsize=10,
        color="black",
    )

    plt.xticks([-1 / 6, 0, 1 / 6, 1 / 3], ["-1/6", "0", "1/6", "1/3"])
    plt.yticks([-1 / 6, 0, 1 / 6, 1 / 3], ["-1/6", "0", "1/6", "1/3"])
    plt.xlim(-1 / 6, 1 / 3)
    plt.ylim(0, 1 / 3)
    plt.xlabel(r"$\xi$")
    plt.ylabel(r"$\eta$")
    plt.tight_layout()
    plt.legend()

    plt.savefig(join(SAVE_PATH, "triangle.pdf"), bbox_inches="tight", pad_inches=0.02)


lumely_triangle()
