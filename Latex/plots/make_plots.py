#!/usr/bin/env python3
"""Turn the CSV output of `sph_experiments` into the paper's figures.

Usage (from the repository root):
    ./build/sph_experiments Latex/data
    python3 Latex/plots/make_plots.py [data_dir] [figure_dir]

Defaults: data_dir = Latex/data, figure_dir = Latex/figures.

Colour roles (validated with the dataviz palette checks):
  * categorical  - blue / orange / aqua for simulation / reference / analytic,
  * ordinal      - one-hue blue ramp for ordered sweep values,
  * sequential   - blue ramp for pressure,
  * diverging    - blue <-> grey <-> red for the density error (red = compressed).
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data"
FIGURES = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "figures"

# --- palette -----------------------------------------------------------------
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"
WALL = "#b5b4ae"
SURFACE = "#ffffff"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RED = "#e34948"
NEUTRAL = "#f0efec"
ORDINAL = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "seq_blue", ["#9ec5f4", "#6da7ec", "#3987e5", "#1c5cab", "#0d366b"])
DIVERGING = LinearSegmentedColormap.from_list("div_blue_red", [BLUE, NEUTRAL, RED])

TEXT_WIDTH = 6.0  # inches, matches the A4 / 2.8 cm margin layout

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.edgecolor": INK_SECONDARY,
    "axes.labelcolor": INK,
    "xtick.color": INK_SECONDARY,
    "ytick.color": INK_SECONDARY,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 1.6,
    "legend.frameon": False,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "pdf.fonttype": 42,
})


# --- data helpers ---------------------------------------------------------------
def read_csv(name):
    """Read a numeric CSV into a dict of numpy arrays (non-numeric columns kept as str)."""
    path = DATA / name
    with path.open() as f:
        rows = [r for r in csv.reader(line for line in f if not line.startswith("#"))]
    header, body = rows[0], rows[1:]
    columns = {}
    for i, key in enumerate(header):
        values = [r[i] for r in body]
        try:
            columns[key] = np.array([float(v) for v in values])
        except ValueError:
            columns[key] = np.array(values)
    return columns


def read_meta():
    meta = {"dam": {}, "hydro": {}}
    with (DATA / "scene_meta.csv").open() as f:
        for row in csv.DictReader(f):
            meta[row["scene"]][row["key"]] = float(row["value"])
    return meta


def value_label(value):
    return f"{value:g}"


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / name)
    plt.close(fig)
    print(f"wrote {FIGURES / name}")


def draw_particles(ax, snap, values, cmap, norm, size=3.0):
    """Scatter fluid particles coloured by `values`; walls in neutral grey."""
    wall = snap["fixed"] > 0.5
    ax.scatter(snap["x"][wall], snap["y"][wall], s=size, c=WALL, linewidths=0, rasterized=True)
    sc = ax.scatter(snap["x"][~wall], snap["y"][~wall], s=size, c=values[~wall],
                    cmap=cmap, norm=norm, linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    return sc


# --- figures -------------------------------------------------------------------
def fig_hydro_pressure(meta):
    m = meta["hydro"]
    H = m["fluid_height"]
    scale = m["rest_density"] * m["gravity"] * H
    prof = read_csv("hydro_profile.csv")
    depth = prof["depth"] / H

    fig, ax = plt.subplots(figsize=(TEXT_WIDTH * 0.55, 2.8))
    ax.scatter(prof["pressure"] / scale, depth, s=3, c="#c9c8c2", linewidths=0,
               label="instantaneous", rasterized=True)
    ax.scatter(prof["pressure_time_averaged"] / scale, depth, s=4, c=BLUE, linewidths=0,
               label="2 s time average", rasterized=True)
    d = np.linspace(0, 1, 2)
    ax.plot(d, d, color=ORANGE, linestyle="--", label=r"$p = \rho_0 g d$")
    ax.set_xlabel(r"pressure $p \,/\, \rho_0 g H$")
    ax.set_ylabel(r"depth $d \,/\, H$")
    ax.set_xlim(0, 1.3)
    ax.set_ylim(1.02, 0)
    ax.legend(loc="upper right")
    save(fig, "hydro_pressure.pdf")


def fig_hydro_settling():
    ts = read_csv("hydro_timeseries.csv")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(TEXT_WIDTH * 0.55, 3.2), sharex=True)
    ax1.plot(ts["t"], 100 * ts["max_density_error"], color=BLUE, label="maximum")
    ax1.plot(ts["t"], 100 * ts["mean_density_error"], color=ORANGE, linestyle="--", label="mean")
    ax1.set_ylabel(r"$\rho/\rho_0 - 1$ [%]")
    ax1.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    ax2.semilogy(ts["t"], ts["mean_kinetic_energy"], color=BLUE)
    ax2.set_ylabel(r"mean $\frac{1}{2}|\mathbf{v}|^2$ [px$^2$/s$^2$]")
    ax2.set_xlabel("time $t$ [s]")
    ax2.set_xlim(0, ts["t"].max())
    fig.align_ylabels()
    save(fig, "hydro_settling.pdf")


def martin_moyce():
    """Experimental reference, if the user has filled in the template."""
    path = DATA / "martin_moyce_n2.csv"
    if not path.exists():
        return None
    data = read_csv("martin_moyce_n2.csv")
    if len(data.get("T", [])) == 0:
        return None
    return data


def fig_dam_front(meta):
    m = meta["dam"]
    n = np.sqrt(m["fluid_height"] / m["fluid_width"])
    tank_limit = m["tank_width"] / m["fluid_width"]
    front = read_csv("dam_front.csv")
    reach = front["T"] <= front["T"][np.argmax(front["Z_over_a"] >= tank_limit - 0.02)] + 0.05

    fig, ax = plt.subplots(figsize=(TEXT_WIDTH * 0.55, 2.8))
    ax.plot(front["T"][reach], front["Z_over_a"][reach], color=BLUE, label="SPH (this work)")
    T = np.linspace(0, 1.2, 50)
    ax.plot(T, 1 + 2 * n * T, color=AQUA, linestyle=":", linewidth=1.8,
            label=r"shallow water, $Z/a = 1 + 2nT$")
    mm = martin_moyce()
    if mm is not None:
        ax.plot(mm["T"], mm["Z_over_a"], linestyle="none", marker="o", markersize=4.5,
                markerfacecolor=ORANGE, markeredgecolor=SURFACE, markeredgewidth=0.8,
                label=r"Martin \& Moyce (1952), $n^2 = 2$")
    ax.axhline(tank_limit, color=INK_SECONDARY, linewidth=0.8)
    ax.text(0.05, tank_limit + 0.05, "far wall", color=INK_SECONDARY, fontsize=8, va="bottom")
    ax.set_xlabel(r"$T = t\sqrt{g/a}$")
    ax.set_ylabel(r"front position $Z/a$")
    ax.set_xlim(0, 2.6)
    ax.set_ylim(0.9, tank_limit + 0.4)
    ax.legend(loc="lower right")
    save(fig, "dam_front.pdf")


def fig_dam_snapshots(meta):
    m = meta["dam"]
    scale = m["rest_density"] * m["gravity"] * m["fluid_height"]
    times = [0, 1, 2, 3]
    fig, axes = plt.subplots(2, 2, figsize=(TEXT_WIDTH, 3.6))
    norm = plt.Normalize(0, 1)
    sc = None
    for ax, T in zip(axes.flat, times):
        snap = read_csv(f"dam_snapshot_T{T}.csv")
        sc = draw_particles(ax, snap, snap["pressure_smoothed"] / scale, SEQUENTIAL, norm, size=2.6)
        ax.set_title(f"$T = {T}$", color=INK)
    cbar = fig.colorbar(sc, ax=axes, shrink=0.8, pad=0.02, extend="max")
    cbar.set_label(r"pressure $p \,/\, \rho_0 g H$")
    cbar.outline.set_visible(False)
    save(fig, "dam_snapshots.pdf")


STIFFNESS_LABELS = {0.0625: "1/16", 0.25: "1/4", 1.0: "1", 4.0: "4"}


def stiffness_label(value):
    return STIFFNESS_LABELS.get(float(value), value_label(value))


def fig_sweep_stiffness():
    sw = read_csv("sweep_k.csv")
    f = sw["value"]  # k / k_ref
    fig, axes = plt.subplots(1, 3, figsize=(TEXT_WIDTH, 2.3))

    ax = axes[0]
    ax.loglog(f, 100 * sw["hydro_max_density_error"], color=BLUE, marker="o", markersize=4,
              label="column, max")
    ax.loglog(f, 100 * sw["dam_mean_max_density_error"], color=ORANGE, marker="s", markersize=4,
              linestyle="--", label="dam break, mean of max")
    ff = np.geomspace(f.min(), f.max(), 50)
    # Hydrostatic compression: rho/rho0 - 1 ~ rho0 g H / (7 k), which is 0.5 % at k_ref.
    ax.loglog(ff, 0.5 / ff, color=AQUA, linestyle=":", linewidth=1.8,
              label=r"$\rho_0 gH/(7k)$")
    ax.set_xlabel(r"$k \,/\, k_\mathrm{ref}$")
    ax.set_ylabel(r"density error [%]")
    ax.set_title("(a) compressibility", loc="left")
    handles, labels = ax.get_legend_handles_labels()

    ax = axes[1]
    ax.plot(f, sw["dam_substeps_per_frame"], color=BLUE, marker="o", markersize=4)
    ax.set_xscale("log")
    ax.set_xlabel(r"$k \,/\, k_\mathrm{ref}$")
    ax.set_ylabel("substeps per frame")
    ax.set_title("(b) cost", loc="left")

    ax = axes[2]
    ax.loglog(f, 100 * sw["hydro_frame_to_frame_change"], color=BLUE, marker="o", markersize=4)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"$k \,/\, k_\mathrm{ref}$")
    ax.set_ylabel(r"$\langle|\Delta p|\rangle$ per frame [%]")
    ax.set_title("(c) pressure noise", loc="left")

    for ax in axes:
        ax.set_xticks(f)
        ax.set_xticklabels([stiffness_label(v) for v in f])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    fig.tight_layout(w_pad=1.2)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3)
    save(fig, "sweep_stiffness.pdf")


def fig_sweep_viscosity():
    sw = read_csv("sweep_visc.csv")
    nu = sw["value"]
    fig, axes = plt.subplots(1, 2, figsize=(TEXT_WIDTH, 2.4))

    ax = axes[0]
    for colour, value in zip(ORDINAL, nu):
        series = read_csv(f"sweep_hydro_visc_{value_label(value)}.csv")
        ax.semilogy(series["t"], series["mean_kinetic_energy"], color=colour,
                    label=rf"$\nu = {value_label(value)}$")
    ax.set_xlabel("time $t$ [s]")
    ax.set_ylabel(r"mean $\frac{1}{2}|\boldsymbol{v}|^2$ [px$^2$/s$^2$]")
    ax.set_xlim(0, series["t"].max())
    ax.set_title("(a) hydrostatic column", loc="left")
    ax.legend(loc="upper left", ncol=2)
    ax.set_ylim(None, ax.get_ylim()[1] * 8)

    ax = axes[1]
    for colour, value in zip(ORDINAL, nu):
        front = read_csv(f"sweep_front_visc_{value_label(value)}.csv")
        ax.plot(front["T"], front["Z_over_a"], color=colour, label=rf"$\nu = {value_label(value)}$")
    ax.set_xlabel(r"$T = t\sqrt{g/a}$")
    ax.set_ylabel(r"front position $Z/a$")
    ax.set_xlim(0, 2.4)
    ax.set_title("(b) dam-break front", loc="left")
    ax.legend(loc="upper left")
    fig.tight_layout(w_pad=1.5)
    save(fig, "sweep_viscosity.pdf")


def fig_boundary():
    """Mirrored vs. zero boundary pressure, and residual wall motion with mirroring."""
    fig = plt.figure(figsize=(TEXT_WIDTH, 3.6))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], wspace=0.3, hspace=0.6)

    variants = [(0, r"$p_{i_b} = 0$", ORANGE), (1, r"mirrored, $p_{i_b} = p_i$", BLUE)]
    for row, (flag, label, colour) in enumerate(variants):
        ax = fig.add_subplot(grid[row, 0])
        layers = read_csv(f"layers_mirror_{flag}.csv")
        ax.bar(layers["height"], 100 * layers["count_fraction"], width=0.9, color=colour)
        ax.set_xlim(0, 60)
        ax.set_title(label, loc="left", fontsize=8, pad=2)
        ax.grid(axis="x", visible=False)
        if row == 0:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(r"height above floor [px] ($\Delta x = 10$ px)")
        ax.set_ylabel("share [%]")
    fig.text(0.02, 0.98, "(a) particle heights in the settled column", fontsize=9, va="top")

    ax = fig.add_subplot(grid[:, 1])
    snap = read_csv("hydro_snapshot.csv")
    drift = read_csv("hydro_drift.csv")
    wall = snap["fixed"] > 0.5
    ax.scatter(snap["x"][wall], snap["y"][wall], s=1.0, c=WALL, linewidths=0, rasterized=True)
    drift_speed = np.hypot(drift["drift_vx"], drift["drift_vy"])
    sc = ax.scatter(drift["x"], drift["y"], s=1.6, c=drift_speed, cmap=SEQUENTIAL,
                    norm=plt.Normalize(0, 12), linewidths=0, rasterized=True)
    # Bin-averaged drift as arrows (30 px cells).
    cell = 30.0
    ix = np.floor(drift["x"] / cell).astype(int)
    iy = np.floor(drift["y"] / cell).astype(int)
    xs, ys, us, vs = [], [], [], []
    for cx in np.unique(ix):
        for cy in np.unique(iy):
            sel = (ix == cx) & (iy == cy)
            if sel.sum() < 4:
                continue
            xs.append(drift["x"][sel].mean())
            ys.append(drift["y"][sel].mean())
            us.append(drift["drift_vx"][sel].mean())
            vs.append(drift["drift_vy"][sel].mean())
    ax.quiver(xs, ys, us, vs, color=INK, angles="xy", scale_units="xy", scale=0.25,
              width=0.006, headwidth=3.5)
    ax.set_aspect("equal")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("(b) residual drift with mirroring", loc="left")
    cbar = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02, extend="max")
    cbar.set_label("mean drift speed [px/s]")
    cbar.outline.set_visible(False)
    save(fig, "boundary.pdf")


def fig_pairing():
    ratios = [1.2, 1.3, 1.6, 2.0]
    sweep = read_csv("sweep_hdx.csv")
    neighbours = dict(zip(sweep["value"], sweep["hydro_mean_neighbors"]))
    fig, axes = plt.subplots(1, 4, figsize=(TEXT_WIDTH, 1.9))
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-0.02, vmax=0.02)
    sc = None
    for ax, ratio in zip(axes, ratios):
        snap = read_csv(f"density_hydro_hdx_{value_label(ratio)}.csv")
        window = (snap["x"] > 60) & (snap["x"] < 240) & (snap["y"] < 150)
        snap = {k: v[window] for k, v in snap.items()}
        sc = draw_particles(ax, snap, snap["density_error"], DIVERGING, norm, size=5)
        n = neighbours.get(ratio, float("nan"))
        ax.set_title(rf"$h/\Delta x = {value_label(ratio)}$" + f"\n({n:.0f} neighbours)", fontsize=8)
    cbar = fig.colorbar(sc, ax=list(axes), shrink=0.8, pad=0.02, extend="both")
    cbar.set_label(r"$\rho/\rho_0 - 1$")
    cbar.ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
    cbar.outline.set_visible(False)
    save(fig, "pairing.pdf")


def fig_density_grid():
    studies = [
        ("k", [0.0625, 0.25, 1, 4], lambda v: rf"$k = {stiffness_label(v)}\,k_\mathrm{{ref}}$"),
        ("visc", [300, 600, 1000, 2000], lambda v: rf"$\nu = {value_label(v)}$"),
        ("hdx", [1.2, 1.3, 1.6, 2.0], lambda v: rf"$h/\Delta x = {value_label(v)}$"),
    ]
    fig, axes = plt.subplots(3, 4, figsize=(TEXT_WIDTH, 2.9))
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-0.03, vmax=0.03)
    sc = None
    for row, (param, values, label) in enumerate(studies):
        for col, value in enumerate(values):
            ax = axes[row, col]
            snap = read_csv(f"density_dam_{param}_{value_label(value)}.csv")
            window = snap["y"] < 230
            snap = {k: v[window] for k, v in snap.items()}
            sc = draw_particles(ax, snap, snap["density_error"], DIVERGING, norm, size=1.2)
            ax.set_title(label(value), fontsize=8)
    cbar = fig.colorbar(sc, ax=axes, shrink=0.75, pad=0.02, extend="both")
    cbar.set_label(r"density error $\rho/\rho_0 - 1$")
    cbar.ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
    cbar.outline.set_visible(False)
    save(fig, "density_grid.pdf")


def main():
    meta = read_meta()
    fig_hydro_pressure(meta)
    fig_hydro_settling()
    fig_dam_front(meta)
    fig_dam_snapshots(meta)
    fig_sweep_stiffness()
    fig_sweep_viscosity()
    fig_pairing()
    fig_boundary()
    fig_density_grid()


if __name__ == "__main__":
    main()
