"""Render the gridded maps for the dashboard from the grids cached by build_data.py.

    python pipeline/make_maps.py

Writes transparent PNGs to docs/maps/ plus docs/data/maps.json (colour scales, so the page can
draw the legends in its own theme). Run build_data.py first.
"""
import argparse
import json
from pathlib import Path

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm, LinearSegmentedColormap, to_hex

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--grids", default=ROOT / "pipeline" / "_grids", type=Path)
p.add_argument("--countries", default=r"F:\Ctree\Firememo_comp\data\boundaries\ne_50m_admin_0_countries.zip", type=Path)
p.add_argument("--out", default=ROOT / "docs", type=Path)
args = p.parse_args()

YEARS = range(2002, 2026)
EARLY, LATE = range(2002, 2012), range(2016, 2026)      # two 10-year windows
GROUPS = ["forest", "deforestation", "other", "savanna", "cropland", "peat"]
KT_PER_TG = 1e3
OUTLINE = (0.537, 0.529, 0.506, 0.55)                    # mid grey, reads on light and dark
LAND = (0.537, 0.529, 0.506, 0.16)
PROJ, PC = ccrs.Robinson(), ccrs.PlateCarree()
W_PX = 1800

(args.out / "maps").mkdir(parents=True, exist_ok=True)
lon = np.arange(-180, 180, 0.25) + 0.125
lat = np.arange(-90, 90, 0.25) + 0.125                   # south-up rows, like the grids
world = gpd.read_file(f"zip://{args.countries}")
world = world[world.ADMIN != "Antarctica"].geometry

grids = {y: np.load(args.grids / f"C_groups_{y}.npy").astype("float64") for y in YEARS}
stack = np.stack([grids[y] for y in YEARS])              # (year, group, lat, lon), Tg C
mean_by_group = stack.mean(axis=0) * KT_PER_TG           # kt C per cell per year
layers = {"all": mean_by_group.sum(axis=0), **{g: mean_by_group[i] for i, g in enumerate(GROUPS)}}

forest_clean = stack[:, 0]
forest_all = stack[:, 0] + stack[:, 1]
total = stack.sum(axis=1)


def window(a, yrs):
    return a[[YEARS.index(y) for y in yrs]].mean(axis=0) * KT_PER_TG


changes = {"all": window(total, LATE) - window(total, EARLY),
           "forest_overall": window(forest_all, LATE) - window(forest_all, EARLY),
           "forest_clean": window(forest_clean, LATE) - window(forest_clean, EARLY)}


def canvas():
    fig = plt.figure(figsize=(W_PX / 200, W_PX / 200 * 0.5), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1], projection=PROJ)
    ax.set_global()
    ax.set_extent([-180, 180, -58, 84], crs=PC)
    ax.spines["geo"].set_visible(False)
    ax.patch.set_alpha(0)
    fig.patch.set_alpha(0)
    ax.add_geometries(world, crs=PC, facecolor=LAND, edgecolor="none", zorder=1)
    return fig, ax


def finish(fig, ax, name):
    ax.add_geometries(world, crs=PC, facecolor="none", edgecolor=OUTLINE, linewidth=0.25, zorder=3)
    fig.savefig(args.out / "maps" / f"{name}.png", transparent=True)
    plt.close(fig)


def stops(cmap, n=9):
    return [to_hex(cmap(i / (n - 1))) for i in range(n)]


scales = {"mean": {}, "change": {}}
seq = plt.get_cmap("inferno_r")
pool = layers["all"][layers["all"] > 0]
vmin = 10 ** np.floor(np.log10(np.percentile(pool, 5)))
vmax = 10 ** np.ceil(np.log10(pool.max()))
for name, arr in layers.items():
    fig, ax = canvas()
    a = np.where(arr > vmin * 0.2, arr, np.nan)
    ax.pcolormesh(lon, lat, a, transform=PC, cmap=seq, norm=LogNorm(vmin, vmax),
                  shading="nearest", rasterized=True, zorder=2)
    finish(fig, ax, f"mean_{name}")
    scales["mean"][name] = {"total_TgC": float(arr.sum() / KT_PER_TG)}
    print("mean", name, f"{arr.sum() / KT_PER_TG:.0f} Tg C/yr")
scales["mean_scale"] = {"type": "log", "vmin": float(vmin), "vmax": float(vmax), "stops": stops(seq),
                        "unit": "kt C per 0.25° cell per year"}

div = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#e9e8e3", "#ef9a8e", "#b3261e"])
lim = float(np.nanpercentile(np.abs(changes["all"][changes["all"] != 0]), 99))
lim = float(10 ** np.ceil(np.log10(lim)))
linthresh = lim / 100
for name, arr in changes.items():
    fig, ax = canvas()
    a = np.where(np.abs(arr) > linthresh / 4, arr, np.nan)       # near-zero stays see-through
    ax.pcolormesh(lon, lat, a, transform=PC, cmap=div,
                  norm=SymLogNorm(linthresh=linthresh, vmin=-lim, vmax=lim, base=10),
                  shading="nearest", rasterized=True, zorder=2)
    finish(fig, ax, f"change_{name}")
    scales["change"][name] = {"net_TgC": float(arr.sum() / KT_PER_TG)}
    print("change", name, f"{arr.sum() / KT_PER_TG:+.0f} Tg C/yr")
scales["change_scale"] = {"type": "symlog", "lim": lim, "linthresh": linthresh, "stops": stops(div),
                          "unit": "kt C per 0.25° cell per year",
                          "early": [EARLY[0], EARLY[-1]], "late": [LATE[0], LATE[-1]]}
scales["period"] = [YEARS[0], YEARS[-1]]
(args.out / "data" / "maps.json").write_text(json.dumps(scales, indent=1), encoding="utf8")
for f in sorted((args.out / "maps").glob("*.png")):
    print(f"{f.stat().st_size / 1e6:5.2f} MB  {f.name}")
