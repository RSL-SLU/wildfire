"""Gridded map layers for the dashboard: yearly 0.25 deg fire carbon and 1 deg trends.

    python pipeline/make_grid_layers.py

Reads the yearly group grids cached by build_data.py (Tg C per 0.25 deg cell, south-up).

docs/data/grid/   two lossless RGB PNGs per year (and for the 2002-2025 mean), 1440 x 720, north-up
    {key}_a.png   R = all fires, G = forest, B = deforestation
    {key}_b.png   R = savanna & grassland, G = cropland, B = peat
    Each byte is 0 for no fire, else 1-255 on one log scale between lo and hi (scale.json).
    Tundra & other = all - the five others, so it needs no channel of its own.

docs/data/trend/  one RGB PNG per fire type, 360 x 180 (1 deg), north-up, 2002-2025
    R = Theil-Sen trend, % of the cell's mean per year: 1..255 maps -CAP..+CAP (0 = no data)
    G = 255 where the Mann-Kendall test gives p < 0.05
    B = the cell's mean fire carbon on the same log scale (for the tooltip)
"""
import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--grids", default=ROOT / "pipeline" / "_grids", type=Path)
p.add_argument("--out", default=ROOT / "docs" / "data", type=Path)
args = p.parse_args()
(args.out / "grid").mkdir(parents=True, exist_ok=True)
(args.out / "trend").mkdir(parents=True, exist_ok=True)

YEARS = list(range(2002, 2026))
# group order in the cached grids
FOREST, DEFOR, OTHER, SAVANNA, CROP, PEAT = range(6)
TYPES = {"all": [0, 1, 2, 3, 4, 5], "forest_overall": [FOREST, DEFOR], "forest_clean": [FOREST],
         "deforestation": [DEFOR], "savanna": [SAVANNA], "cropland": [CROP], "peat": [PEAT], "other": [OTHER]}
LO, HI = 1e-5, 8.0                       # Tg C per 0.25 deg cell per year (10 t .. 8 Mt)
CAP = 15.0                               # trend colour range, % per year
MIN_MEAN_1DEG = 0.005                    # Tg C/yr; below this a 1 deg trend is noise

g = np.stack([np.flipud(np.load(args.grids / f"C_groups_{y}.npy").astype("float64").transpose(1, 2, 0)).transpose(2, 0, 1)
              for y in YEARS])                                   # (year, group, 720, 1440), north-up
print("grids", g.shape, f"max cell {g.sum(axis=1).max():.2f} Tg C")


def q(v):
    """Tg C -> byte on the log scale."""
    t = (np.log(np.clip(v, LO, HI)) - np.log(LO)) / (np.log(HI) - np.log(LO))
    return np.where(v >= LO, np.round(t * 254) + 1, 0).astype("uint8")


def save_rgb(chans, path):
    Image.fromarray(np.dstack(chans), "RGB").save(path, optimize=True)


def write_year(arr, key):                                          # arr: (group, 720, 1440)
    save_rgb([q(arr.sum(axis=0)), q(arr[FOREST]), q(arr[DEFOR])], args.out / "grid" / f"{key}_a.png")
    save_rgb([q(arr[SAVANNA]), q(arr[CROP]), q(arr[PEAT])], args.out / "grid" / f"{key}_b.png")


write_year(g.mean(axis=0), "mean")
for i, y in enumerate(YEARS):
    write_year(g[i], str(y))

# ---- 1 deg trends ---------------------------------------------------------------------------
g1 = g.reshape(len(YEARS), 6, 180, 4, 360, 4).sum(axis=(3, 5))  # (year, group, 180, 360)
yrs = np.array(YEARS, dtype=float)
ii, jj = np.triu_indices(len(YEARS), k=1)
n = len(YEARS)
var_s = n * (n - 1) * (2 * n + 5) / 18
counts = {}
for name, grp in TYPES.items():
    s = g1[:, grp].sum(axis=1).reshape(n, -1)                     # (year, cells)
    m = s.mean(axis=0)
    ok = m >= MIN_MEAN_1DEG
    d = s[jj][:, ok] - s[ii][:, ok]                                # pairwise differences, later minus earlier
    slope = np.median(d / (yrs[jj] - yrs[ii])[:, None], axis=0)
    S = np.sign(d).sum(axis=0)
    z = np.where(S == 0, 0, (S - np.sign(S)) / np.sqrt(var_s))
    pval = 2 * (1 - norm.cdf(np.abs(z)))
    pct = np.full(m.shape, np.nan); pct[ok] = 100 * slope / m[ok]
    sig = np.zeros(m.shape, bool); sig[ok] = pval < 0.05
    r = np.where(ok, np.round((np.clip(np.nan_to_num(pct), -CAP, CAP) + CAP) / (2 * CAP) * 254) + 1, 0).astype("uint8")
    gch = np.where(sig, 255, 0).astype("uint8")
    b = q(m / 16)                                                  # mean per 0.25 deg cell, same scale as the grids
    save_rgb([r.reshape(180, 360), gch.reshape(180, 360), b.reshape(180, 360)], args.out / "trend" / f"{name}.png")
    counts[name] = (int(ok.sum()), int(sig.sum()))

json.dump({"years": YEARS, "lo": LO, "hi": HI, "unit": "Tg C per 0.25 deg cell per year",
           "a": ["all", "forest_clean", "deforestation"], "b": ["savanna", "cropland", "peat"],
           "trend": {"cap": CAP, "period": [YEARS[0], YEARS[-1]], "min_mean_1deg_TgC": MIN_MEAN_1DEG,
                     "cells": {k: {"with_data": v[0], "significant": v[1]} for k, v in counts.items()}}},
          open(args.out / "grid" / "scale.json", "w"), indent=1)
size = lambda d: sum(f.stat().st_size for f in (args.out / d).glob("*.png")) / 1e6
print(f"grid: {size('grid'):.1f} MB, trend: {size('trend'):.2f} MB")
for k, (a, b) in counts.items():
    print(f"  trend {k:15} {a:6,} cells, {b:5,} significant ({100 * b / max(a, 1):.0f}%)")
