"""0.25 deg land-cover images for the map's "Land cover" layer.

    python pipeline/make_landcover_png.py

For every 0.25 deg cell: which land-cover group burned the most carbon, and how much fire carbon
there was in total. One lossless PNG per year plus the 2002-2025 mean, north-up, 1440 x 720:
    R = group index + 1 (0 = no fire), in GROUPS order
    G = amount, 0-255 on a log scale between lo and hi (docs/data/lc/scale.json)
The page decodes the pixels and draws them itself, so nothing is resampled or blurred.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--grids", default=ROOT / "pipeline" / "_grids", type=Path)
p.add_argument("--out", default=ROOT / "docs" / "data" / "lc", type=Path)
args = p.parse_args()
args.out.mkdir(parents=True, exist_ok=True)

YEARS = list(range(2002, 2026))
GROUPS = ["forest", "deforestation", "other", "savanna", "cropland", "peat"]

grids = {y: np.load(args.grids / f"C_groups_{y}.npy").astype("float64") for y in YEARS}
mean = sum(grids.values()) / len(YEARS)
pos = mean.sum(axis=0)
pos = pos[pos > 0]
hi = float(np.percentile(pos, 99.5))
lo = hi / 1000


def encode(g, name):
    total = g.sum(axis=0)
    dom = g.argmax(axis=0)
    show = total >= lo
    r = np.where(show, dom + 1, 0).astype("uint8")
    t = np.clip((np.log(np.maximum(total, lo)) - np.log(lo)) / (np.log(hi) - np.log(lo)), 0, 1)
    gch = np.where(show, np.round(t * 254) + 1, 0).astype("uint8")
    rgb = np.dstack([r, gch, np.zeros_like(r)])
    Image.fromarray(np.flipud(rgb), "RGB").save(args.out / f"{name}.png", optimize=True)
    return int(show.sum())


n = encode(mean, "lc_mean")
for y in YEARS:
    encode(grids[y], f"lc_{y}")
json.dump({"groups": GROUPS, "lo": lo, "hi": hi, "unit": "Tg C per 0.25 deg cell per year",
           "years": YEARS, "width": 1440, "height": 720, "north_up": True},
          open(args.out / "scale.json", "w"), indent=1)
size = sum(f.stat().st_size for f in args.out.glob("*.png"))
print(f"mean layer: {n:,} cells with fire; scale {lo:.2e}-{hi:.3f} Tg C; {size / 1e6:.2f} MB for {len(YEARS) + 1} images")
