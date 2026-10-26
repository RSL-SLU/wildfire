"""Outline the 14 GFED basis regions ("biomes") as polygons for the globe.

    python pipeline/make_regions.py

Traces the basis_regions grid in a GFED5.1 ecosystem file into one (multi)polygon per region
and writes docs/data/regions.geojson.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import xarray as xr
from rasterio.features import shapes
from rasterio.transform import from_origin
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--gfed", default=r"F:\Ctree\GFED5.1", type=Path)
p.add_argument("--out", default=ROOT / "docs" / "data" / "regions.geojson", type=Path)
args = p.parse_args()

CODES = ["BONA", "TENA", "CEAM", "NHSA", "SHSA", "EURO", "MIDE", "NHAF", "SHAF", "BOAS", "CEAS",
         "SEAS", "EQAS", "AUST"]

with xr.open_dataset(args.gfed / "Ecosystem" / "GFED5.1_ecosystem_2022.nc") as d:
    grid = np.flipud(d["basis_regions"].values.astype("int16"))      # north-up
transform = from_origin(-180, 90, 0.25, 0.25)

parts = {i: [] for i in range(1, 15)}
for geom, val in shapes(grid, mask=grid > 0, transform=transform, connectivity=4):
    parts[int(val)].append(shape(geom))

features = []
for i, code in enumerate(CODES, start=1):
    g = unary_union(parts[i]).simplify(0.2, preserve_topology=True)
    g = g.buffer(0)
    # drop specks smaller than ~4 grid cells; they only add triangles to the globe
    if g.geom_type == "MultiPolygon":
        g = unary_union([q for q in g.geoms if q.area > 0.25])
    features.append({"type": "Feature", "properties": {"code": code},
                     "geometry": json.loads(json.dumps(mapping(g)))})
    print(code, g.geom_type, f"{len(getattr(g, 'geoms', [g]))} parts")

fc = {"type": "FeatureCollection", "features": features}
txt = json.dumps(fc, separators=(",", ":"))
txt = __import__("re").sub(r"(\d+\.\d{3})\d+", r"\1", txt)                # 3 decimals is plenty
args.out.write_text(txt, encoding="utf8")
print(f"{len(txt) / 1e6:.2f} MB -> {args.out}")
