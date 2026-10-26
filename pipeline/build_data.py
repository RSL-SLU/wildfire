"""Build the dashboard data from GFED5.1 (+ NRT) files.

Reads the raw GFED netCDF files, adds them up per year by land-cover class for the world,
the 14 GFED basis regions and every country, and writes compact JSON for the web page in
docs/data/. Also writes the same numbers as tidy CSV in docs/data/csv/ for download.

    python pipeline/build_data.py                      # paths below
    python pipeline/build_data.py --gfed F:/GFED5.1 --nrt F:/GFED5.1NRT --countries ne_50m.zip

Sources
    1997-2001  GFED5.1 monthly (total carbon only; 1997-2000 on a 1 deg grid)
    2002-2022  GFED5.1 ecosystem files (16 land-cover classes, carbon + burned area)
    2023-2025  GFED5.1 NRT, reprocessed, combined VIIRS (CMB), same 16 classes
"""
import argparse
import json
from datetime import date
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from rasterio.features import rasterize
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]

p = argparse.ArgumentParser()
p.add_argument("--gfed", default=r"F:\Ctree\GFED5.1", type=Path)
p.add_argument("--nrt", default=r"F:\GFED5.1NRT\Reprocessed\Monthly\Ecosystem\CMB", type=Path)
p.add_argument("--countries", default=r"F:\Ctree\Firememo_comp\data\boundaries\ne_50m_admin_0_countries.zip", type=Path)
p.add_argument("--out", default=ROOT / "docs" / "data", type=Path)
p.add_argument("--grids", default=ROOT / "pipeline" / "_grids", type=Path,
               help="where yearly 0.25 deg grids are cached for make_maps.py (not committed)")
args = p.parse_args()

TOTAL_YEARS = range(1997, 2026)
CLASS_YEARS = range(2002, 2026)
GFED51_LAST = 2022
G_PER_TG = 1e12
M2_PER_KM2 = 1e6
SUB = 5                                     # 0.05 deg subcells for border sharing

CLASSES = ["Tundra", "Sparse boreal forest", "Boreal forest", "Temperate grassland",
           "Temperate shrubland", "Temperate mosaic", "Temperate forest", "Tropical grassland",
           "Tropical shrubland", "Open savanna", "Woody savanna", "Tropical forest", "Other",
           "Cropland", "Peat", "Deforestation"]
# Stack order is also the validated colour order (see docs/assets/app.js).
GROUPS = [
    ("forest", "Forest", [2, 3, 7, 12]),
    ("deforestation", "Deforestation", [16]),
    ("other", "Tundra & other", [1, 13]),
    ("savanna", "Savanna, grassland & shrubland", [4, 5, 6, 8, 9, 10, 11]),
    ("cropland", "Cropland (agricultural)", [14]),
    ("peat", "Peat", [15]),
]
REGIONS = [("BONA", "Boreal North America"), ("TENA", "Temperate North America"),
           ("CEAM", "Central America"), ("NHSA", "Northern Hemisphere South America"),
           ("SHSA", "Southern Hemisphere South America"), ("EURO", "Europe"),
           ("MIDE", "Middle East"), ("NHAF", "Northern Hemisphere Africa"),
           ("SHAF", "Southern Hemisphere Africa"), ("BOAS", "Boreal Asia"),
           ("CEAS", "Central Asia"), ("SEAS", "Southeast Asia"), ("EQAS", "Equatorial Asia"),
           ("AUST", "Australia and New Zealand")]

ECO = args.gfed / "Ecosystem"
MON = args.gfed / "Monthly"
args.out.mkdir(parents=True, exist_ok=True)
(args.out / "csv").mkdir(exist_ok=True)
args.grids.mkdir(parents=True, exist_ok=True)

# ---- regions: GFED basis regions, a true partition of the 0.25 deg grid (south-up) --------
with xr.open_dataset(ECO / "GFED5.1_ecosystem_2022.nc") as d:
    REGION_OF = d["basis_regions"].values.astype(int)            # 0 = ocean, 1..14
REGION_FLAT = REGION_OF.ravel()

# ---- countries: share of each 0.25 deg cell held by each country, painted at 0.05 deg -----
countries = gpd.read_file(f"zip://{args.countries}")[["ADMIN", "ADM0_A3", "CONTINENT", "geometry"]]
countries.columns = ["name", "iso3", "continent", "geometry"]
NY, NX = 720 * SUB, 1440 * SUB
ids = rasterize(((g, i + 1) for i, g in enumerate(countries.geometry)), out_shape=(NY, NX),
                transform=from_origin(-180, 90, 0.25 / SUB, 0.25 / SUB), fill=0, dtype="int32")
ids = np.flipud(ids)                                             # south-up, like GFED
cell = np.repeat(np.arange(720), SUB)[:, None] * 1440 + np.repeat(np.arange(1440), SUB)[None, :]
land = ids > 0
uniq, count = np.unique(cell[land].astype("int64") * 1000 + ids[land], return_counts=True)
CELL_IDX, CTRY_IDX, WEIGHT = uniq // 1000, (uniq % 1000).astype(int) - 1, count / SUB ** 2
del ids, cell, land
print(f"{len(countries)} countries, {len(uniq):,} (cell, country) pairs")


def by_region(grid):
    return np.bincount(REGION_FLAT, weights=grid.ravel(), minlength=15)[1:15]


def by_country(grid):
    return np.bincount(CTRY_IDX, weights=grid.ravel()[CELL_IDX] * WEIGHT, minlength=len(countries))


# ---- 1 deg -> 0.25 deg for 1997-2000 (share each 1 deg cell among its land cells) ---------
land_q = (REGION_OF > 0).astype(float)
land_count = land_q.reshape(180, 4, 360, 4).sum(axis=(1, 3))


def to_quarter_degree(g1):
    n = np.where(land_count > 0, land_count, 16)
    share = np.repeat(np.repeat(g1 / n, 4, axis=0), 4, axis=1)
    has_land = np.repeat(np.repeat(land_count > 0, 4, axis=0), 4, axis=1)
    return np.where(has_land, share * land_q, share)


# ---- readers: every one yields south-up 0.25 deg arrays --------------------------------
def read_classes(year):
    """Yield (class index 0-15, month 0-11, carbon Tg, burned area km2), one 2-D grid each.
    Reads one class (GFED5.1) or one month (NRT) at a time to keep memory small."""
    if year <= GFED51_LAST:
        f = ECO / f"GFED5.1_ecosystem_{year}.nc"
        with xr.open_dataset(f, group="carbon_emissions_partitioning") as c, \
             xr.open_dataset(f, group="burned_area_partitioning") as b:
            assert [n.split("_", 2)[2] for n in c.data_vars] == CLASSES
            for k, (nc_, nb_) in enumerate(zip(c.data_vars, b.data_vars)):
                C = np.nan_to_num(c[nc_].values.astype("float64")) / G_PER_TG
                BA = np.nan_to_num(b[nb_].values.astype("float64")) / M2_PER_KM2
                for m in range(12):
                    yield k, m, C[m], BA[m]
    else:
        files = sorted((args.nrt / str(year)).glob("*.nc"))
        assert len(files) == 12, f"{year}: expected 12 NRT files, found {len(files)}"
        for m, f in enumerate(files):
            with xr.open_dataset(f) as d:
                assert list(d.lct.values) == CLASSES, f"class order changed in {f.name}"
                C = np.nan_to_num(d["EM"].values[0].astype("float64")) / G_PER_TG
                BA = np.nan_to_num(d["BA"].values[0].astype("float64")) / M2_PER_KM2
            for k in range(16):
                yield k, m, C[k], BA[k]


def read_total(year):
    with xr.open_dataset(MON / f"GFED5.1_monthly_{year}.nc") as d:
        g = np.nan_to_num(d["C"].sum("time").values.astype("float64"))
    if g.shape == (180, 360):
        g = to_quarter_degree(g)
    return g / G_PER_TG


# ---- the main loop -----------------------------------------------------------------------
NC, NR, NK = 16, 14, len(countries)
YC = list(CLASS_YEARS)
C_world = np.zeros((NC, len(YC))); BA_world = np.zeros((NC, len(YC)))
C_reg = np.zeros((NR, NC, len(YC))); BA_reg = np.zeros((NR, NC, len(YC)))
C_ctry = np.zeros((NK, NC, len(YC))); BA_ctry = np.zeros((NK, NC, len(YC)))
# monthly carbon by class, summed over years (for the seasonal cycle)
M_world = np.zeros((NC, 12)); M_reg = np.zeros((NR, NC, 12)); M_ctry = np.zeros((NK, NC, 12))
T_world, T_reg, T_ctry = {}, {}, {}
source = {}

for year in TOTAL_YEARS:
    if year in CLASS_YEARS:
        j = YC.index(year)
        Cy = np.zeros((NC, 720, 1440))                           # annual carbon per class
        for k, m, C, BA in read_classes(year):
            Cy[k] += C
            cr, cc = by_region(C), by_country(C)
            C_world[k, j] += C.sum(); BA_world[k, j] += BA.sum()
            C_reg[:, k, j] += cr; BA_reg[:, k, j] += by_region(BA)
            C_ctry[:, k, j] += cc; BA_ctry[:, k, j] += by_country(BA)
            M_world[k, m] += C.sum(); M_reg[:, k, m] += cr; M_ctry[:, k, m] += cc
        np.save(args.grids / f"C_groups_{year}.npy",                # (6, 720, 1440), GROUPS order
                np.stack([Cy[[c - 1 for c in cls]].sum(axis=0) for _, _, cls in GROUPS]).astype("float32"))
        total = Cy.sum(axis=0)
        src = "GFED5.1" if year <= GFED51_LAST else "GFED5.1 NRT (CMB)"
    else:
        total, src = read_total(year), "GFED5.1 (total only)"
    np.save(args.grids / f"C_total_{year}.npy", total.astype("float32"))
    T_world[year], T_reg[year], T_ctry[year] = total.sum(), by_region(total), by_country(total)
    source[year] = src
    print(f"{year}  {src:22}  {total.sum():7.0f} Tg C")

# ---- checks -------------------------------------------------------------------------------
in_ctry = np.array([T_ctry[y].sum() / T_world[y] for y in TOTAL_YEARS])
in_reg = np.array([T_reg[y].sum() / T_world[y] for y in TOTAL_YEARS])
print(f"share inside a country: {100 * in_ctry.min():.1f}-{100 * in_ctry.max():.1f}%   "
      f"inside a region: {100 * in_reg.min():.2f}-{100 * in_reg.max():.2f}%")
assert in_ctry.min() > 0.95 and in_reg.min() > 0.999
known = {2002: 3284, 2015: 3670, 2022: 2881, 2023: 3888}         # notebook 03 totals
for y, v in known.items():
    assert abs(T_world[y] - v) < 1.5, f"{y}: {T_world[y]:.1f} vs {v}"
ba22, ba23 = BA_world[:, YC.index(2022)].sum(), BA_world[:, YC.index(2023)].sum()
print(f"burned area 2022 (GFED5.1) {ba22 / 1e6:.2f} M km2, 2023 (NRT) {ba23 / 1e6:.2f} M km2")
assert 0.5 < ba23 / ba22 < 2.0, "burned-area units differ between GFED5.1 and NRT"

# ---- write JSON ---------------------------------------------------------------------------


def r(a, nd=4):
    """Round for JSON; drop trailing zeros."""
    return [float(f"{v:.{nd}f}") for v in np.asarray(a, dtype=float)]


def entity(C, BA, M, T):
    # Classes with no fire at all are left out, so small countries stay small.
    keep = [k for k in range(NC) if C[k].sum() > 0 or BA[k].sum() > 0]
    return {"total": r(T),
            "C": {str(k + 1): r(C[k]) for k in keep},
            "BA": {str(k + 1): r(BA[k], 2) for k in keep},
            "month": {str(k + 1): r(M[k] / len(YC)) for k in keep}}


series = {
    "world": entity(C_world, BA_world, M_world, [T_world[y] for y in TOTAL_YEARS]),
    "regions": {code: entity(C_reg[i], BA_reg[i], M_reg[i], [T_reg[y][i] for y in TOTAL_YEARS])
                for i, (code, _) in enumerate(REGIONS)},
    "countries": {},
}
for i, row in countries.iterrows():
    T = [T_ctry[y][i] for y in TOTAL_YEARS]
    if max(T) > 0:
        series["countries"][row.iso3] = entity(C_ctry[i], BA_ctry[i], M_ctry[i], T)

meta = {
    "generated": date.today().isoformat(),
    "units": {"C": "Tg C per year", "BA": "km2 per year", "month": "Tg C per month (mean)"},
    "years_total": list(TOTAL_YEARS),
    "years_class": YC,
    "source_by_year": {str(y): s for y, s in source.items()},
    "classes": [{"id": k + 1, "name": n} for k, n in enumerate(CLASSES)],
    "groups": [{"key": k, "name": n, "classes": c} for k, n, c in GROUPS],
    "regions": [{"code": c, "name": n} for c, n in REGIONS],
    "countries": [{"iso3": row.iso3, "name": row["name"], "continent": row.continent}
                  for _, row in countries.iterrows() if row.iso3 in series["countries"]],
}
(args.out / "meta.json").write_text(json.dumps(meta, separators=(",", ":")), encoding="utf8")
(args.out / "series.json").write_text(json.dumps(series, separators=(",", ":")), encoding="utf8")

# simplified country outlines for the choropleth
shapes = countries[countries.iso3.isin(series["countries"])][["iso3", "name", "geometry"]].copy()
shapes["geometry"] = shapes.geometry.simplify(0.08, preserve_topology=True)
(args.out / "countries.geojson").write_text(shapes.to_json(drop_id=True, to_wgs84=True)
                                             .replace(".00000000", ""), encoding="utf8")

# ---- tidy CSV downloads -------------------------------------------------------------------
rows = []
names = {"world": [("WORLD", "World", series["world"])],
         "region": [(c, n, series["regions"][c]) for c, n in REGIONS],
         "country": [(m["iso3"], m["name"], series["countries"][m["iso3"]]) for m in meta["countries"]]}
for level, items in names.items():
    for code, name, e in items:
        for k, vals in e["C"].items():
            for y, v, b in zip(YC, vals, e["BA"][k]):
                rows.append((level, code, name, y, source[y], int(k), CLASSES[int(k) - 1], v, b))
pd.DataFrame(rows, columns=["level", "code", "name", "year", "source", "class_id", "class",
                            "carbon_TgC", "burned_area_km2"]) \
  .to_csv(args.out / "csv" / "gfed_by_class_2002_2025.csv", index=False)

rows = [(lvl, code, name, y, source[y], v)
        for lvl, items in names.items() for code, name, e in items
        for y, v in zip(TOTAL_YEARS, e["total"])]
pd.DataFrame(rows, columns=["level", "code", "name", "year", "source", "carbon_TgC"]) \
  .to_csv(args.out / "csv" / "gfed_total_1997_2025.csv", index=False)

for f in sorted(args.out.rglob("*.*")):
    print(f"{f.stat().st_size / 1e6:6.2f} MB  {f.relative_to(args.out)}")
