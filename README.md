# GFED Fire Carbon Explorer

An interactive dashboard of fire carbon emissions and burned area from **GFED5.1**, 1997–2025, for the
world, the 14 GFED basis regions and every country, split by land cover (forest, deforestation,
savanna/grassland/shrubland, cropland, peat, tundra & other).

**Live site:** https://rsl-slu.github.io/wildfire/ (once GitHub Pages is switched on, see below)

## What is in it

**Dashboard (`index.html`).** Pick an area (World, one of the 14 GFED biome regions, or a country), a
fire type and a period. Chart cards on the left update for that area: key numbers, emissions every
year by land cover, trend by fire type, forest fires with and without deforestation, what burns,
when fires burn, and the largest emitters.

The map on the right has four views, all on the 0.25° GFED grid in a Robinson projection:

| Map view | Shows |
|---|---|
| Amount | Fire carbon per cell, for the 2002–2025 average or any year |
| Unusual year | A year against each cell's 2002–2025 mean (more or less fire than usual) |
| Trend | Theil–Sen trend per 1° cell, 2002–2025; faded where not significant (Mann–Kendall) |
| Land cover | Which land cover burns most in each cell |

Country or biome outlines are clickable, totals can be shown as circles, and a timeline plays the
years. The page follows the browser's light or dark setting.

**Tables (`explore.html`).** Sortable tables, older summary maps, methods, and CSV downloads.

Trends in the charts are the Theil–Sen slope with a Mann–Kendall test, computed in the browser for
whatever period and fire type you pick.

## Data

| Years | Product | Content |
|---|---|---|
| 1997–2001 | GFED5.1 monthly | Total carbon only |
| 2002–2022 | GFED5.1 ecosystem files | Carbon and burned area in 16 land-cover classes |
| 2023–2025 | GFED5.1 NRT, reprocessed, combined VIIRS (CMB) | Same 16 classes |

GFED5.1 is licensed CC BY 4.0. Cite: van der Werf, G. R. et al. (2025), Global Fire Emissions Database
version 5.1 (https://www.globalfiredata.org). Country boundaries: Natural Earth 1:50m (public domain).

## Layout

```
docs/                    the website (GitHub Pages serves this folder)
  index.html             the dashboard: page, styles and code in one file
  explore.html           tables, methods and downloads
  data/                  JSON the pages read, plus CSV downloads in data/csv/
    grid/                yearly 0.25 deg fire carbon by fire type (lossless PNG, decoded in the page)
    trend/               1 deg trend maps, one per fire type
    lc/                  0.25 deg "what burns most" maps, one per year
  maps/                  summary map images for explore.html
pipeline/
  build_data.py          raw GFED netCDF -> docs/data/*.json, CSV, and cached grids
  make_grid_layers.py    cached grids -> docs/data/grid/ and docs/data/trend/
  make_landcover_png.py  cached grids -> docs/data/lc/
  make_regions.py        GFED basis regions -> docs/data/regions.geojson
  make_maps.py           cached grids -> docs/maps/*.png (explore.html)
```

## Rebuild the data

The raw GFED files are large and are not in this repository. Download them from
https://www.globalfiredata.org, then:

```bash
pip install -r pipeline/requirements.txt
python pipeline/build_data.py --gfed /path/to/GFED5.1 --nrt /path/to/GFED5.1NRT/Reprocessed/Monthly/Ecosystem/CMB \
       --countries /path/to/ne_50m_admin_0_countries.zip
python pipeline/make_grid_layers.py
python pipeline/make_landcover_png.py
python pipeline/make_regions.py --gfed /path/to/GFED5.1
python pipeline/make_maps.py --countries /path/to/ne_50m_admin_0_countries.zip
```

`--gfed` must hold `Ecosystem/GFED5.1_ecosystem_YYYY.nc` (2002–2022) and `Monthly/GFED5.1_monthly_YYYY.nc`
(1997–2001). The script stops with an error if totals drift from known values or less than 95% of fire
carbon lands inside a country.

## View it locally

The page loads its data with `fetch`, so it needs a web server, not a double-click:

```bash
python -m http.server -d docs 8000      # then open http://localhost:8000
```

## Publish

Settings → Pages → Source: *Deploy from a branch*, branch `main`, folder `/docs`.
