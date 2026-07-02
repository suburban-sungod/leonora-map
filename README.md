# Leonora map

One-off public deploy of the Leonora–Gwalia gold research map, same UI as [coopers-creek-map](https://github.com/suburban-sungod/coopers-creek-map).

- `index.html` — Leaflet map: Minedex sites, mining tenements, GSWA 1:100K interpreted bedrock, faults/shears. All overlays fetched live from Landgate SLIP public services.
- `site_articles.json` — per-mine Trove articles, merged into Minedex popups by site_code.
- `lidar_tiles/` — 30m Copernicus GLO-30 DEM hillshade (z10-14). No open LiDAR exists for this area (checked ELVIS 2026-07-02); this gives relief context only — drainage, breakaways, the Gwalia pit.

Generated from the gold-research skill template. Source of truth: `~/Documents/Projects/gold-research/leonora/`.
