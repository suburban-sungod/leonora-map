"""
Render hillshade tile pyramid from a Vicmap / ELVIS DEM.

Usage:
    python3 make_tiles.py path/to/dem.tif [--zmin 13] [--zmax 18] [--out lidar_tiles]

Reads a single-band elevation GeoTIFF, computes hillshade, reprojects to
Web Mercator, and writes {z}/{x}/{y}.jpg tiles compatible with the
Leaflet 'LiDAR 50cm (local tiles)' layer in index.html.

Deps: rasterio, numpy, Pillow  (all already installed in your env).
"""

import argparse
import math
import os
import sys

import numpy as np
import rasterio
from PIL import Image
from rasterio.warp import Resampling, calculate_default_transform, reproject


TILE_SIZE = 256
WEBMERC_RES = 156543.03392804097  # m/px at z=0 equator


def deg2tile(lat, lon, z):
    x = (lon + 180.0) / 360.0 * (1 << z)
    y = (1 - math.log(math.tan(math.radians(lat)) + 1 / math.cos(math.radians(lat))) / math.pi) / 2 * (1 << z)
    return int(x), int(y)


def tile_bounds_3857(z, x, y):
    n = 1 << z
    res = WEBMERC_RES / n
    minx = -20037508.342789244 + x * TILE_SIZE * res
    maxy = 20037508.342789244 - y * TILE_SIZE * res
    maxx = minx + TILE_SIZE * res
    miny = maxy - TILE_SIZE * res
    return minx, miny, maxx, maxy, res


def hillshade(elev, res_m, azimuth=315, altitude=45, z_factor=1.0):
    """Lambert-shaded hillshade. Returns uint8 [0..255]."""
    az = math.radians(azimuth - 90)  # cartographic → atan2
    alt = math.radians(altitude)
    # central-difference gradients (m/m)
    dzdx = (np.roll(elev, -1, axis=1) - np.roll(elev, 1, axis=1)) * z_factor / (2 * res_m)
    dzdy = (np.roll(elev, -1, axis=0) - np.roll(elev, 1, axis=0)) * z_factor / (2 * res_m)
    slope = np.arctan(np.hypot(dzdx, dzdy))
    aspect = np.arctan2(dzdy, -dzdx)
    sh = np.cos(alt) * np.cos(slope) + np.sin(alt) * np.sin(slope) * np.cos(az - aspect)
    sh = np.clip(sh, 0, 1)
    out = (sh * 255).astype(np.uint8)
    out[0, :] = out[1, :]
    out[-1, :] = out[-2, :]
    out[:, 0] = out[:, 1]
    out[:, -1] = out[:, -2]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dem", help="Path to elevation GeoTIFF from ELVIS")
    ap.add_argument("--zmin", type=int, default=13)
    ap.add_argument("--zmax", type=int, default=18)
    ap.add_argument("--out", default="lidar_tiles")
    ap.add_argument("--azimuth", type=float, default=315)
    ap.add_argument("--altitude", type=float, default=45)
    ap.add_argument("--z_factor", type=float, default=2.5,
                    help="Vertical exaggeration. >1 makes shallow relief pop.")
    args = ap.parse_args()

    print(f"Opening {args.dem}")
    with rasterio.open(args.dem) as src:
        print(f"  CRS: {src.crs}, size: {src.width}x{src.height}, "
              f"res: {src.res}, bounds: {src.bounds}")

        # Reproject DEM to Web Mercator at native-ish resolution.
        dst_crs = "EPSG:3857"
        transform, w, h = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        dem_merc = np.full((h, w), np.nan, dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1),
            destination=dem_merc,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=dst_crs,
            resampling=Resampling.bilinear,
            src_nodata=src.nodata,
            dst_nodata=np.nan,
        )
        # Fill nodata with neighbour mean so hillshade doesn't get holes
        if np.isnan(dem_merc).any():
            mean = np.nanmean(dem_merc)
            dem_merc = np.where(np.isnan(dem_merc), mean, dem_merc)

        res_m = transform.a
        print(f"  Reprojected to Web Mercator at {res_m:.2f} m/px ({w}x{h})")
        print("  Computing hillshade…")
        shade = hillshade(dem_merc, res_m,
                          azimuth=args.azimuth,
                          altitude=args.altitude,
                          z_factor=args.z_factor)
        print(f"  Hillshade range: {shade.min()}..{shade.max()}")

        # Geographic bounds of the source for tile range calc
        with rasterio.Env():
            geo_bounds = rasterio.warp.transform_bounds(src.crs, "EPSG:4326",
                                                       *src.bounds)
        west, south, east, north = geo_bounds
        print(f"  Geo bounds: lon {west:.4f}…{east:.4f}, "
              f"lat {south:.4f}…{north:.4f}")

        mercator_origin_x = -20037508.342789244
        mercator_origin_y = 20037508.342789244

        total_written = 0
        for z in range(args.zmin, args.zmax + 1):
            x_min, y_max_lat = deg2tile(south, west, z)
            x_max, y_min_lat = deg2tile(north, east, z)
            tile_res = WEBMERC_RES / (1 << z)
            print(f"  z={z}: tiles x {x_min}..{x_max}, y {y_min_lat}..{y_max_lat}, "
                  f"{tile_res:.2f} m/px")
            zwritten = 0
            for tx in range(x_min, x_max + 1):
                for ty in range(y_min_lat, y_max_lat + 1):
                    minx, miny, maxx, maxy, _ = tile_bounds_3857(z, tx, ty)
                    # source-array indices for this tile bbox
                    col0 = (minx - transform.c) / transform.a
                    col1 = (maxx - transform.c) / transform.a
                    row0 = (transform.f - maxy) / -transform.e  # transform.e is negative
                    row1 = (transform.f - miny) / -transform.e
                    col0, col1 = sorted([col0, col1])
                    row0, row1 = sorted([row0, row1])
                    if col1 < 0 or row1 < 0 or col0 >= w or row0 >= h:
                        continue
                    # Resample by nearest-neighbour onto 256x256
                    xs = np.linspace(col0, col1, TILE_SIZE)
                    ys = np.linspace(row0, row1, TILE_SIZE)
                    xi = np.clip(xs.astype(int), 0, w - 1)
                    yi = np.clip(ys.astype(int), 0, h - 1)
                    tile = shade[np.ix_(yi, xi)]
                    out_dir = os.path.join(args.out, str(z), str(tx))
                    os.makedirs(out_dir, exist_ok=True)
                    Image.fromarray(tile, mode="L").convert("RGB").save(
                        os.path.join(out_dir, f"{ty}.jpg"),
                        quality=82, optimize=True)
                    zwritten += 1
            print(f"    wrote {zwritten} tiles")
            total_written += zwritten

        print(f"\nDone — {total_written} tiles in {args.out}/")
        print("Next: git add lidar_tiles && git commit -m 'Add LiDAR tiles' && git push")


if __name__ == "__main__":
    sys.exit(main())
