#!/usr/bin/env python3

import argparse
import geopandas as gpd
import numpy as np
import os
import pandas as pd
import subprocess
from netCDF4 import Dataset
from shapely.geometry import Point

VERSION = "2.0"

CROPS = [
    "Bean",
    "Cassava",
    "Lentil",
    "Maize",
    "Millet",
    "Potato",
    "Rice",
    "Sorghum",
    "Soybean",
    "SweetPotato",
    "Wheat",
]

LUIDS = [
    "10",
    "11",
    "12",
    "20",
    "30",
    "40",
]

V19_DIR = "lookup_1.9"
DATA_DIR = "./data"
SHP36 = "gadm36.shp"

ADMIN_FILES = {
    "CONUS": "conusadm.csv",
    "global": "gadm.csv",
}

MASK_URLS = {
    "CONUS": "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/nldas/NLDAS_masks-veg-soil.nc4",
    "global": "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/gldas/ELEV/GLDASp4_elevation_025d.nc4",
}
MASK_FILES = {
    "CONUS": "NLDAS_masks-veg-soil.nc4",
    "global": "GLDASp4_elevation_025d.nc4",
}
MASK_FIELDS = {
    "CONUS": "CONUS_mask",
    "global": "GLDAS_elevation",
}
LDAS = {
    "CONUS": "NLDAS",
    "global": "GLDAS",
}
LA1S = {
    "CONUS": 25.0625,
    "global": -59.875,
}
LO1S = {
    "CONUS": -124.9375,
    "global": -179.875,
}
DIS = {
    "CONUS": 0.125,
    "global": 0.25,
}
DJS = {
    "CONUS": 0.125,
    "global": 0.25,
}
NIS = {
    "CONUS": 464,
    "global": 1440,
}


def read_ldas_grids(range):
    '''Read in LDAS grid information from mask/elevation file

    Use mask/elevation netCDF file to read in the grids, and create a land mask to filter out open water grids.
    '''
    # Download mask/elevation file
    cmd = [
        "wget",
        "-N",       # Avoid downloading new copies if file already exists
        MASK_URLS[range],
        "-P",
        DATA_DIR,
    ]
    subprocess.run(
        cmd,
    )

    # Read in grids and elevations
    with Dataset(f"{DATA_DIR}/{MASK_FILES[range]}") as nc:
        mask_array = nc[MASK_FIELDS[range]][0]
        mask_array = np.ma.filled(mask_array.astype(float), np.nan)
        lats, lons = np.meshgrid(nc["lat"][:], nc["lon"][:], indexing="ij")

    if range == "CONUS":
        # NLDAS CONUS mask is labeled 0/1. Convert 0 to NAN to be consistent with GLDAS
        mask_array[mask_array == 0] = np.nan

    # Mask sea/non-CONUS grid lat/lon as nan
    lats[np.isnan(mask_array)] = np.nan
    lons[np.isnan(mask_array)] = np.nan

    return [lats, lons], mask_array


def closest_grid(lat, lon, coord):
    '''Find closest grid to an input site
    '''
    lats = coord[0]
    lons = coord[1]
    dist = np.sqrt((lons - lon)**2 + (lats - lat)**2)
    closest = np.unravel_index(np.argmin(dist, axis=None), dist.shape)

    return closest


def find_grid(range, lat, lon, coord, mask_array):
    '''Find closest land grid to an input site

    This function finds the closest unmasked grid and the closest masked grid to the specified site. By comparing the
    two grids, it will determine if the specified grid is a land point.
    '''

    closest = (round((lat - LA1S[range]) / DJS[range]), round((lon - LO1S[range]) / DIS[range]))

    if np.isnan(mask_array[closest]):   # If closest grid is a sea/non-CONUS grid
        closest = closest_grid(lat, lon, coord)
        print(f"Nearest {LDAS[range]} land grid to {lat:.3f}x{lon:.3f} is {coord[0][closest]}x{coord[1][closest]}")

    _lat = coord[0][closest]
    _lon = coord[1][closest]

    grid = f"%.3f%sx%.3f%s" % (abs(_lat), "S" if _lat < 0.0 else "N", abs(_lon), "W" if _lon < 0.0 else "E")

    return grid


def main(params):

    range = params["range"]

    # Read LDAS grid and elevation data
    print(f"Read {LDAS[range]} grids...")
    coord, mask_array = read_ldas_grids(range)

    # Create directory for lookup files
    os.makedirs(f"crop_lookup_{range}_{VERSION}", exist_ok=True)

    # Read projection
    print("Read projection...")
    t = gpd.read_file(f"{DATA_DIR}/{SHP36}", bbox=(0,0,0,0))

    # Read administrative regions
    print(f"Read administrative regions...")
    admin_df = pd.read_csv(
        ADMIN_FILES[range],
        usecols=["UID", "NAME_0", "NAME_1", "NAME_2"],
        index_col="UID",
    )

    # Convert crop lookup tables
    print("Convert crop look-up tables:")
    for crop in CROPS:
        for lu in LUIDS:
            print(" ", crop, lu)
            df = pd.read_csv(
                f"{DATA_DIR}/{V19_DIR}/{crop}_{lu}_v9_Lookup.csv",
                usecols=["Gadm_ID", "Crop_X", "Crop_Y", "avg_soil_file", "major_soil_file"],
                index_col="Gadm_ID",
            )

            df = df.join(admin_df, how="inner")

            if not df.empty:
                print("    Convert coordinate...")
                df["avg_soil_file"] = df["avg_soil_file"].str.replace("_v", "_")    # Fix soil file names
                df["major_soil_file"] = df["major_soil_file"].str.replace("_v", "_")    # Fix soil file names
                df["geometry"] = df.apply(lambda x: Point(x["Crop_X"], x["Crop_Y"]), axis=1)

                ## Convert coordinates
                df = gpd.GeoDataFrame(df, crs=t.crs)
                df = df.to_crs(epsg=4326)
                df["lat"] = df["geometry"].y
                df["lon"] = df["geometry"].x

                ## Find LDAS grids
                print(f"    Project to {LDAS[range]} grids...")
                df[f"{LDAS[range]}_grid"] = df.apply(lambda x: find_grid(range, x["lat"], x["lon"], coord, mask_array), axis=1)

                ## Write to file
                print(f"    Write to file...")
                df = df[["NAME_0", "NAME_1", "NAME_2", "lat", "lon", f"{LDAS[range]}_grid", "avg_soil_file", "major_soil_file"]]
                df.index.name = "GADM_ID"
                df.sort_values(by=['GADM_ID'], inplace=True)

                fn = f"crop_lookup_{range}_{VERSION}/{crop}_{lu}_lookup_{range}_{VERSION}.csv"
                df.to_csv(
                    fn,
                )

def _main():
    parser = argparse.ArgumentParser(description="Convert look-up tables")
    parser.add_argument(
        "--range",
        default="CONUS",
        choices={"global", "CONUS"},
        help="Range",
    )
    args = parser.parse_args()

    main(vars(args))

if __name__ == "__main__":
    _main()
