#!/usr/bin/env python3

import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import os
import subprocess
from netCDF4 import Dataset

DATA_DIR = "../data"
STATE_DIR = "../by_state"
MASK_URL = "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/nldas/NLDAS_masks-veg-soil.nc4"
MASK_FILE = "NLDAS_masks-veg-soil.nc4"
NC_FIELDS = {
    "MASK": "CONUS_mask",
}
SHP = "cb_2018_us_state_500k.shp"
SEVEN_ZIP = "../bin/7zzs"

cmd = [
    "wget",
    "-N",       # Avoid downloading new copies if file already exists
    MASK_URL,
    "-P",
    f"{DATA_DIR}",
]
subprocess.run(
    cmd,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)

with Dataset(f"{DATA_DIR}/{MASK_FILE}") as nc:
    mask_array = nc[NC_FIELDS["MASK"]][0]

    lats, lons = np.meshgrid(nc["lat"][:], nc["lon"][:], indexing="ij")

    lats = lats.flatten()
    lons = lons.flatten()
    masks = mask_array.flatten()

grids = []

for k in range(mask_array.size):
    #if not math.isnan(elev_array[np.unravel_index(k, elev_array.shape)]):
    if masks[k] == 1:
        grids.append(k)

t = gpd.read_file(f"{DATA_DIR}/{SHP}")
t.set_index("STUSPS", inplace=True)
t_wgs84 = t.to_crs(epsg=4326)

os.makedirs(STATE_DIR, exist_ok=True)
for k in grids:
    grid_lat = lats[k]
    grid_lon = lons[k]
    point = Point(grid_lon, grid_lat)
    try:
        state = t_wgs84.index[t_wgs84.contains(point)].to_list()[0]
    except:
        state = t_wgs84.distance(point).sort_values().index[0]

    grid = "%.3f%sx%.3f%s" % (
        abs(grid_lat), "S" if grid_lat < 0.0 else "N", abs(grid_lon), "W" if grid_lon < 0.0 else "E"
    )

    f = f"NLDAS_{grid}.weather"
    cmd = [
        SEVEN_ZIP,
        "a",
        f"{STATE_DIR}/NLDAS_{state}_1979-2021.7z",
        f,
    ]
    subprocess.run(
        cmd,
        #stdout=subprocess.DEVNULL,
        #stderr=subprocess.DEVNULL,
    )
