#!/usr/bin/env python3

import argparse
import numpy as np
import os
import subprocess
import pandas as pd
import sys
import warnings
from netCDF4 import Dataset

warnings.filterwarnings("ignore", category=DeprecationWarning)

DATA_DIR = "./data"
WEATHER_DIR = "./weather"
LOCATION_FILE = "./location.txt"
ADMIN_FILE = "conusadm.csv"
MASK_URL = "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/nldas/NLDAS_masks-veg-soil.nc4"
MASK_FILE = "NLDAS_masks-veg-soil.nc4"
NC_FIELD = "CONUS_mask"
SEVEN_ZIP = "./bin/7zzs"
LA1 = 25.0625
LO1 = -124.9375
DI = 0.125
DJ = 0.125
STATES = {
    "AK": "Alaska",
    "AL": "Alabama",
    "AR": "Arkansas",
    "AZ": "Arizona",
    "CA": "California",
    "CO": "Colorado",
    "CT": "Connecticut",
    "DE": "Delaware",
    "FL": "Florida",
    "GA": "Georgia",
    "HI": "Hawaii",
    "IA": "Iowa",
    "ID": "Idaho",
    "IL": "Illinois",
    "IN": "Indiana",
    "KS": "Kansas",
    "KY": "Kentucky",
    "LA": "Louisiana",
    "MA": "Massachusetts",
    "MD": "Maryland",
    "ME": "Maine",
    "MI": "Michigan",
    "MN": "Minnesota",
    "MO": "Missouri",
    "MS": "Mississippi",
    "MT": "Montana",
    "NC": "North Carolina",
    "ND": "North Dakota",
    "NE": "Nebraska",
    "NH": "New Hampshire",
    "NJ": "New Jersey",
    "NM": "New Mexico",
    "NV": "Nevada",
    "NY": "New York",
    "OH": "Ohio",
    "OK": "Oklahoma",
    "OR": "Oregon",
    "PA": "Pennsylvania",
    "RI": "Rhode Island",
    "SC": "South Carolina",
    "SD": "South Dakota",
    "TN": "Tennessee",
    "TX": "Texas",
    "UT": "Utah",
    "VA": "Virginia",
    "VT": "Vermont",
    "WA": "Washington",
    "WI": "Wisconsin",
    "WV": "West Virginia",
    "WY": "Wyoming",
    "DC": "District of Columbia",
    "AS": "American Samoa",
    "GU": "Guam",
    "MP": "Northern Mariana Islands",
    "PR": "Puerto Rico",
    "UM": "United States Minor Outlying Islands",
    "VI": "U.S. Virgin Islands",
}


def read_location_file():
    '''Read locations from a location file
    '''

    sites = {}

    # Read in all lines with locations
    with open(LOCATION_FILE, "r") as f:
        lines = [line.strip().split() for line in f
                if line.strip() and line.strip()[0] != "#" and line.strip()[0] != "L"]

    if not lines:
        sys.exit(f"Location file {LOCATION_FILE} is empty.")

    # Parse each line
    for line in lines:
        lat = float(line[0])
        lon = float(line[1])

        if len(line) == 3:          # Site name is defined
            name = line[2]
        else:                       # Site name is not defined
            name = "%.3f%sx%.3f%s" % (abs(lat), "S" if lat < 0.0 else "N", abs(lon), "W" if lon < 0.0 else "E")

        sites[name] = (lat, lon)

    return sites


def read_counties(counties):
    '''Read counties from command line arguments
    '''

    sites = {}

    # Read administrative region list
    admin_df = pd.read_csv(
        ADMIN_FILE,
        usecols=["NAME_1", "NAME_2", "lon", "lat"],
    )
    admin_df.set_index(["NAME_1", "NAME_2"], inplace=True)

    # Find lat/lon of each county
    for c in counties:
        sites[f"{c[0]}_{c[1]}"] = (
            admin_df.loc[STATES[c[1]], c[0]]["lat"], admin_df.loc[STATES[c[1]], c[0]]["lon"]
        )

    return sites


def read_ldas_grids():
    '''Read in LDAS grid information from mask/elevation file

    Use mask/elevation netCDF file to read in the grids, and create a land mask to filter out open water grids.
    '''
    # Download mask/elevation file
    os.makedirs(DATA_DIR, exist_ok=True)
    cmd = [
        "wget",
        "-N",       # Avoid downloading new copies if file already exists
        MASK_URL,
        "-P",
        DATA_DIR,
    ]
    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Read in grids and elevations
    with Dataset(f"{DATA_DIR}/{MASK_FILE}") as nc:
        mask_array = nc[NC_FIELD][0]
        mask_array = np.ma.filled(mask_array.astype(float), np.nan)
        lats, lons = np.meshgrid(nc["lat"][:], nc["lon"][:], indexing="ij")

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


def find_grid(lat, lon, coord, mask_array):
    '''Find closest land grid to an input site

    This function finds the closest unmasked grid and the closest masked grid to the specified site. By comparing the
    two grids, it will determine if the specified grid is a land point.
    '''

    closest = (round((lat - LA1) / DJ), round((lon - LO1) / DI))

    if np.isnan(mask_array[closest]):   # If closest grid is a water/non-CONUS grid
        closest = closest_grid(lat, lon, coord)
        print(f"Nearest NLDAS land grid to {lat:.3f}x{lon:.3f} is {coord[0][closest]}x{coord[1][closest]}")

    _lat = coord[0][closest]
    _lon = coord[1][closest]

    grid = f"%.3f%sx%.3f%s" % (abs(_lat), "S" if _lat < 0.0 else "N", abs(_lon), "W" if _lon < 0.0 else "E")

    return grid


def extract_weather_file(site, latlon):
    '''Extract weather files from weather file archive
    '''

    print(f"Generate weather file {WEATHER_DIR}/{site}.weather")
    # Find grid for specified location
    coord, mask = read_ldas_grids()
    grid = find_grid(latlon[0], latlon[1], coord, mask)

    # Write LDAS grid to weather file for reference
    with open(f"{WEATHER_DIR}/{site}.weather", "w") as f:
        f.write(f"# NLDAS grid {grid}\n")

    # Extract from archive
    cmd = f"{SEVEN_ZIP} e NLDAS_weather/NLDAS_CONUS_1979-2021.7z NLDAS_{grid}.weather -so >> {WEATHER_DIR}/{site}.weather"
    subprocess.run(
        cmd,
        shell=True,
    )


def main(params):
    # Get locations from location file or command line arguments
    if not params["county"]:
        print("Use location file")
        sites = read_location_file()
    else:
        sites = read_counties(params["county"])

    # Extract weather files for each location
    os.makedirs(WEATHER_DIR, exist_ok=True)
    [extract_weather_file(s, sites[s]) for s in sites]


def _main():
    parser = argparse.ArgumentParser(description="Extract Cycles weather files from archive")

    parser.add_argument(
        "--county",
        nargs="+",
        type=lambda arg : arg.split(","),
    )

    args = parser.parse_args()

    main(vars(args))


if __name__ == "__main__":
    _main()
