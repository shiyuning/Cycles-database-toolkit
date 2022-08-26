#!/usr/bin/env python3

import geopandas as gpd

SHP36 = "gadm36.shp"
GADM = "gadm.csv"
CONUSADM = "conusadm.csv"
DATA_DIR = "./data"

def main():
    # Read in gadm shapefile
    print("Reading shapefile...")
    t = gpd.read_file(f"{DATA_DIR}/{SHP36}")

    # Calculate centroids
    print("Calculating centroids...")
    t["centroid"] = t.centroid
    t.drop("geometry", inplace=True, axis=1)    # Drop geometry after calculating centroids to save time
    t.set_geometry("centroid", inplace=True)
    t.sort_values(by=['UID'], inplace=True)

    # Convert gadm coordinates to WGS84
    print("Converting coordinates...")
    t_wgs84 = t.to_crs(epsg=4326)
    t_wgs84["lat"] = t_wgs84["centroid"].y
    t_wgs84["lon"] = t_wgs84["centroid"].x

    # Write to csv files
    t_out = t_wgs84[["UID", "NAME_0", "NAME_1", "NAME_2", "lon", "lat"]]
    t_out.to_csv(GADM, index=False)
    ## Filter CONUS areas
    t_out[(t_out["NAME_0"] == "United States") & (t_out["NAME_1"] != "Hawaii") & (t_out["NAME_1"] != "Alaska")].to_csv(CONUSADM, index=False)


if __name__ == "__main__":
    main()
