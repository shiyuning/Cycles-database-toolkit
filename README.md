# Cycles-database-toolkit



- `gadm.py`:
  Read the global administrative region database, and generate `gadm.csv` and `conusadm.csv` with the names and locations (latitudes and longitudes) of all Level-3 (e.g., county level) administrative regions, which can be used to generate crop lookup tables.
- `lookup_table.py`:
  Convert `v1.9` lookup tables to `v2.0`.
  The conversion can be global, or for CONUS only.

  Soil archive name errors in `v1.9` are fixed;
  region names, and their latitudes and longitudes of crop areas (converted from World_Goode_Homolosine_Land projection) are added to the tables;
  and LDAS (GLDAS for global and NLDAS for CONUS) grids are used for weather file lookup.
- `soil_files.py`:
  Convert `v1.9` soil files to `v2.0`.

  Soil archive name errors in `v1.9` are fixed;
  and soil file formats are updated to be compatible with latest Cycles versions.

- `NLDAS-to-Cycles-CONUS.py`:
  Generate Cycles weather files for all CONUS grids

- `archive-by-state.by`:
  Archive CONUS weather files by state.