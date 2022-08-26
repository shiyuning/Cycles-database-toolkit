#!/usr/bin/env python3

import os
import pandas as pd
import subprocess

VERSION = "2.0"
SEVEN_ZIP = "./bin/7zzs"
COMPRESS = "./bin/compress.sh"
TEMP = "temp0"
TEMP_NEW = "temp1"
DATA_DIR = "data"
LOOKUP_V19_DIR = "lookup_1.9"
SOIL_V19_DIR = {
    "avg": "avg_1.9",
    "major": "major_1.9",
}
NAME = "soil_global"
TAB = "\t"
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

TYPES = [
    "avg",
    "major",
]

VARS = [
    "LAYER",
    "THICK",
    "CLAY",
    "SAND",
    "ORGANIC",
    "BD",
    "FC",
    "PWP",
    "SON",
    "NO3",
    "NH4",
    "BYP_H",
    "BYP_V",
]

def convert_soil(f, f_new):
    """Convert soil files in database to latest Cycles format
    """
    NO3 = ["10", "10", "7", "4", "4", "4"]

    # Read original soil file
    soil_str = []
    comments = []
    with open(f"{TEMP}/{f}", "r") as fp:
        for line in fp:
            if line.strip() and line.strip()[0] != "#":
                soil_str.append(line.strip())
            else:
                comments.append(line)

    # Write new soil file
    with open(f"{TEMP_NEW}/{f_new}", "w") as fp:
        ## Write comment lines
        comments[2] = comments[2].replace("v1.9", f"v{VERSION}")
        fp.writelines(comments)

        ## Curve number line
        fp.write("%-15.15s\t" % soil_str[0].split()[0])
        fp.write("%-7.7s\t" % soil_str[0].split()[1])
        fp.write("%s\n" % " ".join(soil_str[0].split()[2:]))
        ## Slope line
        fp.write("%-15.15s\t" % soil_str[1].split()[0])
        try:
            fp.write("%-5.2f\n" % float(soil_str[1].split()[1]))
        except:
            pass
        fp.write(soil_str[2] + "\n")    # Total number of layers line

        ## Header line
        hl = VARS + soil_str[3].split()[10:]
        fp.write(f"{TAB.join(format(n, '7.7s') for n in hl)}")
        fp.write("\n")

        ## Layers
        for kline in range(4, len(soil_str)):
            tmp = soil_str[kline].split()[0:8]
            tmp.append("-999")  # Add SON
            tmp.append(NO3[int(soil_str[kline].split()[0]) - 1])    # Add NO3
            tmp.append("1")  # Add NH4
            tmp.extend(["0.0", "0.0"])  # Add bypass parameters
            tmp.extend(soil_str[kline].split()[10:])
            fp.write(f"{TAB.join(format(c, '7.7s') for c in tmp)}")
            fp.write("\n")

def main():
    # Create new directories for temp files and new archives
    for t in TYPES:
        os.makedirs(f"{NAME}_{t}_{VERSION}", exist_ok=True)

    for crop in CROPS:
        for lu in LUIDS:
            print(crop, lu)
            # Read look-up tables
            df = pd.read_csv(
                f"{DATA_DIR}/{LOOKUP_V19_DIR}/{crop}_{lu}_v9_Lookup.csv",
                usecols=["Gadm_ID", "avg_soil_file", "major_soil_file"],
                index_col="Gadm_ID",
            )

            for t in TYPES:
                os.makedirs(TEMP, exist_ok=True)
                os.makedirs(TEMP_NEW, exist_ok=True)

                archive = f"./{DATA_DIR}/{SOIL_V19_DIR[t]}/{crop}{lu}_v9_{t}.7z"
                cmd = [
                    SEVEN_ZIP,
                    "e",
                    f"-o{TEMP}",
                    "-y",
                    archive,
                ]
                subprocess.run(cmd)

                for _, row in df.iterrows():
                    f = row[f"{t}_soil_file"]

                    if os.path.exists(f"{TEMP}/{f}"):
                        ## Fix soil file name
                        f_new = f.replace("_v", "_")

                        # Convert soil file format
                        convert_soil(f, f_new)

                # Add new soil files to archive
                archive = f"{NAME}_{t}_{VERSION}/{crop}_{lu}_{t}_{VERSION}.7z"
                cmd = [
                    COMPRESS,
                    TEMP_NEW,
                    archive,
                    f"{crop}_{lu}_sg_{t}",
                    SEVEN_ZIP,
                ]
                subprocess.run(
                    cmd,
                    #stdout=subprocess.DEVNULL,
                    #stderr=subprocess.DEVNULL,
                )

                # Delete temp files
                cmd = [
                    "rm",
                    "-rf",
                    TEMP,
                    TEMP_NEW,
                ]
                subprocess.run(cmd)


if __name__ == "__main__":
    main()
