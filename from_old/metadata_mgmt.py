# This script handles different metadata formats for termini projects.
# Note: FID in ArcGIS is represented in python-geopandas as an index.
# Note: Shape in ArcGIS is represented in python-geopandas as a geometry.
# 
# - nsidc-v2 (NSIDC-0723 MEaSUREs, version 2):
#   FID* | Shape* | GlacierID | QualFlag | DATE | GrnlndcNam | Official_n | AltName | ImgSource
#
# - nwgris (NW Greenland long-term)
#   FID* | Shape* | GlacierID | QualFlag | Satellite | Date | ImageID | Author
#
# - sentinel (Sentinel-1 sub-seasonal/everything)
#   FID* | Shape* | GlacierID | QualFlag | TRACK | DATE | ORBIT | SAT
# 
# - icepicks
#   FID* | Shape* | GlacierID | QualFlag | Satellite | Date | ImageID | Author
# 
# Taryn Black, May 2020

import geopandas as gpd

# %% USER-SPECIFIED PARAMETERS
termini_file = "/mnt/d/GreenlandMapping/termini_20200204_20200209.shp"
metadata_standard = "nsidc-v2"
updated_filename = "/mnt/d/GreenlandMapping/termini_1920.shp"

# Required for nsidc_v2 - file including various glacier names
glacierid_file = "/mnt/d/GreenlandMapping/glacier_info/glacier_ids/GlacierIDs.shp"

# %% LOAD FILE
termini = gpd.read_file(termini_file)

print("Analyzing metadata for file   %s" % termini_file.rsplit("/")[-1])
print("Original terminus file attributes:")
print(" | ".join(termini.columns.to_list()))

if glacierid_file:
    glacierid = gpd.read_file(glacierid_file)
    glacierid = glacierid.drop(labels=["POINT_X", "POINT_Y", "geometry"], axis="columns")

# %% DEFINE REQUIRED KEYS FOR EACH METADATA STANDARD
standard_keys = {
    "nsidc-v2" : ["GlacierID", "QualFlag", "DATE", "GrnlndcNam", "Official_n",  "AltName", "ImgSource", "geometry"],
    "nwgris" : ["GlacierID", "QualFlag", "Satellite", "Date", "ImageID", "Author", "geometry"],
    "sentinel" : ["GlacierID", "QualFlag", "TRACK", "DATE", "ORBIT", "SAT", "geometry"],
    "icepicks" : ["GlacierID", "QualFlag", "Satellite", "Date", "ImageID", "Author", "geometry"]
}

print("\nMetadata standard: %s -- with attributes:" % metadata_standard)
print(" | ".join(standard_keys[metadata_standard]))

# %% CHECK WHETHER METADATA STANDARDS ARE PRESENT IN TERMINUS KEYS
print("\nChecking terminus attributes against metadata standards...")
for key in standard_keys[metadata_standard]:
    if key in termini.keys():
        continue
    else:
        print("- Adding " + key)
        if key == "GlacierID":
            print("Error! The GlacierID must be defined during mapping.")
            break
        elif key == "QualFlag":
            print("Error! The QualFlag must be defined during mapping.")
            break
        elif key in ["DATE", "Date"]:
            if "DATE" in standard_keys[metadata_standard]:
                termini.rename(columns={"Date": "DATE"}, inplace=True)
            elif "Date" in standard_keys[metadata_standard]:
                termini.rename(columns={"DATE": "Date"}, inplace=True)
            else:
                print("Error! The Date must be defined during mapping.")
                break
        elif key in ["GrnlndcNam", "Official_n", "AltName"]:
            termini = termini.merge(glacierid, on="GlacierID")
        elif key in ["ImgSource", "Satellite"]:
            satellite = input("--> Which satellite did you use for these termini? (e.g. Landsat-8, Sentinel-1) : ")
            satellite = satellite.title()
            termini[key] = satellite
        elif key == "ImageID":
            image = input("--> Did you trace all of these termini from one image file? (Y/N) : ")
            if image.upper() == "Y":
                image_id = input("--> What is the unique file identifier? ")
                termini[key] = image_id
            elif image.upper() == "N":
                print("That sucks. Should've tracked that during tracing! (error!)")
                break
            else: 
                print("Error! Select Y or N.")
                break
        elif key == "Author":
            author = input("--> Who traced these termini? LASTNAME_FIRSTNAME : ")
            author = "_".join(author.split()).upper()
            termini[key] = author
        elif key in ["TRACK", "DATE", "ORBIT", "SAT"]:
            print("Error! Run sentinel_metadata.py on this file first.")
            break
        elif key == "geometry":
            print("Error! File must have geometry (shapes)!")
            break
        else:
            print("Oh no! %s hasn't been defined as a key." % key)
            break

# %% REMOVE TERMINUS KEYS THAT AREN'T IN CHOSEN METADATA STANDARDS
for key in termini.keys():
    if key not in standard_keys[metadata_standard]:
        termini = termini.drop(labels=key, axis='columns')

print("\nFinal terminus file attributes:")
print(" | ".join(termini.columns.to_list()))

# %% SAVE AS NEW SHAPEFILE
print("\nWriting file: %s" % updated_filename)
termini.to_file(updated_filename)
