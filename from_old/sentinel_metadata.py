# This script adds Sentinel-1 swath metadata to glacier termini, depending which swath the glacier falls in. The swaths must be manipulated to remove overlaps and preserve only the topmost swaths
# Taryn Black, May 2020
 
import geopandas as gpd

# %% USER-SPECIFIED PARAMETERS
# Shapefile containing metadata for each swath of the Sentinel-1 mosaic
sentinel_file = "/mnt/d/GreenlandMapping/GLMosaic_uncalibrated.2020-02-04.2020-02-09.shp"

# Shapefile containing glacier termini associated with Sentinel-1 mosaic
termini_file = "/mnt/d/GreenlandMapping/termini_20200204_20200209.shp"

# %% LOAD DATA
swaths = gpd.read_file(sentinel_file)
swaths = swaths.to_crs(epsg=3413)
termini = gpd.read_file(termini_file)

# %% PREPARE SENTINEL SWATH POLYGONS FOR OVERLAP ERASE
# Reset index to save original order of layers as a column
swaths = swaths.reset_index()
# swaths.plot(ec='black', alpha=0.5)

# Dissolve multiple polygons into single polygons per track
swaths = swaths.dissolve(by="TRACK", as_index=False)
# swaths.plot(ec='black', alpha=0.5)

# %% ERASE LOWER POLYGON WHERE SWATH POLYGONS OVERLAP
# First, reverse sort on original index so that polygons are ordered from bottom to top
swaths = swaths.sort_values('index', ascending=False)

# Loop through each polygon and erase underlying polygons where they overlap
for p in range(0, len(swaths)):
    swaths['geometry'].iloc[p+1:] = swaths['geometry'].iloc[p+1:].difference(swaths['geometry'].iloc[p])
# swaths.plot(ec='black', alpha=0.5)

# %% CLEAN SWATH METADATA IN PREPARATION FOR JOIN
# Remove extraneous attributes from swaths dataframe
swaths = swaths.drop(axis='columns', labels='index')

# Convert swaths columns to non-object dtypes
swaths['TRACK'] = swaths['TRACK'].astype('int64')
swaths['DATE'] = swaths['DATE'].astype('str')
swaths['ORBIT'] = swaths['ORBIT'].astype('int64')
swaths['SAT'] = swaths['SAT'].astype('str')

# %% JOIN SWATH METADATA TO TERMINUS METADATA BASED ON LOCATION
# Spatially join swath metadata to their interesecting termini
termini_with_dates = gpd.sjoin(termini, swaths, how='left')

# ax = swaths.plot(alpha=0.7, ec='white')
# termini.plot(ax=ax, color='red', lw=0.1)

# %% REMOVE DUPLICATED TERMINI
## TODO: make this work for multiple duplicates (usually it's just GID #120, but just in case)
# Find duplicates (termini that cross >1 swath show up twice in dataframe after spatial join)
dupes = termini_with_dates[termini_with_dates.GlacierID.duplicated(keep=False)]

# Intersect the duplicated termini with the swaths to get segments over each swath
dupes_sections = gpd.overlay(dupes, swaths, 'intersection')

# Identify the track of the swath containing the largest section of the terminus
keep_track = dupes_sections[dupes_sections.length == dupes_sections.length.max()].TRACK_2.iloc[0]

# Keep only the terminus record associated with the larger section's swath
keep_dupe = dupes[dupes.TRACK == keep_track]

# Drop duplicates from the original dataframe and add back the one to keep
termini_with_dates = termini_with_dates.drop_duplicates(subset='GlacierID', keep=False)
termini_with_dates = termini_with_dates.append(keep_dupe)

# %% ADDITIONAL METADATA CLEANING FOR TERMINUS DATAFRAME
# Remove index_right column (extraneous attribute from spatial join)
termini_with_dates = termini_with_dates.drop(labels='index_right', axis='columns')

# Set GlacierID as index (other index is extraneous and will be treated as an attribute in ArcGIS)
termini_with_dates = termini_with_dates.set_index('GlacierID')

# Sort by GlacierID
termini_with_dates = termini_with_dates.sort_values('GlacierID')

# %% SAVE THE FINALIZED TERMINUS DATA BACK TO THE ORIGINAL FILE
termini_with_dates.to_file(termini_file)
