# Import glacier terminus shapefiles, calculate areas, plot, etc
# T Black, April 2020

# %% Import Modules
import os
import manage
import geopandas as gpd
import gchange
from shapely.ops import polygonize_full
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gplots as gpt

# %% User Specified Parameters
# Absolute path for shapefile containing glacier reference boxes
# "/home/teblack/GreenlandTermini/glacier_boxes_extended.shp"
box_file = '/home/teblack/glacier_boxes_extended.shp'

# Absolute path for shapefiles containing glacier terminus data
# "/home/teblack/GreenlandTermini/NWGreenland_annual_termini/"
path_base = '/home/teblack/NWGreenland_annual_termini/'

# List of glaciers to omit in analysis
omit_glaciers = [76, 77, 78]  # 76-78 disconnected from ice sheet

# Time range to cover
hydroyear_start = 1972
hydroyear_end = 2019
# TODO: incorporate these into data ingestion

# %% Get terminus shapefiles list and reference box information
# Get list of files with terminus data
termini_files = [
    path_base+f for f in os.listdir(path_base) if f.endswith('.shp') and f[0].isdigit()]

# Get glacier reference boxes and filter to NW Greenland study area
boxes = manage.shp2gdf(box_file)
boxes = boxes.loc[1:90]
# remove disconnected glaciers
boxes = boxes.drop(omit_glaciers, errors='ignore')
GIDS = boxes.GlacierID.values

# %% Initialize dictionary of Glacier objects to store glacier metadata and observations
all_observations = {id: manage.Glacier(id) for id in GIDS}
for id in all_observations:
    all_observations[id].box = boxes.loc[id].geometry
    all_observations[id].officialname = boxes.loc[id].Official_n
    all_observations[id].greenlandicname = boxes.loc[id].GrnlndcNam
    all_observations[id].alternativename = boxes.loc[id].AltName

# %% Extract glacier observation records from shapefiles and compile all
# Loop over shapefiles and create glacier observation records
for f in termini_files:
    # import shapefile to geodataframe and remove disconnected glaciers
    termini = manage.shp2gdf(f)
    termini = termini.drop(omit_glaciers, errors='ignore')
    # Loop through each entry and create an observation object with metadata
    for id in GIDS:
        if id not in termini.GlacierID:
            continue
        box = all_observations[id].box
        obs = manage.TerminusObservation(gid=termini.loc[id].GlacierID,
                                         qflag=termini.loc[id].QualFlag,
                                         satellite=termini.loc[id].Satellite,
                                         date=termini.loc[id].Date,
                                         imageID=termini.loc[id].ImageID,
                                         author=termini.loc[id].Author,
                                         geometry=termini.loc[id].geometry)
        obs.area = gchange.glacierArea(obs, box)
        all_observations[id].add_observation(obs)

# Ensure that all glacier timeseries are sorted by date
for id in GIDS:
    all_observations[id].sort_by_date()

# %% Build dataframe of all observed areas, and interpolate to estimate missing areas
# Create an empty dataframe of id vs hydroyear to fill with glacier areas.
YEARS = list(range(hydroyear_start, hydroyear_end+1))
areas = pd.DataFrame(index=GIDS, columns=YEARS, dtype="float64")
for id in GIDS:
    id_hydroyears = all_observations[id].extract("hydroyear")
    id_areas = all_observations[id].extract("area")
    areas.at[id, id_hydroyears] = id_areas
    
# Linearly interpolate between observations for each glacier (axis=1), and do not extrapolate prior to first observation (limit_area='inside')
interpolated_areas = areas.interpolate(method='linear', axis=1, limit_direction='forward')

# %% Calculate various metrics for plotting
total_annual_area = interpolated_areas.agg("sum", axis="index")
change_annual_area = np.diff(total_annual_area)
count_annual_obs = areas.count(axis="index")
count_annual_all = interpolated_areas.count(axis="index")

# %% Plots!
# identify first year with complete data (observations or interpolations for all glaciers)
i = 0
for yr in YEARS:
    if count_annual_all[yr] == len(GIDS):
        first_complete_year = yr
        first_complete_year_idx = i
        break
    else:
        i += 1
YEARS_full = YEARS[first_complete_year_idx:]

# # %% Plot total glaciated area in each hydrological year
# fig = plt.figure()
# gpt.annualArea(fig, YEARS_full, total_annual_area)

# # %% Plot annual change in area
# fig = plt.figure()
# gpt.annualAreaChange(fig, YEARS_full[1:], change_annual_area[first_complete_year_idx:])

# %% Plot total glaciated area and annual area change (as subplots)
fig = plt.figure(figsize=(10, 4.5))
gpt.annualArea(fig, YEARS_full, total_annual_area[first_complete_year_idx:] - total_annual_area.iloc[first_complete_year_idx], subplots=True, idx=121)
gpt.annualAreaChange(fig, YEARS_full[1:], change_annual_area[first_complete_year_idx:], subplots=True, idx=122)

# %% Number of glaciers observed per year
fig = plt.figure()
gpt.numberObserved(fig, YEARS, count_annual_obs)

# %%
