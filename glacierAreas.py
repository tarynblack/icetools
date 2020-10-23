# Get glacier areas from shapefiles
# THIS IS THE GOOD FILE - Dec 2019

import os
import geopandas as gpd
from shapely.ops import polygonize_full
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import cm

# ////////// USER-SPECIFIED PARAMETERS //////////
# box_file : absolute path for shapefile containing glacier reference boxes
# termini_files : absolute path(s) for shapefile(s) containing glacier termini
#                 data. Attributes must include: Date, QualFlag, geometry.
box_file = "/home/teblack/GreenlandTermini/glacier_boxes_extended.shp"
path_base = "/home/teblack/GreenlandTermini/NWGreenland_annual_termini/"
termini_files = [path_base+f for f in os.listdir(path_base) if
                 f.endswith(".shp") and f[0].isdigit()]
# ///////////////////////////////////////////////


def geomReadReprojectReindex(file, epsg=3574):
    """Reads a shapefile of glacier data and reprojects to a specified EPSG.
    Default EPSG:3574 (WGS 84 / North Pole LAEA Atlantic) (unit=meters)
    Result is a geodataframe containing the shapefile data.
    Reindex the gdf so that index=GlacierID for direct selection by ID."""
    gdf = gpd.read_file(file)
    gdf = gdf.to_crs(epsg=epsg)
    # TODO: check whether glacier has multiple entries
    # (then will have multiple indices of that value, which screws up indexing)
    gdf = gdf.sort_values(by='GlacierID').set_index('GlacierID', drop=False)
    return gdf


def getGlacierInfo(termini_gdf, box_gdf, id):
    """Get terminus and metadata for a given glacier ID in the geodataframe of
    termini data, as well as the glacier's reference box."""
    terminus = termini_gdf.loc[id]
    box = box_gdf.loc[id]
    return terminus, box


def getGlacierArea(terminus, box):
    """Get area of polygon created by intersection of glacier terminus and
    reference box. Default area for EPSG:3574 is m2; return in km2."""
    outline = terminus.geometry.union(box.geometry)
    poly = polygonize_full(outline)[0]
    if poly.is_empty:
        print("%s: Glacier %s trace and box do not overlap" %
              (terminus.Date, terminus.GlacierID))
    area_km = poly.area / 10**6
    return area_km


def getHydrologicalYear(date):
    """Determine hydrological year of a given date. Greenland hydrological year
    is defined as September 1 through August 31. Returns the starting year of
    the hydrological year (aligned with September)."""
    date = pd.to_datetime(date)
    if pd.notnull(date):
        if date.month >= 9:
            hydroyear = date.year
        elif date.month < 9:
            hydroyear = date.year - 1
        return hydroyear


def dayOfHydroyear(date):
    """Convert date to number of days since start of hydrological year (dohy),
    defined as September 1. Analagous to day-of-year."""
    date = pd.to_datetime(date)
    if pd.notnull(date):
        hydroyear = getHydrologicalYear(date)
        start_hydroyear = pd.to_datetime("%s-09-01" % str(hydroyear))
        dohy = (date - start_hydroyear).days
        return dohy


def get_season(date):
    """Input a date and return the climatological season of the date, as int
    (int necessary for applying this function to row/col of a dataframe).
    Season/integer mapping:
        1 : Autumn (Sept. 1 through Nov. 30)
        2 : Winter (Dec. 1 through Feb. 28)
        3 : Spring (Mar. 1 through May 31)
        4 : Summer (Jun. 1 through Aug. 31)"""
    date = pd.to_datetime(date)
    if date.month in [9, 10, 11]:
        return 1 #'AUT'
    elif date.month in [12, 1, 2]:
        return 2 #'WIN'
    elif date.month in [3, 4, 5]:
        return 3 #'SPR'
    elif date.month in [6, 7, 8]:
        return 4 #'SUM'
season_dict = {1 : 'Autumn (SON)',
               2 : 'Winter (DJF)',
               3 : 'Spring (MAM)',
               4 : 'Summer (JJA)'}


def createGlacierRecord(terminus, area):
    """Create a tuple representing a single glacier record, containing:
        (date, area (km^2), qualflag, hydrological_year)"""
# TODO: try making all of this a class, turn unzip into a method
    date = pd.to_datetime(terminus.Date)
    hydroyear = getHydrologicalYear(date)
    record = (date, area, terminus.QualFlag, hydroyear)
    return record


def unzipGlacierRecords(glacier):
    """Take in the list of individual records for a glacier, and return each
    component of the records (date, area, etc.) as its own list."""
    dates = [obs[0] for obs in glacier]
    areas = [obs[1] for obs in glacier]
    qualflags = [obs[2] for obs in glacier]
    hydroyears = [obs[3] for obs in glacier]
    return dates, areas, qualflags, hydroyears


def normGlacierArea(areas):
    """Normalize glacier area such that 0=smallest extent and 1=largest extent
    in the area time series. Input is a list of all areas for the glacier."""
    normalized_area = [(a - np.nanmin(areas)) / (np.nanmax(areas) - np.nanmin(areas))
                       for a in areas]
    return normalized_area


# ////////// DO THE THING //////////

boxes = geomReadReprojectReindex(box_file)
boxes = boxes.loc[1:92]  # Only want NW Greenland glaciers

# Initialize a dictionary to store glacier records
# Key = glacier ID
# Value = list of records, record=(date, area, qualflag)
gIDs = boxes.GlacierID.values
glacier_records = {k: [] for k in gIDs}

# Get a list of all years (start of hydroyear) that were loaded
all_years = [int(y[len(path_base):len(path_base)+4]) for y in termini_files]
all_years.sort()

# Also initialize a dataframe to store just areas (maybe easier to work with)
glacier_area_df = pd.DataFrame(index=gIDs, columns=all_years)

# Create a text file to track which traces and boxes don't overlap
review_traces = open(path_base+"review_traces.txt", "w")

# Loop over all shapefiles, calculate area, and add data to dict and df
for f in termini_files:
    termini = geomReadReprojectReindex(f)
    for id in termini.GlacierID:
        terminus, box = getGlacierInfo(termini, boxes, id)
        area_km = getGlacierArea(terminus, box)
        if area_km == 0:
            review_traces.write("%s %s\n" % (terminus.Date, id))
        record = createGlacierRecord(terminus, area_km)
        # TODO: check whether date record exists for glacier before appending
        glacier_records[id].append(record)
        hydroyear = record[3]
        glacier_area_df.at[id, hydroyear] = area_km
    review_traces.write("\n")
review_traces.close()

# Sort all records by date
for g in glacier_records:
    glacier_records[g].sort()
    
# Convert DF columns from object to float for interpolation
glacier_area_df = glacier_area_df.apply(pd.to_numeric)
# Linearly interpolate areas between observations for each glacier (axis=1).
# Do not extrapolate prior to first observation (limit_area='inside').
glacier_area_df_fill = glacier_area_df.interpolate(method='linear', axis=1,
                                                   limit_area='inside')

LINEWIDTH = 5
MARKERSIZE = 15
TITLESIZE = 24
LABELSIZE = 20
TICKLABELSIZE = 16
LEGENDSIZE = 16

# Sum areas for each hydrological year and plot
total_annual_area = glacier_area_df_fill.agg("sum", axis="index")
plt.figure()
plt.plot(all_years, total_annual_area, '.-', linewidth=4, markersize=16)
plt.xlabel("Hydrological year", fontsize=16)
plt.ylabel("Total area (km$^2$)", fontsize=16)
plt.title("Estimated total NW Greenland glacier area", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()

# AGU 2019!!
# Plot change in area per year
change_annual_area = np.diff(total_annual_area)
plt.figure(figsize=(7, 7.5))
plt.xlabel("Hydrological year", fontsize=LABELSIZE)
plt.bar(all_years[6:], change_annual_area[5:], color='cornflowerblue')
plt.ylabel("$\Delta$Area (km$^2$)", fontsize=LABELSIZE)
plt.title("Total Annual $\Delta$Area", fontsize=TITLESIZE)
plt.xticks(fontsize=TICKLABELSIZE)
plt.yticks(fontsize=TICKLABELSIZE)
plt.setp(plt.gca().get_xticklabels()[::2], visible=False)
plt.grid(color='lightgray')
plt.savefig("/home/teblack/plots/AGU_area_change.png", transparent=True, bbox_inches='tight')

# AGU 2019!
# Total area and annual change in area - subplots
fig = plt.figure(figsize=(10,4.5))
#plt.suptitle("Total Glacier Area Changes", fontsize=TITLESIZE)
ax1 = fig.add_subplot(121)
p1 = plt.plot(all_years[5:], total_annual_area[5:] - total_annual_area.iloc[5], '.-', color='cornflowerblue', linewidth=LINEWIDTH-1, markersize=MARKERSIZE-2)
plt.ylabel("$\Delta$Area (km$^2$)", fontsize=LABELSIZE)
plt.xlabel("Hydrological year", fontsize=LABELSIZE)
plt.title("Cumulative Area Change", fontsize=TITLESIZE)
plt.xticks(fontsize=TICKLABELSIZE)
plt.yticks(fontsize=TICKLABELSIZE)
#plt.setp(plt.gca().get_xticklabels()[::2], visible=False)
plt.grid(color='lightgray')

ax2 = fig.add_subplot(122)
p2 = plt.bar(all_years[6:], change_annual_area[5:], color='cornflowerblue')
#plt.ylabel("$\Delta$Area (km$^2$)", fontsize=LABELSIZE)
plt.xlabel("Hydrological year", fontsize=LABELSIZE)
plt.title("Annual Area Change", fontsize=TITLESIZE)
plt.xticks(fontsize=TICKLABELSIZE)
plt.yticks(fontsize=TICKLABELSIZE)
#plt.setp(plt.gca().get_xticklabels()[::2], visible=False)
plt.grid(color='lightgray')
plt.savefig("/home/teblack/plots/AGU_area_changes.png", transparent=True, bbox_inches='tight')

# Plot summary of glacier observations over time
#   - Bar plot of number of glaciers observed per year
#   - Scatter plot of observations for each glacier (blue=obs, gray=interp)
annual_counts = glacier_area_df.count(axis='index')
plt.figure()
plt.bar(all_years, annual_counts, color='cornflowerblue')
plt.xlabel('Hydrological year', fontsize=16)
plt.ylabel('Number of glaciers observed', fontsize=16)
plt.title('Terminus observations, 1972-2018', fontsize=20)

# AGU 2019!!
plt.figure(figsize=(4.5,15))
for g in gIDs:
    years = unzipGlacierRecords(glacier_records[g])[3]
    no_interp_years = glacier_area_df_fill.loc[g].index[glacier_area_df_fill.loc[g].isna()]
    interp_years = list(set(all_years) - set(years) - set(no_interp_years))
    h1 = plt.scatter(years, [g]*len(years), color='cornflowerblue', s=12)
    h2 = plt.scatter(interp_years, [g]*len(interp_years), edgecolors='cornflowerblue', facecolors='none', s=12)
plt.xlabel('Hydrological year', fontsize=LABELSIZE)
plt.ylabel('Glacier ID', fontsize=LABELSIZE)
plt.title('Terminus observations', fontsize=TITLESIZE)
plt.xticks(fontsize=TICKLABELSIZE)
plt.yticks(fontsize=TICKLABELSIZE)
ax = plt.gca()
#ax.legend((h1, h2), ('Observed', 'Interpolated'), fontsize=LEGENDSIZE, loc='left', bbox_to_anchor=(1, 0.5))
plt.grid(color='lightgray')
plt.savefig("/home/teblack/plots/AGU_terminus_observations.png", transparent=True, bbox_inches='tight')

# Plot normalized area for all glaciers (from observations only)
normalized_area_df = pd.DataFrame(index=gIDs, columns=all_years)
plt.figure()
for g in gIDs:
    _, areas, _, hydroyears = unzipGlacierRecords(glacier_records[g])
    normalized_area = normGlacierArea(areas)
    for h in hydroyears:
        normalized_area_df.at[g, h] = normalized_area[hydroyears.index(h)]
    plt.plot(hydroyears, normalized_area, '.-', color='silver')
mean_normalized_area = normalized_area_df.agg("mean", axis="index")
plt.plot(normalized_area_df.keys(), mean_normalized_area)
plt.xlabel('Hydrological year')
plt.ylabel('Normalized glacier area')
plt.title('Normalized glacier areas (observations only)')

# Plot normalized area for all glaciers (from interpolated observations)
normalized_area_interp_df = pd.DataFrame(index=gIDs, columns=all_years)
plt.figure()
for g in gIDs:
    areas = glacier_area_df_fill.loc[g].values
    normalized_area = normGlacierArea(areas)
    normalized_area_interp_df.loc[g] = normalized_area
    plt.plot(all_years, normalized_area, color='silver')
mean_normalized_area = normalized_area_interp_df.agg("mean", axis="index")
plt.plot(normalized_area_interp_df.keys(), mean_normalized_area,
         '.-', linewidth=LINEWIDTH, markersize=MARKERSIZE)
plt.xlabel('Hydrological year', fontsize=LABELSIZE)
plt.ylabel('Fractional glacier area', fontsize=LABELSIZE)
plt.title('Fractional Glacier Areas (With Interpolated Areas)', fontsize=TITLESIZE)
plt.xticks(fontsize=TICKLABELSIZE)
plt.yticks(fontsize=TICKLABELSIZE)
plt.grid(color='lightgray')

# Dig into timing of glacier observations each year
observation_dates = pd.DataFrame(index=gIDs, columns=all_years)
for g in gIDs:
    dates, _, _, hyears = unzipGlacierRecords(glacier_records[g])
    for idx in range(len(dates)):
        observation_dates.loc[g, hyears[idx]] = dates[idx]

observation_dohy = pd.DataFrame(index=gIDs, columns=all_years)
for year in all_years:
    observation_dohy[year] = observation_dates[year].apply(dayOfHydroyear)

observation_dohy_count = pd.DataFrame(0, index=all_years, columns=range(366))
for year in all_years:
    dohy_counts = observation_dohy[year].value_counts()
    observation_dohy_count.loc[year] = dohy_counts
    observation_dohy_count.fillna(value=0)

observation_season = pd.DataFrame(index=gIDs, columns=all_years)
for year in all_years:
    observation_season[year] = observation_dates[year].apply(get_season)

observation_season_counts = pd.DataFrame(index=all_years, columns=[1, 2, 3, 4])
for year in all_years:
    seasonal_counts = observation_season[year].value_counts()
    observation_season_counts.loc[year] = seasonal_counts
observation_season_counts = observation_season_counts.rename(columns=season_dict)

# Scatter plot of day of hydroyear of observations vs year
plt.figure()
ax = plt.gca()
summer_start_dohy = dayOfHydroyear('%s-06-01' % all_years[0])
summer_end_dohy = dayOfHydroyear('%s-08-31' % all_years[0])
winter_start_dohy = dayOfHydroyear('%s-12-01' % all_years[0])
winter_end_dohy = dayOfHydroyear('%s-02-28' % all_years[0])
ax.add_patch(Rectangle((all_years[0], summer_start_dohy),
                       len(all_years), summer_end_dohy - summer_start_dohy,
                       color='red', alpha=0.1, ec='none'))
ax.add_patch(Rectangle((all_years[0], winter_start_dohy),
                       len(all_years), winter_end_dohy - winter_start_dohy,
                       color='blue', alpha=0.1, ec='none'))
for year in all_years:
    dohys = observation_dohy[year]
    plt.scatter([year]*len(dohys), dohys, color='cornflowerblue')
plt.xlabel('Hydrological year')
plt.ylabel('Days since September 1')
plt.title('Timing of observations')

# Scatter plot of observation day of hydroyear, scaled by #observations
plt.figure()
#ax = plt.gca()
#summer_start_dohy = dayOfHydroyear('%s-06-01' % all_years[0])
#summer_end_dohy = dayOfHydroyear('%s-08-31' % all_years[0])
#winter_start_dohy = dayOfHydroyear('%s-12-01' % all_years[0])
#winter_end_dohy = dayOfHydroyear('%s-02-28' % all_years[0])
#ax.add_patch(Rectangle((all_years[0], summer_start_dohy),
#                       len(all_years), summer_end_dohy - summer_start_dohy,
#                       color='red', alpha=0.1, ec='none'))
#ax.add_patch(Rectangle((all_years[0], winter_start_dohy),
#                       len(all_years), winter_end_dohy - winter_start_dohy,
#                       color='blue', alpha=0.1, ec='none'))
cmap = cm.get_cmap('cool')
for year in all_years:
    obs = observation_dohy_count.loc[year].dropna().sort_values(ascending=False)
    dohys = obs.keys()
    counts = obs.values
    colors = cmap(np.array(counts.tolist())/92)
    plt.scatter([year]*len(dohys), dohys, s=(counts**1.5).tolist(), color=colors)
plt.xlabel('Hydrological year', fontsize=16)
plt.ylabel('Days since September 1', fontsize=16)
plt.title('Timing of observations, scaled to # observations on date', fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()

# Heat map of observation day of hydroyear
plt.figure()


# Stacked bar plot of number observations in each season per year
observation_season_counts.plot.bar(stacked=True, colormap='viridis', alpha=0.7)
plt.xlabel('Hydrological year', fontsize=16)
plt.ylabel('Number of observations', fontsize=16)
plt.title('Seasonal distribution of terminus observations', fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

# Box-whisker plot
observation_dohy.plot.box()
plt.xlabel('Hydrological year')
plt.ylabel('Days since September 1')
plt.title('Spread of observations')
plt.xticks(rotation=90)

# %% RECREATE OLD MATLAB SCRIPT THAT LOOKS AT ADV/RETR AND SIGNIFICANCE
# AGU 2019!
glacier_stats = pd.DataFrame(index=gIDs, columns={'Lost', 'Significant'})
#np.zeros((len(gIDs), 1), dtype=bool)
#sig = np.zeros((len(gIDs), 1), dtype=bool)
for g in gIDs:
    stdev = np.nanstd(glacier_area_df.loc[g])
    data = glacier_area_df.loc[g]
    data = data[~np.isnan(data)]
    change = data.iloc[0] - data.iloc[-1]
    # Did glacier net lose area? (start area is larger than end area)
    if change > 0:
        glacier_stats.loc[g, 'Lost'] = True
    else:
        glacier_stats.loc[g, 'Lost'] = False
    # Was the loss significant? (is the net change bigger than the stdev)
    if abs(change) > stdev:
        glacier_stats.loc[g, 'Significant'] = True
    else:
        glacier_stats.loc[g, 'Significant'] = False

# %% Plot fractional area change, and average FOR EACH GLACIER GROUP
# AGU 2019!
# Get fractional glacier areas for each zone separately
norm_areas_zone1 = normalized_area_interp_df.loc[1:24]
norm_areas_zone2 = normalized_area_interp_df.loc[25:46]
norm_areas_zone3 = normalized_area_interp_df.loc[47:69]
norm_areas_zone4 = normalized_area_interp_df.loc[70:92]

# Calculate annual mean for everything and for individual glacier zones
mean_normalized_area = normalized_area_interp_df.agg("mean", axis="index")
mean_zone1 = norm_areas_zone1.agg("mean", axis="index")
mean_zone2 = norm_areas_zone2.agg("mean", axis="index")
mean_zone3 = norm_areas_zone3.agg("mean", axis="index")
mean_zone4 = norm_areas_zone4.agg("mean", axis="index")

COLORMAP = cm.Blues
color = iter(COLORMAP(np.linspace(0.2,1,4)))

fig = plt.figure(figsize=(10,10))
# Plot individual glaciers
for g in gIDs:
    pi1 = plt.plot(normalized_area_interp_df.keys(),
                   normalized_area_interp_df.loc[g,:].values,
                   color='silver', linewidth=0.7)
plt.setp(pi1, label='individual glaciers')
# Plot zone and total annual averages
pm5 = plt.plot(mean_normalized_area.keys(), mean_normalized_area, '.-',
               color='darkorange', linewidth=LINEWIDTH*3, markersize=MARKERSIZE*0,
               label='Total mean')
gcolor = next(color)
pm1 = plt.plot(mean_zone1.keys(), mean_zone1, '.-', color=gcolor, alpha=1,
               linewidth=LINEWIDTH, markersize=MARKERSIZE, label='Zone 1 mean')
gcolor = next(color)
pm2 = plt.plot(mean_zone2.keys(), mean_zone2, '.-', color=gcolor, alpha=1,
               linewidth=LINEWIDTH, markersize=MARKERSIZE, label='Zone 2 mean')
gcolor = next(color)
pm3 = plt.plot(mean_zone3.keys(), mean_zone3, '.-', color=gcolor, alpha=1,
               linewidth=LINEWIDTH, markersize=MARKERSIZE, label='Zone 3 mean')
gcolor = next(color)
pm4 = plt.plot(mean_zone4.keys(), mean_zone4, '.-', color=gcolor, alpha=1,
               linewidth=LINEWIDTH, markersize=MARKERSIZE, label='Zone 4 mean')

plt.xlabel('Hydrological Year', fontsize=LABELSIZE)
plt.ylabel('Fractional area', fontsize=LABELSIZE)
plt.title('Fractional Glacier Areas', fontsize=TITLESIZE)
plt.xticks(fontsize=TICKLABELSIZE)
plt.yticks(fontsize=TICKLABELSIZE)
plt.grid(color='lightgray')
plt.legend()
plt.savefig("/home/teblack/plots/AGU_frac_areas.png", transparent=True, bbox_inches='tight')
