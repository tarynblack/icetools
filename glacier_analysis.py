#!/usr/bin/env python3
# For Kenai Fjords National Park - Coastal Glaciers project (2020)
# Import glacier geospatial information, calculate changes, and plot
# Taryn Black, August 2020

# %% Import modules
import geopandas as gpd
import pandas as pd
import manage
import glaciermetrics as gm
import matplotlib.pyplot as plt
import gplots as gp

# %% User-specified parameters
# Path for file geodatabase containing glacier information
fgdb = '/mnt/d/KEFJ_CoastalGlaciers/CoastalGlaciers.gdb'

# Path to save analysis output
outpath = '/mnt/d/KEFJ_CoastalGlaciers/analysis/'

# Path/filename to save area/length spreadsheet
data_spreadsheet = '/mnt/d/KEFJ_CoastalGlaciers/analysis/\
    CoastalGlaciers_data.xlsx'

# File geodatabase layer names
outlines_layer = 'Glacier_Outlines'
reflines_layer = 'Glacier_Reference_Lines'
points_layer = 'Glacier_Points'
centerlines_layer = 'Glacier_Centerlines'
boxes_layer = 'Glacier_Boxes'

# %% Load glacier layer information
outlines = gpd.read_file(fgdb, layer=outlines_layer, driver='FileGDB')

reflines = gpd.read_file(fgdb, layer=reflines_layer, driver='FileGDB')
reflines.sort_values(by='Glacier_ID', inplace=True)
reflines.set_index('Glacier_ID', drop=False, inplace=True)

points = gpd.read_file(fgdb, layer=points_layer, driver='FileGDB')
points.sort_values(by='Glacier_ID', inplace=True)
points.set_index('Glacier_ID', drop=False, inplace=True)

centerlines = gpd.read_file(fgdb, layer=centerlines_layer, driver='FileGDB')
centerlines.sort_values(by='Glacier_ID', inplace=True)
centerlines.set_index('Glacier_ID', drop=False, inplace=True)

refboxes = gpd.read_file(fgdb, layer=boxes_layer, driver='FileGDB')
refboxes.sort_values(by='Glacier_ID', inplace=True)
refboxes.set_index('Glacier_ID', drop=False, inplace=True)

# %% Get scope of glaciers and time
GIDS = points.Glacier_ID.values

outlines['Year'] = pd.to_datetime(outlines['Year_'], format='%Y')
del outlines['Year_']
YEAR_START = outlines.Year.min().year
YEAR_END = outlines.Year.max().year
YEARS = range(YEAR_START, YEAR_END+1)

DATE_START = outlines.Source_Date.min()
DATE_END = outlines.Source_Date.max()

START_DECADE = YEAR_START - YEAR_START%10
END_DECADE = YEAR_END + (10 - YEAR_END%10)
RANGE_DECADES = range(START_DECADE, END_DECADE, 10)
DECADES = [pd.to_datetime('{}-01-01'.format(y)).date() for y in RANGE_DECADES]

# %% Initialize dictionary of Glacier objects to store glacier information
all_glaciers = {id: manage.Glacier(id) for id in GIDS}
for id in all_glaciers:
    all_glaciers[id].refline = reflines.loc[id].geometry
    all_glaciers[id].refbox = refboxes.loc[id].geometry
    if centerlines.loc[id].ndim > 1:
        all_glaciers[id].centerline = centerlines.loc[id].geometry.iloc[0]
    else:
        all_glaciers[id].centerline = centerlines.loc[id].geometry
    all_glaciers[id].officialname = points.loc[id].Official_Name
    all_glaciers[id].unofficialname = points.loc[id].Unofficial_Name
    all_glaciers[id].fjordname = points.loc[id].Fjord_Name
    all_glaciers[id].armname = points.loc[id].Arm_Name

# %% Construct an observation time series for each glacier
for id in GIDS:
    
    # Get reference line and all observations for glacier ID
    glacier = outlines.query('Glacier_ID == @id')
    refline = all_glaciers[id].refline
    refbox  = all_glaciers[id].refbox
    cenline = all_glaciers[id].centerline
    
    # Loop through all observations and process data
    for n in range(len(glacier)):
        observation = glacier.iloc[n]

        # Create a Terminus Observation for a row in geodataframe
        obs = manage.TerminusObservation(gid=observation.Glacier_ID,
                                         qflag=observation.Quality_Flag,
                                         termination=observation.Termination_Type,
                                         imageid=observation.Image_ID,
                                         sensor=observation.Sensor,
                                         date=observation.Source_Date,
                                         circadiandate=observation.Circadian_Date,
                                         year=observation.Year,
                                         season=observation.Season,
                                         geometry=observation.geometry)
        
        # Calculate glacier area against reference line
        obs.area = gm.glacierArea(obs, refline)

        # Calculate glacier terminus area (excluding sides) against ref box
        obs.termarea = gm.glacierArea(obs, refbox)

        # Find outline intersection with centerline and resulting length
        inx_point, sublength = gm.centerlineIntersection(obs, cenline)
        obs.centerlineintersection = inx_point
        obs.length = sublength
        
        # Add glacier observation to time series
        all_glaciers[id].add_observation(obs)
    
    # Ensure that all observations are sorted by date
    all_glaciers[id].sort_by_date()

    # Extract time series of area, length, and dates from all observations
    all_glaciers[id].areas = all_glaciers[id].extract('area')
    all_glaciers[id].termareas = all_glaciers[id].extract('termarea')
    all_glaciers[id].lengths = all_glaciers[id].extract('length')
    all_glaciers[id].dates = all_glaciers[id].extract('date')

# %% For each glacier, create plots and calculate other metrics

for id in GIDS:
    glacier = all_glaciers[id]
    print('Analyzing glacier #{}: {}'.format(id, gp.getGlacierName(glacier)))

    # Plot relative area over time
    fig = plt.figure(figsize=(10, 4.5))
    gp.totalRelativeMeasure(fig, glacier, 'area')
    plt.savefig('{}{:02}_relativearea.png'.format(outpath, id), \
        bbox_inches='tight')
    
    # Plot relative length over time
    fig = plt.figure(figsize=(10, 4.5))
    gp.totalRelativeMeasure(fig, glacier, 'length')
    plt.savefig('{}{:02}_relativelength.png'.format(outpath, id), \
        bbox_inches='tight')
    
    # Plot relative area (with/without sides) and length over time, together
    fig = plt.figure(figsize=(10, 4.5))
    gp.totalRelativeMeasureCompare(fig, glacier)
    plt.savefig('{}{:02}_sizecompare.png'.format(outpath, id), \
        bbox_inches='tight')

    # Plot spring and fall area separately
    fig = plt.figure(figsize=(10, 4.5))
    gp.seasonRelativeMeasure(fig, glacier, 'area', spring=True, autumn=True)
    plt.savefig('{}{:02}_seasonarea.png'.format(outpath, id), \
        bbox_inches='tight')
    
    # Plot spring and fall length separately
    fig = plt.figure(figsize=(10, 4.5))
    gp.seasonRelativeMeasure(fig, glacier, 'length', spring=True, autumn=True)
    plt.savefig('{}{:02}_seasonlength.png'.format(outpath, id), \
        bbox_inches='tight')

    # # Plot seasonal area change between each measurement
    # fig = plt.figure(figsize=(10, 4.5))
    # gp.seasonMeasureChange(fig, glacier, 'area', spring=True, autumn=True)
    # plt.savefig('{}{:02}_seasonareachange.png'.format(outpath, id), \
    #     bbox_inches='tight')

    # # Plot seasonal length change between each measurement
    # fig = plt.figure(figsize=(10, 4.5))
    # gp.seasonMeasureChange(fig, glacier, 'length', spring=True, autumn=True)
    # plt.savefig('{}{:02}_seasonlengthchange.png'.format(outpath, id), \
    #     bbox_inches='tight')

    # Plot decadal net area change
    fig = plt.figure(figsize=(10, 4.5))
    gp.decadalMeasureChange(fig, glacier, 'area', DECADES)
    plt.savefig('{}{:02}_decadalareachange.png'.format(outpath, id), \
        bbox_inches='tight')

    # Plot decadal net length change
    fig = plt.figure(figsize=(10, 4.5))
    gp.decadalMeasureChange(fig, glacier, 'length', DECADES)
    plt.savefig('{}{:02}_decadallengthchange.png'.format(outpath, id), \
        bbox_inches='tight')

    # # Plot annual area change for past two decades
    # fig = plt.figure(figsize=(10, 4.5))
    # gp.individualMeasureChange(fig, glacier, 'area', date_start='2000-01-01')
    # plt.savefig('{}{:02}_indivareachange.png'.format(outpath, id), \
    #     bbox_inches='tight')

    # # Plot annual length change for past two decades
    # fig = plt.figure(figsize=(10, 4.5))
    # gp.individualMeasureChange(fig, glacier, 'length', date_start='2000-01-01')
    # plt.savefig('{}{:02}_indivlengthchange.png'.format(outpath, id), \
    #     bbox_inches='tight')

    plt.close('all')
    
    # Calculate output metrics and save to file         
    with open('{}{:02}_metrics.txt'.format(outpath, id), 'w+') as f:
        name = gp.getGlacierName(glacier)
        f.write('ID #{}: {}\n'.format(id, name))
        f.write('{} total observations\n'.format(len(glacier.extract('area'))))
        f.write('Analyzed on: {}\n\n'.format(pd.Timestamp.utcnow().round('s')))

        # Calculate overall net area change and average change rate
        cumul_areachange, change_dates, num_obs = gm.netMeasureChange(
            glacier, 'area')
        areachange_rate_total = gm.rateMeasureChange(glacier, 'area')

        cumul_lenchange, _, _ = gm.netMeasureChange(glacier, 'length')
        lenchange_rate_total = gm.rateMeasureChange(glacier, 'length')

        f.write('Observed date range: {} to {}\n'.format(
            change_dates[0], change_dates[1]))
        f.write('Number of observations: {}\n'.format(num_obs))
        f.write('Net area change: {:.3f} km2\n'.format(
            cumul_areachange.iloc[-1]))
        f.write('Net length change: {:.3f} km\n'.format(
            cumul_lenchange.iloc[-1]))
        f.write('Rate of area change: {:.3f} km2/yr\n'.format(
            areachange_rate_total))
        f.write('Rate of length change: {:.3f} km/yr\n\n'.format(
            lenchange_rate_total))

        # Calculate net area change and average change rate for each decade
        for start_year in DECADES:
            end_year = gm.addDecade(start_year)
            cumul_areachange, change_dates, num_obs = gm.netMeasureChange(
                glacier, 'area', start_year, end_year)
            areachange_rate = gm.rateMeasureChange(
                glacier, 'area', start_year, end_year)

            cumul_lenchange, _, _ = gm.netMeasureChange(
                glacier, 'length', start_year, end_year)
            lenchange_rate = gm.rateMeasureChange(
                glacier, 'length', start_year, end_year)
            
            # Write outputs to text file
            f.write('Subset date range: {} to {}\n'.format(
                start_year, end_year))
            f.write('Observed date range: {} to {}\n'.format(
                change_dates[0], change_dates[1]))
            f.write('Number of observations: {}\n'.format(num_obs))
            f.write('Net area change: {:.3f} km2\n'.format(
                cumul_areachange.iloc[-1]))
            f.write('Net length change: {:.3f} km\n'.format(
                cumul_lenchange.iloc[-1]))
            f.write('Rate of area change: {:.3f} km2/yr\n'.format(
                areachange_rate))
            f.write('Rate of length change: {:.3f} km/yr\n\n'.format(
                lenchange_rate))

# %% Create spreadsheet of glacier length/area data
with pd.ExcelWriter(data_spreadsheet, date_format='YYYY-MM-DD') as writer:
    for id in GIDS:
        glacier = all_glaciers[id]
        ss_data = {'area (km2)': glacier.areas,
                   'terminus area (km2)': glacier.termareas,
                   'length (km)': glacier.lengths,
                   'date': glacier.dates}
        ss_df = pd.DataFrame(ss_data)
        ss_df.to_excel(writer, \
            sheet_name='{} {}'.format(id, gp.getGlacierName(glacier)))

# %% Create summary plots of all glaciers

# Plot all observations per glacier over time (scatter plot)
fig = plt.figure(figsize=(10, 4.5))
gp.glacierObservations(fig, all_glaciers)
plt.savefig('{}observation_timeseries.png'.format(outpath), \
    bbox_inches='tight')

# Plot area change for all glaciers
fig = plt.figure(figsize=(10, 4.5))
gp.measureSummary(fig, all_glaciers, 'area')
plt.savefig('{}summary_areachange.png'.format(outpath), \
    bbox_inches='tight')

# Plot length change for all glaciers
fig = plt.figure(figsize=(10, 4.5))
gp.measureSummary(fig, all_glaciers, 'length')
plt.savefig('{}summary_lengthchange.png'.format(outpath), \
    bbox_inches='tight')

# Plot normalized area change for all glaciers
fig = plt.figure(figsize=(10, 4.5))

plt.savefig('{}summary_normareachange.png'.format(outpath), \
    bbox_inches='tight')

# Plot normalized length change for all glaciers
fig = plt.figure(figsize=(10, 4.5))

plt.savefig('{}summary_normlengthchange.png'.format(outpath), \
    bbox_inches='tight')

print('Done.')