# Taryn Black, November 2019
# Work with sea ice data
# initial data from https://nsidc.org/data/g02202

# %% PROJECT IDEAS
# TODO: get length of ice-free season and icy season (plot)
# TODO: write neighbor pixel section to get individual locations rather than slice from indices
# TODO: check for annual vs hydrological year on plots!!

# Import modules
import os
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import scipy.ndimage
from scipy.spatial.distance import cdist
import seaborn as sns

# %% USER-DEFINED PARAMETERS ==================================================
# Path for location of data
datadir = "/Volumes/insar5/teblack/data/seaice/cdr-monthly/"

datavar = "goddard_nt_seaice_conc_monthly"

# Define region boundaries (lat/lon)
bounds = {'N': 82,
          'S': 68,
          'E': -40,
          'W': -76}

# Pixel distance offshore to filter data (parameters for offshore_pixels())
offshore_distance = 1
outer_distance = 2  # optional, must be greater than offshore_distance
offshore_type = 'ring'  # optional, 'adjacent' or 'ring'

# Median center geographic coordinates of glacier groups (grouped by ~23), to
# use for finding nearby grid cells.
centers = [[70.8, -50.8],
           [74.6, -56.5],
           [76.3, -61.8],
           [77.4, -68.1]]

# %% FUNCTIONS ================================================================


def reduce_datavar(ds, datavar):
    """Reduce full dataset to a specific data variable, for convenience.
    Possible variables for this dataset are:
        projection
        seaice_conc_monthly_cdr
        stdev_of_seaice_conc_monthly_cdr
        melt_onset_day_seaice_conc_monthly_cdr
        qa_of_seaice_conc_monthly_cdr
        goddard_merged_seaice_conc_monthly
        goddard_nt_seaice_conc_monthly
        goddard_bt_seaice_conc_monthly
    See https://nsidc.org/data/g02202#title16 for variable descriptions."""
    ds_datavar = ds[datavar]
    return ds_datavar


def geographic_filter(ds, bounds):
    """Filter full dataset down to a specific region defined by lat/lon bounds.
    'bounds' is a dictionary where the keys are cardinal directions (NESW) and
    the values are the boundary edges in degrees."""
    ds_filtered = ds.where((ds.latitude >= bounds['S']) &
                           (ds.latitude <= bounds['N']) &
                           (ds.longitude >= bounds['W']) &
                           (ds.longitude <= bounds['E']), drop=True)
    return ds_filtered


def get_hydroyear(date):
    """From a date, determine hydrological year. Hydrological year in Greenland
    is defined as September 1 through August 31 (Ettema et al 2009)."""
    pdate = pd.to_datetime(date)
    if pdate.month >= 9:
        hydroyear = pdate.year
    elif pdate.month < 9:
        hydroyear = pdate.year - 1
    return hydroyear


def surface_mask(ds, surface_type):
    """Creates a binary mask for a surface type, where the input surface type
    is True, and all other surfaces are False. Also, mask of 1 (T) and nan (F).
    Possible surface types are: 'pole_hole', 'unused', 'coastal', 'land_mask',
    and 'missing_data'."""
    surface_dict = {'pole_hole': -5,
                    'unused': -4,
                    'coastal': -3,
                    'land_mask': -2,
                    'missing_data': -1,
                    'water': 0}
    surface_flag = surface_dict[surface_type]
    if surface_flag < 0:
        mask = ds.where(ds == surface_flag, other=0).values.astype(np.bool)
    else:
        "Sea ice concentration of zero still needs to be flagged as ocean."""
        mask = ds.where(ds >= 0).values
        mask += 1
        mask[np.isnan(mask)] = 0
        mask = mask.astype(bool)
    if np.ndim(mask) > 2:
        mask = mask.squeeze()
    return mask


def nan_surface_mask(ds, surface_type):
    """Create a binary surface type mask using surface_mask(), and convert to a
    boolean mask to a mask of ones (True) and nans (False). This is for
    plotting convenience."""
    mask = surface_mask(ds, surface_type)
    mask_nan = (1.0 + np.zeros_like(mask))*mask
    mask_nan[mask_nan == 0] = np.nan
    return mask_nan


def neighbor_pixels(mask, distance):
    """Find pixels neighboring the True pixels in a boolean mask, to a given
    distance of pixels away. Output is a binary mask of neighbor and mask
    pixels as True, and all other pixels as False."""
    neighbors = scipy.ndimage.binary_dilation(mask, iterations=distance)
    return neighbors


def offshore_pixels(ds, distance, offshore_type="adjacent", outer_distance=2):
    """Finds ocean pixels within some specified distance from the coast, and
    returns their sea ice concentration values and grid indices.
    If offshore_type is set to "ring", function will return ocean pixels up to
    a given outer_distance away from the coast and excluding the pixels within
    the nearer distance, forming an offshore ring of pixels.
    Function variables are:
        ds : input dataset (see surface_mask() )
        distance : number of pixels away from coast to capture
        offshore_type : can be "adjacent" or "ring"
        outer_distance : number of pixels away from the coast to capture, less
        the pixels within 'distance', if offshore_type is 'ring'"""
    coast_mask = surface_mask(ds, 'coastal')
    coast_neighbors = neighbor_pixels(coast_mask, distance)
    nearshore_vals = (ds.where(coast_neighbors == True).where(ds >= 0).values)
    nearshore_vals = nearshore_vals.squeeze()
    nearshore_idx = list(zip(np.where(np.isfinite(nearshore_vals))[1],
                         np.where(np.isfinite(nearshore_vals))[0]))
    if offshore_type == "adjacent":
        values = nearshore_vals
        indices = nearshore_idx
    elif offshore_type == "ring":
        if not outer_distance > distance:
            print("Hey! outer_distance must be greater than distance\n"
                  "Calculating as type 'adjacent' instead.")
            values = nearshore_vals
            indices = nearshore_idx
        else:
            outer_pixels = neighbor_pixels(coast_mask, outer_distance)
            outer_vals = (ds.where(outer_pixels == True).where(ds >= 0).values)
            values = outer_vals.squeeze()
            outer_idx = list(zip(np.where(np.isfinite(outer_vals))[1],
                                 np.where(np.isfinite(outer_vals))[0]))
            indices = list(set(outer_idx) - set(nearshore_idx))
    return values, indices


def plot_background(ds):
    """Plot masks of land, coast, and ocean pixels from a dataset ds to go
    behind plotted data. Uses nan_surface_mask()"""
    # TODO: better projection
    land = nan_surface_mask(ds, 'land_mask')
    ocean = nan_surface_mask(ds, 'water')
    coast = nan_surface_mask(ds, 'coastal')
    plt.imshow(land, cmap='gray', alpha=0.3)
    plt.imshow(coast, cmap='gray')
    plt.imshow(ocean, cmap='Wistia', alpha=0.2)


def plot_pixel_timeseries(data, x, y):
    """Plot timeseries of sea ice concentration from a single pixel, where x
    and y are pixel grid coordinates (xgrid and ygrid)."""
    px_data = data.isel(xgrid=x, ygrid=y)
    time = px_data.time.values
    vals = px_data.values
    plt.plot(time, vals, '.-')
#    data.isel(xgrid=x, ygrid=y).plot()


def plot_pixel_location(ds, x, y):
    """Plot marker at the location of a single pixel, on a mask background."""
    # TODO: tick labels as lat/lon rather than grid coords?
    plot_background(da_nwg)
    plt.scatter(x, y, marker='.', s=200, c='cyan', edgecolors='gray')


def plot_pixel_info(data, x, y, fig, gs, row, ncols):
    """Plot pixel timeseries and location on a modified grid.
    Can accept data as a list to plot multiple time series on one axes."""
    if not isinstance(data, list):
        data = [data]
    for d in data:
        ax1 = fig.add_subplot(gs[row, 0:ncols-1])
        plot_pixel_timeseries(d, x, y)
    ax1.set_title("Sea ice concentration at %.2f N, %.2f E" %
                  (data[0].isel(xgrid=x, ygrid=y).latitude,
                   data[0].isel(xgrid=x, ygrid=y).longitude))
    ax2 = fig.add_subplot(gs[row, ncols-1])
    plot_pixel_location(d.isel(time=0), x, y)
    return ax1, ax2


def nearest_to_latlon(data, coordinate, tolerance=0.1):
    """Find grid coordinate nearest to a given geographic coordinate
    [latitude, longitude], within a given tolerance."""
    gridpt = data.where((data.latitude > coordinate[0] - tolerance) &
                        (data.latitude < coordinate[0] + tolerance) &
                        (data.longitude > coordinate[1] - tolerance) &
                        (data.longitude < coordinate[1] + tolerance),
                        drop=True)
    gridpt = [gridpt.isel(xgrid=0, ygrid=0).xgrid.values,
              gridpt.isel(xgrid=0, ygrid=0).ygrid.values]
#    gridpt = gridpt.isel(xgrid=0, ygrid=0) # ensure only 1 coord returned
    return gridpt


def period_average(data, period, start='DEC', window=12, season='DJF'):
    """Calculate average sea ice concentration over a defined period.
    Options are:
        - 'calyear'  : annual average, grouped by calendar year
        - 'hydroyear': annual average, grouped by hydrological year (Sep-Aug)
        - 'annual'   : annual average, start month defined by 'start'
        - 'rolling'  : rolling average, window size defined by 'window'
        - 'seasonal' : annual averages for a season (DJF, MAM, JJA, SON)"""
    if period is 'calyear':
        avg_data = data.groupby('time.year').mean(dim='time', skipna=True)
        avg_data = avg_data.rename({'year': 'time'})
        # TODO: figure out how to deal with incomplete years
    elif period is 'hydroyear':
        avg_data = data.groupby('hydroyear').mean(dim='time', skipna=True)
    elif period is 'annual':
        avg_data = data.resample(time='AS-%s' % start).mean()
    elif period is 'rolling':
        avg_data = data.rolling(time=window).mean(dim='time', skipna=True)
    elif period is 'seasonal':
        seasons_avgs = data.resample(time="QS-DEC").mean()
        season_annavg = seasons_avgs.where(seasons_avgs['time.season'] ==
                                           season, drop=True)
        avg_data = season_annavg
    return avg_data


def ice_presence(data, minimum=15, start='SEP'):
    # TODO: write docstring
    ice_concentration = data.where(data >= minimum, other=0)
    ice_boolean = ice_concentration.where(ice_concentration == 0, other=1)
    ice_season_length = (ice_boolean.where(ice_boolean == 1).
                         resample(time='AS-%s' % start).count())
    noice_season_length = (ice_boolean.where(ice_boolean == 0).
                           resample(time='AS-%s' % start).count())
    return ice_concentration, ice_boolean, ice_season_length, noice_season_length


def fullyear_filter(data, start='SEP'):
    # TODO: write docstring
    month_dict = {'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6,
                  'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12}
    filtered_data = data
    while filtered_data[0]['time.month'] != month_dict[start]:
        filtered_data = filtered_data[1:]
    if start == 'JAN':
        while filtered_data[-1]['time.month'] != month_dict['DEC']:
            filtered_data = filtered_data[0:-1]
    else:
        while filtered_data[-1]['time.month'] != month_dict[start]-1:
            filtered_data = filtered_data[0:-1]
    return filtered_data


def point_ice_concentration_data(ice_conc_data, x, y, yearrange, start='SEP'):
    """...""" # TODO: write docstring
    shift = {'JAN': 0, 'FEB': 11, 'MAR': 10, 'APR': 9, 'MAY': 8, 'JUN': 7,
             'JUL': 6, 'AUG': 5, 'SEP': 4, 'OCT': 3, 'NOV': 2, 'DEC': 1}
    months_shift = np.roll(np.arange(1, 13), shift[start])
    ice_df = pd.DataFrame(index=months_shift, columns=yearrange)
    ice_conc_pt = ice_conc_data.isel(xgrid=x, ygrid=y)
    for k in ice_conc_pt:
        month = int(k['time.month'].values)
        hydroyear = int(k.hydroyear.values)
        ice_df.at[month, hydroyear] = np.ndarray.item(k.values.astype(float))
    icydata = np.asarray(ice_df, dtype=float)
    ice_conc_HY = pd.DataFrame(data=icydata,
                               index=months_shift, columns=yearrange)
    # TODO: make this a separate function
    ice_conc_annavg = ice_conc_HY.mean(axis=0, skipna=False)
    ice_conc_annavg = np.asarray(ice_conc_annavg).reshape(1, len(yearrange))
    ice_conc_df = pd.DataFrame(data=ice_conc_annavg,
                               index=np.arange(1), columns=yearrange)
    blank_row = pd.DataFrame(data=np.full((1, len(yearrange)), 100), 
                             index=np.arange(-999,-998), columns=yearrange)
    ice_conc_df = ice_conc_df.append(blank_row)
    ice_conc_df = ice_conc_df.append(ice_conc_HY)
    return ice_conc_df


def plot_concentration_heatmap(ice_conc_data, start='SEP'):
    """...""" # TODO: docstring
    months = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May',
              6: 'June', 7: 'July', 8: 'August', 9: 'September', 10: 'October',
              11: 'November', 12: 'December'}
    shift = {'JAN': 0, 'FEB': 11, 'MAR': 10, 'APR': 9, 'MAY': 8, 'JUN': 7,
             'JUL': 6, 'AUG': 5, 'SEP': 4, 'OCT': 3, 'NOV': 2, 'DEC': 1}
    month_lbl = [months[i] for i in np.roll(np.arange(1, 13), shift[start])]
    ylbls = ['Ann. avg.', '']
    ylbls.extend(month_lbl)
    
    # Annotations for annual average only
    annolbls = np.empty_like(ice_conc_data, dtype='U3')
    annavg = ice_conc_data.loc[0,:].replace(to_replace=np.nan, value=0)
    annavg = [str(round(i)) for i in annavg]
    annolbls[0,:] = annavg
    annolbls[0,:] = ['' if i=='0' else i for i in annolbls[0,:]]
    
#    plt.figure()
    sic = sns.heatmap(ice_conc_data, vmin=0, vmax=100, square=True,
                      cmap='Blues_r', linewidth=0.1, yticklabels=ylbls,
                      cbar_kws={'shrink': 0.4, 'pad': 0.03},
                      annot=annolbls, fmt='')
    sic.set_facecolor('lightgray')
    sic.set_yticklabels(sic.get_yticklabels(), rotation=0)
#    sic.set_xticklabels(sic.get_xticklabels(), rotation=45)
#    sic.set_xlabel('Hydrological year')
#    plt.title(("Monthly and annually averaged sea ice concentration \n"
#               "%.2f N, %.2f E") %
#              (data[0].isel(xgrid=x, ygrid=y).latitude,
#               data[0].isel(xgrid=x, ygrid=y).longitude))
    return sic


# %% LOAD DATA ================================================================

# List all data files in directory
data_files = [datadir+file for file in os.listdir(datadir) if file.endswith(".nc")]

# Concatenate all data into a single dataset
fdata = xr.open_dataset(data_files[0], mask_and_scale=False)
for file in data_files[1:]:
    ds = xr.open_dataset(file, mask_and_scale=False)
    fdata = xr.concat([fdata, ds], dim='time')
fdata = fdata.sortby('time')

# Add coordinate "hydrological year" for easier analysis later
fdata = fdata.assign_coords(hydroyear=('time', [get_hydroyear(d) for
                                                d in fdata.time.values]))

# Reduce and filter dataset to desired variable, location, etc.
rdata = reduce_datavar(fdata, datavar)
data = geographic_filter(rdata, bounds)
da_nwg = data.isel(time=0)

# %% DEFINE OTHER VARIABLES ===================================================
month_dict = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May',
              6: 'June', 7: 'July', 8: 'August', 9: 'September', 10: 'October',
              11: 'November', 12: 'December'}

startyear = data['time.year'].values.min()
endyear = data['time.year'].values.max()
yearrange = np.arange(startyear, endyear+1)

# %% GET PIXELS NEAR COAST
nearshore_vals, nearshore_idx = offshore_pixels(da_nwg, offshore_distance,
                                                offshore_type=offshore_type,
                                                outer_distance=outer_distance)

# !!! ADDITIONAL FILTERING SPECIFIC TO TEST COORDINATES -- REMOVING PIXELS NEAR BAFFIN ISLAND !!!
nearshore_idx = [coord for coord in nearshore_idx if coord[0] > 20]

# %% Plot background - land and coastline - plus sea ice concentration
plt.figure()
plot_background(da_nwg)
plt.imshow(nearshore_vals, cmap='Blues_r')
plt.clim(0, 100)
plt.colorbar().set_label("Sea ice concentration", size=12)
plt.title("Nearshore sea ice concentration")
# TODO: add date to title
# TODO: make this a function - for all area and just nearshore

plt.figure()
plot_background(da_nwg)
plt.scatter([xy[0] for xy in nearshore_idx], [xy[1] for xy in nearshore_idx])
plt.title("Nearshore pixel locations")

# %% Get grid coordinates for four handpicked geographic locations
gpts = [nearest_to_latlon(data, pt, tolerance=0.2) for pt in centers]
nearshore_gpts = [[data.isel(xgrid=x, ygrid=y).xgrid.values,
                   data.isel(xgrid=x, ygrid=y).ygrid.values]
                  for x, y in nearshore_idx]
distances = cdist(np.asarray(gpts), np.asarray(nearshore_gpts), 'euclidean')
dist_idx = np.argmin(distances, axis=1)
nearest_gpts = [nearshore_idx[i] for i in dist_idx]
nearest_gpts.sort(key=lambda coord: coord[0])
nearest_gpts_latlon = [(float(data.isel(xgrid=x, ygrid=y).latitude.values),
                        float(data.isel(xgrid=x, ygrid=y).longitude.values))
                       for x, y in nearest_gpts]

# %% Plot full time series for those four points
nrows = len(nearest_gpts)
ncols = 5
fig = plt.figure()
gs = fig.add_gridspec(nrows, ncols)
for n in range(0, nrows):
    row = n
    x = nearest_gpts[n][0]
    y = nearest_gpts[n][1]
    ax1, ax2 = plot_pixel_info(data, x, y, fig, gs, row, ncols)
    ax1.set_ylim(-5, 105)
    if row == nrows-1:
        plt.setp(ax1.get_xticklabels(), visible=True)
        ax1.set_xlabel("Year")
    else:
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_xlabel(None)
plt.suptitle("Sea ice concentration along Greenland coast, 1978-2018", size=18)

# %% AVERAGE CONCENTRATION PLOTS FOR THE FOUR SELECTED POINTS
nrows = len(nearest_gpts)
ncols = 5

# Average annual concentration
fig = plt.figure()
gs = fig.add_gridspec(nrows, ncols)
for n in range(0, nrows):
    row = n
    x = nearest_gpts[n][0]
    y = nearest_gpts[n][1]
    annavg = period_average(data, 'annual', start='SEP')
    ax1, ax2 = plot_pixel_info(annavg, x, y, fig, gs, row, ncols)
#    ax1.set_ylim(-5, 105)
    if row == nrows-1:
        plt.setp(ax1.get_xticklabels(), visible=True)
        ax1.set_xlabel("Year")
    else:
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_xlabel(None)
plt.suptitle("Sep-Aug avg sea ice concentration along Greenland coast, 1978-2018", size=18)

# Plot both summer and winter on top of each other
fig = plt.figure()
gs = fig.add_gridspec(nrows, ncols)
for n in range(0, nrows):
    row = n
    x = nearest_gpts[n][0]
    y = nearest_gpts[n][1]
    winter = period_average(data, 'seasonal', season='DJF')
    summer = period_average(data, 'seasonal', season='JJA')
    ax1, ax2 = plot_pixel_info([winter, summer], x, y, fig, gs, row, ncols)
    ax1.set_ylim(-5, 105)
    if row == nrows-1:
        plt.setp(ax1.get_xticklabels(), visible=True)
        ax1.set_xlabel("Year")
        ax1.legend('Winter', 'Summer')
    else:
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_xlabel(None)
plt.suptitle("DJF and JJA sea ice concentration along Greenland coast, 1978-2018", size=18)

# Plot winter and summer timeseries separately
fig = plt.figure()
gs = fig.add_gridspec(nrows, ncols)
for n in range(0, nrows):
    row = n
    x = nearest_gpts[n][0]
    y = nearest_gpts[n][1]
    winter = period_average(data, 'seasonal', season='DJF')
    ax1, ax2 = plot_pixel_info(winter, x, y, fig, gs, row, ncols)
#    ax1.set_ylim(-5, 105)
    if row == nrows-1:
        plt.setp(ax1.get_xticklabels(), visible=True)
        ax1.set_xlabel("Year")
    else:
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_xlabel(None)
plt.suptitle("DJF sea ice concentration along Greenland coast, 1978-2018", size=18)

fig = plt.figure()
gs = fig.add_gridspec(nrows, ncols)
for n in range(0, nrows):
    row = n
    x = nearest_gpts[n][0]
    y = nearest_gpts[n][1]
    summer = period_average(data, 'seasonal', season='JJA')
    ax1, ax2 = plot_pixel_info(summer, x, y, fig, gs, row, ncols)
#    ax1.set_ylim(-5, 105)
    if row == nrows-1:
        plt.setp(ax1.get_xticklabels(), visible=True)
        ax1.set_xlabel("Year")
    else:
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_xlabel(None)
plt.suptitle("JJA sea ice concentration along Greenland coast, 1978-2018", size=18)

# %% DURATION OF ICY/ICE-FREE SEASONS
ice_concentration, ice_boolean, ice_season_length, noice_season_length = ice_presence(data, start='SEP')

LINEWIDTH = 5
MARKERSIZE = 15
TITLESIZE = 24
LABELSIZE = 20
TICKLABELSIZE = 16
LEGENDSIZE = 16

n=0
# Plot heatmap of monthly ice concentration
for pt in nearest_gpts:
    x = pt[0]
    y = pt[1]
    df = point_ice_concentration_data(ice_concentration, x, y,
                                      yearrange, start='SEP')
    fig = plt.figure(figsize=(10,7))
    if y == min([pt[1] for pt in nearest_gpts]): # GREATEST latitude
        plt.suptitle("Monthly and Annually Averaged Sea Ice Concentration", size=TITLESIZE)
        plt.title(("%.2f N, %.2f E") % 
                  (data[0].isel(xgrid=x, ygrid=y).latitude,
                   data[0].isel(xgrid=x, ygrid=y).longitude),
                   size=TITLESIZE*0.8)
    else:
        plt.title(("%.2f N, %.2f E") % 
                  (data[0].isel(xgrid=x, ygrid=y).latitude,
                   data[0].isel(xgrid=x, ygrid=y).longitude), size=TITLESIZE*0.8)
    sic = plot_concentration_heatmap(df, start='SEP')
    ax = plt.gca()
    cbar = sic.collections[0].colorbar
    cbar.ax.tick_params(labelsize=TICKLABELSIZE)
    plt.tick_params(axis='both', labelsize=TICKLABELSIZE)
    plt.setp(ax.get_yticklabels()[1::2], visible=False)
    if y == max([pt[1] for pt in nearest_gpts]): # LEAST latitude
#        plt.setp(ax.get_xticklabels(), visible=False)
        [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 5 != 0]
        plt.xlabel("Hydrological Year", size=LABELSIZE)
    else:
        plt.setp(ax.get_xticklabels(), visible=False)
    n += 1
    plt.savefig(("/home/teblack/plots/AGU_SIC%s.png" % n), transparent=False, bbox_inches='tight')
    

# Bar plot - length of ice season
#plt.figure()
#plt.bar(ice_season_length.isel(xgrid=x, ygrid=y)['time.year'],
#        ice_season_length.isel(xgrid=x, ygrid=y).values)
