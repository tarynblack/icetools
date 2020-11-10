# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 10:35:43 2018

@author: teblack
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# CSV file containing glacier area info
path = '/home/teblack/GreenlandTermini/AGU/'
file = 'Annual_Areas_AGU2018.csv'

# Dataset info
all_glaciers = np.arange(1, 93)
all_years = np.arange(1985, 2018)

# Load glacier area tables into dictionary (key = single observation)
# This code assumes one file with all observations
data = pd.read_csv(path+file, parse_dates=['Date_'])
data = data.to_dict('index')

# Initialize dictionary to reorganize data by glacierID. Observations are a
# list for each glacier; each list item is a dictionary with keys for date and
# area. Also convert pandas timestamp to python datetime.
# glacier_areas = {GlacierID: [{'date': date, 'area': area}]}
glacier_areas = {key: [] for key in all_glaciers}
for key in data:
    glacier = data[key]['GlacierID']
    observation = {'date': data[key]['Date_'].date(),
                   'area': data[key]['Shape_Area']}
    glacier_areas[glacier].append(observation)

# Sort each glacier's data from oldest to youngest
for g in glacier_areas:
    glacier_areas[g] = sorted(glacier_areas[g], key=lambda d: d['date'])

# %%
# FUNCTIONS ===================================================================


def assign_hydro_year(d):
    """Assign a hydrological year to a given date. Greenland hydrological year
    runs from September 1 (YYYY-09-01) to August 31 (YYYY-08-31)."""
    if d.month >= 9:
        hydro_year = int(d.year)
    elif d.month < 9:
        hydro_year = int(d.year - 1)
    return hydro_year


def get_glacier_area_data(g):
    """Extract areas and dates for a glacier."""
    areas = np.empty(0)
    dates = np.empty(0)
    hydro_year = np.empty(0)
    for obs in glacier_areas[g]:
        areas = np.append(areas, obs['area'])
        dates = np.append(dates, obs['date'])
        hydro_year = np.append(hydro_year, assign_hydro_year(obs['date']))
    return areas, dates, hydro_year


def calc_glacier_area_change(g):
    """Calculate glacier area change between dates. For a given date, the area
    change given is the change since the last date. The first date in the
    record has an area change of zero."""
    areas, dates, hydro_year = get_glacier_area_data(g)
    area_change = np.diff(areas)
    area_change = np.insert(area_change, 0, 0)
    return area_change, dates, hydro_year


def calc_glacier_area_cumsum(g):
    """Calculate cumulative glacier area change."""
    area_change, dates, hydro_year = calc_glacier_area_change(g)
    area_cumsum = np.cumsum(area_change)
    return area_cumsum, dates, hydro_year


def calc_glacier_area_cumsum_norm(g, normtype):
    """Normalize cumulative area change by the most-retreated area.
    Allowable normtypes are:
    min: normalize by minimum area in record
    minmax: min-max normalization, (val - min)/(max - min)"""
    area_cumsum, dates, hydro_year = calc_glacier_area_cumsum(g)
    if normtype == 'min':
        area_cumsum_norm = area_cumsum/abs(area_cumsum.min())
    elif normtype == 'minmax':
        area_cumsum_norm = (area_cumsum - area_cumsum.min()) / \
                           (area_cumsum.max() - area_cumsum.min())
    return area_cumsum_norm, dates, hydro_year


def nan_helper(data):
    """Helper to handle indices and logical indices of NaNs. Courtesy of @eat
    and @snake_charmer on Stack Overflow.
    Input:
        data - 1D numpy array with possible NaNs
    Output:
        nans - indices of NaNs
        x    - a helper function to convert logical ind of NaNs to real ind
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x = nan_helper(data)
        >>> data[nans] = np.interp(x(nans), x(~nans), data[~nans])
    """
    return np.isnan(data), lambda z: z.nonzero()[0]


def make_hydro_array(data, field):
    """Turn a dictionary structure of data into an array. Each row is a
    glacier, each column is a hydrological year. Elements are observations (or
    derived calculations) for the specified field, which is a key name in the
    dictionary. NaNs where no values exist."""
    HYdata = np.empty([len(all_glaciers), len(all_years)])
    for g in np.arange(len(all_glaciers)):
        agg = all_glaciers[g]
        for y in np.arange(len(all_years)):
            ayy = all_years[y]
            if ayy in data[agg]['hydro_year']:
                ind = np.where(data[agg]['hydro_year'] == ayy)[0][0]
                HYdata[g, y] = data[agg][field][ind]
            elif ayy not in data[agg]['hydro_year']:
                HYdata[g, y] = np.nan
    return HYdata


def interp_hydro_array(HYdata):
    """Linearly interpolate NaN values in array of data per hydrological year.
    This makes per-year comparison math easier than accessing a dictionary."""
    HYdata_interp = np.copy(HYdata)
    for g in np.arange(HYdata_interp.shape[0]):
        nans, x = nan_helper(HYdata_interp[g, ])
        HYdata_interp[g, nans] = np.interp(x(nans), x(~nans),
                                           HYdata_interp[g, ~nans])
    return HYdata_interp


def vert_offset_lines(prev_line, next_line):
    """Add an offset to arrays so that they are vertically offset when plotted
    (y is essentially a relative scale).
    prev_line : previous line (already plotted)
    next_line : the line that needs to be vertically offset
    prev_line and next_line must be the same size (use interpolated arrays if
    needed)"""
    vert_offset = max(prev_line - next_line)
    shift_next_line = next_line + vert_offset
    return vert_offset, shift_next_line


# %%
# AREA CHANGE CALCULATIONS

# Annual area change for each glacier
glacier_area_change = {key: {} for key in all_glaciers}
for g in glacier_area_change:
    area_change, dates, hydro_year = calc_glacier_area_change(g)
    glacier_area_change[g] = {'date': dates,
                              'hydro_year': hydro_year,
                              'area_change': area_change}

# Create array of area change for each glacier for each hydro year
HY_area_change = make_hydro_array(glacier_area_change, 'area_change')
HY_area_change_interp = interp_hydro_array(HY_area_change)

# Cumulative annual area change for each glacier
glacier_area_cumsum = {key: {} for key in all_glaciers}
for g in glacier_area_cumsum:
    area_cumsum, dates, hydro_year = calc_glacier_area_cumsum(g)
    glacier_area_cumsum[g] = {'date': dates,
                              'hydro_year': hydro_year,
                              'area_cumsum': area_cumsum}

# Create array of cumulative area change interpolated for each year
HY_area_cumsum = make_hydro_array(glacier_area_cumsum, 'area_cumsum')
HY_area_cumsum_interp = interp_hydro_array(HY_area_cumsum)

# Normalized cumulative annual area change for each glacier
normtype = 'minmax'  # can be 'min' or 'minmax'
glacier_area_cumsum_norm = {key: {} for key in all_glaciers}
for g in glacier_area_cumsum_norm:
    area_cumsum_norm, dates, hydro_year = calc_glacier_area_cumsum_norm(g, normtype)
    glacier_area_cumsum_norm[g] = {'date': dates,
                                   'hydro_year': hydro_year,
                                   'area_cumsum_norm': area_cumsum_norm}

# Create array of norm cum area change for each glacier for each hydro year
HY_area_csnorm = make_hydro_array(glacier_area_cumsum_norm, 'area_cumsum_norm')
HY_area_csnorm_interp = interp_hydro_array(HY_area_csnorm)

# Calculate total cumulative annual area change (all glaciers, per year)
HY_area_cumsum_interp_total = HY_area_cumsum_interp.sum(axis=0)

# %%
# STATISTICS AND ANALYSIS

# Get total area change at end of record and sort by magnitude
total_area_change = []  # TODO: make this array rather than list
for g in glacier_area_cumsum:
    total_area_change.append((g, glacier_area_cumsum[g]['area_cumsum'][-1]))
total_area_change_sorted = sorted(total_area_change,
                                  key=lambda glacier: glacier[1])

total_area_change_vals = []  # TODO: make this array rather than list
for g in np.arange(len(total_area_change)):
    total_area_change_vals.append(total_area_change[g][1])

# Identify glaciers that have had net advance in record and create a boolean
# mask to hide them from data as needed.
advancing = []
for g in np.arange(len(total_area_change)):
    if total_area_change[g][1] > 0:
        advancing.append(total_area_change[g][0])
advancing_mask = np.full(len(all_glaciers), True, dtype=bool)
advancing_mask[np.subtract(advancing, 1)] = False

# From visual inspection, four glaciers had much greater loss than others.
# Identify these and create a boolean mask to hide them from data as needed.
# Also identify their names for labeling purposes
big_losers = []
for g in np.arange(4):
    big_losers.append(total_area_change_sorted[g][0])
big_losers_mask = np.full(len(all_glaciers), True, dtype=bool)
big_losers_mask[np.subtract(big_losers, 1)] = False

# Calculate average and median change and average normcum change for each year
avg_annual_change = np.mean(HY_area_change_interp, axis=0)
med_annual_change = np.median(HY_area_change_interp, axis=0)
avg_normcum_change = np.mean(HY_area_csnorm_interp, axis=0)

# Get average and median cumulative change at end of record
total_avg_change = np.mean(total_area_change, axis=0)[1]
total_med_change = np.median(total_area_change, axis=0)[1]

# Identify glaciers with very positive norm area changes (and which year)
# (want to identify the spiky glaciers in the early parts of the record).
# Create a boolean mask to hide them from data as needed.
if normtype == 'min':
    thresh = 1
    spiky_ind = np.where(HY_area_csnorm > thresh)
    spiky_gy = np.empty([2, len(spiky_ind[0])])
    for sp in np.arange(len(spiky_ind[0])):
        spiky_gy[0, sp] = all_glaciers[spiky_ind[0][sp]]
        spiky_gy[1, sp] = all_years[spiky_ind[1][sp]]
    spiky_glaciers = np.unique(spiky_gy[0])
elif normtype == 'minmax':
    thresh = 0.5
    spiky_ind = np.where(HY_area_csnorm_interp[:,0] < thresh)
    spiky_glaciers = spiky_ind[0] + 1
spiky_glaciers_mask = np.full(len(all_glaciers), True, dtype=bool)
spiky_glaciers_mask[(spiky_glaciers - 1).tolist()] = False

# Calculate average and median annual change with big-loss glaciers removed
HY_area_change_interp_nobig = HY_area_change_interp[big_losers_mask]
avg_annual_change_nobig = np.mean(HY_area_change_interp_nobig, axis=0)
med_annual_change_nobig = np.median(HY_area_change_interp_nobig, axis=0)

# Assign names for notable glaciers
glacier_names = {3: "Jakobshavn",
                 9: "Store",
                 35: "Alison",
                 42: "Kjer",
                 78: "Berlingske",
                 92: "Humboldt"}

# Turn all_years into a datetime that starts at the beginning of the hydro year
# Also define xlims - beginning of 1985, end of 2017 hydro years
all_years_dt = pd.to_datetime(all_years, format='%Y') + pd.DateOffset(month=9)
hy_lims = pd.to_datetime([1985, 2019], format='%Y') + pd.DateOffset(month=3)
# %%
# PLOTS FOR AGU

# Total median area change per year, with and without outliers
# (outliers are the four glaciers with the biggest losses)
plt.figure(figsize=(8, 5.5), dpi=300)
plt.title('Annual glacier area change', fontsize=24)
plt.xlabel('Year', fontsize=18)
plt.ylabel('$\Delta$area ($km^2$)', fontsize=18)
plt.plot(all_years_dt, med_annual_change/1E6,
         marker='o', mec='none', linewidth=3, markersize=9,
         color='royalblue', label='median')
plt.plot(all_years_dt, med_annual_change_nobig/1E6,
         marker='o', mec='none', linewidth=2, markersize=8,
         color='orange', label='median (outliers removed)')
plt.grid()
plt.xlim(hy_lims)
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/med_area_change.png',
            bbox_inches='tight', pad_inches=0, dpi=600)

# %%
# Histogram of total loss
areas_to_km = np.array(total_area_change_vals)/1E6
plt.figure(figsize=(8, 5.5), dpi=300)
plt.title('Total area change per glacier', fontsize=24)
plt.xlabel('Area change ($km^2$) from 1985 to 2018', fontsize=18)
plt.ylabel('Number of glaciers', fontsize=18)
plt.hist(np.clip(areas_to_km, -50, 10), bins=np.arange(-50, 10, 2.5),
         color='royalblue')
plt.grid()
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/hist_area_change.png',
            bbox_inches='tight', pad_inches=0, dpi=600)

# %%
# Cumulative area change for each glacier EXCEPT four big losers
colors = plt.cm.plasma(np.linspace(0, 1, 92))
plt.figure(figsize=(8, 5.5), dpi=300)
plt.title('Cumulative area change', fontsize=24)
plt.xlabel('Year', fontsize=18)
plt.ylabel('$\Delta$area ($km^2$)', fontsize=18)
for g in glacier_area_cumsum:
    if g in big_losers:
        continue
    else:
        plt.plot(glacier_area_cumsum[g]['date'],
                 glacier_area_cumsum[g]['area_cumsum']/1E6,
                 marker='o', mec='none', markersize=6, color=colors[g])
plt.grid()
plt.xlim(hy_lims)
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/cum_area_change_nobigs.png',
            bbox_inches='tight', pad_inches=0, dpi=600)

# %%
# Cumulative area change for big losers (and others, in gray)
colors = {92: 'royalblue', 3: 'orange',
          35: 'mediumseagreen', 42: 'mediumorchid'}
plt.figure(figsize=(8, 5.5), dpi=300)
plt.title('Cumulative area change', fontsize=24)
plt.xlabel('Year', fontsize=18)
plt.ylabel('$\Delta$area ($km^2$)', fontsize=18)
for g in glacier_area_cumsum:
    if g in big_losers:
        continue
    else:
        plt.plot(glacier_area_cumsum[g]['date'],
                 glacier_area_cumsum[g]['area_cumsum']/1E6,
                 color='silver', marker='.', mec='none')
for g in big_losers:
    plt.plot(glacier_area_cumsum[g]['date'],
             glacier_area_cumsum[g]['area_cumsum']/1E6,
             marker='o', mec='none', linewidth=2, markersize=6,
             color=colors[g], label=glacier_names[g])
plt.grid()
plt.xlim(hy_lims)
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='lower left', fontsize=12)
plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/cum_area_change_bigs.png',
            bbox_inches='tight', pad_inches=0, dpi=600)

# %%
# Cumulative area changes for individual glaciers (single plots)
# Humboldt
plt.figure(figsize=(5.5, 2), dpi=300)
plt.title('Cumul. $\Delta$area: Humboldt', fontsize=24)
plt.xlabel('Year', fontsize=18)
plt.ylabel('$\Delta$area ($km^2$)', fontsize=18)
plt.plot(glacier_area_cumsum[92]['date'],
         glacier_area_cumsum[92]['area_cumsum']/1E6,
         marker='o', mec='none', linewidth=3, markersize=6, color='royalblue')
plt.grid()
plt.xlim(hy_lims)
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/humboldt_area_change.png',
            bbox_inches='tight', pad_inches=0, dpi=600)

# Store
plt.figure(figsize=(5.5, 2), dpi=300)
plt.title('Cumul. $\Delta$area: Store', fontsize=24)
plt.xlabel('Year', fontsize=18)
plt.ylabel('$\Delta$area ($km^2$)', fontsize=18)
plt.plot(glacier_area_cumsum[9]['date'],
         glacier_area_cumsum[9]['area_cumsum']/1E6,
         marker='o', mec='none', linewidth=3, markersize=6, color='royalblue')
plt.grid()
plt.xlim(hy_lims)
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/store_area_change.png',
            bbox_inches='tight', pad_inches=0, dpi=600)

# Berlingske
plt.figure(figsize=(5.5, 2), dpi=300)
plt.title('Cumul. $\Delta$area: Berlingske', fontsize=24)
plt.xlabel('Year', fontsize=18)
plt.ylabel('$\Delta$area ($km^2$)', fontsize=18)
plt.plot(glacier_area_cumsum[78]['date'],
         glacier_area_cumsum[78]['area_cumsum']/1E6,
         marker='o', mec='none', linewidth=3, markersize=6, color='royalblue')
plt.grid()
plt.xlim(hy_lims)
ylim2 = plt.ylim()
plt.ylim(-0.1, ylim2[1])
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/berlingske_area_change.png',
            bbox_inches='tight', pad_inches=0, dpi=600)

# Jakobshavn
plt.figure(figsize=(5.5, 2), dpi=300)
plt.title('Cumul. $\Delta$area: Jakobshavn', fontsize=24)
plt.xlabel('Year', fontsize=18)
plt.ylabel('$\Delta$area ($km^2$)', fontsize=18)
plt.plot(glacier_area_cumsum[3]['date'],
         glacier_area_cumsum[3]['area_cumsum']/1E6,
         marker='o', mec='none', linewidth=3, markersize=6, color='royalblue')
plt.grid()
plt.xlim(hy_lims)
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/jakobshavn_area_change.png',
            bbox_inches='tight', pad_inches=0, dpi=600)

# %%
# Normalized cumulative change in area for each glacier (shows timing well)
# All glaciers in gray, overlay average norm cum area change in darker
# Need to turn all_years into a datetime
plt.figure(figsize=(8, 5.5), dpi=300)
plt.title('Normalized cumul. $\Delta$area', fontsize=24)
plt.xlabel('Year', fontsize=18)
plt.ylabel('min-max norm. $\Delta$area', fontsize=18)
for g in all_glaciers:
    dates = glacier_area_cumsum_norm[g]['date']
    areas_cumsum_norm = glacier_area_cumsum_norm[g]['area_cumsum_norm']
    line_ind, = plt.plot(dates, areas_cumsum_norm, color='silver', marker='o',
                         mec='none', label='individual glaciers',
                         linewidth=0.75, markersize=2)
line_avg, = plt.plot(all_years_dt, avg_normcum_change, color='royalblue',
                     linewidth=6, label='average')
plt.grid()
plt.xlim(hy_lims)
plt.ylim(-0.02, 1.02)
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(handles=[line_ind, line_avg], fontsize=12)
plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/norm_area_change.png',
            bbox_inches='tight', pad_inches=0, dpi=600)

# %%
# Scatter plot of terminus observations
plt.figure(figsize=(6, 4.5), dpi=300)
plt.title('Terminus observations, 1985-2018', fontsize=24)
plt.xlabel('Year', fontsize=18)
plt.ylabel('Glacier ID', fontsize=18)
for g in all_glaciers:
    areas, dates, hydro_year = get_glacier_area_data(g)
    plt.scatter(hydro_year, np.full(len(hydro_year), g),
                marker='o', c='royalblue', edgecolors='face')
plt.grid()
plt.xlim(1984, 2018)
plt.ylim(0, 93)
plt.tick_params(axis='both', direction='out', length=6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.savefig('/home/teblack/GreenlandTermini/AGU/figs/observational_record.png',
#            bbox_inches='tight', pad_inches=0, dpi=600)

# %%
# Plot total cumulative annual area change
plt.figure
plt.title('Total cumulative area change, 1985-2018')
plt.xlabel('Year')
plt.ylabel('Cumulative area change ($km^2$)')
plt.plot(all_years_dt, HY_area_cumsum_interp_total/1E6)
plt.grid()
plt.xlim(hy_lims)
plt.tick_params(axis='both', direction='out', length=6)

# %%
# Average normalized cumulative area change --> combine with all-norms figure
#plt.figure()
#plt.title('Average normalized cumulative area change')
#plt.xlabel('Year')
#plt.ylabel('Normalized cumulative area change')
#plt.plot(all_years, avg_normcum_change, marker='o', mec='none')
#plt.grid()

# %%
# Plot all cumulative area changes with vertical shift (~latitude plot)
# This is tricky because some have big area changes... redo without big_losers
#plt.figure()
#plt.title('Cumulative area change, separated...')
#plt.xlabel('Year')
#plt.ylabel('Cumulative area change ($m^2$)')
#total_offset = 0
#for g in glacier_area_cumsum:
#    if g == 1:
#        plt.plot(glacier_area_cumsum[g]['date'],
#                 glacier_area_cumsum[g]['area_cumsum'],
#                 marker='o', mec='none')
#    else:
#        offset_val, offset_line = vert_offset_lines(HY_area_cumsum_interp[g-2],
#                                                    HY_area_cumsum_interp[g-1])
#        total_offset += offset_val
#        plt.plot(glacier_area_cumsum[g]['date'],
#                 glacier_area_cumsum[g]['area_cumsum'] + total_offset,
#                 marker='o', mec='none')
#plt.grid()

# %%
# Plot norm cumulative area changes with vertical shift
# Because they are all 0 -> 1 can just shift each up by one from last
# Still hard to see variation, so scale by some factor. They'll end up 
# overlapping a bit but maybe it will look more like a seismogram where most
# lines are separated except where there's an event.
#scale = 20
#plt.figure()
#plt.title('Normalized cumulative area change, separated...')
#plt.xlabel('Year')
#plt.ylabel('Normalized cumulative area change')
#for g in glacier_area_cumsum_norm:
#    shift_GACN = scale*glacier_area_cumsum_norm[g]['area_cumsum_norm'] + g - 1
#    plt.plot(glacier_area_cumsum_norm[g]['date'], shift_GACN,
#             marker='o', mec='none')
#plt.grid()
# %%
# OTHER PLOTS (TESTS, ETC)

# Scatter plot of observation instances for each glacier (buggy)
#plt.figure()
#plt.title('Glacier satellite observations, 1985-2018')
#plt.xlabel('Year')
#plt.ylabel('Glacier ID')
#for g in all_glaciers:
#    areas, dates, hydro_year = get_glacier_area_data(g)
#    plt.scatter(dates, np.full(dates.shape, g), marker='o', color='royalblue')
#plt.grid()
#plt.tick_params(axis='both', direction='out', length=6)
##plt.axis('tight')
#plt.yticks(np.arange(0, 92, step=5))
##hy = [h for h in hydro_year]
##hydro_year_dt = pd.to_datetime(np.add(hy, 1), format='%Y')
#plt.xticks(all_years_dt)

#plt.xticks(np.arange(1985, 2018), rotation=90, fontsize=12)
#plt.yticks(np.arange(0, 92, step=5), fontsize=12)
    
# Plot single glacier areas
#areas1, dates1 = get_glacier_area_data(1)
#plt.figure()
#plt.title('Glacier 1 area')
#plt.plot(dates1, areas1, marker='o')
#plt.show()

# Plot all glacier areas on one plot
#plt.figure()
#plt.title('Total glacier areas, 1985-2018')
#plt.xlabel('Year')
#plt.ylabel('Area (sq m)')
#for g in all_glaciers:
#    areas, dates = get_glacier_area_data(g)
#    plt.plot(dates, areas, marker='o')
#plt.show()

# Plot single glacier area change
#area_change1, dates1 = calc_glacier_area_change(1)
#plt.figure()
#plt.title('Glacier 1 area change')
#plt.plot(dates1, area_change1, marker='o')
#plt.show()

## Plot single glacier cumulative area change
#area_cumsum1, dates1 = calc_glacier_area_cumsum(1)
#plt.figure()
#plt.title('Glacier 1 cumulative area change')
#plt.plot(dates1, area_cumsum1, marker='o')

## Plot normalized cumulative area change for a single glacier
#area_cumsum_norm1, dates1 = calc_glacier_area_cumsum_norm(1)
#plt.figure()
#plt.title('Glacier 1 normalized cumulative area change')
#plt.plot(dates1, area_cumsum_norm1, marker='o')