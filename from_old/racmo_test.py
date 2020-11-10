import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

datadir = "~/data/RACMO/data1km_1972-2018/"
datasuf = ".1972-2018.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc"

#precipdata = xr.open_dataset(datadir + "precip" + datasuf, decode_times=False)
smbdata = xr.open_dataset(datadir + "smb_rec" + datasuf, decode_times=False)
#runoffdata = xr.open_dataset(datadir + "runoff" + datasuf, decode_times=False)
#snowmeltdata = xr.open_dataset(datadir + "snowmelt" + datasuf, decode_times=False)
#t2mdata = xr.open_dataset(datadir + "t2m" + datasuf, decode_times=False)
icemask = xr.open_dataset(datadir + "Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc")

# For some reason, most files have lat and lon in polarstereo meters, but SMB
# has x,y as integers, and icemask has x,y as polarstereo meters. Fix.
icemask = icemask.rename({'y': 'lat', 'x': 'lon'})
smbdata = smbdata.rename({'y': 'lat', 'x': 'lon'})
smbdata = smbdata.assign_coords(lon=(smbdata.lon * 1000 + icemask.lon.min()),
                                lat=(smbdata.lat * 1000 + icemask.lat.min()))

# Boundaries for NW Greenland (WGS84 NSIDC polar stereo north; est. from map)
bnd_S = -2400000
bnd_N = -705557
bnd_W = -670279
bnd_E = 52916
bounds_NWG = [bnd_S, bnd_N, bnd_W, bnd_E]

# Time limits
startyr = 1972
endyr = 2018

# %% FUNCTIONS

def time_ind2date(ds):
    """Convert RACMO time units from "Months since 1972-01-15 00:00:00" to datetime64"""
    time = ds.time.values
    dtime = pd.to_timedelta(time, 'M')
    reftime = np.datetime64('1972-01-15')
    newtime = reftime + dtime
    ds['time'] = newtime


def subset_geog_data(ds, bounds):
    """ Bounds are [bottom top left right] """
    data_subset = ds.sel(lon=slice(bounds[2], bounds[3]),
                         lat=slice(bounds[0], bounds[1]))
    return data_subset


def ds_to_da(ds, varname):
    da = ds.to_array(dim=varname, name=varname)
    da[varname].attrs = ds[varname].attrs
    return da


def clean_racmo_to_da(ds, bounds, datavar, icemask, topomask):
    """After reading the data: fix the time units, subset to NW Greenland,
    convert to a dataarray, and mask to ice-only. Inputs are:
        - data     : the dataset as read in by xarray
        - bounds   : subset region boundaries
        - datavar  : the name of the data variable in the dataset
        - icemask  : RACMO icemask, regional subset if necessary
        - topomask : RACMO topography where ice is present"""
    time_ind2date(ds)
    data_NWG = subset_geog_data(ds, bounds)
    data_da = ds_to_da(data_NWG, datavar)
    data_masked = data_da.where(icemask.Promicemask == 3)
    data_masked['topo'] = topomask
    return data_masked


def get_topobin_stats(da, topo):
    binstat = np.empty((4, topo.size-1, da.time.size))
    for t in np.arange(da.time.size):
        bindata = da.isel(time=t).groupby_bins('topo', topo)
        av = bindata.mean().data
        md = bindata.median().data
        mx = bindata.max().data
        mn = bindata.min().data
        binstat[:, :, t] = [av, md, mx, mn]

    stats = xr.Dataset({'average': (['topo_bins', 'time'], binstat[0, :, :]),
                        'medians': (['topo_bins', 'time'], binstat[1, :, :]),
                        'maximum': (['topo_bins', 'time'], binstat[2, :, :]),
                        'minimum': (['topo_bins', 'time'], binstat[3, :, :])},
                        coords={'topo_bins': bindata.last().coords['topo_bins'].data,
                                'time': da.coords['time'].data})
    return stats


def plot_average_timeseries(ds_stats, varname, varunits):
    plt.figure()
    for tb in np.arange(ds_stats.topo_bins.size):
            label = '%s - %s m' % (ds_stats.topo_bins.data[tb].left,
                                   ds_stats.topo_bins.data[tb].right)
            ds_stats.average.isel(topo_bins=tb).plot(label=label)
    plt.legend(title='Elevation', bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.title('%s, binned by elevation' % varname)
    plt.ylabel('Average %s (%s)' % (varname, varunits))


def plot_topobin_stats(ds_stats, bin, varname, varunits):
    ds_stats.average.isel(topo_bins=bin).plot(label='average')
    ds_stats.maximum.isel(topo_bins=bin).plot(label='maximum')
    ds_stats.minimum.isel(topo_bins=bin).plot(label='minimum')
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.title('%s, %s-%s m' % (varname, ds_stats.topo_bins.data[bin].left,
                               ds_stats.topo_bins.data[bin].right))
    plt.ylabel('%s (%s)' % (varname, varunits))
    plt.xlabel('time')


def plot_topobin_stats_subplots(ds_stats, varname, varunits):
    plt.figure()
    # Find actual number of bins (not all bins have data)
    for tb in np.arange(ds_stats.topo_bins.size):
        if ds_stats.maximum.isel(topo_bins=tb).sum() == 0:
            tb_max = tb
            break
    for tb in np.arange(tb_max):
        if tb == 0:
            ax0 = plt.subplot(tb_max, 1, tb+1)
        else:
            plt.subplot(tb_max, 1, tb+1, sharex=ax0)
        plot_topobin_stats(ds_stats, tb, varname, varunits)


def hydroyear_data(da, year):
    """Get data from a hydrological year. Greenland hydrological year is
    defined as September 1 through August 31.
    Input year is the starting year (aligned with September)."""
    da_hy = da.sel(time=slice('%s-09-01' % str(year),
                              '%s-08-31' % str(year+1)))
    return da_hy


def hydroyear_cumsum(da, year):
    """Calculate cumulative sum of a variable over a hydrological year.
    Hydrological year in Greenland is defined as Sept. 1 through Aug. 31.
    Input year is the starting year (aligned with September)."""
    da_hy = hydroyear_data(da, year)
    da_cumsum = da_hy.cumsum(dim='time')
    return da_cumsum


def plot_annual_cumsum(da_stats, topo_bin, varname, varunits):
    """Plot the annual cumulative sum of a variable's average within a 
    topographic bin."""
    annual_colors = plt.cm.plasma(np.linspace(0, 1, endyr-startyr))
    plt.figure()
    for year in np.arange(startyr, endyr):
        hydroyear_cumsum(da_stats, year).average.isel(topo_bins=topo_bin).plot(
                color=annual_colors[year-startyr], label=year)
    plt.legend(bbox_to_anchor=(1.01, 1.1), loc='upper left')
    plt.title('Cumulative annual %s, %s - %s m' % (varname,
            da_stats.topo_bins.data[topo_bin].left, 
            da_stats.topo_bins.data[topo_bin].right))
    plt.ylabel('%s (%s)' % (varname, varunits))
    plt.xlabel('Months since September')


def annual_average(da_stats, year):
    """Calculate the average of a variable in each topo bin over a hydrological 
    year."""
    da_hy = hydroyear_data(da_stats, year)
    ann_avg_eachbin = da_hy.average.mean(dim='time')
    return ann_avg_eachbin


def annual_average_allyears(da_stats, topo_bin, startyr, endyr):
    """Calculate the annual average of a variable for each year."""
    ann_avg = np.zeros(np.arange(startyr, endyr).size)
    for year in np.arange(startyr, endyr):
        ann_avg[year-startyr] = annual_average(da_stats, year).isel(
                topo_bins=topo_bin)
    return ann_avg


def plot_annual_average(da_stats, topo_bin, startyr, endyr, varname, varunits):
    """Plot the annual average of a variable in each topo bin."""
    ann_avg = annual_average_allyears(da_stats, topo_bin, startyr, endyr)
    times = np.arange(str(startyr), str(endyr), dtype='datetime64[Y]')
    h, = plt.plot(times, ann_avg, '.-', linewidth=2)
#    plt.title('Annual average %s, %s - %s m' % (varname,
#            da_stats.topo_bins.data[topo_bin].left, 
#            da_stats.topo_bins.data[topo_bin].right))
    plt.ylabel('%s (%s)' % (varname, varunits))
#    plt.xlabel('Year')
    return h


def plot_annual_average_bar(da_stats, topo_bin, startyr, endyr, varname, varunits):
    """Plot the annual average of a variable in each topo bin."""
    ann_avg = annual_average_allyears(da_stats, topo_bin, startyr, endyr)
    times = np.arange(startyr, endyr)
    plt.figure()
    plt.bar(times, ann_avg)
    plt.title('Annual average %s, %s - %s m' % (varname,
            da_stats.topo_bins.data[topo_bin].left, 
            da_stats.topo_bins.data[topo_bin].right))
    plt.ylabel('%s (%s)' % (varname, varunits))
    plt.xlabel('Year')


def cumulative_anomaly(da_stats, topo_bin, annual=False):
    """Calculate cumulative anomaly of a variable in a topo bin."""
    
    if annual == True:
        ann_avg = annual_average_allyears(da_stats, topo_bin, startyr, endyr)
        series_avg = np.nanmean(ann_avg)
        series_anomaly = ann_avg - series_avg
        cumul_anomaly = np.nancumsum(series_anomaly)
    else:
        series_avg = da_stats.isel(topo_bins=topo_bin).average.mean()
        series_anomaly = da_stats.isel(topo_bins=topo_bin).average - series_avg
        cumul_anomaly = series_anomaly.cumsum(dim='time')
    
    return cumul_anomaly


def plot_cumulative_anomaly(da_stats, topo_bin, startyr, endyr, varname, varunits, annual=False):
    """Plot the cumulative anomaly of a variable in a topo bin."""
    cumul_anom = cumulative_anomaly(da_stats, topo_bin, annual)
#    if type(cumul_anom) == xr.core.dataarray.DataArray:
#        cumul_anom.plot()
#        plt.title('Cumulative anomaly %s, %s - %s m' % (varname,
#                                    da_stats.topo_bins.data[topo_bin].left, 
#                                    da_stats.topo_bins.data[topo_bin].right))
#        plt.ylabel('%s (%s)' % (varname, varunits))
#    elif type(cumul_anom) == np.ndarray:
    times = np.arange(str(startyr), str(endyr), dtype='datetime64[Y]')
    h, = plt.plot(times, cumul_anom, '.-')
#        plt.title('Cumulative anomaly %s, %s - %s m' % (varname,
#                                    da_stats.topo_bins.data[topo_bin].left, 
#                                    da_stats.topo_bins.data[topo_bin].right))
#        plt.ylabel('%s (%s)' % (varname, varunits))
#        plt.xlabel('Year')
    return h
        

def seasonal_average(da_stats, season, topo_bin, startyr, endyr):
    """Calculate seasonal average of data in each year for a topo bin. Seasons
    can be either 'DJF' (winter) or 'JJA' (summer)."""
    season_data = da_stats.where(da_stats.time.dt.season==season)
    season_avg = annual_average_allyears(season_data, topo_bin, startyr, endyr)
    if season == 'DJF':
        times = pd.date_range(str(startyr+1), periods=endyr-startyr, freq='AS-JAN')
    elif season == 'JJA':
        times = pd.date_range(str(startyr+1), periods=endyr-startyr, freq='AS-JUL')
    return season_avg, times


def plot_seasonal_average(da_stats, season, topo_bin, startyr, endyr, varname, varunits):
    """Plot seasonal average of variable in each year for a topo bin."""
    data, time = seasonal_average(da_stats, season, topo_bin, startyr, endyr)
    h, = plt.plot(time, data, '.-')
#    plt.title('Average %s %s, %s - %s m' % (season, varname,
#                                    da_stats.topo_bins.data[topo_bin].left, 
#                                    da_stats.topo_bins.data[topo_bin].right))
#    plt.ylabel('%s (%s)' % (varname, varunits))
#    plt.xlabel('Time')
    return h

# %%
# Topography: subset to NWG, apply ice-sheet mask
icemask_NWG = subset_geog_data(icemask, bounds_NWG)
topo_NWG = icemask_NWG.Topography.where(icemask.Promicemask == 3)
#topo_NWG.name = 'topo'

## Define topographic bins, in 500m intervals.
#topo_range = np.arange(0, 4500, 500)
#
## Process and plot data! Clean up the RACMO data, get topographic stats, plot:
## - time series of average data for each topographic band, one plot.
## - time series of max/avg/min for each topo band, on subplots.
#
#print("Cleaning SMB...")
#smb_masked = clean_racmo_to_da(smbdata, bounds_NWG, 'SMB_rec', icemask, topo_NWG)
#print("Getting SMB stats...")
#smb_stats = get_topobin_stats(smb_masked, topo_range)
#print("Plotting SMB...")
#smbname = 'SMB'
#smbunits = smb_masked.SMB_rec.units
#plot_average_timeseries(smb_stats, smbname, smbunits)
#plot_topobin_stats_subplots(smb_stats, smbname, smbunits)
##plot_annual_cumsum(smb_stats, 0, smbname, smbunits)
#plot_annual_average(smb_stats, 0, startyr, endyr, smbname, smbunits)
##plot_annual_average_bar(smb_stats, 0, startyr, endyr, smbname, smbunits)
## Also plot cumulative SMB (which doesn't make as much sense for other vars)
#plt.figure()
#smb_cumsum = smb_stats.average.isel(topo_bins=0).cumsum(dim='time').plot()
#plt.title('Cumulative SMB 1972-2018, %s - %s m' % (
#        smb_stats.topo_bins.data[0].left, smb_stats.topo_bins.data[0].right))
#plt.ylabel('SMB (%s)' % smb_masked.SMB_rec.units)
##plot_cumulative_anomaly(smb_stats, 0, startyr, endyr, smbname, smbunits)
##plot_cumulative_anomaly(smb_stats, 0, startyr, endyr, smbname, smbunits, annual=True)
#print("Done")

# %% SMB in four latitude bands for four glacier groups; 0-500m elevation
# 12/5/2019 - AGU
# Define topographic bins, in 500m intervals.
topo_range = np.arange(0, 501, 500)

# Define latitude bands ([S, N]) and new region boundaries
band1 = [-2400000, -1810000]
bounds1 = band1 + bounds_NWG[2:4]
band2 = [-1810000, -1515000]
bounds2 = band2 + bounds_NWG[2:4]
band3 = [-1515000, -1395000]
bounds3 = band3 + bounds_NWG[2:4]
band4 = [-1395000, -1030000]
bounds4 = band4 + bounds_NWG[2:4]

print("Cleaning SMB...")
smb_masked1 = clean_racmo_to_da(smbdata, bounds1, 'SMB_rec', icemask, topo_NWG)
smb_masked2 = clean_racmo_to_da(smbdata, bounds2, 'SMB_rec', icemask, topo_NWG)
smb_masked3 = clean_racmo_to_da(smbdata, bounds3, 'SMB_rec', icemask, topo_NWG)
smb_masked4 = clean_racmo_to_da(smbdata, bounds4, 'SMB_rec', icemask, topo_NWG)
print("Getting SMB stats...")
smb_stats1 = get_topobin_stats(smb_masked1, topo_range)
smb_stats2 = get_topobin_stats(smb_masked2, topo_range)
smb_stats3 = get_topobin_stats(smb_masked3, topo_range)
smb_stats4 = get_topobin_stats(smb_masked4, topo_range)
print("Plotting SMB...")
smbname = 'SMB'
smbunits = smb_masked1.SMB_rec.units
print("Done")

# %% Plot SMB
# 12/5/2019 - AGU

COLORMAP = cm.Blues

LINEWIDTH = 5
MARKERSIZE = 15
TITLESIZE = 24
LABELSIZE = 20
TICKLABELSIZE = 16
LEGENDSIZE = 16

# Plot cumulative anomaly SMB for each zone
fig, ax = plt.subplots(figsize=(10,5))
color = iter(COLORMAP(np.linspace(0.2,1,4)))
h1 = plot_cumulative_anomaly(smb_stats1, 0, startyr, endyr, smbname, smbunits, annual=True)
h1.set_color(next(color))
h1.set_linewidth(LINEWIDTH)
h1.set_markersize(MARKERSIZE)
h2 = plot_cumulative_anomaly(smb_stats2, 0, startyr, endyr, smbname, smbunits, annual=True)
h2.set_color(next(color))
h2.set_linewidth(LINEWIDTH)
h2.set_markersize(MARKERSIZE)
h3 = plot_cumulative_anomaly(smb_stats3, 0, startyr, endyr, smbname, smbunits, annual=True)
h3.set_color(next(color))
h3.set_linewidth(LINEWIDTH)
h3.set_markersize(MARKERSIZE)
h4 = plot_cumulative_anomaly(smb_stats4, 0, startyr, endyr, smbname, smbunits, annual=True)
h4.set_color(next(color))
h4.set_linewidth(LINEWIDTH)
h4.set_markersize(MARKERSIZE)
ax.legend((h4, h3, h2, h1), ('Zone 4', 'Zone 3', 'Zone 2', 'Zone 1'), fontsize=LEGENDSIZE)
plt.tick_params(axis='both', labelsize=TICKLABELSIZE)
plt.grid(color='lightgray')
plt.xlabel("Year", size=LABELSIZE)
plt.ylabel('%s (%s)' % (smbname, smbunits), size=LABELSIZE)
plt.title("Cumulative Anomaly of Surface Mass Balance\n0-500 m Elevation", size=TITLESIZE)
plt.savefig("/home/teblack/plots/AGU_SMB_cumanom.png", transparent=True, bbox_inches='tight')

# %% PLOT ANNUAL AND SEASONAL AVERAGE SMB FOR EACH ZONE AS SUBPLOTS
fig = plt.figure(figsize=(10,10))
plt.suptitle("Average Surface Mass Balance\n0-500 m Elevation", size=TITLESIZE)

# Plot annual average for each zone
#fig, ax = plt.subplots()
ax1 = fig.add_subplot(311)
color = iter(COLORMAP(np.linspace(0.2,1,4)))
a1 = plot_annual_average(smb_stats1, 0, startyr, endyr, smbname, smbunits)
a1.set_color(next(color))
a1.set_linewidth(LINEWIDTH)
a1.set_markersize(MARKERSIZE)
a2 = plot_annual_average(smb_stats2, 0, startyr, endyr, smbname, smbunits)
a2.set_color(next(color))
a2.set_linewidth(LINEWIDTH)
a2.set_markersize(MARKERSIZE)
a3 = plot_annual_average(smb_stats3, 0, startyr, endyr, smbname, smbunits)
a3.set_color(next(color))
a3.set_linewidth(LINEWIDTH)
a3.set_markersize(MARKERSIZE)
a4 = plot_annual_average(smb_stats4, 0, startyr, endyr, smbname, smbunits)
a4.set_color(next(color))
a4.set_linewidth(LINEWIDTH)
a4.set_markersize(MARKERSIZE)
ax1.legend((a4, a3, a2, a1), ('Zone 4', 'Zone 3', 'Zone 2', 'Zone 1'), fontsize=LEGENDSIZE*.9)
plt.tick_params(axis='both', labelsize=TICKLABELSIZE)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.grid(color='lightgray')
#plt.xlabel("Year", size=LABELSIZE)
plt.ylabel("")
plt.title("Annual", size=TITLESIZE*0.8)

# Add winter (DJF) average
ax2 = fig.add_subplot(312, sharex=ax1)
color = iter(cm.Purples(np.linspace(0.2,1,4)))
h1 = plot_seasonal_average(smb_stats1, 'DJF', 0, startyr, endyr, smbname, smbunits)
h1.set_color(next(color))
h1.set_linewidth(LINEWIDTH)
h1.set_markersize(MARKERSIZE)
h2 = plot_seasonal_average(smb_stats2, 'DJF', 0, startyr, endyr, smbname, smbunits)
h2.set_color(next(color))
h2.set_linewidth(LINEWIDTH)
h2.set_markersize(MARKERSIZE)
h3 = plot_seasonal_average(smb_stats3, 'DJF', 0, startyr, endyr, smbname, smbunits)
h3.set_color(next(color))
h3.set_linewidth(LINEWIDTH)
h3.set_markersize(MARKERSIZE)
h4 = plot_seasonal_average(smb_stats4, 'DJF', 0, startyr, endyr, smbname, smbunits)
h4.set_color(next(color))
h4.set_linewidth(LINEWIDTH)
h4.set_markersize(MARKERSIZE)
ax2.legend((h4, h3, h2, h1), ('Zone 4', 'Zone 3', 'Zone 2', 'Zone 1'), fontsize=LEGENDSIZE*.9)
plt.title('Winter (DJF)', size=TITLESIZE*0.8)
plt.tick_params(axis='both', labelsize=TICKLABELSIZE)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.grid(color='lightgray')
plt.ylabel('%s (%s)' % (smbname, smbunits), size=LABELSIZE)

# Add summer (JJA) average
ax3 = fig.add_subplot(313, sharex=ax1)
color = iter(cm.Oranges(np.linspace(0.2,1,4)))
hh1 = plot_seasonal_average(smb_stats1, 'JJA', 0, startyr, endyr, smbname, smbunits)
hh1.set_color(next(color))
hh1.set_linewidth(LINEWIDTH)
hh1.set_markersize(MARKERSIZE)
hh2 = plot_seasonal_average(smb_stats2, 'JJA', 0, startyr, endyr, smbname, smbunits)
hh2.set_color(next(color))
hh2.set_linewidth(LINEWIDTH)
hh2.set_markersize(MARKERSIZE)
hh3 = plot_seasonal_average(smb_stats3, 'JJA', 0, startyr, endyr, smbname, smbunits)
hh3.set_color(next(color))
hh3.set_linewidth(LINEWIDTH)
hh3.set_markersize(MARKERSIZE)
hh4 = plot_seasonal_average(smb_stats4, 'JJA', 0, startyr, endyr, smbname, smbunits)
hh4.set_color(next(color))
hh4.set_linewidth(LINEWIDTH)
hh4.set_markersize(MARKERSIZE)
ax3.legend((hh4, hh3, hh2, hh1), ('Zone 4', 'Zone 3', 'Zone 2', 'Zone 1'), fontsize=LEGENDSIZE*.9)
plt.title('Summer (JJA)', size=TITLESIZE*0.8)
plt.tick_params(axis='both', labelsize=TICKLABELSIZE)
plt.grid(color='lightgray')
plt.xlabel("Year", size=LABELSIZE)
#plt.ylabel('%s (%s)' % (smbname, smbunits), size=LABELSIZE)

plt.savefig("/home/teblack/plots/AGU_SMB_annavg.png", transparent=True, bbox_inches='tight')

# %%

print("Cleaning precipitation...")
precip_masked = clean_racmo_to_da(precipdata, bounds_NWG, 'precipcorr', icemask, topo_NWG)
print("Getting precipitation stats...")
precip_stats = get_topobin_stats(precip_masked, topo_range)
print("Plotting precipitation...")
precipname = 'Precipitation'
precipunits = precip_masked.precipcorr.units
#plot_average_timeseries(precip_stats, precipname, precipunits)
#plot_topobin_stats_subplots(precip_stats, precipname, precipunits)
#plot_annual_cumsum(precip_stats, 0, precipname, precipunits)
#plot_annual_average(precip_stats, 0, startyr, endyr, precipname, precipunits)
#plot_annual_average_bar(precip_stats, 0, startyr, endyr, precipname, precipunits)
print("Done")

print("Cleaning runoff...")
runoff_masked = clean_racmo_to_da(runoffdata, bounds_NWG, 'runoffcorr', icemask, topo_NWG)
print("Getting runoff stats...")
runoff_stats = get_topobin_stats(runoff_masked, topo_range)
print("Plotting runoff...")
runoffname = 'Runoff'
runoffunits = runoff_masked.runoffcorr.units
#plot_average_timeseries(runoff_stats, runoffname, runoffunits)
#plot_topobin_stats_subplots(runoff_stats, runoffname, runoffunits)
#plot_annual_cumsum(runoff_stats, 0, runoffname, runoffunits)
#plot_annual_average(runoff_stats, 0, startyr, endyr, runoffname, runoffunits)
#plot_annual_average_bar(runoff_stats, 0, startyr, endyr, runoffname, runoffunits)
print("Done")

print("Cleaning snowmelt...")
snowmelt_masked = clean_racmo_to_da(snowmeltdata, bounds_NWG, 'snowmeltcorr', icemask, topo_NWG)
print("Getting snowmelt stats...")
snowmelt_stats = get_topobin_stats(snowmelt_masked, topo_range)
print("Plotting snowmelt...")
snowmeltname = 'Snowmelt'
snowmeltunits = snowmelt_masked.snowmeltcorr.units
#plot_average_timeseries(snowmelt_stats, snowmeltname, snowmeltunits)
#plot_topobin_stats_subplots(snowmelt_stats, snowmeltname, snowmeltunits)
#plot_annual_cumsum(snowmelt_stats, 0, snowmeltname, snowmeltunits)
#plot_annual_average(snowmelt_stats, 0, startyr, endyr, snowmeltname, snowmeltunits)
#plot_annual_average_bar(snowmelt_stats, 0, startyr, endyr, snowmeltname, snowmeltunits)
print("Done")

print("Cleaning temperature...")
t2m_masked = clean_racmo_to_da(t2mdata, bounds_NWG, 't2mcorr', icemask, topo_NWG)
print("Getting temperature stats...")
t2m_stats = get_topobin_stats(t2m_masked, topo_range)
print("Plotting temperature...")
t2mname = '2m temperature'
t2munits = t2m_masked.t2mcorr.units
#plot_average_timeseries(t2m_stats, t2mname, t2munits)
#plot_topobin_stats_subplots(t2m_stats, t2mname, t2munits)
#plot_annual_average(t2m_stats, 0, startyr, endyr, t2mname, t2munits)
print("Done")

