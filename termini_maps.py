#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 17:53:05 2019

@author: teblack
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import rasterio as rio
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.pyplot import cm

mask_file = "/Volumes/insar5/teblack/data/GIMP_mask/GimpOceanMask_90m.shp"
ice_tif = "/Volumes/insar5/teblack/data/GIMP_mask/GimpIceMask_90m_v1.1.tif"
NWGdir = "/home/teblack/GreenlandTermini/NWGreenland_annual_termini"
vel_file = "/Volumes/insar5/teblack/data/MEASURES_annual_vel_mosaics/2018/greenland_vel_mosaic500_2018_vv_v01.tif"
gid_file = "/home/teblack/GreenlandTermini/glacier_info/glacier_ids/GlacierIDs.shp"

seaice_pts = pd.DataFrame({'id':[4, 3, 2, 1],
                           'label':['A', 'B', 'C', 'D'],
                           'Lat_m':[-1212474.517, -1462469.016, -1637465.192, -2062456.796],
                           'Lon_m':[-612487.261, -462490.139, -387491.693, -362492.549],
                           'Lat_deg':[77.508, 75.909, 74.557, 70.842],
                           'Lon_deg':[-71.801, -62.549, -58.314, -54.968]})

# %% Read files

# Read ocean mask file (~Greenland outline)
gimp = gpd.read_file(mask_file)
gimp.to_crs(epsg=3574)

# Read all terminus shapefiles into one geodataframe
files = os.listdir(NWGdir)
path = [os.path.join(NWGdir, i) for i in files
        if ".shp" in i and ".xml" not in i and "template" not in i]
termini = gpd.GeoDataFrame(pd.concat([gpd.read_file(i) for i in path],
                                     sort=False),
                           crs=gpd.read_file(path[0]).crs)
termini.to_crs(epsg=3574)
hydroyears = [d.year if d.month>=9 else d.year+1 for d in pd.to_datetime(termini.Date)]
termini['hydroyear'] = hydroyears

# Read glacier ID points
gids = gpd.GeoDataFrame(gpd.read_file(gid_file), crs=gpd.read_file(gid_file).crs)
gids = gids.where(gids.GlacierID <= 92)

# Read ice mask tif
iceR = rio.open(ice_tif)
ice = iceR.read_band(1)

# Read velocity mosaic
speedR = rio.open(vel_file)
speed = speedR.read_band(1)
#speed[speed < 10] = speed.max()

# %% PLOT MAP OF GREENLAND, TERMINI, SEA ICE POINTS, GLACIER IDS
fig, ax = plt.subplots(figsize=(46, 24))
fig.set_figheight(46)
ax.set_aspect('equal')

# Plot Greenland outline
gimp.plot(ax=ax, color='silver', linewidth=3)

# Plot speed?
#extent=[gimp.geometry.envelope.bounds.minx.min()-50000,
#        gimp.geometry.envelope.bounds.maxx.max()+50000,
#        gimp.geometry.envelope.bounds.miny.min()-50000,
#        gimp.geometry.envelope.bounds.maxy.max()+50000]
#plt.imshow(speed, extent=extent, cmap='Greys', norm=LogNorm(vmin=10, vmax=20000, clip=True), alpha=0.4)

# Plot termini on Greenland outline, colored by observation date
termini.plot(ax=ax, column='hydroyear', cmap='rainbow_r', linewidth=4, legend=True)

# Set map boundaries based on observed area
bounds = termini.geometry.bounds
bound_buffer = 100000
plt.xlim([bounds.minx.min()-bound_buffer, bounds.maxx.max()+bound_buffer])
plt.ylim([bounds.miny.min()-bound_buffer, bounds.maxy.max()+bound_buffer])
ax.axis('off')

# Plot sea ice observation points
COLORMAP=cm.Blues
colors = COLORMAP(np.linspace(1, 0.2, 4))
plt.scatter(seaice_pts.Lon_m.values, seaice_pts.Lat_m.values, s=2000, color=colors)

# Plot glacier IDs as labels, colored by group
#gids.plot(ax=ax)

# Plot separating lines and labels for the four glacier zones
plt.hlines(y=-1810000, xmin=-332000, xmax=-100000, color='gray', linestyles='dashed', linewidth=6)
plt.hlines(y=-1515000, xmin=-360000, xmax=-100000, color='gray', linestyles='dashed', linewidth=6)
plt.hlines(y=-1395000, xmin=-520000, xmax=-100000, color='gray', linestyles='dashed', linewidth=6)
plt.hlines(y=-1030000, xmin=-360000, xmax=-100000, color='gray', linestyles='dashed', linewidth=6)
#band1 = [-2400000, -1810000]
#band2 = [-1810000, -1515000]
#band3 = [-1515000, -1395000]
#band4 = [-1395000, -1030000]

plt.savefig('/home/teblack/plots/AGU_map.png', transparent=True, bbox_inches='tight')

# %% Plot sea ice pixel locations
#seaice_lon = [pt[0] for pt in seaice_pts]
#seaice_lat = [pt[1] for pt in seaice_pts]
#fig, ax = plt.subplots()
#ax.set_aspect('equal')
#gimp.plot(ax=ax, color='silver', linewidth=0.5)
#plt.scatter(seaice_lon, seaice_lat)

## %% Plot speed?
#fig, ax = plt.subplots()
#plt.imshow(speed, cmap='Greys_r', norm=LogNorm(vmin=0.1, vmax=13000))
