# Functions to handle Sentinel-1 mosaic attributes
# Taryn E. Black, UW ESS/APL-PSC, teblack@uw.edu
# Created 27 October 2020

from os.path import exists
from subprocess import run

import geopandas as gpd
import pandas as pd

# Need termini file, and out path
f_termini = '/mnt/e/greenland-termini/Greenland_Glaciers.gdb'
termini_layer = 'glacier_termini'
mos_in = '/mnt/s'
mos_out = '/mnt/e/greenland-termini/swath-attributes'


def preProcessTermini(f_termini):
    """Load termini and remove duplicated rows."""
    termini = gpd.read_file(f_termini, layer=termini_layer).to_crs(epsg=3413)
    # Remove duplicated termini rows
    termini.drop_duplicates(keep='first', \
        subset=['Glacier_ID', 'Image_ID', 'Create_User'], inplace=True)
    return termini

def getMosaicName(imageid):
    mosaic_name = '{}.shp'.format(imageid)
    return mosaic_name

def getMosaicFile(imageid, fin, fout):
    mosaic_name = getMosaicName(imageid)
    f_mosaic = '{}/{}'.format(fout, mosaic_name)
    if not exists(f_mosaic):
        command = 'python3 getSentinelSwathAttributes.py --basedir=\'{}\' --imageid=\'{}\' --version=2 --out={}'.format(
                fin, imageid, fout)
        run(command)
    if exists(f_mosaic):
        mosaic = gpd.read_file(f_mosaic)
    return mosaic

def joinToTermini(mosaic, termini):
    """Join mosaic metadata to glacier terminus metadata based on location. Reset index in case duplicate termini are created (which duplicates the index)."""
    termini = gpd.sjoin(termini, mosaic, how='left')
    termini.reset_index(inplace=True)
    termini.rename(columns={'Image_ID_left': 'Image_ID'}, inplace=True)
    return termini

def removeDuplicateTermini(mosaic, termini):
    """Termini that cross >1 swath show up twice in dataframe after spatial join. Find these and remove duplicates."""
    # Find duplicates
    duplicates = termini[termini.Glacier_ID.duplicated(keep=False)]
    if not duplicates.empty:
        for id in duplicates.Glacier_ID.unique():
            # Intersect duplicates with mosaic to get segments over each swath
            duplicates_id = duplicates[duplicates.Glacier_ID == id]
            duplicates_sections = gpd.overlay(duplicates_id, mosaic, 'intersection')
            # Identify the swath track containing the largest section of the terminus
            keep_track = duplicates_sections[duplicates_sections.length == \
                duplicates_sections.length.max()].Image_Tile_2.iloc[0]
            # Keep only the terminus associated with the larger section's track.
            drop_index = duplicates_id[duplicates_id.Image_Tile != keep_track].index
            termini.drop(drop_index, inplace=True)
        # Drop duplicates from original dataframe and add back the one to keep
        # termini.drop_duplicates(subset='Glacier_ID', keep=False, inplace=True)
        # termini = termini.append(keep_duplicate)
    return termini

def setTerminiAttributes(termini):
    """Add new attributes based on mosaic join with duplicates removed."""
    termini.Image_Tile_Coordinate = termini.Image_Tile
    termini.Source_Date = termini.DATE
    termini.Image_Reference_System = 'STO'
    return termini

def cleanResults(termini, matching_termini):
    matching_termini.set_index('index', inplace=True)
    for i in matching_termini.index:
        termini.at[i, 'Image_Tile_Coordinate'] = matching_termini.at[i, \
            'Image_Tile_Coordinate']
        termini.at[i, 'Image_Reference_System'] = matching_termini.at[i, \
            'Image_Reference_System']
        termini.at[i, 'Source_Date'] = matching_termini.at[i, 'Source_Date']
        termini = termini.astype({
            'Glacier_ID'             : 'uint32',
            'Quality_Flag'           : 'string',
            'Image_ID'               : 'string',
            'Sensor'                 : 'string',
            'Image_Tile_Coordinate'  : 'string',
            'Image_Reference_System' : 'string',
            'Source_Date'            : 'string',
            'Create_User'            : 'string',
            'Create_Date'            : 'datetime64',
            'Digitization_Method'    : 'string',
            'Edit_User'              : 'string',
            'Edit_Date'              : 'datetime64',
            'Comment'                : 'string',
            'SHAPE_Length'           : 'float64',
            'geometry'               : 'geometry'})
    return termini


termini = preProcessTermini(f_termini)
list_imageids = termini.Image_ID.unique()
list_mosaics = [m for m in list_imageids if m is not None and 'SEN1' in m]
for mos_id in list_mosaics:
    mosaic = getMosaicFile(mos_id, mos_in, mos_out)
    
    matching_termini = termini[termini.Image_ID == mos_id]
    matching_termini = joinToTermini(mosaic, matching_termini)
    matching_termini = removeDuplicateTermini(mosaic, matching_termini)
    matching_termini = setTerminiAttributes(matching_termini)
    termini = cleanResults(termini, matching_termini)

termini.to_file(f_termini, driver='FileGDB', layer=termini_layer)
