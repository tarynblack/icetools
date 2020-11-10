# Take in shapefile of overlapping polygons with Sentinel swath attributes and return a new file with the polygon overlaps removed
# Taryn E. Black, UW ESS/APL-PSC, teblack@uw.edu
# Created 28 October 2020

import argparse
import os
import geopandas as gpd
import pandas as pd

def get_args():
    """Usage statement."""
    parser = argparse.ArgumentParser(description='Remove overlaps from Sentinel swath polygons in swath mosaic shapefile.')
    parser.add_argument('--basedir', help='Path to look for swath mosaic directories.')
    parser.add_argument('--mosaic', help='Input Sentinel swath mosaic file')
    parser.add_argument('--imageid', help='Sentinel-1 image ID, i.e. \'SEN1_NSIDC_0723_VX_YYYYMMDD_YYYYMMDD\'')
    parser.add_argument('--version', help='Sentinel mosaic NSIDC version number', default=2)
    parser.add_argument('--out', help='Path to store output file', default=str(os.getcwd()))
    args = parser.parse_args()
    return args

def getSentinelImageID(f_mosaic, version):
    dates = f_mosaic.split('/')[-1].split('.')[0].replace('-','')
    startdate = dates[4:8] + dates[0:4]
    enddate = dates[12:16] + dates[8:12]
    daterange = '{}_{}'.format(startdate, enddate)
    sentinelID = 'SEN1_NSIDC_0723_V{}_{}'.format(version, daterange)
    return sentinelID

def getMosaicFromImageID(image_id, basedir=None):
    startdate, enddate = image_id.split('_')[-2:]
    dirnameYMD = '{}-{}-{}'.format(
        startdate[0:4], startdate[4:6], startdate[6:8])
    startdateMDY = '{}-{}-{}'.format(
        startdate[4:6], startdate[6:8], startdate[0:4])
    enddateMDY = '{}-{}-{}'.format(
        enddate[4:6], enddate[6:8], enddate[0:4])
    f_mosaic = '{}/{}-{}.shp'.format(
        dirnameYMD, startdateMDY, enddateMDY)
    if basedir is not None:
        f_mosaic = '{}/{}'.format(basedir, f_mosaic)
    return f_mosaic

def getSentinelTileCoordinates(swath):
    tile_coordinate = "{}_{:0>3}_{:0>6}".format(swath.SAT, swath.TRACK, swath.ORBIT)
    return tile_coordinate

def preProcessMosaic(f_mosaic, version):
    """Load mosaic and pre-process attributes."""
    mosaic = gpd.read_file(f_mosaic).to_crs(epsg=3413)
    # Save original order of mosaic layers as a column
    mosaic.reset_index(inplace=True)
    # Convert attributes to non-object dtypes
    mosaic['TRACK'] = mosaic['TRACK'].astype('int64')
    mosaic['DATE'] = mosaic['DATE'].astype('str')
    mosaic['ORBIT'] = mosaic['ORBIT'].astype('int64')
    mosaic['SAT'] = mosaic['SAT'].astype('str')
    # Add columns for Image ID and Tile Coordinate
    mosaic['Image_ID'] = getSentinelImageID(f_mosaic, version)
    mosaic['Image_Tile_Coordinate'] = pd.Series(dtype='str')
    for i in range(len(mosaic)):
        mosaic['Image_Tile_Coordinate'].iloc[i] = getSentinelTileCoordinates(
            mosaic.iloc[i])
    # Drop columns for TRACK, ORBIT, SAT (no longer needed)
    mosaic.drop(['TRACK', 'ORBIT', 'SAT'], axis=1, inplace=True)
    return mosaic

def dissolveSwaths(mosaic):
    """Dissolve multiple polygons into single polygons based on TRACK."""
    mosaic = mosaic.dissolve(by="Image_Tile_Coordinate", as_index=False)
    return mosaic

def eraseLowerSwaths(mosaic):
    # Reverse sort on original index so that polygons are ordered bottom to top
    mosaic.sort_values('index', ascending=False, inplace=True)

    # Loop through each polygon and erase underlying polygons where they overlap
    for p in range(0, len(mosaic)):
        mosaic['geometry'].iloc[p+1:] = mosaic['geometry'].iloc[p+1:].difference(
            mosaic['geometry'].iloc[p])
    mosaic.drop('index', axis=1, inplace=True)
    return mosaic

def main():
    args = get_args()

    if args.imageid:
        if args.basedir:
            f_mosaic = getMosaicFromImageID(
                args.imageid, args.basedir)
        else:
            f_mosaic = getMosaicFromImageID(args.imageid)
        if not os.path.exists(f_mosaic):
            print(
                'Oh no! Swath mosaic file {} does not exist'.format(
                    f_mosaic))
            mosaic = None
        else:
            mosaic = preProcessMosaic(f_mosaic, args.version)
    elif args.mosaic:
        if not os.path.exists(args.mosaic):
            print(
                'Oh no! Swath mosaic file {} does not exist'.format(
                    f_mosaic))
            mosaic = None
        else:
            mosaic = preProcessMosaic(args.mosaic, args.version)
    else:
        print('Must specify either imageid or mosaic.')
        mosaic = None
    
    if mosaic is not None:        
        mosaic = dissolveSwaths(mosaic)
        mosaic = eraseLowerSwaths(mosaic)
        
        sentinelid = getSentinelImageID(args.mosaic, args.version)
        fout = '{}/{}.shp'.format(args.out, sentinelid)
        mosaic.to_file(fout)
    # else: mosaic=None
    
    # return mosaic

main()