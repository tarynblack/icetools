# Take in shapefile of overlapping polygons with Radarsat swath attributes and return a new file with the polygon overlaps removed
# Taryn E. Black, UW ESS/APL-PSC, teblack@uw.edu
# Created 9 November 2020

import argparse
import os
import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPolygon


def get_args():
    """Usage statement."""
    parser = argparse.ArgumentParser(description='Remove overlaps from Radarsat swath polygons in swath mosaic shapefile.')
    parser.add_argument('--basedir', help='Path to look for swath mosaic directories.')
    parser.add_argument('--mosaic', help='Input Radarsat swath mosaic file')
    parser.add_argument('--imageid', help='Radarsat-1 image ID, i.e. \'RAD1_NSIDC_0633_VX_RRR_YYYY_YYYY\'')
    parser.add_argument('--out', help='Path to store output file', default=str(os.getcwd()))
    args = parser.parse_args()
    return args

def getRadarsatImageID(f_mosaic):
    f_mosaic = f_mosaic.split('/')[-1]
    version = f_mosaic.split('_')[-1].split('.')[0].upper()
    resolution = f_mosaic.split('_')[0].split('m')[0]
    years = f_mosaic.split('_')[0].split('band')[-1]
    radarsatID = 'RAD1_NSIDC_0633_{}_{}_{}'.format(version, resolution, years)
    return radarsatID

def getRadarsatMosaicFromImageID(image_id, basedir=None):
    version = image_id.split('_')[-3].lower()
    resolution = image_id.split('_')[-2]
    years = image_id.split('_')[-1]
    startyear = '20{}'.format(years[0:2])
    endyear = '20{}'.format(years[2:4])
    dirname = '{}_{}'.format(startyear, endyear)
    f_mosaic = '{}/{}mCband{}_{}.1.shp'.format(
        dirname, resolution, years, version)
    if basedir is not None:
        f_mosaic = '{}/{}'.format(basedir, f_mosaic)
    return f_mosaic

def getRadarsatTileCoordinates(swath):
    tile_coordinate = "RAD1_{:0>3}_{:0>5}".format(swath.Track, swath.Orbit)
    return tile_coordinate

def preProcessMosaic(f_mosaic):
    """Load mosaic and pre-process attributes."""
    mosaic = gpd.read_file(f_mosaic).to_crs(epsg=3413)
    # Save original order of mosaic layers as a column
    mosaic.reset_index(inplace=True)
    # Convert attributes to non-object dtypes
    mosaic['Track'] = mosaic['Track'].astype('int64')
    mosaic['Date'] = mosaic['Date'].astype('str')
    mosaic['Orbit'] = mosaic['Orbit'].astype('int64')
    # Add columns for Image ID and Tile Coordinate
    mosaic['Image_ID'] = getRadarsatImageID(f_mosaic)
    mosaic['Image_Tile_Coordinate'] = pd.Series(dtype='str')
    for i in range(len(mosaic)):
        mosaic['Image_Tile_Coordinate'].iloc[i] = getRadarsatTileCoordinates(
            mosaic.iloc[i])
    # Drop columns for Track, Orbit, SAT (no longer needed)
    mosaic.drop(['Track', 'Orbit', 'Sensor'], axis=1, inplace=True)
    return mosaic

def dissolveSwaths(mosaic):
    """Dissolve multiple polygons into single polygons based on Track."""
    mosaic = mosaic.dissolve(by="Image_Tile_Coordinate", as_index=False)
    return mosaic

def eraseLowerSwaths(mosaic):
    # Reverse sort on original index so that polygons are ordered bottom to top
    mosaic.sort_values('index', ascending=False, inplace=True)
    mosaic['geometry'] = mosaic.buffer(0.01)

    # Loop through each polygon and erase underlying polygons where they overlap
    for p in range(0, len(mosaic)):
        mosaic['geometry'].iloc[p+1:] = mosaic['geometry'].iloc[p+1:].difference(
            mosaic['geometry'].iloc[p])
    mosaic.drop('index', axis=1, inplace=True)
    return mosaic

def eraseLineStrings(mosaic):
    for shape in mosaic.index:
        if mosaic.loc[shape].geometry.type == 'GeometryCollection':
            polys = [g for g in mosaic.loc[shape].geometry if g.type != 'LineString']
            polygons = MultiPolygon(polys)
            mosaic['geometry'].loc[shape] = polygons
    return mosaic


def main():
    args = get_args()

    if args.imageid:
        radarsatid = args.imageid
        if args.basedir:
            f_mosaic = getRadarsatMosaicFromImageID(
                args.imageid, args.basedir)
        else:
            f_mosaic = getRadarsatMosaicFromImageID(args.imageid)
        if not os.path.exists(f_mosaic):
            print(
                'Oh no! Swath mosaic file {} does not exist'.format(
                    f_mosaic))
            mosaic = None
        else:
            mosaic = preProcessMosaic(f_mosaic)
    elif args.mosaic:
        radarsatid = getRadarsatImageID(args.mosaic)
        if not os.path.exists(args.mosaic):
            print(
                'Oh no! Swath mosaic file {} does not exist'.format(
                    args.mosaic))
            mosaic = None
        else:
            f_mosaic = args.mosaic.split('/')[-1]
            mosaic = preProcessMosaic(args.mosaic)
    else:
        print('Must specify either imageid or mosaic.')
        mosaic = None
    
    if mosaic is not None:        
        mosaic = dissolveSwaths(mosaic)
        mosaic = eraseLowerSwaths(mosaic)
        mosaic = eraseLineStrings(mosaic)
        
        fout = '{}/{}.shp'.format(args.out, radarsatid)
        mosaic.to_file(fout)
    # else: mosaic=None
    
    # return mosaic

main()