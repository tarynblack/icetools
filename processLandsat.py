#!/usr/bin/env python3
# Identify, download, reproject, and colorize Landsat scenes for a given set \
# of Landsat scene parameters
# Taryn Black, 29 July 2020

import argparse
import os
from subprocess import run, check_output


class SmartFormatter(argparse.HelpFormatter):
    """From https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text"""
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def pidDeconstruct(parser, args):
    """Helper function for known mode"""
    if 'pid=' in args.pid:
        pid = args.pid.replace('\'','').split('=')[1]
    args.pid = pid
    args.sensor = pid[0:4]
    args.level = pid[5:9]
    args.path = pid[10:13]
    args.row = pid[13:16]
    args.year = pid[17:21]
    args.month = pid[21:23]
    args.day = pid[23:25]
    args.collection = pid[35:37]
    args.tier = pid[38:40]
    return args


def pidConstruct(parser, args):
    """Helper function for search mode"""
    arg_dict = vars(args)
    for arg in arg_dict:
        argval = arg_dict[arg]
        if type(argval)==list:
            if len(argval)>2:
                default_val = parser.get_default(arg)[0]
                vars(args)[arg].remove(default_val)
                newargval = vars(args)[arg]
                vars(args)[arg] = str(newargval).translate(str.maketrans({
                    '[':'{', ']':'}', '\'':'', ' ':''
                }))
            elif len(argval)==2:
                default_val = parser.get_default(arg)[0]
                vars(args)[arg].remove(default_val)
                newargval = vars(args)[arg]
                vars(args)[arg] = str(newargval).translate(str.maketrans({
                    '[':'', ']':'', '\'':'', ' ':''
                }))
            elif len(argval)==1:
                vars(args)[arg] = str(argval).strip('[]').strip('\'')
            
    pid = '{}_{}_{}{}_{}{}{}_*_{}_{}'.format(
        args.sensor, args.level, args.path, args.row, args.year, args.month, \
        args.day, args.collection, args.tier)
    args.pid = pid
    return args


def strZPad2(input):
    output = '{:0>2}'.format(input)
    return output


def strZPad3(input):
    output = '{:0>3}'.format(input)
    return output


def get_args():
    """Get Landsat image search and processing parameters."""
    # TODO: figure out how to format range of options e.g. month -> {06..09}

    parser = argparse.ArgumentParser(description='Identify, download, and \
        reproject Landsat scenes from Google Cloud Platform.',
        formatter_class=SmartFormatter)
    
    parser.add_argument('--epsg',
        required=True,
        help='EPSG code for spatial reference system, see epsg.io for codes')
    parser.add_argument('--basedir', 
        default=str(os.getcwd()), 
        help='Path under which to store subdirectories of downloaded images.')
    
    subparsers = parser.add_subparsers(title='Processing mode',
        description='Define a single known product ID or search for images that\
             match a set of parameters',
        help='')

    parser_known = subparsers.add_parser('known', 
        description='Define a Landsat Product ID to download and process.\
            A typical PID looks like: \
            LC08_L1TP_069018_20190722_20190901_01_T1 \
            [SENSOR]_[LEVEL]_[PATHROW]_[DATE]_[PROCESSDATE]_[COLL]_[TIER]',
        help='known product ID mode')
    parser_known.set_defaults(func=pidDeconstruct)
    parser_known.add_argument('pid', 
        type=str, 
        help='USGS Product ID for scene')
    
    parser_search = subparsers.add_parser('search',
        description='Define search parameters to find matching Landsat images. \
            Leave blank or use \'?\' wildcard in place of characters to \
            include all options for a parameter. \
            A typical Landsat product ID looks like: \
            LC08_L1TP_069018_20190722_20190901_01_T1 \
            [SENSOR]_[LEVEL]_[PATHROW]_[DATE]_[PROCESSDATE]_[COLL]_[TIER]',
        formatter_class=SmartFormatter,
        help='parameter search mode')
    parser_search.set_defaults(func=pidConstruct)
    sensor_help = '''R|Four-character Landsat sensor identifier\n\
        LC08 : Landsat-8 Operational Land Imager (OLI) \n\
               and Thermal Infrared Sensor (TIRS)\n\
        LE07 : Landsat-7 Enhanced Thematic Mapper Plus (ETM+)\n\
        LT05 : Landsat-5 Thematic Mapper (TM)\n\
        LM05 : Landsat-5 Multispectral Scanner (MSS)\n\
        LT04 : Landsat-4 Thematic Mapper (TM)\n\
        LM04 : Landsat-4 Multispectral Scanner (MSS)\n\
        LM03 : Landsat-3 Multispectral Scanner (MSS)\n\
        LM02 : Landsat-2 Multispectral Scanner (MSS)\n\
        LM01 : Landsat-1 Multispectral Scanner (MSS)'''
    parser_search.add_argument('--sensor',
        metavar='[LC08, LE07, LT05, LT04, LM05, LM04, LM03, LM02, LM01]',
        choices=['????', 'LC08', 'LE07', 'LT05', 'LT04',
            'LM05', 'LM04', 'LM03', 'LM02', 'LM01'], 
        action='append',
        default=['????'],
        help=sensor_help)
    level_help = '''R|Landsat processing level, for georegistration accuracy\n\
        L1TP : Level-1 Precision and Terrain Correction\n\
        L1GT : Level-1 Systematic Terrain Correction\n\
        L1GS : Level-1 Systematic Correction\n\
        For additional information, see https://www.usgs.gov/land-resources/nli/landsat/landsat-levels-processing'''
    parser_search.add_argument('--level',
        metavar='[L1TP, L1GT, L1GS]',
        choices=['????', 'L1TP', 'L1GT', 'L1GS'],
        action='append',
        default=['????'],
        help=level_help)
    parser_search.add_argument('--path',
        metavar='[001-251]',
        type=strZPad3,
        choices=['???']+
                [range(1,251+1)]+
                ['{:0>3}'.format(p) for p in [*range(1,251+1)]],
        action='append',
        default=['???'],
        help='Three-character Landsat orbital path number')
    parser_search.add_argument('--row',
        metavar='[001-248]',
        type=strZPad3,
        choices=['???']+
                [range(1,248+1)]+
                ['{:0>3}'.format(r) for r in [*range(1,248+1)]],
        action='append',
        default=['???'],
        help='Three-character Landsat orbital row number')
    parser_search.add_argument('--year',
        action='append',
        default=['????'],
        help='Four-digit calendar year in which image was collected')
    parser_search.add_argument('--month',
        metavar='[01-12]',
        type=strZPad2,
        choices=['??']+
                [range(1,12+1)]+
                ['{:0>2}'.format(m) for m in [*range(1,12+1)]],
        action='append',
        default=['??'],
        help='Two-digit calendar month in which image was collected')
    parser_search.add_argument('--day',
        metavar='[01-31]',
        type=strZPad2,
        choices=['??']+
                [range(1,31+1)]+
                ['{:0>2}'.format(d) for d in [*range(1,31+1)]],
        action='append',
        default=['??'],
        help='Two-digit calendar day on which image was collected')
    parser_search.add_argument('--collection', 
        metavar='[01]',
        choices=['??', '01'],
        action='append',
        default=['??'],
        help='Management structure for Landsat archive')
    tier_help = '''R|Category of Landsat Level-1 products, based on data quality and processing level\n\
        T1 : Tier-1 (highest quality, suitable for time-series analysis)\n\
        T2 : Tier-2 (lower quality, do not meet geometric standards for T1)\n\
        RT : Real-Time (available until processed to T1 or T2)
        For additional information, see https://www.usgs.gov/land-resources/nli/landsat/landsat-collection-1'''
    parser_search.add_argument('--tier',
        metavar='[T1, T2, RT]',
        choices=['??', 'T1', 'T2', 'RT'],
        action='append',
        default=['??'],
        help=tier_help)
    
    args = parser.parse_args()
    args.func(parser_search, args)

    return args


def getGoogleCloudSearchPath(args):
    """Use the given parameters to create a product ID and Google Cloud bucket path, including wildcards, in order to search GCP for matching directories. If a product ID is given as a parameter, use that and ignore other Landsat parameters."""
    search_pid = args.pid
    bucket = 'gs://gcp-public-data-landsat/{}/{}/{}/{}/{}'.format(
        args.sensor, args.collection, args.path, args.row, search_pid)
    return search_pid, bucket


def getMatchingGCPDirectories(args):
    """List Google Cloud Platform directories that match search parameters and save to file."""
    search_pid, bucket = getGoogleCloudSearchPath(args)
    list_file = 'landsat_scenelist_{}.txt'.format(search_pid)
    list_file = list_file.replace(' ','')
    cmdstr = 'gsutil ls -d {}>{}/{}'.format(bucket, args.basedir, list_file)
    print('\nUsing gsutil to search directories that match:\n{}'.format(bucket))
    run(cmdstr, shell=True)
    return list_file


def chooseSatelliteBands(sensor):
    """Choose which bands to download, depending on the sensor."""   
    rgb_bands = {
        'LC08': [4, 3, 2, 8], # Landsat-8: 4-red, 3-green, 2-blue, 8-pan
        'LE07': [3, 2, 1, 8], # Landsat-7: 3-red, 2-green, 1-blue, 8-pan
        'LT05': [3, 2, 1],    # Landsat-4,5 TM: 3-red, 2-green, 1-blue
        'LT04': [3, 2, 1],    # Landsat-4,5 TM: 3-red, 2-green, 1-blue
        'LM05': [3, 2, 1],    # Landsat-4,5 MSS: 3-NIR, 2-red, 1-green
        'LM04': [3, 2, 1],    # Landsat-4,5 MSS: 3-NIR, 2-red, 1-green
        'LM03': [6, 5, 4],    # Landsat-1,2,3 MSS: 6-NIR, 5-red, 4-green
        'LM02': [6, 5, 4],    # Landsat-1,2,3 MSS: 6-NIR, 5-red, 4-green
        'LM01': [6, 5, 4]     # Landsat-1,2,3 MSS: 6-NIR, 5-red, 4-green
        }
    bands = rgb_bands.get(sensor)
    return bands


def checkDirectoryExists(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)


def reprojectScene(epsg, base_dir, band_file, dest_dir):
    """Use GDAL to reproject Landsat scene.

    These commands are copied from Ian Joughin (getgoogle.py).
    """
    proj = str(check_output(
        'gdalsrsinfo -o proj4 "EPSG:{}"'.format(epsg), shell=True), 'utf-8'
               ).strip()
    proj_wkt = str(check_output(
        'gdalsrsinfo -o wkt "EPSG:{}"'.format(epsg), shell=True), 'utf-8'
                  ).strip()
    proj_name = proj_wkt.split(',')[0].split('[')[1].strip('"')
    print('--> Reprojecting with gdalwarp to EPSG:{} - {}'.format(
        epsg, proj_name))
    checkDirectoryExists('{}/temp'.format(base_dir))
    cmdstr = 'gdalwarp \
        -srcnodata 0 \
        -tr 15.0 15.0 \
        -dstnodata 0 \
        -co "COMPRESS=DEFLATE" \
        -co "TILED=YES" \
        -co "ZLEVEL=9" \
        -co "PREDICTOR=2" \
        -co "BLOCKXSIZE=512" \
        -co "BLOCKYSIZE=512" \
        -r cubic \
        -t_srs \'{}\' \
        {}/temp/{} {}/{}'.format(proj, base_dir, band_file, dest_dir, band_file)
    run(cmdstr, shell=True)
    small_file = 'small_' + band_file
    checkDirectoryExists('{}/small'.format(dest_dir))
    cmdstr = 'gdalwarp \
        -srcnodata 0 \
        -tr 100.0 100.0 \
        -dstnodata 0 \
        -co "COMPRESS=DEFLATE" \
        -co "TILED=YES" \
        -co "ZLEVEL=9" \
        -co "PREDICTOR=2" \
        -co "BLOCKXSIZE=512" \
        -co "BLOCKYSIZE=512" \
        -r near \
        -t_srs \'{}\' \
        {}/temp/{} \
        {}/small/{} ; \
        gdaladdo \
        -r average \
        {}/small/{} \
        2 4 8 16 ; '.format(
            proj, base_dir, band_file, dest_dir, small_file, \
            dest_dir, small_file)
    run(cmdstr, shell=True)


def downloadScene(bands, pid, epsg, base_dir, gcp_path, dest_dir):
    """Download Landsat scene bands and metadata from Google Cloud Platform.

    The scene is downloaded to a temp directory, then reprojected and saved to the destination directory, which is named for the image product ID.
    """
    checkDirectoryExists('{}/temp'.format(base_dir))
    for band in bands:
        print('--- Band {}:'.format(band))
        band_file = '{}_B{}.TIF'.format(pid, band)
        command = 'gsutil cp {}{} {}/temp/'.format(
            gcp_path, band_file, base_dir)
        run(command, shell=True)
        reprojectScene(epsg, base_dir, band_file, dest_dir)
    command = 'gsutil cp {}{}_MTL.txt {} ; \
        rm {}/temp/*'.format(gcp_path, pid, dest_dir, base_dir)
    run(command, shell=True)


def getBandFilePaths(base_dir, pid, bands):
    """Get file paths for red, green, blue, and panchromatic band files"""
    red   = '{}/{}/{}_B{}.TIF'.format(base_dir, pid, pid, bands[0])
    green = '{}/{}/{}_B{}.TIF'.format(base_dir, pid, pid, bands[1])
    blue  = '{}/{}/{}_B{}.TIF'.format(base_dir, pid, pid, bands[2])
    if len(bands) == 4:
        pan = '{}/{}/{}_B{}.TIF'.format(base_dir, pid, pid, bands[3])
    else:
        pan = ''
    return red, green, blue, pan


def createColorComposite(base_dir, pid, bands):
    """Use RGB bands to create a true-color or near-true-color composite image, 
    and pansharpen if the panchromatic band is available"""
    
    red, green, blue, pan = getBandFilePaths(base_dir, pid, bands)
    
    if len(bands) == 4:
        composite_TIF = '{}/{}_pan-composite.TIF'.format(base_dir, pid)
        print('Panchromatic band present. \
            \nCompositing bands {}, {}, {}; pansharpening band {}...'.format(
                bands[0], bands[1], bands[2], bands[3]))
        cmd = 'gdal_pansharpen.py {} {} {} {} {}'.format(
            pan, red, green, blue, composite_TIF)
        # call(cmd, shell=True)
    else:
        composite_TIF = '{}/composite/{}_composite.TIF'.format(base_dir, pid)
        composite_VRT = '{}/composite/{}_composite.vrt'.format(base_dir, pid)
        print('No panchromatic band present. \
            \nCompositing bands {}, {}, and {}...'.format(
                bands[0], bands[1], bands[2]))
        cmd = 'gdalbuildvrt -separate {} {} {} {}; \
            gdal_translate {} {}'.format(
            composite_VRT, red, green, blue, composite_VRT, composite_TIF)
    run(cmd, shell=True)
    
    if os.path.exists(composite_TIF):
        print('-> Created {}'.format(composite_TIF))
    else:
        print('-> Failed to create {}'.format(composite_TIF))


def main():
    args = get_args()
    list_file = getMatchingGCPDirectories(args)
    if os.path.exists('{}/{}'.format(args.basedir, list_file)):
        print('-> Scene list saved to {}/{}'.format(args.basedir, list_file))
    else:
        print('-> Failed to create {}/{}'.format(args.basedir, list_file))
    
    with open('{}/{}'.format(args.basedir, list_file)) as scene_list:
        total_images = sum(1 for _ in scene_list)

    scene_list = open('{}/{}'.format(args.basedir, list_file))
    for counter, gcp_path in enumerate(scene_list, start=1):
        gcp_path = gcp_path.strip()
        pid = gcp_path.split('/')[-2]
        print('\nAccessing PID:  {}'.format(pid))
        print('Image {} of {}'.format(counter, total_images))
        
        dest_dir = '{}/{}'.format(args.basedir, pid)
        checkDirectoryExists(dest_dir)
        
        sensor = pid[0:4]
        bands = chooseSatelliteBands(sensor)
        
        print('--> Downloading scene to:  {}'.format(dest_dir))
        downloadScene(bands, pid, args.epsg, args.basedir, gcp_path, \
            dest_dir)
        
        print('--> Creating color composite in:  {}'.format(args.basedir))
        createColorComposite(args.basedir, pid, bands)
    scene_list.close()

main()
