# Functions to evaluate geodatabase attributes based on Landsat Product ID
# T Black, 24 July 2020

from datetime import datetime

null_value = str(-999)

landsat_ids = ['LC08', 'LE07', 'LT05', 'LT04', \
               'LM05', 'LM04', 'LM03', 'LM02', 'LM01']

def getSensor(id):
    sensor = id[0:4]
    return sensor

def getTileCoordinate(id):
    sensor = getSensor(id)
    if sensor in landsat_ids:
        tile_coordinate = id[10:16]
    else:
        tile_coordinate = null_value
    return tile_coordinate

def getReferenceSystem(id):
    sensor = getSensor(id)
    if sensor in landsat_ids:
        sensor_number = int(sensor[2:4])
        if sensor_number > 3:
            reference_system = 'WRS-2'
        elif sensor_number <= 3:
            reference_system = 'WRS-1'
    else:
        reference_system = null_value
    return reference_system

def getSourceDate(id):
    sensor = getSensor(id)
    if sensor in landsat_ids:
        year = id[17:21]
        month = id[21:23]
        day = id[23:25]
        source_date = year + '-' + month + '-' + day
    else:
        source_date = null_value
    return source_date

def getCircadianDate(id):
    date = getSourceDate(id)
    if date:
        year = int(date[0:4])
        month = int(date[5:7])
        day = int(date[8:10])
        day_of_year = datetime(year, month, day).timetuple().tm_yday
    else:
        day_of_year = null_value
    return day_of_year

def getYear(id):
    date = getSourceDate(id)
    if date:
        year = date[0:4]
    else:
        year = null_value
    return year

def getSeason(id):
    date = getSourceDate(id)
    if date:
        month = int(date[5:7])
        day = int(date[8:10])
        if month in [4, 5, 6]:
            season = 'SPR'
        elif month in [7]:
            season = 'SUM'
        elif month in [8]:
            if day < 17:
                season = 'SUM'
            else:
                season = 'AUT'
        elif month in [9]:
            season = 'AUT'
        elif month in [10]:
            if day < 17:
                season = 'AUT'
            else:
                season = 'WIN'
        elif month in [11, 12, 1, 2, 3]:
            season = 'WIN'
    else:
        season = null_value
    return season