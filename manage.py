# Classes and functions for handling intake and management of glacier terminus data.

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, LineString


class Glacier:
    def __init__(self, gid):
        self.gid = gid
        self.refline = LineString()
        self.refbox = LineString()
        self.centerline = LineString()
        self.officialname = ""
        self.greenlandicname = ""
        self.alternativename = ""
        self.unofficialname = ""
        self.obsseries = []
        self.areas = []
        self.widths = []
        self.termareas = []
        self.lengths = []
        self.dates = []
        self.seasons = []
        self.hydroyears = []
        self.missingyears = []
        self.interpyears = []
        self.interpareas = []
        self.interptermareas = []
        self.interplengths = []

    def add_observation(self, observation):
        if observation.gid != self.gid:
            print('Cannot add glacier %s observation to glacier %s observation series' % (observation.gid, self.gid))
        self.obsseries.append(observation)

    def sort_by_date(self):
        self.obsseries = sorted(self.obsseries, key=lambda k: k.date)
    
    def extract(self, attr):
        """Extract list of a given attribute for each TerminusObservation in the Glacier observation series"""
        data_list = eval("[obs.%s for obs in self.obsseries]" % attr)
        data_list = pd.Series(data_list)
        return data_list
      
    def getMissingYears(self, year_list):
        """Identify years with missing data for an attribute"""
        observed_hydroyears = self.extract('hydroyear')
        missing_years = list(set(year_list) - set(observed_hydroyears))
        return missing_years
    
    def interpolateMeasurements(self, attr, year_list):
        """Interpolate values between observations of area or length. Do not interpolate prior to first observation."""
        observations = pd.DataFrame(data={attr: self.extract(attr).values},\
                                    index=self.extract('hydroyear'))
        missing_years = self.getMissingYears(year_list)
        missing_data = pd.DataFrame(data={attr: None}, \
                                    index=missing_years)
        interpolated_data = observations.append(missing_data, sort=True).sort_index()
        interpolated_data = interpolated_data.interpolate(method='linear', \
            limit_direction='forward')
        return interpolated_data

    def getInterpolatedYears(self, attr, year_list):
        """Identify years in which data have been interpolated."""
        interpolated_data = self.interpolateMeasurements(attr, year_list)
        interpolated_years_index = interpolated_data.dropna().index
        observed_years = self.hydroyears
        interpolated_years = list(
            set(interpolated_years_index) - set(observed_years))
        return interpolated_years

class TerminusObservation:
    def __init__(self, gid, qflag, termination, imageid, sensor, date, circadiandate, year, geometry):
        # Attributes that must be defined on instantiation
        self.gid = gid
        self.qflag = qflag
        self.termination = termination
        self.imageid = imageid
        self.sensor = sensor
        self.date = pd.to_datetime(date)
        self.circadiandate = circadiandate
        self.year = year
        self.season = getSeason(self.date)
        self.geometry = geometry
        # Attributes that are determined from initial instance attributes
        self.hydroyear = hydrologicalYear(self.date)
        self.dayofhydroyear = dayOfHydroyear(self.date)
        # Attributes that are calculated elsewhere...
        self.area = 0.0
        self.width = 0.0
        self.length = 0.0
        self.termarea = 0.0
        self.centerlineintersection = Point()


def shp2gdf(file, epsg=3574):
    """Reads a shapefile of glacier data and reprojects to a specified EPSG (default is EPSG:3574 - WGS 84 / North Pole LAEA Atlantic, unit=meters)
    Result is a geodataframe containing the shapefile data.
    Also reindex the gdf so that index=GlacierID for direct selection by ID."""
    gdf = gpd.read_file(file)
    gdf = gdf.to_crs(epsg=epsg)
    # TODO: check whether glacier has multiple entries
    # (then will have multiple indices of that value, which screws up indexing)
    gdf = gdf.sort_values(by='GlacierID').set_index('GlacierID', drop=False)
    return gdf


def glacierInfo(termini_gdf, box_gdf, gid):
    """Get terminus and metadata for a given glacier ID in the geodataframe of
    termini data, as well as the glacier's reference box."""
    terminus = termini_gdf.loc[gid]
    box = box_gdf.loc[gid]
    return terminus, box


def hydrologicalYear(date):
    """Determine hydrological year of a given date. Greenland hydrological year
    is defined as September 1 through August 31. Returns the starting year of
    the hydrological year (aligned with September)."""
    date = pd.to_datetime(date)
    if pd.notnull(date):
        if date.month >= 9:
            hydroyear = date.year
        elif date.month < 9:
            hydroyear = date.year - 1
        return hydroyear


def dayOfHydroyear(date):
    """Convert date to number of days since start of hydrological year,
    defined as September 1. Analogous to day-of-year."""
    date = pd.to_datetime(date)
    if pd.notnull(date):
        hydroyear = hydrologicalYear(date)
        start_hydroyear = pd.to_datetime('%s-09-01' % str(hydroyear))
        day_of_hydroyear = (date - start_hydroyear).days
        return day_of_hydroyear


def getSeason(date):
    """Determine Northern Hemisphere season based on date."""
    date = pd.to_datetime(date)
    month = date.month
    if month in [3, 4, 5]:
        season = 'spring'
    elif month in [6, 7, 8]:
        season = 'summer'
    elif month in [9, 10, 11]:
        season = 'autumn'
    elif month in [12, 1, 2]:
        season = 'winter'
    return season

