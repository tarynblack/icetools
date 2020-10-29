#!/usr/bin/env python3
# Calculate metrics for glacier change
# Taryn Black, August 2020

import pandas as pd
from shapely import ops
import gplots as gp


def glacierArea(terminus, box):
    """Calculate the area of the polygon created by the intersection of a glacier terminus and the glacier's reference box. The default area for EPSG:3574 is in m2; this script returns in km2."""
    outline = terminus.geometry.union(box)
    gpoly = ops.polygonize_full(outline)[0]
    if gpoly.is_empty:
        print("%s: Glacier %s trace and box do not overlap" % (terminus.date, terminus.gid))
    area_km = gpoly.area / 10**6
    return area_km


def centerlineIntersection(terminus, centerline):
    """Locate intersection between a glacier outline observation and glacier centerline geometry. Split the centerline at the intersection and calculate the length of the substring. Return intersection point and substring length."""
    intersection_point = centerline.intersection(terminus.geometry)
    split_centerline = ops.split(centerline, terminus.geometry)
    substring_length = split_centerline[0].length
    return intersection_point, substring_length


def terminusWidth(terminus, box):
    """Compute box width at points where terminus intersects box, and return average width in km."""
    terminus = terminus.geometry
    if box.geom_type == 'MultiLineString':
        box = ops.linemerge(ops.MultiLineString(box))
    else:
        box = ops.LineString(box)
    # Split box into two halves
    half1 = ops.substring(box, 0, 0.5, normalized=True)
    half2 = ops.substring(box, 0.5, 1, normalized=True)
    # Get terminus intersections with box halves
    tx1 = terminus.intersection(half1)
    tx2 = terminus.intersection(half2)
    # Get distance from intersection point to other half
    tx1_dist = tx1.distance(half2)
    tx2_dist = tx2.distance(half1)
    # Get average distance across box, i.e. average terminus width
    average_width = (tx1_dist + tx2_dist) / 2 / 10**3
    return average_width


def addDecade(start_date):
    start_date = pd.to_datetime(start_date).date()
    end_date = pd.date_range(start_date, periods=2, freq='9Y')[-1]
    end_date = end_date.date()
    return end_date


def filterByDates(glacier, measure, date_start, date_end):
    measures = glacier.extract(measure)
    dates = glacier.extract('date')

    # Filter to data between selected dates
    if date_start:
        date_start = pd.to_datetime(date_start)
        measures = measures.where(dates >= date_start).dropna()
        dates = dates.where(dates >= date_start).dropna()
    if date_end: 
        date_end = pd.to_datetime(date_end)
        measures = measures.where(dates <= date_end).dropna()
        dates = dates.where(dates <= date_end).dropna()
    return measures, dates


def netMeasureChange(glacier, measure, date_start=None, date_end=None):
    
    measures, dates = filterByDates(glacier, measure, date_start, date_end)
    num_obs = len(measures)

    if num_obs == 0:
        print('No observations for {} (#{}) between {} and {}'.format(
            gp.getGlacierName(glacier), glacier.gid, date_start, date_end))
        measures = pd.Series(0.0)
        net_change_dates = (date_start, date_end)
    
    else:
        if measures.index[0] != 0:
        # In order to calculate area difference in time subsets, there needs \
        # to be one area measurement prior to the beginning of the subset. \
        # Then e.g. if there is only one measurement in the subset, the area \
        # change is relative to the previous measurement.
            newindex = measures.index[0] - 1
            prev_measure = glacier.extract(measure).loc[newindex]
            measures.loc[newindex] = prev_measure
            measures.sort_index(inplace=True)
            prev_date = glacier.extract('date').loc[newindex]
            dates.loc[newindex] = prev_date
            dates.sort_index(inplace=True)
        net_change_dates = (dates.iloc[0].date(), dates.iloc[-1].date())

    # Calculate inter-measurement, cumulative, and net change
    measure_change = measures.diff()
    cumul_measure_change = measure_change.cumsum()
    cumul_measure_change.iloc[0] = 0.0
    return cumul_measure_change, net_change_dates, num_obs


def rateMeasureChange(glacier, measure, date_start=None, date_end=None):
    # Get net area/length change and the dates over which it's calculated
    cumul_measure_change, net_change_dates, _ = netMeasureChange(
        glacier, measure, date_start, date_end)
    net_measure_change = cumul_measure_change.iloc[-1]
    date_first = net_change_dates[0]
    date_last = net_change_dates[1]

    # Get length of time range in years (SMS=semi-month start frequency \
    # (1st and 15th), slightly more precise than just months)
    years_diff = len(pd.date_range(date_first, date_last, freq='SMS'))/24

    # Average rate of area change per year = net change / years
    rate_measure_change = net_measure_change / years_diff
    return rate_measure_change
    

def normMeasureChange(glacier, measure, date_start=None, date_end=None):
    # Calculate normalized area/length change over a time period.
    # "Normalized" such that:
    # 0 = the smallest extent of the glacier in the time period
    # 1 = the largest extent of the glacier in the time period
    # all other values scaled linearly between those

    cumul_measure_change, _, _ = netMeasureChange(glacier, measure, \
        date_start, date_end)
    
    max_measure = cumul_measure_change.max()
    min_measure = cumul_measure_change.min()
    scaled_measure = (cumul_measure_change - min_measure) / \
        (max_measure - min_measure)
    return scaled_measure


def firstFullYear(all_glaciers):
    """Identify the first year in which all glaciers have either an observed or an interpolated data point."""
    first_full_year = 0
    for g in all_glaciers:
        glacier = all_glaciers[g]
        g_firstyear = glacier.interpareas.dropna().index[0]
        if g_firstyear > first_full_year:
            first_full_year = g_firstyear
    return first_full_year
