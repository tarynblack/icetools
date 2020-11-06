#!/usr/bin/env python3
# Calculate metrics for glacier change
# Taryn Black, August 2020

import pandas as pd
from shapely import ops

def addDecade(start_date):
    start_date = pd.to_datetime(start_date).date()
    end_date = pd.date_range(start_date, periods=2, freq='9Y')[-1]
    end_date = end_date.date()
    return end_date


# def filterByDates(glacier, measure, date_start, date_end):
#     measures = glacier.extract(measure)

#     # Filter to data between selected dates
#     if date_start:
#         date_start = pd.to_datetime(date_start)
#         measures = measures.where(glacier.dates >= date_start).dropna()
#         dates = dates.where(glacier.dates >= date_start).dropna()
#     if date_end: 
#         date_end = pd.to_datetime(date_end)
#         measures = measures.where(glacier.dates <= date_end).dropna()
#         dates = dates.where(glacier.dates <= date_end).dropna()
#     return measures, dates


def firstFullYear(all_glaciers):
    """Identify the first year in which all glaciers have either an observed or an interpolated data point."""
    first_full_year = 0
    for g in all_glaciers:
        g_firstyear = all_glaciers[g].datayears[0]
        if g_firstyear > first_full_year:
            first_full_year = g_firstyear
    return first_full_year


def finalNetChange(glaciers, attr, startdate=None):
    final_net_change = pd.Series(index=glaciers.keys())
    for g in glaciers:
        glacier = glaciers[g]
        cumul_change, _, _ = glacier.cumulativeChange(attr, startdate)
        final_net_change.loc[g] = cumul_change.iloc[-1].values
    return final_net_change


def stdevChange(glaciers, attr, startdate=None):
    """Test whether glacier net change is outside one standard deviation of its variability."""
    change_stdev = pd.Series(index=glaciers.keys())
    for g in glaciers:
        glacier = glaciers[g]
        change_stdev.at[g] = getattr(glacier, attr).std(skipna=True)
    return change_stdev


def significantChange(glaciers, attr, startdate=None):
    """Identify glaciers that have experienced significant change, defined as a net change greater than one standard deviation of the total variability."""
    final_net_change = finalNetChange(glaciers, attr, startdate)
    change_stdev = stdevChange(glaciers, attr, startdate)
    significance = abs(final_net_change) > change_stdev
    return significance


def stableGlaciers(glaciers, attr, startdate=None):
    """Identify glaciers that have been stable in an attribute since startdate."""
    final_net_change = finalNetChange(glaciers, attr, startdate)
    significance = significantChange(glaciers, attr, startdate)
    stable_glaciers = final_net_change.where(significance==False).dropna()
    return stable_glaciers
    

def advancingGlaciers(glaciers, startdate=None):
    """Identify glaciers which have had significant net increase in length since startdate."""
    net_length_change = finalNetChange(glaciers, 'interplengths', startdate)
    significance = significantChange(glaciers, 'interplengths', startdate)
    advancing_glaciers = net_length_change.where(net_length_change > 0).dropna()
    advancing_glaciers = advancing_glaciers.where(significance==True).dropna()
    return advancing_glaciers


def growingGlaciers(glaciers, startdate=None):
    """Identify glaciers which have had significant net increase in area since startdate."""
    net_area_change = finalNetChange(glaciers, 'interpareas', startdate)
    significance = significantChange(glaciers, 'interpareas', startdate)
    growing_glaciers = net_area_change.where(net_area_change > 0).dropna()
    growing_glaciers = growing_glaciers.where(significance==True).dropna()
    return growing_glaciers


def retreatingGlaciers(glaciers, startdate=None):
    """Identify glaciers which have had significant net decrease in length since startdate."""
    net_length_change = finalNetChange(glaciers, 'interplengths', startdate)
    significance = significantChange(glaciers, 'interplengths', startdate)
    retreating_glaciers = net_length_change.where(net_length_change < 0).dropna()
    retreating_glaciers = retreating_glaciers.where(significance==True).dropna()
    return retreating_glaciers


def shrinkingGlaciers(glaciers, startdate=None):
    """Identify glaciers which have had significant net decrease in area since startdate."""
    net_area_change = finalNetChange(glaciers, 'interpareas', startdate)
    significance = significantChange(glaciers, 'interpareas', startdate)
    shrinking_glaciers = net_area_change.where(net_area_change < 0).dropna()
    shrinking_glaciers = shrinking_glaciers.where(significance==True).dropna()
    return shrinking_glaciers


def dominantGlaciers(glaciers, attr, startdate=None):
    """Identify glaciers that dominate in an attribute since startdate. Dominance defined as greater than two standard deviations of the population mean."""
    final_net_change = finalNetChange(glaciers, attr, startdate)
    mean_net_change = final_net_change.mean()
    std_net_change = final_net_change.std()
    dominant_threshold_pos = mean_net_change + 2*std_net_change
    dominant_threshold_neg = mean_net_change - 2*std_net_change
    dominant_glacier_pos = final_net_change.where(
        final_net_change > dominant_threshold_pos).dropna()
    dominant_glacier_neg = final_net_change.where(
        final_net_change < dominant_threshold_neg).dropna()
    dominant_glaciers = dominant_glacier_pos.append(dominant_glacier_neg)
    return dominant_glaciers

def filterGlaciers(glaciers, ids, idtype='remove'):
    """Filter glacier dataset by glacier IDs to keep or remove."""
    glacier_dict_copy = glaciers.copy()
    glacier_ids = glaciers.keys()
    if idtype == 'remove':
        remove_ids = ids
    elif idtype == 'keep':
        remove_ids = glacier_ids - ids
    [glacier_dict_copy.pop(id) for id in remove_ids]
    return glacier_dict_copy

