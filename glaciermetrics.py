#!/usr/bin/env python3
# Calculate metrics for glacier change
# Taryn Black, August 2020

import pandas as pd
from shapely import ops
import gplots as gp


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
