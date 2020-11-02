# Import and plot ocean temperature data
# T Black, August 2019

from os import listdir
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

datadir = "/Volumes/insar5/teblack/data/oceandata/"
resolution = 'b'  # 'b' = low resolution; 'c' = high resolution (CTDs)
years = np.arange(1970, 2021)

header250 = ["Date", "avgDepth", "avgTemp", "lat", "lon"]
headeravg = ["Year", "avgTemp", "numCasts"]


def quadData(path, header):
    data = pd.read_table(path, sep='\s+', header=None, names=header)
    return data


def quadBoundaries(filedir):
    bounds = filedir.split("/")[0].split("_")[1:5]
    edges = {"N": bounds[1],
             "S": bounds[0],
             "E": bounds[2],
             "W": bounds[3]}
    return edges


def assembleQuadDict(data, edges):
    quadrangle = {"edges": edges,
                  "data": data}
    return quadrangle


def getAllFiles(datadir, resolution):
    """load all files...?
    datadir : main directory containing subdirectories of quadrangle data
    resolution : which file type to load; 'b' low res, 'c' high res (CTD)"""
    quaddirs = listdir(path=datadir)
    allQuads = {}
    for d in quaddirs:
        prefix = d.split("_")[0]
        file250 = "%s%s.250" % (prefix, resolution)
        fileavg = "%s%s.avg" % (prefix, resolution)
        quad250 = assembleQuadDict(quadData(datadir+d+'/'+file250, header250),
                                   quadBoundaries(d+'/'+file250))
        quadavg = assembleQuadDict(quadData(datadir+d+'/'+fileavg, headeravg),
                                   quadBoundaries(d+'/'+fileavg))
        allQuads[d] = {'250': quad250, 'avg': quadavg}
    return quaddirs, allQuads


def getQuadDataFromDict(quadsdict, quad, datatype):
    """Extract data from one quadrangle (quad) from the dictionary of all quads
    (quadsdict). Datatype is either '250' (for all quad measurements around
    250m deep), or 'avg' (for annual average measurements in a quad)."""
    data = quadsdict[quad][datatype]['data']
    return data


def getQuadEdgesFromDict(quadsdict, quad, datatype):
    """Extract boundaries of one quadrangle (quad) from the dictionary of all
    quads (quadsdict). Datatype is either '250' (for all quad measurements
    around 250m deep), or 'avg' (for annual average measurements in a quad)."""
    edges = quadsdict[quad][datatype]['edges']
    return edges


def checkQuadWithinBounds(quadedges, NSEW, bound, direction):
    """Check whether a quadrangle falls within prescribed geographic bounds.
    Returns True/False.
    quadedges : dictionary of edges of quadrangle (from getQuadEdgesFromDict)
    NSEW      : which edge to evaluate ('N', 'S', 'E', 'W')
    bound     : geographic boundary value to check against (degrees)
    direction : 'above' or 'below', whether quad should be above/below bound"""
    if direction == 'above':
        return int(quadedges[NSEW]) >= bound
    elif direction == 'below':
        return int(quadedges[NSEW]) <= bound
    else:
        print('Invalid direction, must be \'above\' or \'below\'')


def subsetQuadsGeographic(quadsdict, datatype, bounds):
    """Go through all quadrangles and evaluate whether they fall within
    prescribed geographic boundaries. Return a subset of accepted quads.
    quadsdict : dictionary of all quadrangles and their data
    datatype  : '250' or 'avg'
    bounds    : nested list of boundary conditions, matching
                checkQuadWithinBounds arguments [NSEW, bound, direction], e.g.
                [['N', 70, 'above'],['E', 56, 'below']]. Must be nested even if
                only one bound."""
    subsetQuads = {}
    for d in quadsdict:
        edges = getQuadEdgesFromDict(quadsdict, d, datatype)
        allBounds = []
        for b in bounds:
            allBounds.append(checkQuadWithinBounds(edges, b[0], b[1], b[2]))
        if not all(allBounds):
            continue
        else:
            subsetQuads[d] = allQuads[d]
    return subsetQuads


def summarizeQuadsAnnual(quadsdict, years):
    """Summarize annual average temperature data in a group of quadrangles."""
    quads = list(quadsdict.keys())
    annualData = []
    for y in years:
        yeardata = []
        for q in quads:
            data = getQuadDataFromDict(quadsdict, q, 'avg')
            if y in list(data.Year):
                i = list(data.Year).index(y)
                tempXcast = data.avgTemp[i]*data.numCasts[i]
                qcast = data.numCasts[i]
                yeardata.append([tempXcast, qcast])
        if not yeardata:
            annualData.append([y, np.nan, np.nan])
        elif yeardata:
            yeardata = np.array(yeardata)
            sumCasts = yeardata[:,1].sum()
            annAvgTemp = yeardata[:,0].sum() / sumCasts
            annualData.append([y, annAvgTemp, sumCasts])
    annualData = np.array(annualData)
    return annualData


def annualStationsQuad(quad):
    """Plot number of stations per year (one station can have multiple
    measurements on a single cast). Same as plotted at ocean.ices.dk but as a
    more useful bar plot rather than line plot."""
    years = quad["data"].Year
    casts = quad["data"].numCasts
    plt.figure()
    plt.bar(years, casts)
    plt.xlabel('Year')
    plt.ylabel('# casts')
    plt.title("Annual stations at ~250m in quadrangle %s-%s N, %s-%s W" %
              (quad["edges"]["S"], quad["edges"]["N"],
               quad["edges"]["E"], quad["edges"]["W"]))


def annualStationsSummary(summary, title):
    """Plot number of stations per year for a summary of quadrangles."""
    years = summary[:,0]
    casts = summary[:,2]
    plt.figure()
    plt.bar(years, casts)
    plt.xlabel('Year')
    plt.ylabel('# casts')
    plt.title(title)


def plotTemperatureRecords(quad):
    """Plot individual temperature measurements at lat/lon, depth."""
    lat = quad["data"].lat
    lon = quad["data"].lon
    depth = quad["data"].avgDepth
    temp = quad["data"].avgTemp
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(lon, lat, -depth, c=temp)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Depth (m)')


def plotAnnualAvgRegionalTemp(summary, quadsdict, title):
    years = summary[:,0]
    temps = summary[:,1]
#    plt.figure()
    h, = plt.plot(years, temps, '.-', linewidth=2, markersize=10, c='cornflowerblue')
    for q in quadsdict:
        y = getQuadDataFromDict(quadsdict, q, 'avg').Year
        t = getQuadDataFromDict(quadsdict, q, 'avg').avgTemp
#        y = quadsdict[q]['avg']['data'].Year
#        t = quadsdict[q]['avg']['data'].avgTemp
        plt.scatter(y, t, s=10, c='gray')
#    plt.xlabel('Year', fontsize=16)
#    plt.ylabel('Temperature ($^\circ$C)', fontsize=16)
    plt.title(title, fontsize=20)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
#    plt.grid()
    plt.ylim(-1.2, 5)
    return h


#quad250 = assembleQuadDict(quadData(datadir+file250, header250),
#                             quadBoundaries(file250))
#quadavg = assembleQuadDict(quadData(datadir+fileavg, headeravg),
#                             quadBoundaries(fileavg))

#annualStations(quadavg)
#plotMeasurements(quad250)

quadDirs, allQuads = getAllFiles(datadir, resolution)

#plotTemperatureRecords(allQuads[quadDirs[20]]['250'])

# Subset quadrangles into those north and south of 72N
northof72 = subsetQuadsGeographic(allQuads, '250', [['S', 72, 'above']])
southof72 = subsetQuadsGeographic(allQuads, '250', [['N', 72, 'below']])

# Summarize quadrangle data in geographic subsets
northof72_summary = summarizeQuadsAnnual(northof72, years)
southof72_summary = summarizeQuadsAnnual(southof72, years)

# Plot number of stations per year in each region
annualStationsSummary(northof72_summary, 'Annual stations north of 72$^\circ$N')
annualStationsSummary(southof72_summary, 'Annual stations south of 72$^\circ$N')

# Plot average temperature in each region
#plt.figure()
#plt.subplot(1,2, 1)
#plotAnnualAvgRegionalTemp(northof72_summary, northof72, 'Avg annual 250m ocean temperature N of 72$^\circ$N')
#plt.subplot(1, 2, 2)
#plotAnnualAvgRegionalTemp(southof72_summary, southof72, 'Avg annual 250m ocean temperature S of 72$^\circ$N')

#COLORMAP = cm.Blues

# %%
LINEWIDTH = 5
MARKERSIZE = 15
TITLESIZE = 24
LABELSIZE = 20
TICKLABELSIZE = 16
LEGENDSIZE = 16

fig = plt.figure(figsize=(10, 10))
plt.suptitle("Annual Average Ocean Temperature, 250m Depth", size=TITLESIZE)
#color = iter(COLORMAP(np.linspace(0.2,1,4)))
ax1 = fig.add_subplot(211)
plt.grid(color='lightgray', which='both')
h1 = plotAnnualAvgRegionalTemp(northof72_summary, northof72, 'North of 72$^\circ$N')
#h1.set_color(next(color))
h1.set_linewidth(LINEWIDTH)
h1.set_markersize(MARKERSIZE)
plt.tick_params(axis='both', labelsize=TICKLABELSIZE)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.ylabel('Temperature ($^\circ$C)', size=LABELSIZE)

ax2 = fig.add_subplot(212, sharex=ax1)
plt.grid(color='lightgray', which='both')
h2 = plotAnnualAvgRegionalTemp(southof72_summary, southof72, 'South of 72$^\circ$N')
#h2.set_color(next(color))
h2.set_linewidth(LINEWIDTH)
h2.set_markersize(MARKERSIZE)
plt.tick_params(axis='both', labelsize=TICKLABELSIZE)
plt.xlabel("Year", size=LABELSIZE)
plt.ylabel('Temperature ($^\circ$C)', size=LABELSIZE)

plt.savefig("/home/teblack/plots/AGU_ocntemp.png", transparent=True, bbox_inches='tight')
