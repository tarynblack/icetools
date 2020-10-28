# Functions for common plots for glacier terminus/area analysis

import matplotlib
import matplotlib.pyplot as plt
import glaciermetrics as gm
import pandas as pd


season_c = {'spring': 'mediumseagreen',
            'summer': 'darkorchid',
            'autumn': 'darkorange',
            'winter': 'dodgerblue'}

measure_units = {'length': 'km', \
                 'area': '$km^2$'}

# Pick a color scheme for each termination type, and for each color use different marker shapes for each glacier
# lake = dodgerblue
# land = darkorchid
# tidewater = mediumseagreen
# mixed = darkorange
glacier_design = {
    1  : {'c' : 'dodgerblue',\
          's' : 'o-'},              # Bear, lake-terminating
    2  : {'c' : 'mediumseagreen',\
          's' : 'o-'},              # Aialik, tidewater
    3  : {'c' : 'dodgerblue',\
          's' : '^-'},              # Pedersen, lake-terminating
    4  : {'c' : 'mediumseagreen',\
          's' : '^-'},              # Holgate, tidewater
    5  : {'c' : 'darkorchid',\
          's' : 'o-'},              # South Holgate - West, land-terminating
    6  : {'c' : 'darkorange',\
          's' : 'o-'},              # South Holgate - East, mixed
    7  : {'c' : 'darkorchid',\
          's' : '^-'},              # Northeastern, land-terminating
    8  : {'c' : 'mediumseagreen',\
          's' : 's-'},              # Northwestern, tidewater
    9  : {'c' : 'mediumseagreen',\
          's' : 'D-'},              # Ogive, tidewater
    10 : {'c' : 'mediumseagreen',\
          's' : 'v-'},              # Anchor, tidewater
    11 : {'c' : 'darkorange',\
          's' : '^-'},              # Reconstitution, mixed
    12 : {'c' : 'darkorange',\
          's' : 's-'},              # Southwestern, mixed
    13 : {'c' : 'darkorchid',\
          's' : 's-'},              # Sunrise, land-terminating 
    14 : {'c' : 'darkorchid',\
          's' : 'D-'},              # Paguna, land-terminating
    15 : {'c' : 'mediumseagreen',\
          's' : 'p-'},              # McCarty, tidewater
    16 : {'c' : 'darkorchid',\
          's' : 'v-'},              # Dinglestadt, land-terminating
    17 : {'c' : 'darkorchid',\
          's' : 'p-'},              # Split, land-terminating
    18 : {'c' : 'dodgerblue',\
          's' : 's-'},              # Yalik, lake-terminating
    19 : {'c' : 'darkorange',\
          's' : 'D-'}               # Petrof, mixed
}


def figureProperties(fig, ax, graph):
    linewidth = 3
    markersize = 8
    titlesize = 24
    labelsize = 20
    ticklabelsize = 16
    if ax.get_legend() != None:
        legendsize = 16
    
    ax.axes.title.set_fontsize(titlesize)
    ax.xaxis.label.set_fontsize(labelsize)
    ax.yaxis.label.set_fontsize(labelsize)
    ax.tick_params(labelsize=ticklabelsize)
    ax.grid(color='lightgray')
    if ax.get_legend() != None:
        ax.legend(fontsize=legendsize)
    ax.axhline(linewidth=1.0, color='black')
    ax.set_axisbelow(True)

    if type(graph) == matplotlib.lines.Line2D:
        graph.set_zorder(5)
        graph.set_linewidth(linewidth)
        graph.set_markersize(markersize)


def manageSubplots(fig, subplot_bool, idx):
    if subplot_bool == False:
        ax = fig.add_axes([0., 0., 1., 1.])
    elif subplot_bool == True:
        ax = fig.add_subplot(idx)
    return ax


def getGlacierName(GlacierClass):
    if GlacierClass.officialname:
        name = GlacierClass.officialname
    else:
        name = GlacierClass.unofficialname
    return name


def getSeasonMeasures(dates, measures, seasons, seasonstr):
    dates_season = dates.where(seasons==seasonstr).dropna()
    measures_season = measures.where(seasons==seasonstr).dropna()
    return dates_season, measures_season


def check_measure(measure):
    measure_types = ['length', 'area']
    if measure not in measure_types:
        raise ValueError("Invalid measure type. Expected one of: {}".format(
            measure_types))


def align_yscale(ax1, ax2):
    ax1min, ax1max = ax1.get_ylim()
    ax2min, ax2max = ax2.get_ylim()
    axmin = min(ax1min, ax2min)
    axmax = max(ax1max, ax2max)
    ax1.set_ylim(axmin, axmax)
    ax2.set_ylim(axmin, axmax)


def annotate_bars(ax, anno, x, y):
    for i in range(len(anno)):
        sign = y[i]/abs(y[i])
        ax.annotate(anno[i], 
                    xy=(x[i], 0),#y[i]), 
                    xytext = (0, -sign * 8),
                    textcoords = 'offset points',
                    ha='center', va='center')


def totalRelativeMeasure(fig, glacier, measure, subplots=False, idx=111):
    check_measure(measure)

    ax = manageSubplots(fig, subplots, idx)
    dates = glacier.extract('date') 
    cumul_measures, _, _ = gm.netMeasureChange(glacier, measure)
    name = getGlacierName(glacier)

    graph, = ax.plot(dates, cumul_measures, 'o-', color='mediumblue')
    ax.set_title('{}: Observed {} Change'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('Cumulative {} Change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    figureProperties(fig, ax, graph)


def totalRelativeMeasureCompare(fig, glacier, subplots=False, idx=111):
    dates = glacier.extract('date') 
    name = getGlacierName(glacier)
    cumul_area, _, _ = gm.netMeasureChange(glacier, 'area')
    cumul_termarea, _, _ = gm.netMeasureChange(glacier, 'termarea')
    cumul_length, _, _ = gm.netMeasureChange(glacier, 'length')
    
    # Plot areas on one axis
    ax_area = manageSubplots(fig, subplots, idx)
    gr_area, = ax_area.plot(dates, cumul_area, 'o-', color='darkblue', \
        label='Terminus and lateral area')
    gr_termarea, = ax_area.plot(dates, cumul_termarea, 'o-', color='blue', \
        fillstyle='none', label='Terminus area only')
    figureProperties(fig, ax_area, gr_termarea)
    gr_length, = ax_area.plot([], [], 's-', color='salmon', \
        label='Centerline length')
    ax_area.set_title('{}: Observed Size Changes'.format(name))
    ax_area.set_xlabel('Date')
    ax_area.set_ylabel('Cumulative Area Change ({})'.format(
        measure_units['area']))
    ax_area.legend()
    
    # Plot length on other axis
    ax_length = ax_area.twinx()
    gr_length, = ax_length.plot(dates, cumul_length, 's-', color='salmon', \
        alpha=0.7, label='Centerline length')
    ax_length.set_ylabel('Cumulative Length Change ({})'.format(
        measure_units['length']))
    
    figureProperties(fig, ax_area, gr_area)
    figureProperties(fig, ax_area, gr_termarea)
    figureProperties(fig, ax_length, gr_length)
    ax_length.grid(b=None)
    align_yscale(ax_area, ax_length)


def seasonRelativeMeasure(fig, glacier, measure, \
    spring=False, summer=False, autumn=False, winter=False, \
    subplots=False, idx=111):
    check_measure(measure)

    ax = manageSubplots(fig, subplots, idx)
    dates = glacier.extract('date')
    measures = glacier.extract(measure)
    seasons = glacier.extract('season')
    name = getGlacierName(glacier)

    if spring:
        dates_spr, measures_spr = getSeasonMeasures(dates, measures, \
            seasons, 'SPR')
        cumul_measures_spr = measures_spr.diff().cumsum()
        cumul_measures_spr.iloc[0] = 0.0
        graph, = ax.plot(dates_spr, cumul_measures_spr, \
            'o-', color=season_c['spring'], label='Spring')
        figureProperties(fig, ax, graph)
    if summer:
        dates_sum, measures_sum = getSeasonMeasures(dates, measures, \
            seasons, 'SUM')
        cumul_measures_sum = measures_sum.diff().cumsum()
        cumul_measures_sum.iloc[0] = 0.0
        graph, = ax.plot(dates_sum, cumul_measures_sum, \
            '^-', color=season_c['summer'], label='Summer')
        figureProperties(fig, ax, graph)
    if autumn:
        dates_aut, measures_aut = getSeasonMeasures(dates, measures, \
            seasons, 'AUT')
        cumul_measures_aut = measures_aut.diff().cumsum()
        cumul_measures_aut.iloc[0] = 0.0
        graph, = ax.plot(dates_aut, cumul_measures_aut, \
            's-', color=season_c['autumn'], label='Autumn')
        figureProperties(fig, ax, graph)
    if winter:
        dates_win, measures_win = getSeasonMeasures(dates, measures, \
            seasons, 'WIN')
        cumul_measures_win = measures_win.diff().cumsum()
        cumul_measures_win.iloc[0] = 0.0
        graph, = ax.plot(dates_win, cumul_measures_win, \
            'D-', color=season_c['winter'], label='Winter')
        figureProperties(fig, ax, graph)
    
    ax.set_title('{}: Observed Seasonal {} Change'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('Cumulative {} Change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    ax.legend()
    figureProperties(fig, ax, graph)


def individualMeasureChange(fig, glacier, measure, \
    date_start=None, date_end=None, subplots=False, idx=111):
    check_measure(measure)
    
    ax = manageSubplots(fig, subplots, idx)
    measures, dates = gm.filterByDates(glacier, measure, date_start, date_end)
    measure_change = measures.diff()
    name = getGlacierName(glacier)

    graph = ax.bar(dates, measure_change, width=75, color='mediumblue')
    ax.set_title('{}: {} Change Between Observations'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('{} change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    figureProperties(fig, ax, graph)


def decadalMeasureChange(fig, glacier, measure, decade_startyears, \
    subplots=False, idx=111):
    check_measure(measure)
    
    ax = manageSubplots(fig, subplots, idx)
    name = getGlacierName(glacier)

    decade_startyears = pd.to_datetime(decade_startyears)
    decadal_net_changes = pd.Series(dtype='float64')
    decade_labels = []
    bar_annotations = []
    for startyear in decade_startyears:
        endyear = gm.addDecade(startyear)
        cumul_decadal_measure_change, _, num_obs = gm.netMeasureChange(
            glacier, measure, startyear,endyear)
        net_decadal_measure_change = cumul_decadal_measure_change.iloc[-1]
        midyear = endyear.year - 4
        decadal_net_changes.loc[midyear] = net_decadal_measure_change
        decade_labels.append('{}-{}'.format(startyear.year, startyear.year+9))
        bar_annotations.append('{} obsv'.format(num_obs))
    
    rects = ax.bar(decadal_net_changes.index.values, 
                   decadal_net_changes.values, 
                   color='mediumblue', width=5)
    annotate_bars(ax, bar_annotations, 
                  decadal_net_changes.index.values, decadal_net_changes.values)
    ax.set_title('{}: Observed Decadal {} Change'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Decade')
    ax.set_ylabel('Net {} Change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    xtlocs = ax.get_xticks()
    ax.set_xticks(xtlocs+5)
    ax.set_xticklabels(decade_labels)
    if max(ax.get_ylim()) == 0.0:
        yrange, _ = ax.get_ylim()
        ax.set_ylim(top=abs(0.05*yrange))
    figureProperties(fig, ax, rects)
    ax.set_axisbelow(True)
    ax.grid(axis='x')


def seasonMeasureChange(fig, glacier, measure, \
    spring=False, summer=False, autumn=False, winter=False, \
    subplots=False, idx=111):
    check_measure(measure)
    
    ax = manageSubplots(fig, subplots, idx)
    dates = glacier.extract('date')
    measures = glacier.extract(measure)
    seasons = glacier.extract('season')
    name = getGlacierName(glacier)

    if spring:
        dates_spr, measures_spr = getSeasonMeasures(
            dates, measures, seasons, 'SPR')
        dates_spr_change = dates_spr[dates_spr.notna()]
        measures_spr_change = measures_spr[measures_spr.notna()].diff()
        graph = ax.bar(dates_spr_change, measures_spr_change, \
            width=75, color=season_c['spring'], label='Spring')
        figureProperties(fig, ax, graph)
    if summer:
        dates_sum, measures_sum = getSeasonMeasures(
            dates, measures, seasons, 'SUM')
        dates_sum_change = dates_sum[dates_sum.notna()]
        measures_sum_change = measures_sum[measures_sum.notna()].diff()
        graph = ax.bar(dates_sum_change, measures_sum_change, \
            width=75, color=season_c['summer'], label='Summer')
        figureProperties(fig, ax, graph)
    if autumn:
        dates_aut, measures_aut = getSeasonMeasures(
            dates, measures, seasons, 'AUT')
        dates_aut_change = dates_aut[dates_aut.notna()]
        measures_aut_change = measures_aut[measures_aut.notna()].diff()
        graph = ax.bar(dates_aut_change, measures_aut_change, \
            width=75, color=season_c['autumn'], label='Autumn')
        figureProperties(fig, ax, graph)
    if winter:
        dates_win, measures_win = getSeasonMeasures(
            dates, measures, seasons, 'WIN')
        dates_win_change = dates_win[dates_win.notna()]
        measures_win_change = measures_win[measures_win.notna()].diff()
        graph = ax.bar(dates_win_change, measures_win_change, \
            width=75, color=season_c['winter'], label='Winter')
        figureProperties(fig, ax, graph)
    
    ax.set_title('{}: Observed Seasonal {} Change'.format(
        name, measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('{} change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    ax.legend()
    figureProperties(fig, ax, graph)


def annualObservations(fig, all_glaciers, years_list, show_firstyear=True, \
    subplots=False, idx=111):
    ax = manageSubplots(fig, subplots, idx)
    # xticklabs = []

    for g in all_glaciers:
        glacier = all_glaciers[g]
        # name = getGlacierName(glacier)
        # xticklabs.append(name)
        gid = glacier.gid
        obs_years = glacier.extract('hydroyear')
        interp_years = glacier.interpyears

        graph1 = ax.scatter([gid]*len(obs_years), obs_years, \
            color='mediumblue')
        graph2 = ax.scatter([gid]*len(interp_years), interp_years, \
            edgecolors='gray', facecolors='none')

        figureProperties(fig, ax, graph1)
        figureProperties(fig, ax, graph2)
    
    if show_firstyear:
        first_full_year = gm.firstFullYear(all_glaciers)
        ax.axhline(y=first_full_year, linewidth=3.0, color='red', zorder=0.5)
    ax.set_title('Observation Time Series for Each Glacier')
    ax.set_ylabel('Hydrological Year')
    ax.set_xlabel('Glacier ID')
    # ax.legend()
    # xticklocs = range(1, len(all_glaciers)+1)
    # plt.xticks(xticklocs, xticklabs)
    # plt.xlim(left=0)
    plt.ylim(bottom=years_list[0]-1, top=years_list[-1]+1)

    figureProperties(fig, ax, graph1)
    figureProperties(fig, ax, graph2)



def measureSummary(fig, all_glaciers, measure, subplots=False, idx=111):
    check_measure(measure)

    ax = manageSubplots(fig, subplots, idx)

    for g in all_glaciers:
        glacier = all_glaciers[g]
        dates = glacier.extract('date')
        cumul_measures, _, _ = gm.netMeasureChange(glacier, measure)
        name = getGlacierName(glacier)
        graph, = ax.plot(dates, cumul_measures,
            glacier_design[g]['s'], color=glacier_design[g]['c'],
            label=name)
        figureProperties(fig, ax, graph)
    
    ax.set_title('Summary of Observed {} Changes'.format(measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('Cumulative {} Change ({})'.format(
        measure.capitalize(), measure_units[measure]))
    ax.legend(loc='center left', bbox_to_anchor=(1.04, 0.5))
    figureProperties(fig, ax, graph)


def normMeasureSummary(fig, all_glaciers, measure, subplots=False, idx=111):
    check_measure(measure)

    ax = manageSubplots(fig, subplots, idx)

    for g in all_glaciers:
        glacier = all_glaciers[g]
        dates = glacier.extract('date')
        scaled_measure = gm.normMeasureChange(glacier, measure)
        name = getGlacierName(glacier)
        graph, = ax.plot(dates, scaled_measure,
            glacier_design[g]['s'], color=glacier_design[g]['c'],
            label=name)
        figureProperties(fig, ax, graph)
    
    ax.set_title('Summary of Observed {} Changes, Normalized'.format(
        measure.capitalize()))
    ax.set_xlabel('Date')
    ax.set_ylabel('Scaled Cumulative {} Change'.format(measure.capitalize()))
    plt.ylim(-0.02, 1.02)
    # ax.legend(loc='center left', bbox_to_anchor=(1.04, 0.5))
    figureProperties(fig, ax, graph)
