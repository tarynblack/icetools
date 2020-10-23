# Functions for common plots for glacier terminus/area analysis

import matplotlib
import matplotlib.pyplot as plt


def figureProperties(fig, ax, graph):
    linewidth = 5
    markersize = 15
    titlesize = 24
    labelsize = 20
    ticklabelsize = 16
    # legendsize = 16
    
    ax.axes.title.set_fontsize(titlesize)
    ax.xaxis.label.set_fontsize(labelsize)
    ax.yaxis.label.set_fontsize(labelsize)
    ax.tick_params(labelsize=ticklabelsize)
    ax.grid(color='lightgray')

    if type(graph) == matplotlib.lines.Line2D:
        graph.set_linewidth(linewidth)
        graph.set_markersize(markersize)


def manageSubplots(fig, subplot_bool, idx):
    if subplot_bool == False:
        ax = fig.add_axes([0., 0., 1., 1.])
    elif subplot_bool == True:
        ax = fig.add_subplot(idx)
    return ax


def annualArea(fig, time, annual_area, subplots=False, idx=111):
    ax = manageSubplots(fig, subplots, idx)
    graph, = ax.plot(time, annual_area, '.-')
    ax.set_title("Estimated total area")
    ax.set_xlabel("Hydrological year")
    ax.set_ylabel(r"Cumulative area change (km$^2$)")
    figureProperties(fig, ax, graph)


def annualAreaChange(fig, time, area_change, subplots=False, idx=111):
    ax = manageSubplots(fig, subplots, idx)
    graph = ax.bar(time, area_change, color='cornflowerblue')
    ax.set_title(r"Total Annual $\Delta$Area")
    ax.set_xlabel("Hydrological year")
    ax.set_ylabel(r"$\Delta$Area (km$^2$)")
    figureProperties(fig, ax, graph)


def numberObserved(fig, time, counts, subplots=False, idx=111):
    ax = manageSubplots(fig, subplots, idx)
    graph = ax.bar(time, counts, color='cornflowerblue')
    ax.set_title("Terminus observations, %s-%s" % (time[0], time[-1]))
    ax.set_xlabel("Hydrological year")
    ax.set_ylabel("Number of glaciers observed")
    figureProperties(fig, ax, graph)
    