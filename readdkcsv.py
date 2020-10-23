#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:33:54 2019

@author: ian
"""
import numpy as np
import argparse
from datetime import datetime, timedelta  
import matplotlib.pyplot as plt

minTemp = -5
minDepth = 225
midDepth = 250
maxDepth = 275

class oceanData() :
    def __init__(self,cruise,station,mType,lat,lon,date,bottom) :
        self.date=date
        self.cruise=cruise
        self.station=station
        self.lat=lat
        self.lon=lon
        self.mType=mType
        self.depth=np.array([])
        self.T=np.array([])
        self.bottom=bottom
        
    def addPoint(self,d,T) :
        self.depth=np.append(self.depth,d)
        self.T=np.append(self.T,T)
        
    def extractAverageDepth(self,minD,maxD) :
        # find valid points in range (minD,maxD)
        pts=np.logical_and(np.logical_and(self.depth >= minD,self.depth <= maxD),self.T > minTemp)
        #
        if sum(pts) > 0 :
            avgTemp= np.average(self.T[pts])
            avgDepth=np.average(self.depth[pts])
            print(self.date,avgTemp,avgDepth)
        else :
            avgDepth=np.nan
            avgTemp=np.nan
        return self.date,avgDepth,avgTemp
        
        
def getMyArgs() :
    parser=argparse.ArgumentParser(description='\033[1mRead CSV CTD Data\033[0m',epilog=' ')
    parser.add_argument('csvFile',type=str,nargs=1,help='csv file with data')
    #parser.add_argument('--plotnum',type=int,nargs=1,default=0,help='plot number for predfined plots')
    args=parser.parse_args()
    return args.csvFile[0]
    #

def parseFloatArg(myArg) :
    try :
        a=float(myArg)
    except :
        a=-999
    return a

#
# pulls date for individual l ine 
def parseTempLine(line) :
    if 'yyyy' not in line and len(line) > 20 :
        parts=line.split(',')
        cruise,station,myType=parts[0:3]
        myDate=datetime.strptime(parts[3],"%Y-%m-%dT%H:%M")
        myLat,myLon=parseFloatArg(parts[4]),parseFloatArg(parts[5])
        myBottom,myDepth=parseFloatArg(parts[6]),parseFloatArg(parts[7])
        myTemp=parseFloatArg(parts[8])
        
        return {'cruise' : cruise, 'station' : station, 'type' : myType, 'date' : myDate, 'lat' : myLat, 'lon' : myLon, 'bottom' : myBottom,
                'depth' : myDepth, 'temperature' : myTemp}
    else :
        return None

    
def consolodateRecords(myRecords) :
    #
    myData=[]
    for myRecord in myRecords :
        fRec=myRecord[0]
        # make a new ctd record
        newRec=oceanData(fRec['cruise'],fRec['station'],fRec['type'],fRec['lat'],fRec['lon'],fRec['date'],fRec['bottom'])
        # append to the list of records
        myData.append(newRec)
        # loop through the individual lines
        for ptMeas in myRecord :
            newRec.addPoint(ptMeas['depth'],ptMeas['temperature'])
            #
    # return a list of ocean data records, with one per ctd cast
    return myData

def parseTemps(csvFile) :
    fp =open(csvFile,'r')
    # this is just a fake last record to start things off 
    # records=[ [{ctd lines},{} ], [{}, {}]]
    lastRecord={'cruise' : '', 'station' : '', 'type' : '', 'date' : datetime(1900,1,1), 'lat' : 0., 'lon' : 0., 'bottom' : 0.,
                'depth' : 0., 'temperature' : 0.}
    records,currentRecord=[],[]
    #
    for line in fp :
        tempRecord=parseTempLine(line)
        if tempRecord == None :
            continue
        #
        # if date and location are the same, append this line
        if tempRecord['date'] == lastRecord['date'] and tempRecord['lat'] == lastRecord['lat'] and tempRecord['lon'] == lastRecord['lon'] :
            currentRecord.append(tempRecord)
        # otherwise start a new record
        else :
            # since first record, start a new list of lines for that ctd cast.
            currentRecord=[tempRecord]
            # add this new ctd record, in to the list of other ctd recrods
            records.append(currentRecord)
        lastRecord=tempRecord
    fp.close()
    return records

def main():
    csvFile=getMyArgs()
    #print(csvFile)
    myRecords=parseTemps(csvFile)
    myData = consolodateRecords(myRecords)
    print("First cast in record:")
    print("Depth: ", myData[0].depth)
    print("Temp:  ", myData[0].T)
    print("\nCast date,  Avg. temp [C],  Avg. depth [m]")
    #
    # setup storage
    myDate,myDepth,myTemp250=[],[],[]
    years=np.arange(1900,2020)
    # blank dictionaries to keep track of average (sum) and count
    Tcount={year : 0. for year in years}
    TYear={ year : 0. for year in years}
    #
    # this file has 250 m averages for each cast (225 to 275)
    fpOutAll=open(csvFile.replace('csv', str(midDepth)),'w')
    # loop over all ctd casts 
    for myDatum in myData :
        # my datum is an an ocean data record, so apply the average
        date,d250,T250=myDatum.extractAverageDepth(minDepth, maxDepth)
        # write the individual ctd averages to the .250 file
        if np.isfinite(T250) :
            print(f'{date.strftime("%Y-%m-%d")}  {d250:10.5f}  {T250:10.5f} {myDatum.lat:10.5f} {myDatum.lon:10.5f} ',file=fpOutAll)
        # make list of ponits
        myDate.append(date)
        myDepth.append(d250)
        myTemp250.append(T250)
        # update sums for average
        if T250 > minTemp and np.isfinite(T250) : 
            TYear[date.year] += T250
            Tcount[date.year]+=1.0
    fpOutAll.close()
    #
    # Now compute average for all years with data, use nan for no date
    #
    for year in years :
        if Tcount[year] > 0 :
            TYear[year]=TYear[year]/Tcount[year]
        else :
            TYear[year]=np.nan
    # basically make an array from the dictionaries
    TT=np.array([TYear[year] for year in years])
    #
    # output annual averages to .avg file
    fpOut=open(csvFile.replace('csv','avg'),'w')
    for year in years :
        # only output years with data
        if np.isfinite(TYear[year]) :
            print(year,TYear[year],Tcount[year],file=fpOut)
        
    plt.plot(years,TT,'r*')
    plt.xlabel('Year', fontsize=16)
    plt.ylabel('Average annual temperature ($^\circ$C)', fontsize=16)
    plt.show()  
    
main()
