import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import linspace
from numpy import meshgrid
import numpy as np
from netCDF4 import Dataset
import json
import glob


def extractProfiles(fnames):
    lons = []
    lats=[]
    profileDict = {}
    for fname in fnames:
        f = open(fname, 'r')
        for line in f:
            line = line.split()
            try:
                eyed = int(line[0])
                depth = float(line[9])
                lon = float(line[6])
                lat = float(line[7])
                time = str(line[5])
                pres = float(line[8])
                temp = float(line[11])
                sal = float(line[13])
                if eyed not in profileDict.keys():
                    profileDict[eyed] = {"time":time,"lat":lat,
                        "lon":lon, "pres":[],"sal":[],"temp":[]}
                profileDict[eyed]["sal"].append(sal)
                profileDict[eyed]["temp"].append(temp)
                profileDict[eyed]["pres"].append(pres)
            except:
                print("something went wrong: ",line)

    return profileDict


with open('data/profiles.json', 'w') as outfile:
    json.dump(extractProfiles(glob.glob("data/*.txt")), outfile)



