import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import linspace
from numpy import meshgrid
import numpy as np
from netCDF4 import Dataset
import json
from profile import Profile



def extractProfiles(fname):
    ##Load JSON data into profile objects
    ##and return index of deepest one
    json_file = open(fname) 
    data = json.load(json_file)
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    for p in data.keys():
        profile = Profile(p,data[p])
        if len(profile.ipres)>0:
            profiles.append(profile)
            if data[p]["pres"][-1] > deepestdepth:
                deepestindex = len(profiles)-1
                deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex

 
##how is it interpolating?
def graphInterpolation(ps,i):
    fig,(ax1,ax2)= plt.subplots(1,2)
    print(ps[i].sals)
    print(ps[i].pres)
    print(ps[i].eyed)
    ax1.plot(ps[i].temps,ps[i].pres)
    ax1.plot(ps[i].itemps,ps[i].ipres)
    ax1.plot(ps[i].sals,ps[i].pres)
    ax1.plot(ps[i].isals,ps[i].ipres)
    ax2.plot(ps[i].idensities-1000,ps[i].ipres)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    plt.show()
##graph two profiles
def graphProfiles(ps,i,j):
    fig,(ax1,ax2,ax3)= plt.subplots(1,3)
    ax1.plot(ps[j].itemps,ps[j].ipres)
    ax1.plot(ps[i].itemps,ps[i].ipres)
    ax2.plot(ps[j].isals,ps[j].ipres)
    ax2.plot(ps[i].isals,ps[i].ipres)
    ax3.plot(ps[i].idensities-1000,ps[i].ipres)
    ax3.plot(ps[j].idensities-1000,ps[j].ipres)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()
    plt.show()

profiles,deepestindex = extractProfiles('data/profiles.json')
graphProfiles(profiles,4,5)
