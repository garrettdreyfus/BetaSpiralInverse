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

def cruiseCount(fname):
    json_file = open(fname) 
    data = json.load(json_file)
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    cruises ={"null"}
    for p in data.keys():
        profile = Profile(p,data[p])
        if profile.cruise not in cruises:
            print(len(cruises))
        cruises.add(profile.cruise)
    print(len(cruises))
    print(cruises)


def extractProfilesMonths(fname,months):
    ##Load JSON data into profile objects
    ##and return index of deepest one
    json_file = open(fname) 
    data = json.load(json_file)
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    for p in data.keys():
        profile = Profile(p,data[p])
        if profile.time.month in months :
            if len(profile.ipres)>0:
                profiles.append(profile)
                if data[p]["pres"][-1] > deepestdepth:
                    deepestindex = len(profiles)-1
                    deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex



def search(profiles,deepestindex):
    #Lets look for neutral surfaces every 200 dbar below 1000 dbar
    deeprange = range(200,max(profiles[deepestindex].ipres),200)
    #A dictionary mapping neutral surface pressures to pressures at lat,lon
    surfaces = {}
    for r in deeprange:
        #note: an annoying format to store this information but easily graphable
        surfaces[r]=[[],[],[]]
    #Ever profile
    for j in range(len(profiles)):
        #every neutral surface
        for r in deeprange:
            ns = profiles[deepestindex].neutralDepth(profiles[j],r) 
            if ns != None:
                surfaces[r][0].append(profiles[j].lon)
                surfaces[r][1].append(profiles[j].lat)
                surfaces[r][2].append(ns)
    return surfaces

def graphSurfaces(surfaces,contour=False,profiles=None,deepestindex=None):
    for i in surfaces.keys():
        print(i,len(surfaces[i][1]))
        if len(surfaces[i][0])>3:
            fig,ax = plt.subplots(1,1)
            mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
            mapy.drawmapboundary(fill_color='aqua')
            mapy.fillcontinents(color='coral',lake_color='aqua')
            mapy.drawcoastlines()
            x,y = mapy(surfaces[i][0],surfaces[i][1])
            #Plot the surface 
            if contour:
                plt.tricontourf(x,y,np.asarray(surfaces[i][2]),cmap="plasma")
            else:
                plt.scatter(x,y,c=np.asarray(surfaces[i][2]),cmap="plasma")
            mapy.colorbar()
            #map the reference profile
            if profiles and deepestindex:
                x,y = mapy(profiles[deepestindex].lon,profiles[deepestindex].lat)
                mapy.scatter(x,y,c="red")
            fig.suptitle("NS: "+str(i))
            plt.show()

#
cruiseCount('data/2008profiles.json')

#profiles,deepestindex = extractProfilesMonths('data/profiles.json',range(13))
#print("DONE WITH EXTRACTING PROFILES")
#surfaces = search(profiles,deepestindex)
#print("DONE FINDING SURFACES")
#with open('data/surfaces.json', 'w') as outfile:
    #json.dump(surfaces, outfile)
##json_file = open("data/surfaces.json") 
##surfaces = json.load(json_file)
#print("NOW GRAPHING")
#graphSurfaces(surfaces)
