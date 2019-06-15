import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import linspace
from numpy import meshgrid
import numpy as np
from netCDF4 import Dataset
import json
from profile import Profile



json_file = open('data/profiles.json') 
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
#graphInterpolation(profiles,2)
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

print(len(profiles))
#graphProfiles(profiles,5,9)
deeprange = range(min(profiles[deepestindex].ipres),max(profiles[deepestindex].ipres),200)
deeprange = range(1000,max(profiles[deepestindex].ipres),200)
surfaces = {}
for r in deeprange:
    surfaces[r]=[[],[],[]]
for j in range(len(profiles)):
    flag = False
    for r in deeprange:
        ns = profiles[deepestindex].neutralDepth(profiles[j],r) 
        if ns != None:
            surfaces[r][0].append(profiles[j].lon)
            surfaces[r][1].append(profiles[j].lat)
            surfaces[r][2].append(ns)
    if flag:
        graphProfiles(profiles,deepestindex,j)
for i in surfaces.keys():
    print(i,len(surfaces[i][1]))
    if len(surfaces[i][0])>3:
        fig,ax = plt.subplots(1,1)
        mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
        mapy.drawmapboundary(fill_color='aqua')
        mapy.fillcontinents(color='coral',lake_color='aqua')
        mapy.drawcoastlines()
        x,y = mapy(surfaces[i][0],surfaces[i][1])
        #plt.tricontourf(x,y,surfaces[i][2],cmap="plasma")
        mapy.scatter(x,y,c=(np.asarray(surfaces[i][2])-i))
        mapy.colorbar()
        x,y = mapy(profiles[deepestindex].lon,profiles[deepestindex].lat)
        mapy.scatter(x,y,c="red")
        fig.suptitle("NS: "+str(i))
        plt.show()

        
