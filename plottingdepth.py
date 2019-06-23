import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import linspace
from numpy import meshgrid
import numpy as np
from netCDF4 import Dataset

#get the max depth at each profile
def extractProfileDepth(fname):
    f = open(fname, 'r')
    lons = []
    lats=[]
    maxDepthDict = {}
    for line in f:
        line = line.split()
        try:
            eyed = int(line[0])
            depth = float(line[9])
            lon = float(line[6])
            lat = float(line[7])
            if depth < 0:
                print(eyed)
                break
            if eyed in maxDepthDict.keys():
                if depth > maxDepthDict[eyed][0]:
                    maxDepthDict[eyed][0] = depth
            elif depth>500:
                maxDepthDict[eyed] = [depth,lon,lat]
            lons.append(lon)
            lats.append(lat)
        except:
            print(line)
    return maxDepthDict

#plot stations
def mapLocations(lats,lons):
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    x,y = mapy(lons,lats)
    mapy.scatter(x,y)

#plot a heatmap of maximum depths
def depthHeatMap(depthmap):
    fig, ax1 = plt.subplots(1,1)
    lons = []
    lats = []
    depths = []
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    for i in depthmap.values():
        depths.append(i[0])
        lons.append(i[1])
        lats.append(i[2])
    x,y = mapy(lons,lats)
    mapy.scatter(x,y,c=depths,cmap="plasma")
    mapy.colorbar()
    fig.suptitle("Profiles deeper than 500m 2008 UDASH")
    plt.show()

#return seafloor depth at a certain latitude and longitude
def searchBath(bathDataset,lat,lon):
    spacing = bathDataset.variables["spacing"][0]
    startlon = bathDataset.variables["x_range"][0]
    ret = bathDataset.variables["dimension"][:][0]
    z = bathDataset.variables["z"]
    i = int((90-lat)/spacing)
    j = int((lon-startlon)/spacing)
    if np.isnan(z[j+i*ret]):
        return 0
    else:
        return float(z[j+i*ret])

#map the difference between the deepest profile and bathymetry
def diffHeatMap(depthmap):
    fig, ax1 = plt.subplots(1,1)
    lons = []
    lats = []
    depths = []
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    d = Dataset("data/ver1_netcdf_geo.nc")
    for i in depthmap.values():
        realdepth = searchBath(d,i[2],i[1])
        depthdelta = abs(realdepth) - abs(i[0])
        if abs(depthdelta)<100:
            depths.append(depthdelta)
            lons.append(i[1])
            lats.append(i[2])
    x,y = mapy(lons,lats)
    mapy.scatter(x,y,c=depths,cmap="plasma")
    fig.suptitle("2008 Profiles Within 100 meters of the Seafloor")
    c = mapy.colorbar()
    c.set_label("Depth of Seafloor - Depth of Deepest CTD Reading in meters")
    fig, ax1 = plt.subplots(1,1)
    fig.suptitle("2008")
    ax1.set_xlabel("Depth of Seafloor - Depth of Deepest CTD Reading in meters")
    ax1.set_ylabel("Number of Profiles")
    ax1.hist(depths,bins=50)
    plt.show()

#plot bathymetry on map
def bathCheck():
    fig, ax1 = plt.subplots(1,1)
    lons = []
    lats = []
    depths = []
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    d = Dataset("data/ver1_netcdf_geo.nc")
    for lat in range(65,90):
        for lon in range(-180,180):
            lats.append(lat)
            lons.append(lon)
            depths.append(searchBath(d,lat,lon))
    print(depths)
    x,y = mapy(lons,lats)
    plt.tricontourf(x,y,depths,cmap="plasma")
    mapy.colorbar()
    fig.suptitle("Arctic Bathymetry")
    fig, ax1 = plt.subplots(1,1)
    plt.show()
 
            

maxDepthDict = extractProfileDepth("data/2008.txt")
diffHeatMap(maxDepthDict)
#bathCheck()

