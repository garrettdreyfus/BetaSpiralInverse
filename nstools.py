import csv
from mpl_toolkits.basemap import Basemap
from numpy import linspace
import matplotlib.pyplot as plt
from numpy import meshgrid
import numpy as np
from netCDF4 import Dataset
import json
from geopy.distance import geodesic
from profile import Profile
import random

def extractProfiles(fnames):
    ##Load JSON data into profile objects
    ##and return index of deepest one
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    print(fnames)
    for fname in fnames:
        json_file = open(fname,'r') 
        data = json.load(json_file)
        for p in data.keys():
            profile = Profile(p,data[p])
            if len(profile.ipres)>0:
                profiles.append(profile)
                if data[p]["pres"][-1] > deepestdepth:
                    deepestindex = len(profiles)-1
                    deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex

def cruiseCount(profiles):
    cruises ={"null"}
    for profile in profiles:
        if profile.cruise not in cruises:
            cruises.add(profile.cruise)
    return cruises

def cruiseSearch(profiles,cruisename,year=None):
    results =[]
    for p in profiles:
        if p.cruise == cruisename:
            if year:
                if year == p.time.year:
                    results.append(p)
            else:
                results.append(p)
    return results

def plotCruise(profiles,cruisename,fig=None,ax=None,show=True):
    lats, lons, depths=[],[],[]
    for p in profiles:
        lats.append(p.lat)
        lons.append(p.lon)
        depths.append(np.max(p.pres))

    if not fig and not ax:
        fig,ax = plt.subplots(1,1)

    fig.suptitle(cruisename)
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    x,y = mapy(lons,lats)
    plt.scatter(x,y,c=depths,cmap="plasma")
    mapy.colorbar()
    if show:
        plt.show()

def plotProfiles(profiles,title,specialprofile=None,fig=None,ax=None,show=True):
    lats, lons, depths=[],[],[]
    for p in profiles:
        lats.append(p.lat)
        lons.append(p.lon)
        depths.append(np.max(p.pres))

    if not fig and not ax:
        fig,ax = plt.subplots(1,1)

    fig.suptitle(title)
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    x,y = mapy(lons,lats)
    plt.scatter(x,y,c=depths,cmap="plasma")
    mapy.colorbar()
    if specialprofile:
        x,y = mapy(specialprofile.lon,specialprofile.lat)
        plt.scatter(x,y,c="red")
    if show:
        plt.show()

def plotCruiseAndRef(cruises,refcruises,show=True):
    fig,ax = plt.subplots(1,1)
    lats, lons, depths=[],[],[]
    for p in refcruises:
        lats.append(p.lat)
        lons.append(p.lon)
        depths.append(np.max(p.pres))
    #fig.suptitle(cruisename)
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    x,y = mapy(lons,lats)
    plt.scatter(x,y)
    lats, lons, depths=[],[],[]
    for p in cruises:
        lats.append(p.lat)
        lons.append(p.lon)
        depths.append(np.max(p.pres))
    #fig.suptitle(cruisename)
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    x,y = mapy(lons,lats)
    plt.scatter(x,y,c=depths,cmap="plasma")
    mapy.colorbar()
    if show:
        plt.show()



def transArcticSearch(profiles):
    results = []
    for name in cruiseCount(profiles):
        westernprofile = False
        easternprofile = False
        for profile in cruiseSearch(profiles,name):
            if profile.lat > 70 and abs(profile.lon)<20:
                easternprofile =True
            if profile.lat > 70 and abs(profile.lon)>160:
                westernprofile = True
        if westernprofile and easternprofile:
            results.append(name)
    return results

        


def extractProfilesMonths(fnames,months):
    ##Load JSON data into profile objects
    ##and return index of deepest one
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    for fname in fnames:
        json_file = open(fname) 
        data = json.load(json_file)
        for p in data.keys():
            profile = Profile(p,data[p])
            if profile.time.month in months :
                if len(profile.ipres)>0:
                    profiles.append(profile)
                    if data[p]["pres"][-1] > deepestdepth:
                        deepestindex = len(profiles)-1
                        deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex

def extractProfilesBox(fnames,lonleft,lonright,latbot,lattop):
    ##Load JSON data into profile objects
    ##and return index of deepest one
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    for fname in fnames:
        json_file = open(fname) 
        data = json.load(json_file)
        for p in data.keys():
            profile = Profile(p,data[p])
            print(profile.lat,profile.lon)
            if latbot < profile.lat < lattop and lonleft < profile.lon < lonright:
                if len(profile.ipres)>0:
                    profiles.append(profile)
                    if data[p]["pres"][-1] > deepestdepth:
                        deepestindex = len(profiles)-1
                        deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex

def closestIdentifiedNS(profiles,queryprofile,depth,radius):
    minimumdistance = radius
    minprofile = None
    for p in profiles:
        if depth in p.neutraldepth.keys():
            dist = geodesic((queryprofile.lat,queryprofile.lon),(p.lat,p.lon)).km
            if dist<minimumdistance:
                minimumdistance = dist
                minprofile = p
    #print(minprofile.neutraldepth)
    return minprofile

def profileInBox(profiles,lonleft,lonright,latbot,lattop):
    results=[]
    for p in profiles:
        if lonleft< p.lon <lonright and latbot < p.lat < lattop:
            results.append(p)
    return results

def peerSearch(profiles,deepestindex,depth,profilechoice,radius=500):
    surfaces = {}
    surfaces[depth]=[[],[],[],[]]
    profilechoice.neutraldepth[depth] = depth
    print(profilechoice.eyed)
    references=[]
    references.append(profilechoice)
    print(len(profiles))
    while len(profiles)>0:
        foundcounter = 0
        for p in profiles.copy():
            closest = closestIdentifiedNS(references,p,depth,radius)
            if closest:
                #print(closest.neutraldepth)
                ns = closest.neutralDepth(p,closest.neutraldepth[depth],depthname=depth,searchrange=100) 
                if ns != None:
                    surfaces[depth][0].append(p.lon)
                    surfaces[depth][1].append(p.lat)
                    surfaces[depth][2].append(ns)
                    surfaces[depth][3].append(p.eyed)
                    profiles.remove(p)
                    foundcounter +=1
                    references.append(p)
        if foundcounter ==0:
            break

        print("found: ,",foundcounter,"left: ", len(profiles))
        #plotProfiles(references,"ITS SPREADING",specialprofile = profilechoice)
    return surfaces

def search(profiles,deepestindex):
    #Lets look for neutral surfaces every 200 dbar below 1000 dbar
    deeprange = range(1000,max(profiles[deepestindex].ipres),200)
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
                print(x)
                print(y)
                print(np.asarray(surfaces[i][2]))
                plt.scatter(x,y,c=np.asarray(surfaces[i][2]),cmap="plasma")
            mapy.colorbar()
            #map the reference profile
            if profiles and deepestindex:
                x,y = mapy(profiles[deepestindex].lon,profiles[deepestindex].lat)
                mapy.scatter(x,y,c="red")

            fig.suptitle("NS: "+str(i))
            #plt.show()
            mng = plt.get_current_fig_manager()
            mng.resize(*mng.window.maxsize())
            plt.savefig("refpics/RUN1/PRES/ns"+str(i)+".png")


def getProfileById(profiles,eyed):
    for p in profiles:
        if p.eyed == eyed:
            return p

def filterCruises(profiles,cruisenames):
    finalprofiles = []
    for p in profiles:
        if p.cruise in cruisenames:
            finalprofiles.append(p)
    return finalprofiles
        

#def mapPV(profiles,depth,contour=False):
    #ms = []
    #stds = []
    #for year in range(2008,2016):
        #lat = []
        #lon = []
        #pv = []
        #for p in profiles:
            #if p.time.year ==year: 
                #x = p.potentialVorticity(depth)
                #if x:
                    #lat.append(p.lat)
                    #lon.append(p.lon)
                    #pv.append(x)
        ##latfinal = [] 
        ##lonfinal = []
        ##pvfinal = []
        #m = np.mean(pv)
        #ms.append(m)
        #stdev = np.std(pv)
        #stds.append(stdev)
            ##for i in range(len(pv)):
            ##if abs(pv[i]-m) < stdev:
                ##latfinal.append(lat[i])
                ##lonfinal.append(lon[i])
                ##pvfinal.append(pv[i])
        ##lat = latfinal
        ##lon =lonfinal
        ##pv = pvfinal
                
        #mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
        #mapy.drawmapboundary(fill_color='aqua')
        #mapy.fillcontinents(color='coral',lake_color='aqua')
        #mapy.drawcoastlines()
        #x,y = mapy(lon,lat)
        #if contour:
            #plt.tricontourf(x,y,pv,cmap="RdYlGn")
        #else:
            #plt.scatter(x,y,c=pv,cmap="RdYlGn")
        #mapy.colorbar()
        #plt.title(year)
        #mng = plt.get_current_fig_manager()
        #mng.resize(*mng.window.maxsize())
        #plt.show()
    ##plt.scatter(range(1980,2015),ms,label="means")
    ##plt.scatter(range(1980,2015),stds,label="stdevs")
    ##plt.legend()
    ##plt.show()


def mapPV(profiles,depth,contour=False):
    lat = []
    lon = []
    pv = []
    for p in profiles:
        x = p.potentialVorticity(depth)
        if x:
            lat.append(p.lat)
            lon.append(p.lon)
            pv.append(x)
            print(x)
    latfinal = [] 
    lonfinal = []
    pvfinal = []
    m = np.mean(pv)
    #ms.append(m)
    stdev = np.std(pv)
    for i in range(len(pv)):
        if abs(pv[i]-m) < 2*stdev:
            latfinal.append(lat[i])
            lonfinal.append(lon[i])
            pvfinal.append(pv[i])
    lat = latfinal
    lon =lonfinal
    pv = pvfinal
            
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=0)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    x,y = mapy(lon,lat)
    if contour:
        plt.tricontourf(x,y,pv,cmap="RdYlGn")
    else:
        plt.scatter(x,y,c=pv,cmap="RdYlGn")
    mapy.colorbar()
    #plt.title(year)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

