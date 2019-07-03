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
from scipy.interpolate import Rbf
import pyproj

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



def deepestProfile(profiles):
    deepestindex=-1
    deepestdepth=-1
    for p in range(len(profiles)):
        if profiles[p].ipres[-1] >deepestdepth:
            deepestdepth=profiles[p].ipres[-1]
            deepestindex=p
    return deepestindex


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
            if latbot <= profile.lat <= lattop and lonleft <= profile.lon <= lonright:
                if len(profile.ipres)>0:
                    profiles.append(profile)
                    if data[p]["pres"][-1] > deepestdepth:
                        deepestindex = len(profiles)-1
                        deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex

def removeNorwegianSea(profiles):
    finalprofiles = []
    deepestindex = -1
    deepestdepth = 0
    for pindex in range(len(profiles)):
        p = profiles[pindex]
        if not (p.lat<81 and -23 < p.lon < 20):
            finalprofiles.append(p)
            if p.ipres[-1] > deepestdepth:
                deepestindex = pindex
                deepestdepth = p.ipres[-1]
    return finalprofiles,deepestindex 

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

def addDataToSurfaces(profiles,surfaces,stdevs):
    tempSurfs = {}
    for k in surfaces.keys():
        tempSurf = np.array([[],[],[[],[],[],[]],[]])
        for l in range(len(surfaces[k][0])):
            p = getProfileById(profiles,surfaces[k][3][l])
            t,s = p.atPres(surfaces[k][2][l])
            pv = p.potentialVorticity(surfaces[k][2][l])
            if t and s and pv and p:
                tempSurf[0].append(surfaces[k][0][l])
                tempSurf[1].append(surfaces[k][1][l])
                tempSurf[2][0].append(-surfaces[k][2][l])
                tempSurf[2][1].append(t)
                tempSurf[2][2].append(s)
                tempSurf[2][3].append(pv)
                tempSurf[3].append(surfaces[k][3][l])
        for j in range(len(tempSurf[2])):
            m = np.median(tempSurf[2][j])
            s = np.std(tempSurf[2][j])
            a = np.asarray(np.where(abs(tempSurf[2][j] -m)>stdevs*s))
            np.asarray(tempSurf[2][j])[a[0]] == np.nan

        if len(tempSurf[0])>5:
            tempSurfs[k] = tempSurf

    return tempSurfs



def graphSurfaces(surfaces,quantindex,contour=False,profiles=None,deepestindex=None):
    for i in surfaces.keys():
        if len(surfaces[i][0])>3:
            fig,ax = plt.subplots(1,1)
            mapy = Basemap(projection='ortho', lat_0=90,lon_0=-60)
            mapy.drawmapboundary(fill_color='aqua')
            mapy.fillcontinents(color='coral',lake_color='aqua')
            mapy.drawcoastlines()
            x,y = mapy(surfaces[i][0],surfaces[i][1])
            #Plot the surface 
            if contour:
                plt.tricontourf(x,y,np.asarray(surfaces[i][2][quantindex]),cmap="plasma")
            else:
                plt.scatter(x,y,c=np.asarray(surfaces[i][2][quantindex]),cmap="plasma")
                m = np.median(np.asarray(surfaces[i][2][quantindex]))
                s = np.std(np.asarray(surfaces[i][2][quantindex]))
                plt.clim(m-2*s,m+2*s)
            mapy.colorbar()
            #map the reference profile
            if profiles and deepestindex:
                x,y = mapy(profiles[deepestindex].lon,profiles[deepestindex].lat)
                mapy.scatter(x,y,c="red")

            fig.suptitle("NS: "+str(i))
            mng = plt.get_current_fig_manager()
            mng.resize(*mng.window.maxsize())
            plt.show()
            #plt.savefig("refpics/RUN2/PV/ns"+str(i)+".png")


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

def plotProfile(p):
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.plot(p.itemps,p.ipres)
    ax2.plot(p.isals,p.ipres)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    fig.suptitle(str(p.eyed)+" lat: "+str(p.lat)+" lon: "+ str(p.lon))
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()

def surfacesToXYZ(surfaces):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum = 'WGS84',preserve_units=True)
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum = 'WGS84',preserve_units=True)
    newsurfaces = {}
    for k in surfaces.keys():
        tempSurf = np.array([[],[],[[],[],[],[]],[]])
        for l in range(len(surfaces[k][0])):
            x,y,z = pyproj.transform(lla,ecef,surfaces[k][0][l],surfaces[k][1][l],surfaces[k][2][0][l],radians=False)
            tempSurf[0].append(x)
            tempSurf[1].append(y)
            tempSurf[2][0].append(z)
            tempSurf[2][1].append(surfaces[k][2][1][l])
            tempSurf[2][2].append(surfaces[k][2][2][l])
            tempSurf[2][3].append(surfaces[k][2][3][l])
            tempSurf[3].append(surfaces[k][3][l])
        newsurfaces[k]=tempSurf
    return newsurfaces

def deduplicateXYZ(x,y,z):
    return zip(*set(list(zip(x,y,z))))

def removeDiscontinuities(x,y,z,radius=40):
    x=np.asarray(x)
    y=np.asarray(y)
    z=np.asarray(z)
    final = np.zeros(x.shape)
    for i in range(len(x)):
        if final[i] == False:
            r = np.sqrt((x- x[i])**2 + (y - y[i])**2)
            inside = r<radius*1000
            inside[i]=False
            s = 0
            counter = 0
            for t in range(len(inside)):
                if t:
                    s+=z[i]
                    counter+=1
            if counter >0:
                z[i] = s/counter
                    
            #print(r)
            if np.count_nonzero(final) == 0  :
                final =inside
            else:
                final = np.logical_or(final,inside)
    return x[np.invert(final)],y[np.invert(final)],z[np.invert(final)] 


def createMesh(n,xvals,yvals):
    return np.meshgrid(np.linspace(np.min(xvals),np.max(xvals),n), np.linspace(np.min(yvals),np.max(yvals),n),indexing="xy")

def generateMaskedMesh(x,y):
    xi,yi = createMesh(50,x,y)
    final = np.zeros(xi.shape)
    for i in range(len(x)):
        r = np.sqrt((xi- x[i])**2 + (yi - y[i])**2)
        inside = r<50*1000
        if np.count_nonzero(final) == 0  :
            final =inside
        else:
            final = np.logical_or(final,inside)
    return xi[final],yi[final]


def interpolateSurface(x,y,z,d=None):
    if d:
        print(len(d))
        print(len(x))
        print(len(y))
        print(len(z))
        qrbfi = Rbf(x,y,z,d,function="thin_plate",smooth=10.0)
        rbfi = Rbf(x,y,z,function="thin_plate",smooth=10.0)
        xi,yi = generateMaskedMesh(x,y)
        zi = rbfi(xi,yi)
        ti = qrbfi(xi,yi,zi)
        return xi,yi,zi,ti
    else:
        rbfi = Rbf(x,y,z,function="thin_plate",smooth=10.0)
        xi,yi = generateMaskedMesh(x,y)
        zi = rbfi(xi,yi)
        return xi,yi,zi

def xyzToSurface(x,y,z,d):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum = 'WGS84',preserve_units=True)
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum = 'WGS84',preserve_units=True)
    surface={}
    surface[d] = np.array([[],[],[[],[],[],[]],[]])
    lon, lat, a = pyproj.transform(ecef,lla,x,y,z,radians=False)
    surface[d][0]=lon
    surface[d][1]=lat
    surface[d][2]=[a]
    return surface
        



#def interpolateSurface(surface):
    #for i in surfa
