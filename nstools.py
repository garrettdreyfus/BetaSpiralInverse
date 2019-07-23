import csv
from mpl_toolkits.basemap import Basemap
from numpy import linspace
import matplotlib.pyplot as plt
from numpy import meshgrid
import numpy as np
from netCDF4 import Dataset
import json
from geopy.distance import geodesic
from geopy.distance import great_circle
from profile import Profile
import random
from scipy.interpolate import Rbf
import pyproj
import bathtools
import copy
import gsw
import pygam
import pickle
#from gswmatlab.pyinterface import geo_strf_isopycnal
import scipy.io as sio
import itertools
from scipy.linalg import svd
from sklearn.decomposition import TruncatedSVD
from progress.bar import Bar
import parametertools as ptools


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

def emptySurface():
    return {"lats":[],"lons":[],"ids":[],"data":{"pres":[],"t":[],"s":[]}}

def peerSearch(profiles,deepestindex,depth,profilechoice,radius=500):
    surfaces = {}
    surfaces[depth]= emptySurface()
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
                    surfaces[depth]["lons"].append(p.lon)
                    surfaces[depth]["lats"].append(p.lat)
                    surfaces[depth]["data"]["pres"].append(ns)
                    surfaces[depth]["ids"].append(p.eyed)
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
        surfaces[r]=emptySurface()
    #Ever profile
    for j in range(len(profiles)):
        #every neutral surface
        for r in deeprange:
            ns = profiles[deepestindex].neutralDepth(profiles[j],r) 
            if ns != None:
                surfaces[r]["lons"].append(profiles[j].lon)
                surfaces[r]["lats"].append(profiles[j].lat)
                surfaces[r]["data"]["pres"].append(ns)
    return surfaces

def addDataToSurfaces(profiles,surfaces,stdevs,debug=True):
    if debug:
        print("adding data to surfaces")
    tempSurfs = {}
    negativecount = 0
    for k in surfaces.keys():
        tempSurf = emptySurface()
        for l in range(len(surfaces[k]["lons"])):
            p = getProfileById(profiles,surfaces[k]["ids"][l])
            p.interpolate()
            t,s = p.atPres(surfaces[k]["data"]["pres"])
            pv = p.potentialVorticity(surfaces[k]["data"]["pres"][l],debug=False)
            if pv and pv<0:
                negativecount +=1 
            if t and s and pv and p and pv != np.Inf and pv != np.nan and not np.isnan(t):
                tempSurf["lons"].append(surfaces[k]["lons"][l])
                tempSurf["lats"].append(surfaces[k]["lats"][l])
                tempSurf["data"]["pres"].append(-surfaces[k]["data"]["pres"][l])
                tempSurf["data"]["t"].append(t)
                tempSurf["data"]["s"].append(s)
                tempSurf["data"]["pv"].append(pv)
                tempSurf["ids"].append(surfaces[k][3][l])
        if len(tempSurf[0])>5:
            tempSurfs[k] = tempSurf

        print("ns: ",k," negative count: ",negativecount)
    return tempSurfs


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

def filterSurfacesByLine(surfaces,lon,radius=20):
    print("filtering by crosssection")
    for k in surfaces.keys():
        distFilter = np.zeros(len(surfaces[k]["lons"]), dtype=bool)
        for index in range(len(surfaces[k]["lons"])):
            p = (surfaces[k]["lats"][index],surfaces[k]["lons"][index])
            if ((np.cos(p[0]*(np.pi/180)))*abs(((p[1]+180)%180)-lon)*111)<radius:
                #print(((np.cos(lpoint[0]*(np.pi/180)))*abs(p[1]-lpoint[1])*111))
                distFilter[index] = True
        surfaces[k]["lons"] = np.asarray(surfaces[k]["lons"])[distFilter]
        surfaces[k]["lats"] = np.asarray(surfaces[k]["lats"])[distFilter]
        for j in surfaces[k]["data"].keys():
            surfaces[k]["data"][j] = np.asarray(surfaces[k]["data"][j])[distFilter]
        if len(surfaces[k]["ids"])>0:
            surfaces[k]["ids"] = np.asarray(surfaces[k]["ids"])[distFilter]
    return surfaces
        
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


def addXYToSurfaces(surfaces,debug=True):
    if debug:
        print("converting surfaces to xyz")
    newsurfaces = {}
    for k in surfaces.keys():
        x,y = homemadeXY(surfaces[k]["lons"],surfaces[k]["lats"])
        surfaces[k]["x"]=x
        surfaces[k]["y"]=y
    return surfaces

def deduplicateXYZ(x,y,z):
    return zip(*set(list(zip(x,y,z))))

def removeDiscontinuities(surface,radius=10,debug=True):
    if debug:
        print("removing discontinuities")
    x=np.asarray(surface["x"])
    y=np.asarray(surface["y"])
    z=np.asarray(surface["data"]["pres"])
    final = np.zeros(x.shape)
    #print(x)
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
    final = np.invert(final)
    for k in surface.keys():
        if k == "data":
            for d in surface[k]:
                surface[k][d] = np.asarray(surface[k][d])[final]
        else:
            surface[k] = np.asarray(surface[k])[final]
    return surface

def createMesh(n,xvals,yvals,custom=False):
    if custom:
        return np.meshgrid(np.linspace(np.min(xvals),np.max(xvals),n), np.linspace(np.min(yvals),np.max(yvals),n),indexing="xy")
    else:
        xmin=-1793163
        xmax=971927
        ymin=-1455096
        ymax=1200385
        return np.meshgrid(np.linspace(xmin,xmax,n), np.linspace(ymin,ymax,n),indexing="xy")


def generateMaskedMesh(x,y,radius=150):
    xi,yi = createMesh(125,x,y)
    final = np.zeros(xi.shape)
    neighbors = []
    for i in range(len(x)):
        r = np.sqrt((xi- x[i])**2 + (yi - y[i])**2)
        inside = r<radius*1000
        if np.count_nonzero(final) == 0  :
            final =inside
        else:
            final = final+inside
    for i in range(len(final[0])):
        if final[0][i] > 2:
            final[0][i]=True
        else:
            final[0][i] = False

    indexcount = np.full(final.shape,np.nan)
    count = 0
    for i in range(final.shape[0]):
        for j in range(final.shape[1]):
            if final[i][j]:
                indexcount[i][j] = count
                count+=1

    finalxi=[]
    finalyi=[]
    finalneighbors = []
    finalids=[]
    for l in range(final.shape[0]):
        for k in range(final.shape[1]):
            if final[l][k]:
                finalxi.append(xi[l][k])
                finalyi.append(yi[l][k])
                finalids.append(l*125+k)
            if l != final.shape[0]-1 and k != final.shape[1]-1:
                s = []
                if final[l][k] and final[l+1][k+1] and final[l+1][k] and final[l][k+1]:
                    s.append(int(indexcount[l][k]))
                    s.append(int(indexcount[l][k+1]))
                    s.append(int(indexcount[l+1][k]))
                    s.append(int(indexcount[l+1][k+1]))
                    finalneighbors.append(tuple(s))

    return np.asarray(finalxi),np.asarray(finalyi),finalneighbors,finalids

def removeOutlierSurfaces(surfaces,stdevs=2):
    for k in surfaces.keys():
        m = np.median(surfaces[k]["data"]["pres"])
        s = np.std(surfaces[k]["data"]["pres"])
        filt = np.where(np.abs(surfaces[k]["data"]["pres"]-m)<s*stdevs)
        for field in surfaces[k].keys():
            if field=="data":
                for datafield in surfaces[k][field].keys():
                    surfaces[k][field][datafield]=surfaces[k][field][datafield][filt]
                surfaces[k][field]=surfaces[k][field][filt]
            else:
                surfaces[k][field]=surfaces[k][field][filt]
    return surfaces
    
def interpolateSurface(surface,debug=True):
    #print("######")
    interpsurf={}
    X = np.zeros((len(surface["x"]),2))
    X[:,0]=surface["x"]
    X[:,1]=surface["y"]
    xi,yi,neighbors,finalids = generateMaskedMesh(surface["x"],surface["y"])
    interpdata={}
    interpsurf["x"] =xi
    interpsurf["y"] =yi
    interpsurf["ids"] =finalids
    if len(xi) != len(finalids):
        print("OH NOOOOOO")
    for k in surface["data"].keys():
        notnan = ~np.isnan(surface["data"][k])
        if len(notnan)>10:
            gam = pygam.LinearGAM(pygam.te(0,1)).fit(X,surface["data"][k][notnan])
            Xgrid = np.zeros((yi.shape[0],2))
            Xgrid[:,0] = xi
            Xgrid[:,1] = yi
            interpdata[k] = gam.predict(Xgrid)
        else:
            interpdata[k] = np.asarray([np.nan]*len(xi))
    interpsurf["data"] = interpdata
    interpsurf["data"]["ids"] = finalids
    interpsurf = addLatLonToSurface(interpsurf)
    return interpsurf,neighbors

def homemadeXY(lon,lat):
    x=[]
    y=[]
    for i in range(len(lat)):
        theta = np.deg2rad(lon[i])
        r = ((90-lat[i]) *111*1000)
        x.append(r*np.cos(theta))
        y.append(r*np.sin(theta))
    return np.asarray(x),np.asarray(y)

def addLatLonToSurface(surface,debug = True):
    if debug:
        print("converting xyz back to lat lon a")
    lat = 90-(np.sqrt((surface["x"]**2+surface["y"]**2))/111000.0)
    lon = np.degrees(np.arctan2(surface["y"],surface["x"]))
    #print("lat: ",lat, " lon: ",lon)
    surface["lons"]=lon
    surface["lats"]=lat
    return surface
     
        
def findNeighboringPoints(profiles,lat,lon,radius=30):
    lats =[]
    lons = []
    eyeds= []
    for p in profiles:
        if geodesic((p.lat,p.lon),(lat,lon)).km<radius:
            lats.append(p.lat)
            lons.append(p.lon)
            eyeds.append(p.eyed)
    plt.scatter(lons,lats)
    plt.show()



def nanCopySurfaces(surfaces):
    nancopy = {}
    for k in surfaces.keys():
        nancopy[k]=emptySurface()
        nancopy[k]["lons"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["ids"] = surfaces[k]["ids"]
        nancopy[k]["lats"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["x"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["y"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["u"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["ids"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["v"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["hx"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["h"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["CKVO"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["CKVB"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["hy"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["t"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["s"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["pv"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["pres"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["curl"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["uabs"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["vabs"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["uprime"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["vprime"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["dsdx"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["data"]["dsdy"] = np.full(len(surfaces[k]["lons"]),np.nan)
    return nancopy

def addHToSurfaces(surfaces):
    for k in surfaces.keys():
        surfaces[k]["data"]["h"] = np.full(len(surfaces[k]["lons"]),np.nan)
    minimum = int(np.min(list(surfaces.keys())))
    maximum = int(np.max(list(surfaces.keys())))
    depths = range(minimum,maximum+1,200)
    for j in range(len(depths))[1:-1]:
        for index in range(len(surfaces[depths[j]]["x"])):
            if index%100==0:
                print("adding thickness: ",index,"/",len(surfaces[depths[j]]["x"]),"on NS: ",depths[j])
            eyed = int(surfaces[depths[j]]["ids"][index])
            foundbelow = np.where(np.asarray(surfaces[depths[j-1]]["ids"])==eyed)
            found = index
            foundabove = np.where(np.asarray(surfaces[depths[j+1]]["ids"])==eyed)
            if len(foundbelow)!=0 and len(foundbelow[0]) != 0 and len(foundabove)!=0 and len(foundabove[0]) != 0:
                foundbelow = foundbelow[0][0]
                foundabove = foundabove[0][0]
                tophalf = abs(surfaces[depths[j+1]]["data"]["pres"][foundabove]-surfaces[depths[j]]["data"]["pres"][found])/2.0
                bothalf = abs(surfaces[depths[j]]["data"]["pres"][found]-surfaces[depths[j-1]]["data"]["pres"][foundbelow])/2.0
                surfaces[depths[j]]["data"]["h"][found] = tophalf + bothalf

    return surfaces

def addK(surfaces,cachename=None):
    for k in surfaces.keys():
        surfaces[k]["data"]["CKVB"] = np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["CKVO"] = np.full(len(surfaces[k]["lons"]),np.nan)
        for i in Bar("K "+str(k)+': ').iter(range(len(surfaces[k]["lons"]))):
            lat = surfaces[k]["lats"][i]
            lon = surfaces[k]["lons"][i]
            if not (np.isnan(lat) or np.isnan(lon)):
                pv = surfaces[k]["data"]["pv"][i]
                pres = surfaces[k]["data"]["pres"][i]
                KV = ptools.Kv(lat,lon,pv,pres,cachename)
                surfaces[k]["data"]["CKVO"][i] = KV[0]
                surfaces[k]["data"]["CKVB"][i] = KV[1]
    return surfaces
            



                
 
    
def addPrimes(surfaces,neighbors,distances,debug=False):
    staggered = nanCopySurfaces(surfaces)
    alldxs = []

    for k in Bar('Adding Primes: ').iter(neighbors.keys()):
        for s in neighbors[k]:
            if surfaces[k]["x"][s[0]]>surfaces[k]["x"][s[1]]:
                print("This should never happen")
            if surfaces[k]["y"][s[0]]>surfaces[k]["y"][s[2]]:
                print("This should never happen")

            s=np.asarray(s)
            staggered= averageOverNeighbors(staggered,surfaces,k,s)
            staggered = addGradients(staggered,surfaces,k,s,distances)

    return staggered

def addGradients(staggered,surfaces,k,s,distances):
    #NS thickness slope
    dhdx = []
    dhdy = []
    dhdx.append((surfaces[k]["data"]["h"][s[1]]-surfaces[k]["data"]["h"][s[0]])/distances[k][(s[0],s[1])])
    dhdx.append((surfaces[k]["data"]["h"][s[3]]-surfaces[k]["data"]["h"][s[2]])/distances[k][(s[2],s[3])])
    dhdy.append((surfaces[k]["data"]["h"][s[2]]-surfaces[k]["data"]["h"][s[0]])/distances[k][(s[0],s[2])])
    dhdy.append((surfaces[k]["data"]["h"][s[3]]-surfaces[k]["data"]["h"][s[1]])/distances[k][(s[1],s[3])])
    #stream function gradient
    dpsidx = []
    dpsidy = []
    dpsidx.append((surfaces[k]["data"]["psi"][s[1]]-surfaces[k]["data"]["psi"][s[0]])/distances[k][(s[0],s[1])])
    dpsidx.append((surfaces[k]["data"]["psi"][s[3]]-surfaces[k]["data"]["psi"][s[2]])/distances[k][(s[2],s[3])])
    dpsidy.append((surfaces[k]["data"]["psi"][s[2]]-surfaces[k]["data"]["psi"][s[0]])/distances[k][(s[0],s[2])])
    dpsidy.append((surfaces[k]["data"]["psi"][s[3]]-surfaces[k]["data"]["psi"][s[1]])/distances[k][(s[1],s[3])])
    #salt gradient 
    dsdx = []
    dsdy = []
    dsdx.append((surfaces[k]["data"]["s"][s[1]]-surfaces[k]["data"]["s"][s[0]])/distances[k][(s[0],s[1])])
    dsdx.append((surfaces[k]["data"]["s"][s[3]]-surfaces[k]["data"]["s"][s[2]])/distances[k][(s[2],s[3])])
    dsdy.append((surfaces[k]["data"]["s"][s[2]]-surfaces[k]["data"]["s"][s[0]])/distances[k][(s[0],s[2])])
    dsdy.append((surfaces[k]["data"]["s"][s[3]]-surfaces[k]["data"]["s"][s[1]])/distances[k][(s[1],s[3])])
    staggered[k]["data"]["hx"][s[0]] = np.mean(dhdx)
    staggered[k]["data"]["hy"][s[0]] = np.mean(dhdy)
    staggered[k]["data"]["u"][s[0]] = (1/gsw.f(surfaces[k]["lats"][s[0]]))*np.mean(dpsidy)
    staggered[k]["data"]["v"][s[0]] = (-1/gsw.f(surfaces[k]["lats"][s[0]]))*np.mean(dpsidx)
    staggered[k]["data"]["dsdx"][s[0]] = np.mean(dsdx)
    staggered[k]["data"]["dsdy"][s[0]] =  np.mean(dsdy)
    return staggered


def averageOverNeighbors(staggered,surfaces,k,s):
    lon = np.mean(np.abs(surfaces[k]["lons"][s]))*np.sign(surfaces[k]["lons"][s[0]])
    staggered[k]["lons"][s[0]] = lon
    staggered[k]["data"]["pres"][s[0]] = np.mean(surfaces[k]["data"]["pres"][s])
    staggered[k]["data"]["ids"][s[0]] = surfaces[k]["ids"][s[0]]
    staggered[k]["lats"][s[0]] = np.mean(surfaces[k]["lats"][s])
    x = np.mean(surfaces[k]["x"][s])
    staggered[k]["x"][s[0]] = x
    staggered[k]["data"]["h"][s[0]] = np.mean(surfaces[k]["data"]["h"][s])
    staggered[k]["data"]["pv"][s[0]] = np.mean(surfaces[k]["data"]["pv"][s])
    y = np.mean(surfaces[k]["y"][s])
    staggered[k]["y"][s[0]] = y
    return staggered



def trueDistanceLookup(surface,neighbors):
    lookup = {}
    for square in neighbors:
        for edge in itertools.combinations(square,2):
            p = tuple(sorted(edge))
            if p not in lookup.keys():
                lookup[p] = geodesic((surface["lats"][p[0]],surface["lons"][p[0]]),(surface["lats"][p[1]],surface["lons"][p[1]])).m
    return lookup

def addStreamFunc(surfaces,profiles):
    neutraldepths={}
    for k in surfaces.keys():
        surfaces[k]["data"]["psi"]=np.empty(np.size(surfaces[k]["ids"]))
        surfaces[k]["data"]["psi"][:] =np.nan
        for i in range(len(surfaces[k]["ids"])):
            if surfaces[k]["ids"][i] not in neutraldepths:
                neutraldepths[surfaces[k]["ids"][i]] =[[],[]]
            if abs(k)<3700:
                neutraldepths[surfaces[k]["ids"][i]][0].append(k)
                neutraldepths[surfaces[k]["ids"][i]][1].append(abs(surfaces[k]["data"]["pres"][i]))
    refns = []
    refns_p = []
    for k in neutraldepths.keys():
        if len(neutraldepths[k][0])==18:
            p = getProfileById(profiles,k)
            for j in range(len(neutraldepths[k][0])):
                refns.append(neutraldepths[k][0][j]) 
                refns_p.append(p.sigma2(neutraldepths[k][1][j]))    
            break

    refns = np.asarray(refns)
    refns_d = np.asarray(refns_p)
    s = []
    t = []
    ip = []
    p_ref = []
    ns = []
    ns_p = []
    ks = []
    count =0
    for k in neutraldepths.keys():
        count+=1
        if count > 2000:
            break
        p = getProfileById(profiles,k)
        s.append(p.isals)
        t.append(p.itemps)
        ip.append(np.abs(p.ipres))
        #p_ref.append([10.1235]*len(p.ipres))
        nslabels = neutraldepths[k][0]
        ns.append(np.abs(nslabels[::-1]))
        nsa = np.abs(neutraldepths[k][1])
        #print(nslabels,nsa)
        ns_p.append(nsa)
        #print("refns_d", refns_d)
        #print("refns", refns)
        #print("nslabels", nslabels)
        #print("isin", np.isin(refns,nslabels))
        #print("where isin", np.where(np.isin(refns,nslabels)))
        nsdensref = refns_d[np.where(np.isin(refns,nslabels))[0]]
        #print(len(nsdensref),len(nsa))
        #print("nsdensref",nsdensref)
        psi = p.geoIsopycnal(nsa,nsdensref)
        for depth in range(len(nslabels)):
            #print(psi)
            #print(surfaces[nslabels[depth]]["ids"])
            targets = np.where(np.asarray(surfaces[nslabels[depth]]["ids"]) ==str(k) )
            #print(targets)
            surfaces[nslabels[depth]]["data"]["psi"][targets] = psi[depth]
        p_ref.append([0]*len(p.ipres))
        ks.append(k)
     
    #results = geo_strf_isopycnal(s,t,ip,p_ref,ns,ns_p,ks)
    
    #with open('data/geoisopycnal.pickle', 'wb') as outfile:
        #pickle.dump([results,ks], outfile)
    #print(results)
    return surfaces

def addStreamFuncFromFile(surfaces,profiles,isopycnalfile,referencefile):
    psi = np.asarray(sio.loadmat(isopycnalfile)["geoisopycnal"]).transpose()
    ids = np.asarray(sio.loadmat(referencefile)["ks"])
    tags = np.asarray(sio.loadmat(referencefile)["ns"]).transpose()
    for k in surfaces.keys():
        surfaces[k][2]["psi"] = np.full_like(surfaces[k]["data"]["pres"])
        for i in range(len(surfaces[k]["lons"])):
            if surfaces[k]["ids"][i] in ids:
                col = np.where(ids == surfaces[k]["ids"][i])[0][0]
                row = np.where(k == tags[col])[0]
                if len(row) >0:
                    print(psi[col][row])
                    surfaces[k]["data"]["psi"].append(psi[col][row[0]])
                else:
                    surfaces[k]["data"]["psi"].append(np.nan)
                
            else:
                surfaces[k]["data"]["psi"].append(np.nan)
    return surfaces

def convertOldSurfaces(surfaces):
    newsurfaces = {}
    for k in surfaces.keys():
        surface = {}
        surface["lons"] = surfaces[k][0]
        surface["lats"] = surfaces[k][1]
        surface["data"] = {"pres":surfaces[k][2][0],"t":surfaces[k][2][1],"s":surfaces[k][2][0],"pv":surfaces[k][2][3]}
        surface["ids"] = surfaces[k][3]
        newsurfaces[k] = surface
    return newsurfaces
        


