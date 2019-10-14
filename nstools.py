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
import graph
import copy
import gsw
import pickle
#from gswmatlab.pyinterface import geo_strf_isopycnal
import scipy.io as sio
import itertools
from scipy.linalg import svd
from sklearn.decomposition import TruncatedSVD
from progress.bar import Bar
import parametertools as ptools
from prettytable import PrettyTable
import pdb
import interptools

#From a list of filenames extract a bunch of profile objects
#also returns the profile with the deepestindex because that may be useful
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


##finds deepest profile from given profiels
def deepestProfile(profiles):
    deepestindex=-1
    deepestdepth=-1
    for p in range(len(profiles)):
        if profiles[p].ipres[-1] >deepestdepth:
            deepestdepth=profiles[p].ipres[-1]
            deepestindex=p
    return deepestindex

##finds number of unqie cruises
def cruiseCount(profiles):
    cruises ={"null"}
    for profile in profiles:
        if profile.cruise not in cruises:
            cruises.add(profile.cruise)
    return cruises

##returns all profiles from a cruise
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

## finds all profiles in a small little slice of lat and lon
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

        

##extract profiles only from certain months
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

## extract profiles within a certain box
def extractProfilesBox(fnames,lonleft,lonright,latbot,lattop):
    ##Load JSON data into profile objects
    ##and return index of deepest one
    profiles = []
    deepestindex =[] 
    deepestdepth = 0 
    for fname in fnames:
        json_file = open(fname) 
        data = json.load(json_file)
        for p in Bar('Loading profiles in box:   ').iter(data.keys()):
            if len(data[p]["pres"])>10:
                profile = Profile(p,data[p])
                if latbot <= profile.lat <= lattop and lonleft <= profile.lon <= lonright:
                    if len(profile.ipres)>0:
                        profiles.append(profile)
                        if data[p]["pres"][-1] > deepestdepth:
                            deepestindex = len(profiles)-1
                            deepestdepth=data[p]["pres"][-1]
    return profiles, deepestindex

#I don't like the norwegian sea and neither should you !!!
## this removes profiles that lie within it
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

##finds the closest profile from a list of profiles in which
# a neutral surface has already been found
def closestIdentifiedNS(profiles,queryprofile,depth,radius,speedy=False):
    if not hasattr(closestIdentifiedNS,"cache"):
        closestIdentifiedNS.cache = {}
    minimumdistance = radius
    minprofile = None
    for p in profiles:
        if depth in p.neutraldepth.keys():
            eyeds = tuple(sorted((queryprofile.eyed,p.eyed)))
            if eyeds not in closestIdentifiedNS.cache.keys():
                dist = great_circle((queryprofile.lat,queryprofile.lon),(p.lat,p.lon)).km
                closestIdentifiedNS.cache[eyeds] = dist
            else:
                dist = closestIdentifiedNS.cache[eyeds]
            if dist<minimumdistance:
                minimumdistance = dist
                minprofile = p
            #elif minimumdistance != radius:
                #print(dist)
    #print(minprofile.neutraldepth)
    return minprofile

## find profiles in a certain location
def profileInBox(profiles,lonleft,lonright,latbot,lattop,depth):
    results=[]
    for p in profiles:
        if lonleft< p.lon <lonright and latbot < p.lat < lattop\
            and np.max(abs(p.pres))>depth:
            results.append(p)
    return results

##create an empty surface
def emptySurface():
    return {"lats":[],"lons":[],"ids":[],\
            "data":{"pres":[],"t":[],"s":[],"pv":[],"n^2":[],"alpha":[],\
            "beta":[]}}

##given a seed profile, finds neutral surfaces by comparing each profile 
## to the nearest profile with an identified neutral surface. Only quits
## when all remaining profiles either are not close enough to a profile with a 
## known NS or a NS is not able to be determined for them
def peerSearch(profiles,depth,profilechoice,radius=500):
    surfaces = {}
    surfaces[depth]= emptySurface()
    profilechoice.neutraldepth[depth] = depth
    references=[]
    references.append(profilechoice)
    while len(profiles)>0:
        foundcounter = 0
        closestcounter = 0
        for p in profiles.copy():
            closest = closestIdentifiedNS(references,p,depth,radius)
            if closest:
                closestcounter+=1
                #print(closest.neutraldepth)
                ns = closest.neutralDepth(p,closest.neutraldepth[depth],depthname=depth,searchrange=500) 
                if ns != None:
                    surfaces[depth]["lons"].append(p.lon)
                    surfaces[depth]["lats"].append(p.lat)
                    surfaces[depth]["data"]["pres"].append(ns)
                    surfaces[depth]["ids"].append(p.eyed)
                    profiles.remove(p)
                    foundcounter +=1
                    #references.append(p)
        if foundcounter ==0:
            break

        print("found: ,",foundcounter,"left: ", len(profiles))
        print("close enough: ",closestcounter)
        #plotProfiles(references,"ITS SPREADING",specialprofile = profilechoice)
    return surfaces

## runs the peer search on every neutral surface essentially
def runPeerSearch(profiles,shallowlimit,deeplimit,surfacestep,profilechoice,radius):
    surfaces = {}
    for d in range(shallowlimit,deeplimit,surfacestep):
        print("NSearching: ",d)
        surfaces.update(peerSearch(profiles.copy(),d,profilechoice,radius))
    return surfaces

## faultyway of finding neutral surfaces in the arctic
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

##find the index of the point that lies directly above you on surface
def findAboveIndex(surfaces,k,l):
    if k-200 in surfaces.keys() and k+200 in surfaces.keys():
        eyed = surfaces[k]["ids"][l]
        above = np.where(np.asarray(surfaces[k-200]["ids"]) == eyed)[0]
        below = np.where(np.asarray(surfaces[k+200]["ids"]) == eyed)[0]
        if len(above)>0 and len(below)>0:
            middle = surfaces[k]["data"]["pres"][l]
            above = surfaces[k-200]["data"]["pres"][above[0]]
            below = surfaces[k+200]["data"]["pres"][below[0]]
            return abs(above+middle)/2,abs(below+middle)/2
    return None,None

## calculates data throughout neutral surfaces, and some that require calculation 
## shoots that all into a new surfaces object and returns it
def addDataToSurfaces(profiles,surfaces,stdevs,debug=True):
    tempSurfs = {}
    for k in Bar("Adding data to: ").iter(surfaces.keys()):
        negativecount = 0
        tempSurf = emptySurface()
        monthcount = [0]*13
        monthnegativecount = [0]*13
        for l in range(len(surfaces[k]["lons"])):
            p = getProfileById(profiles,surfaces[k]["ids"][l])
            if p and not np.isnan(p.isals).any() :
                if (surfaces[k]["data"]["pres"][l]+35 in p.ipres and surfaces[k]["data"]["pres"][l]-35 in p.ipres):
                    pv = p.potentialVorticityAtHautala(surfaces[k]["data"]["pres"][l])
                    t,s = p.betweenPres(surfaces[k]["data"]["pres"][l]+35,surfaces[k]["data"]["pres"][l]-35)
                else:
                    pv = np.nan
                    t,s = p.atPres(surfaces[k]["data"]["pres"][l])
                #monthcount[p.time.month]=monthcount[p.time.month]+1
                if pv and pv<0:
                    pv=0
                    #monthnegativecount[p.time.month]=monthnegativecount[p.time.month]+1
                    negativecount +=1 
                if type(pv) == type(None):
                    pv = np.nan
                if ~np.isnan(t) and ~np.isnan(s) and ~np.isnan(pv) and p and pv != np.Inf:
                    tempSurf["lons"].append(surfaces[k]["lons"][l])
                    tempSurf["lats"].append(surfaces[k]["lats"][l])
                    tempSurf["data"]["pres"].append(abs(surfaces[k]["data"]["pres"][l]))
                    tempSurf["data"]["t"].append(t)
                    tempSurf["data"]["s"].append(s)
                    tempSurf["data"]["pv"].append(pv)
                    tempSurf["ids"].append(surfaces[k]["ids"][l])
                    tempSurf["data"]["n^2"].append(pv*(9.8/gsw.f(p.lat)))
                    tempSurf["data"]["alpha"].append(gsw.alpha(s,t,surfaces[k]["data"]["pres"][l]))
                    tempSurf["data"]["beta"].append(gsw.beta(s,t,surfaces[k]["data"]["pres"][l]))


                    
        if len(tempSurf["lats"])>5:
            tempSurfs[k] = tempSurf
        print("\n###########"+str(k)+"#################")
        print(monthcount)
        print(monthnegativecount)
        print("############################")
        if len(surfaces[k]["lons"])>0:
            print("ns: ",k," negative count: ",negativecount/len(surfaces[k]["lons"]),"pv mean:" ,np.mean(tempSurf["data"]["pv"]))
    tempSurfs = addStreamFunc(tempSurfs,profiles)
    return tempSurfs


##return a profile given by an id
def getProfileById(profiles,eyed):
    for p in profiles:
        if p.eyed == eyed:
            return p

##remove cruises that are not mentioned in cruisenames
def filterCruises(profiles,cruisenames):
    finalprofiles = []
    for p in profiles:
        if p.cruise in cruisenames:
            finalprofiles.append(p)
    return finalprofiles

##find points within certain distance from line
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
        
##plot profile t and s
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
    
#just prints out the percentages of each quantity that are nan along a surface 
def surfaceDiagnostic(surfaces):
    diagnostics ={}
    t = PrettyTable(['Property', 'nan%'])
    for k in surfaces.keys():
        for d in surfaces[2000]["data"].keys():
            if d not in diagnostics.keys():
                diagnostics[d]=[0,0]
            diagnostics[d][0] = diagnostics[d][0] + np.count_nonzero(np.isnan(surfaces[k]["data"][d]))
            diagnostics[d][1] = diagnostics[d][1] + np.count_nonzero(~np.isnan(surfaces[k]["data"][d]))
    for d in diagnostics.keys():
        if sum(diagnostics[d])!=0:
            percentage = int(round(diagnostics[d][0]/sum(diagnostics[d]),3)*100)
            t.add_row([d,percentage])
    print(t)


##create a nan copy of a surface
def nanCopySurfaces(surfaces):
    nancopy = {}
    for k in surfaces.keys():
        nancopy[k]=emptySurface()
        nancopy[k]["lons"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["ids"] = surfaces[k]["ids"]
        nancopy[k]["lats"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["x"] = np.full(len(surfaces[k]["lons"]),np.nan)
        nancopy[k]["y"] = np.full(len(surfaces[k]["lons"]),np.nan)
        datafields = ["u","v","hx","h","CKVB","hy","t","s","pv","pres",\
                     "curl","uabs","vabs","uprime","vprime","dsdx","dsdz","dsdy",\
                    "d2sdx2","d2sdy2","dtdx","dtdy","dpdx","dpdy","n^2",\
                    "dqnotdx","dqnotdy","d2thetads2","dalphadtheta",\
                    "alpha","beta","dalphads","dbetads","dalphadp",\
                    "dbetadp","psi","psinew","dqdz","dqdx","dqdy",\
                    "d2qdz2","d2qdx2","d2qdy2","khp","khpterm","khpdz","toph",\
                    "both","dpsidx","dpsidy","ssh","ids","z","knownu"\
                    ,"knownv","diffkr","kapgm","kapredi","nsdiff"]
        for d in datafields:
            nancopy[k]["data"][d] = np.full(len(surfaces[k]["lons"]),np.nan)
    return nancopy
#perform a vertical gradient of a quantity of a surface
#WITH RESPECT TO PRESSURE NOT Z
# out is a surface you are adding to 
#data is a surface with data
#depths is a dictionary of indexs to depths
## above center and below are the indexs of a point at three different levels
## attr is the quantity being derived, outattr is where you should store the result
def vertGrad(out,data,depths,k,above,center,below,attr,outattr,factor=1):
    dattr = data[depths[k+1]]["data"][attr][below]-data[depths[k-1]]["data"][attr][above]
    dz = data[depths[k+1]]["data"]["pres"][below]-data[depths[k-1]]["data"]["pres"][above]
    out[depths[k]]["data"][outattr][center] = factor * dattr/dz
    return out

##calculate the height of each neutral surface by the difference in pressures
## of each neutral surface
def addHeight(surfaces):    
    minimum = int(np.min(list(surfaces.keys())))
    maximum = int(np.max(list(surfaces.keys())))
    depths = range(minimum,maximum+1,200)
    for j in Bar('Adding Heights:   ').iter(range(len(depths))[1:-1]):
        for index in range(len(surfaces[depths[j]]["x"])):
            eyed = int(surfaces[depths[j]]["ids"][index])
            foundbelow = np.where(np.asarray(surfaces[depths[j+1]]["ids"])==eyed)
            found = index
            foundabove = np.where(np.asarray(surfaces[depths[j-1]]["ids"])==eyed)
            if len(foundbelow)!=0 and len(foundbelow[0]) != 0 and len(foundabove)!=0 and len(foundabove[0]) != 0:
                foundbelow = foundbelow[0][0]
                foundabove = foundabove[0][0]
                tophalf = abs(surfaces[depths[j-1]]["data"]["pres"][foundabove]-surfaces[depths[j]]["data"]["pres"][found])/2.0
                bothalf = abs(surfaces[depths[j]]["data"]["pres"][found]-surfaces[depths[j+1]]["data"]["pres"][foundbelow])/2.0
                surfaces[depths[j]]["data"]["h"][found] = tophalf + bothalf
                surfaces[depths[j]]["data"]["toph"][found] = tophalf
                surfaces[depths[j]]["data"]["both"][found] = bothalf
                surfaces[depths[j]]["data"]["dbetadp"][found] = (surfaces[depths[j-1]]["data"]["beta"][foundabove] - surfaces[depths[j+1]]["data"]["beta"][foundbelow])/((tophalf + bothalf)*2)
    return surfaces
#calculates and stores all the necessary vertical gradients
#then proceeds to calculate all the double vertical gradients
def addVerticalGrad(surfaces): 
    minimum = int(np.min(list(surfaces.keys())))
    maximum = int(np.max(list(surfaces.keys())))
    depths = range(minimum,maximum+1,200)
    for j in Bar('Adding Vertical Gradients:   ').iter(range(len(depths))[1:-1]):
        for index in range(len(surfaces[depths[j]]["x"])):
            eyed = int(surfaces[depths[j]]["ids"][index])
            foundbelow = np.where(np.asarray(surfaces[depths[j+1]]["ids"])==eyed)
            found = index
            foundabove = np.where(np.asarray(surfaces[depths[j-1]]["ids"])==eyed)
            if eyed != -999 and len(foundbelow)!=0 and len(foundbelow[0]) != 0 and len(foundabove)!=0 and len(foundabove[0]) != 0:
                foundbelow = foundbelow[0][0]
                foundabove = foundabove[0][0]
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"s","dsdz",factor=-1)
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"pv","dqdz",factor=-1)
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"alpha","dalphadp")
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"beta","dbetadp")
                surfaces[depths[j]]["data"]["khp"][found] = ptools.calculateKHP(surfaces,depths[j],found)
                surfaces[depths[j]]["data"]["khpterm"][found] = ptools.calculateKHP(surfaces,depths[j],found) * surfaces[depths[j]]["data"]["kapredi"][found] 
                
    for j in Bar('Adding Vertical Double Gradients:   ').iter(range(len(depths))[1:-1]):
        for index in range(len(surfaces[depths[j]]["x"])):
            eyed = int(surfaces[depths[j]]["ids"][index])
            foundbelow = np.where(np.asarray(surfaces[depths[j+1]]["ids"])==eyed)
            found = index
            foundabove = np.where(np.asarray(surfaces[depths[j-1]]["ids"])==eyed)
            if eyed != -999 and len(foundbelow)!=0 and len(foundbelow[0]) != 0 and len(foundabove)!=0 and len(foundabove[0]) != 0:
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"dqdz","d2qdz2",factor=-1)
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"khpterm","khpdz",factor=-1)
    return surfaces

## when we stagger our grid in order to have our gradients 
## "be in the right places" if you will we need to average
## the other quantities like position and t,s,v,h, and pres
def averageOverNeighbors(staggered,surfaces,k,s):
    x = np.mean(surfaces[k]["x"][s])
    y = np.mean(surfaces[k]["y"][s])
    staggered[k]["x"][s[0]] = x
    staggered[k]["y"][s[0]] = y
    staggered[k]["lats"][s[0]] = 90-np.sqrt(x**2 + y**2)/111000.0
    staggered[k]["lons"][s[0]] = np.degrees(np.arctan2(y,x))
    for d in surfaces[k]["data"].keys():
        if d in ["t","s","pv","h","pres","knownu","knownv",\
                "psi","n^2","alpha","beta","toph",\
                "both","kapredi","kapgm","diffkr"]:
            staggered[k]["data"][d][s[0]] = np.mean(surfaces[k]["data"][d][s])
        else:
            staggered[k]["data"][d][s[0]] = staggered[k]["data"][d][s[0]]
    return staggered

##terrible name but add the bathymetric variability coeffient of KVB
def addK(surfaces,cachename=None):
    for k in Bar("adding CKVB: ").iter(surfaces.keys()):
        surfaces[k]["data"]["CKVB"] = np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["CKVO"] = np.full(len(surfaces[k]["lons"]),np.nan)
        for i in range(len(surfaces[k]["lons"])):
            lat = surfaces[k]["lats"][i]
            lon = surfaces[k]["lons"][i]
            if not (np.isnan(lat) or np.isnan(lon)):
                pv = surfaces[k]["data"]["pv"][i]
                pres = surfaces[k]["data"]["pres"][i]
                surfaces[k]["data"]["CKVB"][i] = ptools.Kv(lat,lon,pv,pres,cachename)
    return surfaces
            
#add all the horizontal gradients, and then the double gradients 
def addHorizontalGrad(surfaces,neighbors,distances,debug=False):
    alldxs = []
    staggered = nanCopySurfaces(surfaces)

    for k in Bar('Adding Horizontal Gradients: ').iter(neighbors.keys()):
        for s in neighbors[k]:
            s=np.asarray(s)
            if not np.isnan(s).any():
                staggered = addGradients(staggered,surfaces,k,s,distances)
        for s in neighbors[k]:
            s=np.asarray(s)
            staggered= averageOverNeighbors(staggered,surfaces,k,s)
            staggered = addDoubleGradients(staggered,k,s,distances)
    return staggered

##calculate gradient based on neighboring points
def spatialGrad(surfaces,k,distances,s,attr,factorx=1,factory=1):
    dx = []
    dy = []
    dx.append((surfaces[k]["data"][attr][s[1]]-surfaces[k]["data"][attr][s[0]])/distances[k][(s[0],s[1])])
    dx.append((surfaces[k]["data"][attr][s[3]]-surfaces[k]["data"][attr][s[2]])/distances[k][(s[2],s[3])])
    dy.append((surfaces[k]["data"][attr][s[2]]-surfaces[k]["data"][attr][s[0]])/distances[k][(s[0],s[2])])
    dy.append((surfaces[k]["data"][attr][s[3]]-surfaces[k]["data"][attr][s[1]])/distances[k][(s[1],s[3])])
    dx = np.mean(dx)*factorx
    dy = np.mean(dy)*factory
    return dx,dy

##double derivative of conservative temperature with respect to salinity
def d2thetads2(surfaces,k,s):
    temps = surfaces[k]["data"]["t"][s[0:4]]
    salts = surfaces[k]["data"]["s"][s[0:4]]
    l = np.argsort(salts)
    temps = temps[l]
    salts = salts[l]
    grad1 = np.gradient(np.asarray([temps,salts]))[1][0]
    grad2 = np.mean(np.gradient(np.asarray([grad1,salts]))[1][0])
    return grad2

## set the horizontal gradient 
def setSpatialGrad(out,data,k,s,distances,attr,attrx,attry,factorx=1,factory=1,mode="modify"):
    dx,dy = spatialGrad(data,k,distances,s,attr)
    out[k]["data"][attrx][s[0]] = dx*factorx
    out[k]["data"][attry][s[0]] = dy*factory
    return out

#get the graient of two attributes among neighbors
## so like dt/ds
def attrGrad(out,data,k,s,attry,attrx,outattr):
    xs = data[k]["data"][attrx][s]
    ys = data[k]["data"][attry][s]
    sort = np.argsort(xs)
    xs = xs[sort]
    ys = ys[sort]
    out[k]["data"][outattr][s[0]] = np.mean(np.gradient([ys,xs])[1][0])
    return out

##add all non-vertical gradients
def addGradients(staggered,surfaces,k,s,distances):
    #NS thickness slope
    staggered = setSpatialGrad(staggered,surfaces,k,s,distances,"h","hx","hy")
    staggered = setSpatialGrad(staggered,surfaces,k,s,distances,"psi","v","u",(-1/gsw.f(surfaces[k]["lats"][s[0]])),(1/gsw.f(surfaces[k]["lats"][s[0]])))
    staggered = setSpatialGrad(staggered,surfaces,k,s,distances,"psi","dpsidx","dpsidy")
    staggered = setSpatialGrad(staggered,surfaces,k,s,distances,"s","dsdx","dsdy")
    staggered = setSpatialGrad(staggered,surfaces,k,s,distances,"t","dtdx","dtdy")
    staggered = setSpatialGrad(staggered,surfaces,k,s,distances,"pres","dpdx","dpdy")
    staggered = setSpatialGrad(staggered,surfaces,k,s,distances,"n^2","dqnotdx","dqnotdy",gsw.f(surfaces[k]["lats"][s[0]])/9.8,gsw.f(surfaces[k]["lats"][s[0]])/9.8)
    staggered = setSpatialGrad(staggered,surfaces,k,s,distances,"pv","dqdx","dqdy")
    staggered = attrGrad(staggered,surfaces,k,s,"alpha","t","dalphadtheta")
    staggered = attrGrad(staggered,surfaces,k,s,"alpha","s","dalphads")
    staggered = attrGrad(staggered,surfaces,k,s,"beta","s","dbetads")
    staggered[k]["data"]["d2thetads2"][s[0]] =  d2thetads2(surfaces,k,s)

    return staggered

##add all non-vertical double gradients
def addDoubleGradients(staggered,k,s,distances):
    #NS thickness slope
    d2sdx2,bop = spatialGrad(staggered,k,distances,s,"dsdx")
    bop,d2sdy2 = spatialGrad(staggered,k,distances,s,"dsdy")

    d2qdx2,bop = spatialGrad(staggered,k,distances,s,"dqdx")
    bop,d2qdy2 = spatialGrad(staggered,k,distances,s,"dqdy")

    staggered[k]["data"]["d2sdx2"][s[0]] = d2sdx2
    staggered[k]["data"]["d2sdy2"][s[0]] =  d2sdy2

    staggered[k]["data"]["d2qdx2"][s[0]] = d2qdx2
    staggered[k]["data"]["d2qdy2"][s[0]] =  d2qdy2

    return staggered


## this thing is a little bit wild
## so in the matlab version of gsw they have a function that 
## gives the geostrophic stream function on neutral surfaces
## this does not exist in the python/c version of the gsw.
## so I created a port of it to python, annoying part is
## it requires some other functions also not in the python/c gsw
## but luckily they are in the python only gsw
## so I stole those functions and made a frankenstein that is like
## 100000 times faster than the matlab geostrophic function thing ;)
def addStreamFunc(surfaces,profiles):
    ##re organize surfaces so instead of keys of surface levels
    ##keys of id. Essentially assemble list of ids and pressures where neutral
    ## surfaces are defined
    neutraldepths={}
    for k in surfaces.keys():
        surfaces[k]["data"]["psi"]=np.empty(np.size(surfaces[k]["ids"]))
        surfaces[k]["data"]["psi"][:] =np.nan
        for i in range(len(surfaces[k]["ids"])):
            if surfaces[k]["ids"][i] not in neutraldepths:
                neutraldepths[surfaces[k]["ids"][i]] =[[],[]]
            if abs(k)<5000 and int(abs(surfaces[k]["data"]["pres"][i])) not in neutraldepths[surfaces[k]["ids"][i]][1]:
                neutraldepths[surfaces[k]["ids"][i]][0].append(k)
                neutraldepths[surfaces[k]["ids"][i]][1].append(int(abs(surfaces[k]["data"]["pres"][i])))
    refns = []
    refns_p = []
    for k in neutraldepths.keys():
        if len(neutraldepths[k][0])==len(surfaces.keys()):
            p = getProfileById(profiles,k)
            for j in range(len(neutraldepths[k][0])):
                if ~np.isnan(p.sigma2(neutraldepths[k][1][j])) and neutraldepths[k][0][j] not in refns:
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
    for k in Bar("Adding Stream Func to :").iter(neutraldepths.keys()):
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
            if type(surfaces[nslabels[depth]]["ids"][0]) == type(1):
                targets = np.where(np.asarray(surfaces[nslabels[depth]]["ids"]) ==int(k) )
            elif type(surfaces[nslabels[depth]]["ids"][0]) == type("str"):
                targets = np.where(np.asarray(surfaces[nslabels[depth]]["ids"]) ==str(k) )
            else:
                print("What is the type of the ids dude")

            surfaces[nslabels[depth]]["data"]["psi"][targets] = psi[depth]
        p_ref.append([0]*len(p.ipres))
        ks.append(k)
     
    #results = geo_strf_isopycnal(s,t,ip,p_ref,ns,ns_p,ks)
    
    #with open('data/geoisopycnal.pickle', 'wb') as outfile:
        #pickle.dump([results,ks], outfile)
    #print(results)
    return surfaces

##given fields of psi adds velocity fields by taking the spatial derivative
def streamFuncToUV(surfaces,neighbors,distances):
    for k in Bar('Adding uv: ').iter(neighbors.keys()):
        for s in neighbors[k]:
            s=np.asarray(s)
            if not np.isnan(s).any():
                surfaces = setSpatialGrad(surfaces,surfaces,k,s,distances,"psinew","vabs","uabs",(-1/gsw.f(surfaces[k]["lats"][s[0]])),(1/gsw.f(surfaces[k]["lats"][s[0]])))
                surfaces = setSpatialGrad(surfaces,surfaces,k,s,distances,"psiref","vref","uref",(-1/gsw.f(surfaces[k]["lats"][s[0]])),(1/gsw.f(surfaces[k]["lats"][s[0]])))
                surfaces = setSpatialGrad(surfaces,surfaces,k,s,distances,"psisol","vsol","usol",(-1/gsw.f(surfaces[k]["lats"][s[0]])),(1/gsw.f(surfaces[k]["lats"][s[0]])))
    return surfaces

##Add bathymetry and remove points below
def addBathAndMask(surfaces,neighbors,region):
    surfaces = bathtools.addBathToSurfaces(surfaces,region)
    for k in Bar("Bath Masking").iter(surfaces.keys()):
        for l in range(len(surfaces[k]["lats"])):
            if abs(surfaces[k]["data"]["pres"][l]) > abs(surfaces[k]["data"]["z"][l]):
                for d in surfaces[k]["data"].keys():
                    if d != "ids":
                        surfaces[k]["data"][d][l] = np.nan
                surfaces[k]["x"][l] = np.nan
                surfaces[k]["y"][l] = np.nan
                surfaces[k]["lats"][l] = np.nan
                surfaces[k]["lons"][l] = np.nan
                surfaces[k]["ids"][l] = -999
        final = np.zeros(surfaces[k]["x"].shape)
        for s in neighbors[k]:
            keepem = True
            for l in s:
                if np.isnan(surfaces[k]["data"]["t"][l]):
                    keepem = False
            if keepem:
                for l in s:
                    final[l] = True
        for l in range(len(final)):
            if not final[l]:
                for d in surfaces[k]["data"].keys():
                    if d != "ids":
                        surfaces[k]["data"][d][l] = np.nan
                surfaces[k]["x"][l] = np.nan
                surfaces[k]["y"][l] = np.nan
                surfaces[k]["lats"][l] = np.nan
                surfaces[k]["lons"][l] = np.nan
                surfaces[k]["ids"][l] = -999

    return surfaces
   
def addParametersToSurfaces(surfaces,neighbors,distances):
    surfaces = addHeight(surfaces)
    surfaces = addHorizontalGrad(surfaces,neighbors,distances)
    surfaces = addBathAndMask(surfaces,neighbors,"nepb")
    surfaces = addVerticalGrad(surfaces)
    ptools.saveBathVarTermCache(surfaces,"data/bathVarecco.pickle","nepb")
    surfaces = addK(surfaces,"data/bathVarecco.pickle")
    return surfaces

def artificialPSIRef(surfaces,reflevel = 1700):
    for k in Bar("artifical ref").iter(surfaces.keys()):
        surfaces[k]["data"]["psiref"] = np.full_like(surfaces[k]["data"]["pres"],np.nan)
        for l in range(len(surfaces[k]["ids"])):
            if surfaces[k]["ids"][l] in surfaces[reflevel]["ids"]:
                at = np.where(np.asarray(surfaces[reflevel]["ids"]) == surfaces[k]["ids"][l])[0][0]
                surfaces[k]["data"]["psiref"][l] = surfaces[k]["data"]["psi"][l] - surfaces[reflevel]["data"]["psi"][at]
    return surfaces

def normalizePSI(surfaces,reflevel = 1700):
    for k in Bar("artifical ref").iter(surfaces.keys()):
        psimin = np.nanmin(surfaces[k]["data"]["psi"])
        for l in range(len(surfaces[k]["ids"])):
            surfaces[k]["data"]["psi"][l] = surfaces[k]["data"]["psi"][l] - psimin
    return surfaces


    
##this is spicy. it takes a couple mat files and works out the stream function
## on neutral surfaces from it. Just to check the frankenstein
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

def surfaceSubtract(s1,s2,method="dist"):
    surfaces = {}
    if method == "id":
        for k in s1.keys():
            if k in s2.keys():
                tempSurf = emptySurface()
                for l in range(len(s1[k]["ids"])):
                    if s1[k]["ids"][l] in s2[k]["ids"]:
                        j = np.where(np.asarray(s2[k]["ids"]) == s1[k]["ids"][l])[0][0]
                        tempSurf["lats"].append(s1[k]["lats"][l])
                        tempSurf["lons"].append(s1[k]["lons"][l])
                        tempSurf["ids"].append(s1[k]["ids"][l])
                        for d in s1[k]["data"].keys():
                            if d in s2[k]["data"].keys():
                                if d not in tempSurf["data"]:
                                    tempSurf["data"][d] = []
                                tempSurf["data"][d].append((s1[k]["data"][d][l] - s2[k]["data"][d][j])/s2[k]["data"][d][j])
                surfaces[k] = tempSurf
    if method == "dist":
        s1 = interptools.addXYToSurfaces(s1)
        s2 = interptools.addXYToSurfaces(s2)
        
        for k in s1.keys():
            if k in s2.keys():
                tempSurf = emptySurface()
                for l in range(len(s1[k]["ids"])):
                    j = closestPointSurface(s1[k],s2[k],l)
                    tempSurf["lats"].append(s1[k]["lats"][l])
                    tempSurf["lons"].append(s1[k]["lons"][l])
                    tempSurf["ids"].append(s1[k]["ids"][l])
                    for d in s1[k]["data"].keys():
                        if d in s2[k]["data"].keys():
                            if d not in tempSurf["data"]:
                                tempSurf["data"][d] = []
                            tempSurf["data"][d].append((s1[k]["data"][d][l] - s2[k]["data"][d][j])/s2[k]["data"][d][j])
                surfaces[k] = tempSurf
 
    return surfaces

def closestPointSurface(s1,s2,index):
    dists = (s2["x"]-s1["x"][index])**2+(s2["y"]-s1["y"][index])**2
    return np.argmin(dists)

                    
                



