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
from prettytable import PrettyTable


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
            dist = great_circle((queryprofile.lat,queryprofile.lon),(p.lat,p.lon)).km
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
    return {"lats":[],"lons":[],"ids":[],\
        "data":{"pres":[],"t":[],"s":[],"pv":[],"n^2":[],"alpha":[],"beta":[]}}

def peerSearch(profiles,deepestindex,depth,profilechoice,radius=500):
    surfaces = {}
    surfaces[depth]= emptySurface()
    profilechoice.neutraldepth[depth] = depth
    references=[]
    references.append(profilechoice)
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


def runPeerSearch(profiles,deepestindex,shallowlimit,deeplimit,surfacestep,profilechoice,radius):
    surfaces = {}
    for d in range(shallowlimit,deeplimit,surfacestep)[::-1]:
        print("NSearching: ",d)
        surfaces.update(peerSearch(profiles.copy(),deepestindex,d,profilechoice,1000))
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


def addDataToSurfaces(profiles,surfaces,stdevs,debug=True):
    tempSurfs = {}
    for k in Bar("Adding data to: ").iter(surfaces.keys()):
        negativecount = 0
        tempSurf = emptySurface()
        monthcount = [0]*13
        monthnegativecount = [0]*13
        for l in range(len(surfaces[k]["lons"])):
            p = getProfileById(profiles,surfaces[k]["ids"][l])
            above,below = findAboveIndex(surfaces,k,l)
            if p and not np.isnan(p.isals).any() :
                if (above or below):
                    pv = p.potentialVorticityBetween(above,below)
                    t,s = p.betweenPres(above,below)
                else:
                    pv = p.potentialVorticityAt(surfaces[k]["data"]["pres"][l])
                    t,s = p.atPres(surfaces[k]["data"]["pres"][l])
                monthcount[p.time.month]=monthcount[p.time.month]+1
                if pv and pv<0:
                    monthnegativecount[p.time.month]=monthnegativecount[p.time.month]+1
                    negativecount +=1 
                if not pv:
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

        print("ns: ",k," negative count: ",negativecount/len(surfaces[k]["lons"]),"pv mean:" ,np.mean(tempSurf["data"]["pv"]))
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
    x=np.asarray(surface["x"])
    y=np.asarray(surface["y"])
    z=np.asarray(surface["data"]["pres"])
    final = np.zeros(x.shape)
    #print(x)
    for i in Bar("Removing discontinuities: ").iter(range(len(x))):
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
                if final[l][k] and final[l][k+1] and final[l+1][k+1] and final[l+1][k]:
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
    for k in Bar("Interpolating: ").iter(surface["data"].keys()):
        notnan = ~np.isnan(surface["data"][k])
        if np.count_nonzero(notnan)>10:
            gam = pygam.LinearGAM(pygam.te(0,1)).fit(X[notnan],surface["data"][k][notnan])
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

def interpolateSurfaces(surfaces,debug=True):
    interpolatedsurfaces = {}
    neighbors={}
    lookups={}
    for k in surfaces.keys():
        if ~np.isnan(surfaces[k]["data"]["pres"]).any():
            surfaces[k] = removeDiscontinuities(surfaces[k],radius=0.1)
            interpolatedsurfaces[k],neighbors[k] = interpolateSurface(surfaces[k])
            lookups[k] = trueDistanceLookup(interpolatedsurfaces[k],neighbors[k])
    return interpolatedsurfaces,neighbors,lookups
#

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
        percentage = int(round(diagnostics[d][0]/sum(diagnostics[d]),3)*100)
        t.add_row([d,percentage])
    print(t)


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
                    "d2qdz2","d2qdx2","d2qdy2","khp","khpdz"]
        for d in datafields:
            nancopy[k]["data"][d] = np.full(len(surfaces[k]["lons"]),np.nan)
    return nancopy

def fillOutEmptyFields(surfaces):
    for k in surfaces.keys():
        datafields = ["u","v","hx","h","CKVB","hy","t","s","pv","pres",\
                     "curl","uabs","vabs","uprime","vprime","dsdx","dsdz","dsdy",\
                    "d2sdx2","d2sdy2","dtdx","dtdy","dpdx","dpdy","n^2",\
                    "dqnotdx","dqnotdy","d2thetads2","dalphadtheta",\
                    "alpha","beta","dalphads","dbetads","dalphadp",\
                    "dbetadp","psi","dqdz","dqdx","dqdy",\
                    "d2qdz2","d2qdx2","d2qdy2","khp","khpdz","dpsidx","dpsidy"]
        for d in datafields:
            if d not in surfaces[k]["data"].keys():
                surfaces[k]["data"][d] = np.full(len(surfaces[k]["lons"]),np.nan)
    return surfaces


def vertGrad(out,data,depths,k,above,center,below,attr,outattr,factor=1):
    dattr = data[depths[k-1]]["data"][attr][above]-data[depths[k+1]]["data"][attr][below]
    dz = data[depths[k-1]]["data"]["pres"][above]-data[depths[k+1]]["data"]["pres"][below]
    out[depths[k]]["data"][outattr][center] = factor * dattr/dz
    return out

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
                surfaces[depths[j]]["data"]["dbetadp"][found] = (surfaces[depths[j-1]]["data"]["beta"][foundabove] - surfaces[depths[j+1]]["data"]["beta"][foundbelow])/((tophalf + bothalf)*2)
    return surfaces

def calculateKHP(staggered,k,index):
    dalphadtheta = staggered[k]["data"]["dalphadtheta"][index]
    dalphadp = staggered[k]["data"]["dalphadp"][index]
    dalphads = staggered[k]["data"]["dalphads"][index]
    dbetads = staggered[k]["data"]["dbetads"][index]
    alpha = staggered[k]["data"]["alpha"][index]
    beta = staggered[k]["data"]["beta"][index]
    dbetadp = staggered[k]["data"]["dbetadp"][index]
    betaTherm = staggered[k]["data"]["dbetadp"][index]
    alphat = dalphadtheta+2*(alpha/betaTherm)*dalphads-(alpha**2/betaTherm**2)*dbetads
    alphap = dalphadp -(alpha/betaTherm)*dbetadp
    magct = staggered[k]["data"]["dtdx"][index]**2 + staggered[k]["data"]["dtdx"][index]**2
    cdotp = staggered[k]["data"]["dtdx"][index]*staggered[k]["data"]["dpdx"][index]+staggered[k]["data"]["dtdy"][index]*staggered[k]["data"]["dpdy"][index]
    return alphat*magct+alphap*cdotp


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
            if len(foundbelow)!=0 and len(foundbelow[0]) != 0 and len(foundabove)!=0 and len(foundabove[0]) != 0:
                foundbelow = foundbelow[0][0]
                foundabove = foundabove[0][0]
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"s","dsdz",factor=-1)
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"pv","dqdz",factor=-1)
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"alpha","dalphadp")
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"beta","dbetadp")
                surfaces[depths[j]]["data"]["khp"][found] = calculateKHP(surfaces,depths[j],found)
                
    for j in Bar('Adding Vertical Double Gradients:   ').iter(range(len(depths))[1:-1]):
        for index in range(len(surfaces[depths[j]]["x"])):
            eyed = int(surfaces[depths[j]]["ids"][index])
            foundbelow = np.where(np.asarray(surfaces[depths[j+1]]["ids"])==eyed)
            found = index
            foundabove = np.where(np.asarray(surfaces[depths[j-1]]["ids"])==eyed)
            if len(foundbelow)!=0 and len(foundbelow[0]) != 0 and len(foundabove)!=0 and len(foundabove[0]) != 0:
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"dqdz","d2qdz2",factor=-1)
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"khp","khpdz",factor=-1)
    return surfaces

def averageOverNeighbors(staggered,surfaces,k,s):
    staggered[k]["lons"][s[0]] = np.mean(np.abs(surfaces[k]["lons"][s]))*np.sign(surfaces[k]["lons"][s[0]])
    staggered[k]["lats"][s[0]] = np.mean(surfaces[k]["lats"][s])
    staggered[k]["x"][s[0]] = np.mean(surfaces[k]["x"][s])
    staggered[k]["y"][s[0]] = np.mean(surfaces[k]["y"][s])
    for d in surfaces[k]["data"].keys():
        if d != "ids":
            staggered[k]["data"][d][s[0]] = np.mean(surfaces[k]["data"][d][s])
    return staggered

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
            
    
def addHorizontalGrad(surfaces,neighbors,distances,debug=False):
    alldxs = []
    staggered = nanCopySurfaces(surfaces)

    for k in Bar('Adding Horizontal Gradients: ').iter(neighbors.keys()):
        for s in neighbors[k]:
            s=np.asarray(s)
            if not np.isnan(s).any():
                staggered = addGradients(surfaces,surfaces,k,s,distances)
        for s in neighbors[k]:
            s=np.asarray(s)
            if not np.isnan(s).any():
                staggered= averageOverNeighbors(staggered,surfaces,k,s)
                staggered = addDoubleGradients(surfaces,k,s,distances)

    return staggered

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

def d2thetads2(surfaces,k,s):
    temps = surfaces[k]["data"]["t"][s[0:3]]
    salts = surfaces[k]["data"]["s"][s[0:3]]
    l = np.argsort(salts)
    temps = temps[l]
    salts = salts[l]
    d1 = (temps[1] - temps[0])/(salts[1] - salts[0])
    d2 = (temps[2] - temps[1])/(salts[2] - salts[1])
    return (d2-d1)/(salts[2] - salts[0])

def setSpatialGrad(out,data,k,s,distances,attr,attrx,attry,factorx=1,factory=1,mode="modify"):
    dx,dy = spatialGrad(data,k,distances,s,attr)
    out[k]["data"][attrx][s[0]] = np.mean(dx)*factorx
    out[k]["data"][attry][s[0]] = np.mean(dy)*factory
    return out

def attrGrad(out,data,k,s,attry,attrx,outattr):
    xs = data[k]["data"][attrx][s]
    ys = data[k]["data"][attry][s]
    sort = np.argsort(xs)
    xs = xs[sort]
    ys = ys[sort]
    grad = []
    for i in range(1,len(ys)-1):
        grad.append((ys[i+1]-ys[i-1])/(xs[i+1]-xs[i-1]))
    out[k]["data"][outattr][s[0]] = np.mean(grad)
    return out

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


def trueDistanceLookup(surface,neighbors):
    lookup = {}
    for square in Bar("distance calc: ").iter(neighbors):
        for edge in itertools.combinations(square,2):
            p = tuple(sorted(edge))
            if p not in lookup.keys():
                lookup[p] = geodesic((surface["lats"][p[0]],surface["lons"][p[0]]),(surface["lats"][p[1]],surface["lons"][p[1]])).m
    return lookup

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
            if abs(k)<3700 and int(abs(surfaces[k]["data"]["pres"][i])) not in neutraldepths[surfaces[k]["ids"][i]][1]:
                neutraldepths[surfaces[k]["ids"][i]][0].append(k)
                neutraldepths[surfaces[k]["ids"][i]][1].append(int(abs(surfaces[k]["data"]["pres"][i])))
    refns = []
    refns_p = []
    for k in neutraldepths.keys():
        if len(neutraldepths[k][0])==18:
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

def streamFuncToUV(surfaces,neighbors,distances):
    for k in Bar('Adding uv: ').iter(neighbors.keys()):
        for s in neighbors[k]:
            s=np.asarray(s)
            if not np.isnan(s).any():
                surfaces = setSpatialGrad(surfaces,surfaces,k,s,distances,"psinew","vabs","uabs",(-1/gsw.f(surfaces[k]["lats"][s[0]])),(1/gsw.f(surfaces[k]["lats"][s[0]])))
                surfaces = setSpatialGrad(surfaces,surfaces,k,s,distances,"psiref","v","u",(-1/gsw.f(surfaces[k]["lats"][s[0]])),(1/gsw.f(surfaces[k]["lats"][s[0]])))
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

