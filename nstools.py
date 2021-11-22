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
from scipy.interpolate import Rbf, griddata
import pyproj
import bathtools
import graph
import copy
import gsw
import pickle
#from gswmatlab.pyinterface import geo_strf_isopycnal
import scipy.signal
import scipy.io as sio
import itertools
from scipy.linalg import svd
from scipy.linalg import svd
from sklearn.decomposition import TruncatedSVD
from progress.bar import Bar
import parametertools as ptools
from prettytable import PrettyTable
import pdb
import interptools
from scipy import interpolate,integrate
import xarray as xr

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

def timeFilter(profiles,mins,r):
    ps = []
    for p in profiles:
        if abs(p.time-mins)<r:
            ps.append(p)
    return ps

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

##finds the closest profile from a list of profiles in which
# a neutral surface has already been found
def closestIdentifiedNS(profiles,queryprofile,depth,radius,speedy=False):
    if not hasattr(closestIdentifiedNS,"cache"):
        closestIdentifiedNS.cache = {}
    minimumdistance = radius
    minprofile = None
    for p in profiles:
        if depth in p.neutraldepth.keys() and p.presInRange(depth):
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
        #if not hasattr(p,"lon"):
            #pdb.set_trace()
        if lonleft< p.lon <lonright and latbot < p.maplat < lattop\
            and np.max(abs(p.pres))>depth:
            results.append(p)
    return results

##create an empty surface
def emptySurface():
    return {"lats":[],"maplats":[],"lons":[],"ids":[],"cruise":[],\
            "data":{"pres":[],"t":[],"s":[],"pv":[],"n^2":[],"alpha":[],\
            "beta":[],"dalphadp":[],"dalphadtheta":[],"dthetads":[],\
            "psi":[],"drhodz":[],"knownu":[],"knownv":[],"kapredi":[],"kapgm":[],"diffkr":[]}}

##given a seed profile, finds neutral surfaces by comparing each profile 
## to the nearest profile with an identified neutral surface. Only quits
## when all remaining profiles either are not close enough to a profile with a 
## known NS or a NS is not able to be determined for them
def peerSearch(profiles,depth,profilechoice,radius=500,peer=True):
    surfaces = {}
    surfaces[depth]= emptySurface()
    profilechoice.neutraldepth[depth] = depth
    references=[]
    references.append(profilechoice)
    while len(profiles)>0:
        foundcounter = 0
        closestcounter = 0
        for p in Bar("searching...").iter(profiles.copy()):
            closest = closestIdentifiedNS(references,p,depth,radius)
            if closest:
                closestcounter+=1
                #print(closest.neutraldepth)
                ns = closest.neutralDepth(p,closest.neutraldepth[depth],depthname=depth,searchrange=500) 
                #ns = closest.neutralDepthInSituAnom(p,closest.neutraldepth[depth],depthname=depth,searchrange=500)
                if ns != None:
                    dx = closest.x - p.x
                    dy = closest.y - p.y
                    pt,ps,_,_,_,_= p.atPres(ns,interp=True)
                    drhods,drhodt,_ = gsw.rho_first_derivatives(ps,pt,ns)
                    ct,cs,_,_,_,_ = closest.atPres(closest.neutraldepth[depth],interp=True)
                    dsdx,dsdy = (ps-cs)/dx,(ps-cs)/dy
                    dtdx,dtdy = (pt-ct)/dx,(pt-ct)/dy
                    nserror = np.linalg.norm(((drhods*dsdx + drhodt*dtdx),(drhods*dsdy + drhodt*dtdy)))
                    surfaces[depth]["lons"].append(p.lon)
                    surfaces[depth]["lats"].append(p.lat)
                    surfaces[depth]["maplats"].append(p.maplat)

                    surfaces[depth]["data"]["pres"].append(ns)
                    surfaces[depth]["ids"].append(p.eyed)
                    surfaces[depth]["cruise"].append(p.cruise)
                    if p.cruise == "A20":
                        print(p.cruise)
                    profiles.remove(p)
                    foundcounter +=1
                    if peer:
                        references.append(p)
        if foundcounter ==0:
            break

        print("found: ,",foundcounter,"left: ", len(profiles))
        print("close enough: ",closestcounter)
        #plotProfiles(references,"ITS SPREADING",specialprofile = profilechoice)
    return surfaces

def gammaSearch(profiles,depth,gamma,profilechoice,radius=500,peer=True):
    surfaces = {}
    depth = np.asarray(depth)
    for k in depth:
        surfaces[k]= emptySurface()
        surfaces[k]["data"]["gamma"]=[]
    for p in profiles:
        gammas, errl,errh = gamma_n.gamma_n(p.isals,p.itemps,p.ipres,p.lon.copy(),p.lat.copy())
        #plt.plot(gammas,np.asarray(p.ipres))
        #plt.show()
        ssol,tsol,psol = gamma_n.neutralsurfaces(p.isals,p.itemps,p.ipres,gammas,np.asarray(gamma))
        for d in zip(psol,gamma,depth):
            if ~np.isnan(d[0]):
                surfaces[d[2]]["lons"].append(p.lon)
                surfaces[d[2]]["lats"].append(p.lat)
                surfaces[d[2]]["data"]["pres"].append(d[0])
                p.neutraldepth[d[2]]=d[0]
                surfaces[d[2]]["ids"].append(p.eyed)
    return surfaces

## runs the peer search on every neutral surface essentially
def runPeerSearch(profiles,ns,profilechoice,peer=False,radius=10**10,gammas=None):
    surfaces = {}
    if not gammas:
        for d in ns:
            print("NSearching: ",d)
            surfaces.update(peerSearch(profiles.copy(),d,profilechoice,radius,peer))
        return surfaces
    else:
        return gammaSearch(profiles,ns,gammas,profilechoice,radius,peer)

## calculates data throughout neutral surfaces, and some that require calculation 
## shoots that all into a new surfaces object and returns it
def addDataToSurfaces(region,profiles,surfaces,debug=True,noise=0):
    tempSurfs = {}
    for k in Bar("Adding data to: ").iter(surfaces.keys()):
        negativecount = 0
        tempSurf = emptySurface()
        for l in range(len(surfaces[k]["lons"])):
            surfaces[k]["data"]["pres"][l] = abs(surfaces[k]["data"]["pres"][l])
            p = getProfileById(profiles,surfaces[k]["ids"][l])
            lon = surfaces[k]["lons"][l]
            lat = surfaces[k]["lats"][l]
            maplat = surfaces[k]["maplats"][l]
            if p and not np.isnan(p.isals).any() and region.geoFilter(lon,maplat):
                if noise != 0:
                    nspres = int(surfaces[k]["data"]["pres"][l] + np.random.normal(0,noise))
                    p.neutraldepth[k] = nspres
                    if nspres > np.max(p.ipres):
                        nspres = p.ipres[-1]
                    if nspres < np.min(p.ipres[0]):
                        nspres = p.ipres[0]
                else:
                    nspres = int(surfaces[k]["data"]["pres"][l])

                if (nspres+35 in p.ipres and nspres-35 in p.ipres):
                    pv, drhodz = p.potentialVorticityAtHautala(nspres)
                    if pv and pv<0:
                        print(pv)
                    dthetads = p.dthetads(nspres)
                    t,s,alpha,beta,dalphadtheta,dalphadp = p.atPres(abs(surfaces[k]["data"]["pres"][l]),full=True,interp=True)
                    u,v = p.velAtPres(nspres)
                    kapredi,kapgm,diffkr = p.mixAtPres(nspres)
                else:
                    pv = np.nan
                    dsdz = np.nan
                    dthetads = np.nan
                    t,s,alpha,beta,dalphadtheta,dalphadp = p.atPres(abs(surfaces[k]["data"]["pres"][l]),full=True,interp=True)
                    u,v = p.velAtPres(nspres)
                    kapredi,kapgm,diffkr = p.mixAtPres(nspres)
                if (pv and pv<0) or pv==0:
                    #print("second block: ",pv,drhodz)
                    #pv=0
                    negativecount +=1 
                if not isinstance(pv,np.float64):
                    pv = np.nan
                if ~np.isnan(t) and ~np.isnan(s) and ~np.isnan(pv) and p and pv != np.Inf and ~np.isnan(dthetads):
                    tempSurf["lons"].append(lon)
                    tempSurf["lats"].append(lat)
                    tempSurf["maplats"].append(maplat)
                    tempSurf["data"]["pres"].append(abs(nspres))
                    tempSurf["data"]["t"].append(t)
                    tempSurf["data"]["s"].append(s)
                    tempSurf["data"]["pv"].append(pv)
                    tempSurf["data"]["drhodz"].append(drhodz)
                    tempSurf["data"]["dalphadtheta"].append(dalphadtheta)
                    tempSurf["data"]["dalphadp"].append(dalphadp)
                    tempSurf["data"]["dthetads"].append(dthetads)
                    tempSurf["ids"].append(surfaces[k]["ids"][l])
                    tempSurf["data"]["n^2"].append(pv*(9.8/gsw.f(p.lat)))
                    tempSurf["data"]["alpha"].append(alpha)
                    tempSurf["data"]["beta"].append(beta)
                    tempSurf["data"]["knownu"].append(u)
                    tempSurf["data"]["knownv"].append(v)
                    tempSurf["data"]["kapredi"].append(kapredi)
                    tempSurf["data"]["kapgm"].append(kapgm)
                    tempSurf["data"]["diffkr"].append(diffkr)
        if len(tempSurf["lats"])>5:
            tempSurfs[k] = tempSurf
        if len(surfaces[k]["lons"])>0:
            print("ns: ",k," negative count: ",negativecount/len(surfaces[k]["lons"]),"pv mean:" ,np.mean(tempSurf["data"]["pv"]))
    surfaceDiagnostic(tempSurfs)
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

##plot profile t and s
def plotProfile(p):
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.plot(p.itemps,p.ipres)
    ax2.plot(p.isals,p.ipres)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    fig.suptitle(str(p.eyed)+" lat: "+str(p.maplat)+" lon: "+ str(p.lon))
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()
    
#just prints out the percentages of each quantity that are nan along a surface 
def surfaceDiagnostic(surfaces,doprint=True):
    diagnostics ={}
    t = PrettyTable(['Property', 'nan%'])
    for k in surfaces.keys():
        for d in surfaces[k]["data"].keys():
            if d not in diagnostics.keys():
                diagnostics[d]=[0,0]
            diagnostics[d][0] = diagnostics[d][0] + np.count_nonzero(np.isnan(surfaces[k]["data"][d]))
            diagnostics[d][1] = diagnostics[d][1] + np.count_nonzero(~np.isnan(surfaces[k]["data"][d]))
    for d in diagnostics.keys():
        if sum(diagnostics[d])!=0:
            percentage = int(round(diagnostics[d][0]/sum(diagnostics[d]),3)*100)
            t.add_row([d,percentage])
    if doprint:
        print(t)
    return diagnostics

##create a nan copy of a surface
def nanCopySurfaces(surfaces,simple=False):
    nancopy = {}
    for k in surfaces.keys():
        nancopy[k]=emptySurface()
        nancopy[k]["lons"] = surfaces[k]["lons"]
        nancopy[k]["ids"] = surfaces[k]["ids"]
        nancopy[k]["lats"] = surfaces[k]["lats"]
        nancopy[k]["maplats"] = surfaces[k]["maplats"]
        nancopy[k]["x"] = np.full_like(surfaces[k]["lons"],np.nan)
        nancopy[k]["y"] = np.full_like(surfaces[k]["lons"],np.nan)
        datafields = ["u","v","hx","h","CKVB","hy","t","s","pv","pres",\
                     "curl","uabs","vabs","uprime","vprime","dsdx","dsdz","dsdy",\
                    "d2sdx2","d2sdy2","dtdx","dtdy","dpdx","dpdy","n^2",\
                    "dqnotdx","dqnotdy","d2thetads2","dalphadtheta",\
                    "alpha","beta","dalphads","dalphadp","drhodz",\
                      "psi","psinew","dqdz","j","dqdx","dqdy","ddiffkrdz","d2diffkrdz2",\
                    "d2qdz2","d2qdx2","d2qdy2","khp","khpterm","khpdz","toph",\
                    "both","dpsidx","dpsidy","ssh","ids","z","knownu"\
                    ,"knownv","diffkr","kapgm","kapredi",\
                    "nsdiff","dthetads","bathvar","depthmean","kvbdiagnostic"]
        if simple:
            for d in surfaces[k]["data"].keys():
                nancopy[k]["data"][d] = np.full_like(surfaces[k]["data"][d],np.nan,dtype=np.float)
        else:
            for d in datafields:
                nancopy[k]["data"][d] = np.full_like(surfaces[k]["lons"],np.nan,dtype=np.float)
    return nancopy
#perform a vertical gradient of a quantity of a surface
#WITH RESPECT TO PRESSURE NOT Z
# out is a surface you are adding to 
#data is a surface with data
#depths is a dictionary of indexs to depths
## above center and below are the indexs of a point at three different levels
## attr is the quantity being derived, outattr is where you should store the result
def vertGrad(out,data,depths,k,above,center,below,attr,outattr,factor=1):
    return quantVertGrad(out,data,depths,k,above,center,below,attr,"pres",outattr,factor=factor)

def quantVertGrad(out,data,depths,k,above,center,below,attr,quantattr,outattr,factor=1):
    dattr = data[depths[k+1]]["data"][attr][below]-data[depths[k-1]]["data"][attr][above]
    dz = data[depths[k+1]]["data"][quantattr][below]-data[depths[k-1]]["data"][quantattr][above]
    out[depths[k]]["data"][outattr][center] = factor * dattr/dz
    return out



##calculate the height of each neutral surface by the difference in pressures
## of each neutral surface
def addHeight(surfaces):    
    depths = sorted(list(surfaces.keys()))
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
                #surfaces[depths[j]]["data"]["dbetadp"][found] = (surfaces[depths[j-1]]["data"]["beta"][foundabove] - surfaces[depths[j+1]]["data"]["beta"][foundbelow])/((tophalf + bothalf)*2)
    return surfaces
#calculates and stores all the necessary vertical gradients
#then proceeds to calculate all the double vertical gradients
def addVerticalGrad(surfaces): 
    depths = sorted(list(surfaces.keys()))
    for j in Bar('Adding Vertical Gradients:   ').iter(range(len(depths))[1:-1]):
        for index in range(len(surfaces[depths[j]]["x"])):
            eyed = int(surfaces[depths[j]]["ids"][index])
            foundbelow = np.where(np.asarray(surfaces[depths[j+1]]["ids"])==eyed)
            found = index
            foundabove = np.where(np.asarray(surfaces[depths[j-1]]["ids"])==eyed)
            if eyed != -999 and len(foundbelow)!=0 and len(foundbelow[0]) != 0 and len(foundabove)!=0 and len(foundabove[0]) != 0:
                foundbelow = foundbelow[0][0]
                foundabove = foundabove[0][0]
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"pv","dqdz",factor=-1)
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"diffkr","ddiffkrdz",factor=-1)
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"s","dsdz",factor=-1)
                surfaces = quantVertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"dthetads","s","d2thetads2")
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
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"ddiffkrdz","d2diffkrdz2",factor=-1)
                surfaces = vertGrad(surfaces,surfaces,depths,j,foundabove,found,foundbelow,"khpterm","khpdz",factor=-1)
    return surfaces

## when we stagger our grid in order to have our gradients 
## "be in the right places" if you will we need to average
## the other quantities like position and t,s,v,h, and pres
def averageOverNeighbors(staggered,surfaces,k,s):
    x = np.mean(surfaces[k]["x"][s[0]])
    y = np.mean(surfaces[k]["y"][s[0]])
    sign = np.sign(np.sum(surfaces[k]["maplats"][s[0]]))
    staggered[k]["x"][s[0]] = x
    staggered[k]["y"][s[0]] = y
    staggered[k]["lats"][s[0]] = y 
    staggered[k]["maplats"][s[0]] = -y
    staggered[k]["lons"][s[0]] = x
    for d in surfaces[k]["data"].keys():
        if d in ["t","s","pv","h","pres","knownu","knownv",\
                "psi","n^2","alpha","beta","toph","dthetads",\
                "both","kapredi","kapgm","diffkr","dalphadp",\
                "dalphadtheta","drhodz"]:
            #staggered[k]["data"][d][s[0]] = np.mean(surfaces[k]["data"][d][s])
            staggered[k]["data"][d][s[0]] = surfaces[k]["data"][d][s[0]]
        elif d in staggered[k]["data"].keys():
            staggered[k]["data"][d][s[0]] = staggered[k]["data"][d][s[0]]
    return staggered

##terrible name but add the bathymetric variability coeffient of KVB
def addK(surfaces,cachename=None,H_0=1000):
    for k in Bar("adding CKVB: ").iter(surfaces.keys()):
        surfaces[k]["data"]["CKVB"] = np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["bathvar"] = np.full(len(surfaces[k]["lons"]),np.nan)
        surfaces[k]["data"]["depthmean"] = np.full(len(surfaces[k]["lons"]),np.nan)
        for i in range(len(surfaces[k]["lons"])):
            lat = surfaces[k]["maplats"][i]
            lon = surfaces[k]["lons"][i]
            if not (np.isnan(lat) or np.isnan(lon)):
                pv = surfaces[k]["data"]["pv"][i]
                pres = surfaces[k]["data"]["pres"][i]
                CKVB, bathvar,depthmean = ptools.Kv(lat,lon,pv,pres,cachename,H_0=H_0)
                surfaces[k]["data"]["CKVB"][i] = CKVB
                f = gsw.f(surfaces[k]["lats"][i])
                j,_,_ = ptools.jAndDerivatives(np.sqrt(surfaces[k]["data"]["n^2"][i]),f,pv,surfaces[k]["data"]["dqdz"][i])
                surfaces[k]["data"]["j"][i] = j
                surfaces[k]["data"]["bathvar"][i] = bathvar
                surfaces[k]["data"]["depthmean"][i] = depthmean

    #for k in Bar("adding CKVB: ").iter(surfaces.keys()):
        #bathvarmean = np.nanmean(surfaces[k]["data"]["bathvar"])
        #depthmeanmean = np.nanmean(surfaces[k]["data"]["depthmean"])
        #for i in range(len(surfaces[k]["lons"])):
            #surfaces[k]["data"]["kvbdiagnostic"][i] = (depthmeanmean)* (surfaces[k]["data"]["bathvar"][i]-bathvarmean)
            ##surfaces[k]["data"]["kvbdiagnostic"][i] = (depthmeanmean)* (bathvarmean)

    return surfaces
            
#add all the horizontal gradients, and then the double gradients 
def addHorizontalGrad(surfaces,neighbors,distances,ignore,debug=False):
    alldxs = []
    staggered = nanCopySurfaces(surfaces)

    for k in Bar('Adding Horizontal Gradients: ').iter(neighbors.keys()):
        surfaces[k]["data"]["n^2"] = np.full_like(surfaces[k]["data"]["pv"],np.nan)
        for j in range(len(surfaces[k]["lats"])):
            surfaces[k]["data"]["n^2"][j] = surfaces[k]["data"]["pv"][j]*9.8/(gsw.f(surfaces[k]["lats"][j]))
        for s in neighbors[k]:
            s=np.asarray(s)
            if not np.isnan(s).any():
                staggered = addGradients(staggered,surfaces,k,s,distances,ignore)
        for s in neighbors[k]:
            s=np.asarray(s)
            staggered= averageOverNeighbors(staggered,surfaces,k,s)
            staggered = addDoubleGradients(staggered,neighbors,k,s,distances)
    return staggered

def centeredGrad(surfaces,neighbors,k,distances,s,attr,factorx=1,factory=1,single=False):

    if s[1] >= len(neighbors[k]) :
        dx=np.nan
    else:
        xneighbors=False
        for p in neighbors[k]:
            if p[1] == s[0]:
                xneighbors = p
        if xneighbors:
            dx = surfaces[k]["data"][attr][xneighbors[0]]-surfaces[k]["data"][attr][s[1]]
            xdist = distances[k][(s[0],s[1])]+distances[k][(xneighbors[0],xneighbors[1])]
            #print("######")
            #print("point 1: ",surfaces[k]["lats"][xneighbors[0]],surfaces[k]["lons"][xneighbors[0]])
            #print(s[0],s[1])
            #print("point 2: ",surfaces[k]["lats"][s[1]],surfaces[k]["lons"][s[1]])
            #print("dx: ",surfaces[k]["data"][attr][xneighbors[0]]-surfaces[k]["data"][attr][s[1]])
            #print("xdist: ",xdist)
            #pdb.set_trace()
            dx=-dx/xdist
        else:
            dx=np.nan

    if s[2] >= len(neighbors[k]):
        dy=np.nan
    else:
        yneighbors=False
        for p in neighbors[k]:
            if p[2] == s[0]:
                yneighbors = p
        if yneighbors:
            dy = surfaces[k]["data"][attr][s[2]]-surfaces[k]["data"][attr][yneighbors[0]]
            ydist = distances[k][(s[0],s[2])]+distances[k][(yneighbors[0],yneighbors[2])]
            dy = 2*dy/ydist
        else:
            dy=np.nan

    return dx,dy 
    #lets find this center derivative boys

##calculate gradient based on neighboring points
def spatialGrad(surfaces,k,distances,s,attr,factorx=1,factory=1,single=False):
    dx = []
    dy = []

    hautaladist = False
    if hautaladist:
        dx.append((surfaces[k]["data"][attr][s[1]]-surfaces[k]["data"][attr][s[0]])/(distances[k][(s[0],s[1])]))
        if not single:
            dx.append((surfaces[k]["data"][attr][s[3]]-surfaces[k]["data"][attr][s[2]])/(distances[k][(s[0],s[1])]))
    else:
        dx.append((surfaces[k]["data"][attr][s[1]]-surfaces[k]["data"][attr][s[0]])/(distances[k][(s[0],s[1])]))
        if not single:
            dx.append((surfaces[k]["data"][attr][s[3]]-surfaces[k]["data"][attr][s[2]])/(distances[k][(s[2],s[3])]))
    if hautaladist:
        dy.append((surfaces[k]["data"][attr][s[2]]-surfaces[k]["data"][attr][s[0]])/distances[k][(s[0],s[2])])
        if not single:
            dy.append((surfaces[k]["data"][attr][s[3]]-surfaces[k]["data"][attr][s[1]])/distances[k][(s[0],s[2])])
    else:
        dy.append((surfaces[k]["data"][attr][s[2]]-surfaces[k]["data"][attr][s[0]])/distances[k][(s[0],s[2])])
        if not single:
            dy.append((surfaces[k]["data"][attr][s[3]]-surfaces[k]["data"][attr][s[1]])/distances[k][(s[1],s[3])])

    dx = np.nanmean(dx)*factorx
    dy = np.nanmean(dy)*factory
    return dx,dy

def uvError(surfaces,k,distances,s,attr,factorx=1,factory=1,single=False):
    dx = []
    dy = []

    hautaladist = False
    if hautaladist:
        dx.append(np.sqrt(surfaces[k]["data"][attr][s[1]]**2+surfaces[k]["data"][attr][s[0]]**2))
    else:
        dx.append(np.sqrt(surfaces[k]["data"][attr][s[1]]**2+surfaces[k]["data"][attr][s[0]]**2))
    if hautaladist:
        dy.append(np.sqrt(surfaces[k]["data"][attr][s[2]]**2+surfaces[k]["data"][attr][s[0]]**2))
    else:
        dy.append(np.sqrt(surfaces[k]["data"][attr][s[2]]**2+surfaces[k]["data"][attr][s[0]]**2))

    dx = np.nanmean(dx)*factorx
    dy = np.nanmean(dy)*factory
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
def setSpatialGrad(out,data,k,s,distances,attr,attrx,\
        attry,factorx=1,factory=1,mode="modify",centered=False):
    if not centered:
        dx,dy = spatialGrad(data,k,distances,s,attr)
        out[k]["data"][attrx][s[0]] = dx*factorx
        out[k]["data"][attry][s[0]] = dy*factory
    if centered:
        print(centered)
        dx,dy = centeredGrad(data,k,distances,s,attr)
        out[k]["data"][attrx][s[0]] = dx*factorx
        out[k]["data"][attry][s[0]] = dy*factory
    return out

## set the horizontal gradient 
def setUVError(out,data,k,s,distances,attr,attrx,\
        attry,factorx=1,factory=1,mode="modify",centered=False):
    if not centered:
        dx,dy = uvError(data,k,distances,s,attr)
        out[k]["data"][attrx][s[0]] = dx*factorx
        out[k]["data"][attry][s[0]] = dy*factory
    if centered:
        print(centered)
        dx,dy = uvError(data,k,distances,s,attr)
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
def addGradients(staggered,surfaces,k,s,distances,ignore):
    #NS thickness slope
    grads = {"h":["hx","hy",1,1],\
            "psi":["v","u",(1/gsw.f(surfaces[k]["lats"][s[0]])),(-1/gsw.f(surfaces[k]["lats"][s[0]]))],\
            "s":["dsdx","dsdy",1,1],\
            "pres":["dpdx","dpdy",1,1],\
            "n^2":["dqnotdx","dqnotdy",gsw.f(surfaces[k]["lats"][s[0]])/9.8,gsw.f(surfaces[k]["lats"][s[0]])/9.8],\
            "t":["dtdx","dtdy",1,1],\
            "pv":["dqdx","dqdy",1,1]}
    for g in grads.keys():
        if g not in ignore:
            staggered = setSpatialGrad(staggered,surfaces,k,s,distances,g,grads[g][0],grads[g][1],grads[g][2],grads[g][3])
        else:
            staggered[k]["data"][grads[g][0]][s[0]] = surfaces[k]["data"][grads[g][0]][s[0]] 
            staggered[k]["data"][grads[g][1]][s[0]] = surfaces[k]["data"][grads[g][1]][s[0]] 
    #######alpha_wrt_CT_t_exact
    staggered = attrGrad(staggered,surfaces,k,s,"alpha","s","dalphads")

    return staggered

##add all non-vertical double gradients
def addDoubleGradients(staggered,neighbors,k,s,distances):
    #NS thickness slope
    d2sdx2,bop = spatialGrad(staggered,k,distances,s,"dsdx")
    bop,d2sdy2 = spatialGrad(staggered,k,distances,s,"dsdy")

    d2qdx2,bop = spatialGrad(staggered,k,distances,s,"dqdx")
    bop,d2qdy2 = spatialGrad(staggered,k,distances,s,"dqdy")

    staggered[k]["data"]["d2sdx2"][s[0]] = d2sdx2
    staggered[k]["data"]["d2sdy2"][s[0]] = d2sdy2

    staggered[k]["data"]["d2qdx2"][s[0]] = d2qdx2
    staggered[k]["data"]["d2qdy2"][s[0]] = d2qdy2

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
        print(k)
        surfaces[k]["data"]["psi"]=np.empty(np.size(surfaces[k]["ids"]))
        surfaces[k]["data"]["psi"][:] =np.nan
        for i in range(len(surfaces[k]["ids"])):
            if surfaces[k]["ids"][i] not in neutraldepths:
                neutraldepths[surfaces[k]["ids"][i]] =[[],[]]
            if int(abs(surfaces[k]["data"]["pres"][i])) not in neutraldepths[surfaces[k]["ids"][i]][1]:
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
                surfaces = setSpatialGrad(surfaces,surfaces,k,s,distances,"psinew","vabs","uabs",(1/gsw.f(surfaces[k]["maplats"][s[0]])),(-1/gsw.f(surfaces[k]["maplats"][s[0]])))
                surfaces = setSpatialGrad(surfaces,surfaces,k,s,distances,"psiref","vref","uref",(1/gsw.f(surfaces[k]["maplats"][s[0]])),(-1/gsw.f(surfaces[k]["maplats"][s[0]])))
                surfaces = setSpatialGrad(surfaces,surfaces,k,s,distances,"psisol","vsol","usol",(1/gsw.f(surfaces[k]["maplats"][s[0]])),(-1/gsw.f(surfaces[k]["maplats"][s[0]])))
                surfaces = setUVError(surfaces,surfaces,k,s,distances,"psierror","verror","uerror",(surfaces[k]["data"]["vsol"][s[0]]),(surfaces[k]["data"]["usol"][s[0]]))
    return surfaces

def findNeighborSet(neighbors,k,index):
    for l in range(len(neighbors[k])):
        if neighbors[k][l][0] == index:
            return l
    
# Bridge gaps
def neighborCrawl(surfaces,neighbors,distances,k,start_set,dimension,length=15):
    dist=distances[k][(start_set[0],start_set[dimension])]
    new_set = start_set
    depth = 0
    while depth<length:
        if np.isnan(surfaces[k]["x"][new_set[dimension]]):
            old_set = new_set
            new_set_index = findNeighborSet(neighbors,k,new_set[dimension])
            if new_set_index:
                new_set = neighbors[k][new_set_index]
                dist+=distances[k][tuple([old_set[0],new_set[0]])]
            else:
                break
        else:
            return new_set[dimension],dist
        depth+=1
    return start_set[dimension],distances[k][(start_set[0],start_set[dimension])]


##Add bathymetry and remove points below
def addBathAndMask(surfaces,neighbors,distances,region):
    surfaces = bathtools.addBathToSurfaces(surfaces,region)
    for k in Bar("Bath Masking").iter(surfaces.keys()):
        for l in range(len(surfaces[k]["lats"])):
            if -surfaces[k]["data"]["pres"][l] < surfaces[k]["data"]["z"][l]:
                for d in surfaces[k]["data"].keys():
                    if d != "ids":
                        surfaces[k]["data"][d][l] = np.nan
                surfaces[k]["x"][l] = np.nan
                surfaces[k]["y"][l] = np.nan
                #surfaces[k]["lats"][l] = np.nan
                #surfaces[k]["lons"][l] = np.nan
                surfaces[k]["ids"][l] = -999

        tempneighbors = []
        for neighbor_index in range(len(neighbors[k])):
            neighbor_set = neighbors[k][neighbor_index]
            new_set = [neighbor_set[0]]
            for l in range(1,len(neighbor_set)):
                neighbor, dist = neighborCrawl(surfaces,neighbors,distances,k,neighbor_set,l)
                new_set.append(neighbor)
                if not (new_set[0],neighbor) in distances[k]:
                    distances[(new_set[0],neighbor)]=dist

            tempneighbors.append(new_set)
        neighbors[k]=tempneighbors.copy()
        final = np.zeros(surfaces[k]["x"].shape)
        # Iterate through every set of neighbors
        for s in neighbors[k]:
            # if any neighbor is removed set keepem to false
            if not (np.isnan(np.asarray(surfaces[k]["data"]["t"])[np.asarray(s)]).any()):
                for l in s:
                    final[l] = True
        distances[k] = interptools.trueDistanceLookup(surfaces[k],neighbors[k])
        for l in range(len(final)):
            if not final[l]:
                for d in surfaces[k]["data"].keys():
                    if d != "ids":
                        surfaces[k]["data"][d][l] = np.nan
                surfaces[k]["x"][l] = np.nan
                surfaces[k]["y"][l] = np.nan
                surfaces[k]["lats"][l] = np.nan
                surfaces[k]["maplats"][l] = np.nan
                surfaces[k]["lons"][l] = np.nan
                surfaces[k]["ids"][l] = -999
    return surfaces,neighbors,distances
   
def addParametersToSurfaces(region,surfaces,neighbors,distances,ignore=[],H_0=1000):
    #surfaceDiagnostic(surfaces)
    surfaces = addHeight(surfaces)
    #print("after height")
    #surfaceDiagnostic(surfaces)
    surfaces = addHorizontalGrad(surfaces,neighbors,distances,ignore)
    #print("after horizontal grad")
    #surfaceDiagnostic(surfaces)
    surfaces, neighbors, distances = addBathAndMask(surfaces,neighbors,distances,region)
    #print("after bath mask")
    #surfaceDiagnostic(surfaces)
    surfaces = addVerticalGrad(surfaces)
    #print("after vertical grad")
    #surfaceDiagnostic(surfaces)
    ptools.saveBathVarTermCache(surfaces,"data/bathVar.pickle",region)
    #print("after bath var term")
    #surfaceDiagnostic(surfaces)
    surfaces = addK(surfaces,"data/bathVar.pickle",H_0=H_0)
    return surfaces, neighbors, distances

def artificialPSIRef(surfaces,reflevel = 1700):
    for k in Bar("artifical ref").iter(surfaces.keys()):
        surfaces[k]["data"]["psiref"] = np.full_like(surfaces[k]["data"]["pres"],np.nan,dtype=np.double)
        surfaces[reflevel]["ids"] = np.asarray(surfaces[reflevel]["ids"])
        for l in range(len(surfaces[k]["ids"])):
            if surfaces[k]["ids"][l] in surfaces[reflevel]["ids"]:
                at = np.where(surfaces[reflevel]["ids"] == surfaces[k]["ids"][l])[0][0]
                if ~np.isnan(surfaces[reflevel]["data"]["psi"][at]) and ~np.isnan(surfaces[k]["data"]["psi"][l]):
                    surfaces[k]["data"]["psiref"][l] = (surfaces[k]["data"]["psi"][l] - surfaces[reflevel]["data"]["psi"][at])
    for k in Bar("setting ref").iter(surfaces.keys()):
        surfaces[k]["data"]["psi"] = surfaces[k]["data"]["psiref"]
    return surfaces

def normalizePSI(surfaces):
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

def surfaceSubtract(s1,s2,method="dist",metric="%",offset=0):
    surfaces = {}
    if method == "id":
        for k in s1.keys():
            if k in s2.keys():
                tempSurf = emptySurface()
                for l in range(len(s1[k]["ids"])):
                    if s1[k]["ids"][l] in s2[k]["ids"]:
                        j = np.where(np.asarray(s2[k]["ids"]) == s1[k]["ids"][l])[0][0]
                        tempSurf["lats"].append(s1[k]["lats"][l])
                        tempSurf["maplats"].append(s1[k]["maplats"][l])
                        tempSurf["lons"].append(s1[k]["lons"][l])
                        tempSurf["ids"].append(s1[k]["ids"][l])
                        for d in s1[k]["data"].keys():
                            if d in s2[k]["data"].keys() and len(s2[k]["data"][d])>j and len(s1[k]["data"][d])>l:
                                if d not in tempSurf["data"]:
                                    tempSurf["data"][d] = []
                                tempSurf["data"][d].append((s1[k]["data"][d][l] - s2[k]["data"][d][j])/s2[k]["data"][d][j])
                surfaces[k] = tempSurf
    if method == "dist":
        s1 = interptools.addXYToSurfaces(s1)
        s2 = interptools.addXYToSurfaces(s2)
        for k in Bar("Subtracting").iter(s1.keys()):
            if k in s2.keys():
                tempSurf = emptySurface()
                for l in range(len(s1[k]["ids"])):
                    tempSurf["lats"].append(s1[k]["lats"][l])
                    tempSurf["maplats"].append(s1[k]["maplats"][l])
                    tempSurf["lons"].append(s1[k]["lons"][l])
                    tempSurf["ids"].append(s1[k]["ids"][l])
                    for d in s1[k]["data"].keys():
                        if d in s2[k]["data"].keys() and (~np.isnan(s2[k]["data"][d])).any():
                            j = closestPointSurface(s1[k],s2[k],d,l,number=1)
                            if d in s2[k]["data"].keys() and np.max(j) <len(s2[k]["data"][d]) and l <len(s1[k]["data"][d]) :
                                s2value = s2[k]["data"][d][j]
                                if d not in tempSurf["data"]:
                                    tempSurf["data"][d] = []
                                if metric == "-":
                                    tempSurf["data"][d].append((s1[k]["data"][d][l] -offset- s2value ))
                                if metric == "%":
                                    tempSurf["data"][d].append((s1[k]["data"][d][l] -offset- s2value)/s2value)
                                if metric == "/":
                                    tempSurf["data"][d].append((s1[k]["data"][d][l]-offset)/s2value)
                        #elif d not in s2[k]["data"].keys() :
                            #print("missing",d)
                surfaces[k] = tempSurf
    surfaces = interptools.addXYToSurfaces(surfaces)
    return surfaces

def closestPointSurface(s1,s2,d,index,number=1):
    dists = (s2["x"][~np.isnan(s2["data"][d])] - s1["x"][index])**2+(s2["y"][~np.isnan(s2["data"][d])]-s1["y"][index])**2
    #if np.nanmin(dists)>100:
        #print(s1["lons"][index],s1["lats"][index],np.nanmin(dists))
    idx = np.argmin(dists)
    return idx

def inverseReady(surfaces):
    isitnanstr = np.asarray(["alpha","dsdz","hx","hy","dsdx","dsdy","pres","d2sdx2","d2sdy2",\
    "dalphadtheta","dalphads","dalphadp","dtdx","dtdy",\
    "dqnotdx","dqnotdy","dpdx","dpdy","pv","CKVB",\
    "beta","d2qdx2","d2qdy2","d2qdz2","khp","toph","both"])
    diagnostics = surfaceDiagnostic(surfaces)
    fine = True
    missing = []
    for d in isitnanstr:
        if d in diagnostics:
            percentage = int(round(diagnostics[d][0]/sum(diagnostics[d]),3)*100)
            if percentage == 100:
                missing.append(d)
                fine = False
        else:
            missing.append(d)
            fine = False
    if fine:
        print("Hey this is ready!")
    else:
        print("Hey I think you forgot a few things")
        print(missing)

def neutralityError(surfaces):
    for k in surfaces.keys():
        surfaces[k]["data"]["nserror"] = np.full_like(surfaces[k]["data"]["t"],np.nan)
    for k in surfaces.keys():
        for l in range(len(surfaces[k]["lats"])):
            s,t,p = surfaces[k]["data"]["s"][l],surfaces[k]["data"]["t"][l],surfaces[k]["data"]["pres"][l]
            dsdx,dsdy = surfaces[k]["data"]["dsdx"][l],surfaces[k]["data"]["dsdy"][l]
            dtdx,dtdy = surfaces[k]["data"]["dtdx"][l],surfaces[k]["data"]["dtdy"][l]
            ps,pt,pp = gsw.rho_first_derivatives(s,t,p)
            nserror = np.linalg.norm(((ps*dsdx + pt*dtdx),(ps*dsdy + pt*dtdy)))
            surfaces[k]["data"]["nserror"][l] = nserror
    return surfaces


def addOldUnits(surfaces):
    for k in surfaces.keys():
        surfaces[k]["data"]["pottemp"] = gsw.pt_from_CT(surfaces[k]["data"]["s"],surfaces[k]["data"]["t"])
        surfaces[k]["data"]["psu"] = gsw.SP_from_SA(surfaces[k]["data"]["s"],surfaces[k]["data"]["pres"],surfaces[k]["lons"],surfaces[k]["maplats"])
    return surfaces


def meridionalTransport(surfaces,lat,startlon,endlon,startpres,endpres):
    transportsums = {}
    for k in surfaces.keys():
        if  startpres < int(k) < endpres:
            lons = []
            vabs = []
            height = []
            rhos = []
            for l in range(len(surfaces[k]["lats"])):
                if np.abs(surfaces[k]["maplats"][l] - lat)<0.01 and  startlon< surfaces[k]["lons"][l] <endlon:
                    lons.append(surfaces[k]["lons"][l])
                    vabs.append(surfaces[k]["data"]["vabs"][l])
                    height.append(surfaces[k]["data"]["h"][l])
                    rho = gsw.rho(surfaces[k]["data"]["s"][l],surfaces[k]["data"]["t"][l],surfaces[k]["data"]["pres"][l])
                    rhos.append(rho)
            lons=np.asarray(lons)
            vabs=np.asarray(vabs)
            height=np.asarray(height)
            s = np.argsort(lons)
            lons = lons[s]
            vabs = vabs[s]
            height = height[s]
            width = np.nanmin(np.gradient(np.asarray(lons))*np.cos(np.deg2rad(lat))*111*10**3)
            transport = (width*np.asarray(height)*np.asarray(vabs))*rho
            for i in range(len(lons)):
                transportsums.setdefault(lons[i],0)
                if ~np.isnan(transport[i]):
                    transportsums[lons[i]] = transportsums[lons[i]]+transport[i]
    lons = []
    finaltransport = []
    for k in transportsums.keys():
        lons.append(k)
        finaltransport.append(transportsums[k])
    lons = np.asarray(lons)
    finaltransport = np.asarray(finaltransport)
    s = np.argsort(lons)
    totaltransport = np.nansum(finaltransport)
    return totaltransport

def transportDiagnostics(surfaces,vector=["uabs","vabs"]):
    print("vema")
    vema = transportAcross(surfaces,-29.5,-180,-33,3600,100000,"maplats","lons",vector[1])
    print("hunter")
    hunter = transportAcross(surfaces,-29.5,-33,180,3600,100000,"maplats","lons",vector[1])
    print("northern")
    northern = transportAcross(surfaces,-35,-40,-30,3600,100000,"lons","maplats",vector[0])
    print("southern")
    southern = transportAcross(surfaces,-35,-30,-20,3600,100000,"lons","maplats",vector[0])

    two = transportAcross(surfaces,-29.5,-33,180,0,100000,"maplats","lons",vector[1],tempcriteria=2)
    onepointsix = transportAcross(surfaces,-29.5,-33,180,0,100000,"maplats","lons",vector[1],tempcriteria=1.6)
    onepointtwo = transportAcross(surfaces,-29.5,-33,180,0,100000,"maplats","lons",vector[1],tempcriteria=1.2)
    zeropointeight = transportAcross(surfaces,-29.5,-33,180,0,100000,"maplats","lons",vector[1],tempcriteria=0.8)


    curl =0# regionCurl(surfaces,3600,-39,-36,-38,-27)
    print(curl)

    results = {"1.6":onepointsix,\
               "1.2":onepointtwo,\
               "0.8":zeropointeight,\
               "2":two}
    return results


def transportAcross(surfaces,lat,startlon,endlon,startpres,endpres,normalcoord,alongcoord,normalvelocity,tempcriteria=100,unit="volume"):
    transportsums = {}
    surfaces = addOldUnits(surfaces)

    lat = surfaces[list(surfaces.keys())[0]]["maplats"][np.nanargmin(np.abs(surfaces[list(surfaces.keys())[0]]["maplats"] - lat))]
    for k in surfaces.keys():
        if  startpres < int(k) < endpres:
            lons = []
            vabs = []
            height = []
            rhos = []
            for l in range(len(surfaces[k]["lats"])):
                if np.abs(surfaces[k][normalcoord][l] - lat)<0.01 and  startlon< surfaces[k][alongcoord][l] <endlon\
                   and (surfaces[k]["data"]["pottemp"][l]<tempcriteria):
                    lons.append(surfaces[k][alongcoord][l])
                    vabs.append(surfaces[k]["data"][normalvelocity][l])
                    height.append(surfaces[k]["data"]["h"][l])
                    rho = gsw.rho(surfaces[k]["data"]["s"][l],surfaces[k]["data"]["t"][l],surfaces[k]["data"]["pres"][l])
                    rhos.append(rho)
            lons=np.asarray(lons)
            vabs=np.asarray(vabs)
            height=np.asarray(height)
            rhos = np.asarray(rhos)
            s = np.argsort(lons)
            lons = lons[s]
            vabs = vabs[s]
            height = height[s]
            rhos = rhos[s]
            if len(lons)>1:
                width = np.nanmin(np.gradient(np.asarray(lons))*np.cos(np.deg2rad(lat))*111*10**3)
                if unit=="mass":
                    transport = (width*np.asarray(height)*np.asarray(vabs))*rhos
                if unit == "volume":
                    transport = (width*np.asarray(height)*np.asarray(vabs))
                for i in range(len(lons)):
                    transportsums.setdefault(lons[i],0)
                    if ~np.isnan(transport[i]):
                        transportsums[lons[i]] = transportsums[lons[i]]+transport[i]
    lons = []
    finaltransport = []
    for k in transportsums.keys():
        lons.append(k)
        finaltransport.append(transportsums[k])
    lons = np.asarray(lons)
    finaltransport = np.asarray(finaltransport)
    s = np.argsort(lons)
    totaltransport = np.nansum(finaltransport)
    return totaltransport


def externalReference(surfaces,referencefile):
    ref_vs = xr.open_dataset(referencefile)
    levels = np.asarray((sorted(surfaces.keys())))
    below = levels[levels>1000][0]
    above = levels[levels<1000][-1]
    print(below,above)
    ref_1000 = {}
    for below_i in range(len(surfaces[below]["ids"])):
        eyed = surfaces[below]["ids"][below_i]
        lon = surfaces[below]["lons"][below_i]
        lat = surfaces[below]["maplats"][below_i]
        if eyed in surfaces[above]["ids"] and ~np.isnan(lon):
            above_i = surfaces[above]["ids"].index(eyed)
            u_c = -np.interp(1000,[below,above],[surfaces[below]["data"]["u"][below_i],surfaces[below]["data"]["u"][above_i]])
            v_c = -np.interp(1000,[below,above],[surfaces[below]["data"]["v"][below_i],surfaces[below]["data"]["v"][above_i]])
            u0 = float(ref_vs.interp(LONGITUDE=365+lon,LATITUDE=lat).U)
            v0 = float(ref_vs.interp(LONGITUDE=365+lon,LATITUDE=lat).V)
            ref_1000[eyed] = [u_c-u0,v_c-v0]


    for k in Bar("referencing to yomaha from model uv").iter(list(surfaces.keys())):
        surfaces[k]["data"]["yomahau"] =np.full_like(surfaces[k]["ids"],np.nan,dtype=np.double)
        surfaces[k]["data"]["yomahav"] =np.full_like(surfaces[k]["ids"],np.nan,dtype=np.double)
        for l in range(len(surfaces[k]["ids"])):
            eyed = surfaces[k]["ids"][l]
            if eyed in ref_1000.keys():
                surfaces[k]["data"]["yomahau"][l]= -surfaces[k]["data"]["u"][l]-ref_1000[eyed][0]
                surfaces[k]["data"]["yomahav"][l]= -surfaces[k]["data"]["v"][l]-ref_1000[eyed][1]
    return surfaces


def twoCReference(surfaces):
    surfaces = addOldUnits(surfaces)
    levels = np.asarray((sorted(surfaces.keys())))
    below = levels[int(len(levels)/2)]
    ref_2C = {}
    for eyed in surfaces[below]["ids"]:
        for k_i in range(1,len(levels)):
            k = levels[k_i]
            k_prev = levels[k_i-1]
            lon = surfaces[k]["lons"][k_i]
            if eyed in surfaces[k]["ids"] and ~np.isnan(lon) and eyed in surfaces[k_prev]["ids"]:
                cold_i = surfaces[k]["ids"].index(eyed)
                warm_i = surfaces[k_prev]["ids"].index(eyed)
                if surfaces[k]["data"]["pottemp"][cold_i] < 2:
                    warmertemp = surfaces[k_prev]["data"]["pottemp"][warm_i]
                    coldertemp = surfaces[k]["data"]["pottemp"][cold_i]
                    u_c = -np.interp(2,[coldertemp,warmertemp],[surfaces[k]["data"]["u"][cold_i],surfaces[k_prev]["data"]["u"][warm_i]])
                    v_c = -np.interp(2,[coldertemp,warmertemp],[surfaces[k]["data"]["v"][cold_i],surfaces[k_prev]["data"]["v"][warm_i]])
                    ref_2C[eyed] = [u_c,v_c]
                    break


    for k in Bar("referencing to 2C from model uv").iter(list(surfaces.keys())):
        surfaces[k]["data"]["2CU"] =np.full_like(surfaces[k]["ids"],np.nan,dtype=np.double)
        surfaces[k]["data"]["2CV"] =np.full_like(surfaces[k]["ids"],np.nan,dtype=np.double)
        for l in range(len(surfaces[k]["ids"])):
            eyed = surfaces[k]["ids"][l]
            if eyed in ref_2C.keys():
                surfaces[k]["data"]["2CU"][l]= -surfaces[k]["data"]["u"][l]-ref_2C[eyed][0]
                surfaces[k]["data"]["2CV"][l]= -surfaces[k]["data"]["v"][l]-ref_2C[eyed][1]
    return surfaces


def twoCProjection(surfaces,quants):
    surfaces = addOldUnits(surfaces)
    levels = np.asarray((sorted(surfaces.keys())))
    surf_2C = surfaces[levels[0]].copy()
    for d in surf_2C["data"].keys():
        for l in range(len(surf_2C["data"]["t"])):
            surf_2C["data"][d][l]=np.nan
    for eyed_i in range(len(surfaces[levels[0]]["ids"])):
        eyed = surfaces[levels[0]]["ids"][eyed_i]
        for k_i in range(1,len(levels)):
            k = levels[k_i]
            k_prev = levels[k_i-1]
            lon = surfaces[k]["lons"][k_i]
            if eyed in surfaces[k]["ids"] and ~np.isnan(lon) and eyed in surfaces[k_prev]["ids"]:
                cold_i = surfaces[k]["ids"].index(eyed)
                warm_i = surfaces[k_prev]["ids"].index(eyed)
                if surfaces[k]["data"]["pottemp"][cold_i] < 2:
                    warmertemp = surfaces[k_prev]["data"]["pottemp"][warm_i]
                    coldertemp = surfaces[k]["data"]["pottemp"][cold_i]
                    for q in quants:
                        surf_2C["data"][q][eyed_i] = -np.interp(2,[coldertemp,warmertemp],[surfaces[k]["data"][q][cold_i],surfaces[k_prev]["data"][q][warm_i]])
    return {2:surf_2C}


### morris estimate [2.45, 1.59, 0.47]
### SH/GF estimate [-0.51, -0.31, -0.05]
def morrisMixing(hunter_transport=[-0.51, -0.31, -0.05]):
    isotherms = np.array([1.6,1.2,0.8])
    surfacearea = np.array([8*(10**12),7.1*(10**12),5.6*(10**12)])
    thetaz = np.array([2.2*(10**-3),2.1*(10**-3),1.7*(10**-3)])
    vema_transport = np.array([3.99,4.14,4.02])
    vema_temp = np.array([0.06,0.08,0.03])
    hunter_temp=np.array([1.06,0.90,0.60])
    romanche_transport = np.array([-0.77,-0.44,-0.03])
    romanche_temp = np.array([1.16,0.97,0.78])
    equatorial_transport = np.array([-1.89,-1.55,-0.76])
    equatorial_temp = np.array([0.92,0.83,0.7])
    # twt = temperature weighted transport
    twt = vema_transport * (vema_temp-isotherms)
    twt += hunter_transport * (hunter_temp-isotherms)
    twt += romanche_transport * (romanche_temp-isotherms)
    twt += equatorial_transport * (equatorial_temp-isotherms)
    mixing = (twt*(10**6))/np.multiply(surfacearea,-thetaz)
    print(mixing)
    return(mixing)

def zenk1993Blend(surfaces):
    ref = graph.transportRefIsotherm(inv,2,-29.5,-30,-9,1000,6000,ax=None)
    transports = np.asarray([0,0.15,0.5])-ref
    print(transports)
    print(morrisMixing)

def zenk1999Blend(transports):
    pointeight = 0.47
    morrisMixing([transports[0]+pointeight,transports[1]+pointeight,pointeight])





