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
from gswmatlab.pyinterface import geo_strf_isopycnal

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

def addDataToSurfaces(profiles,surfaces,stdevs,debug=True):
    if debug:
        print("adding data to surfaces")
    tempSurfs = {}
    negativecount = 0
    for k in surfaces.keys():
        tempSurf = np.array([[],[],[[],[],[],[]],[]])
        for l in range(len(surfaces[k][0])):
            p = getProfileById(profiles,surfaces[k][3][l])
            p.interpolate()
            t,s = p.atPres(surfaces[k][2][l])
            pv = p.potentialVorticity(surfaces[k][2][l],debug=False)
            if pv and pv<0:
                negativecount +=1 
            if t and s and pv and p and pv != np.Inf and pv != np.nan and not np.isnan(t):
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
            #np.asarray(tempSurf[2][j])[a[0]] == np.nan
            

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
        distFilter = np.zeros(len(surfaces[k][0]), dtype=bool)
        for index in range(len(surfaces[k][0])):
            p = (surfaces[k][1][index],surfaces[k][0][index])
            #print(lpoint,p)
            #if geodesic(lpoint,p).km<radius:
                #distFilter[index] = True
            if ((np.cos(p[0]*(np.pi/180)))*abs(((p[1]+180)%180)-lon)*111)<radius:
                #print(((np.cos(lpoint[0]*(np.pi/180)))*abs(p[1]-lpoint[1])*111))
                distFilter[index] = True
        surfaces[k][0] = np.asarray(surfaces[k][0])[distFilter]
        surfaces[k][1] = np.asarray(surfaces[k][1])[distFilter]
        surfaces[k][2] = [np.asarray(surfaces[k][2][0])[distFilter],np.asarray(surfaces[k][2][1])[distFilter],np.asarray(surfaces[k][2][2])[distFilter],np.asarray(surfaces[k][2][3])[distFilter]]
        if len(surfaces[k][3])>0:
            surfaces[k][3] = np.asarray(surfaces[k][3])[distFilter]
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


def surfacesToXYZPolar(surfaces,debug=True):
    if debug:
        print("converting surfaces to xyz")
    newsurfaces = {}
    for k in surfaces.keys():
        tempSurf = np.array([[],[],[[],[],[],[]],[]])
        x,y = homemadeXY(surfaces[k][0],surfaces[k][1])
        tempSurf[0] = x
        tempSurf[1] = y
        tempSurf[2][0] = surfaces[k][2][0]
        tempSurf[2][1] = surfaces[k][2][1]
        tempSurf[2][2] = surfaces[k][2][2]
        tempSurf[2][3] = surfaces[k][2][3]
        tempSurf[3] = surfaces[k][3]
        newsurfaces[k]=tempSurf
    return newsurfaces



def surfacesToXYZ(surfaces,debug=True):
    if debug:
        print("converting surfaces to xyz")
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

def removeDiscontinuities(x,y,z,auxdata=[],radius=10,debug=True):
    if debug:
        print("removing discontinuities")
    x=np.asarray(x)
    y=np.asarray(y)
    z=np.asarray(z)
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
    if auxdata:
        newaux=[]
        for a in auxdata:
            newaux.append(np.asarray(a)[np.invert(final)])
        return x[np.invert(final)],y[np.invert(final)],z[np.invert(final)],newaux
    else:
        return x[np.invert(final)],y[np.invert(final)],z[np.invert(final)]


def createMesh(n,xvals,yvals):
    return np.meshgrid(np.linspace(np.min(xvals),np.max(xvals),n), np.linspace(np.min(yvals),np.max(yvals),n),indexing="xy")

def generateMaskedMesh(x,y,radius=100):
    xi,yi = createMesh(100,x,y)
    final = np.zeros(xi.shape)
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
    return xi[final],yi[final]

def removeOutlierSurfaces(surfaces,stdevs=2):
    for k in surfaces.keys():
        m = np.median(surfaces[k][2][0])
        s = np.std(surfaces[k][2][0])
        filt = np.where(np.abs(surfaces[k][2][0]-m)<s*stdevs)
        surfaces[k][0]=surfaces[k][0][filt]
        surfaces[k][1]=surfaces[k][1][filt]
        surfaces[k][2][0] =surfaces[k][2][0][filt]
        surfaces[k][2][1] =surfaces[k][2][1][filt]
        surfaces[k][2][2] =surfaces[k][2][2][filt]
        surfaces[k][2][3] =surfaces[k][2][3][filt]
        if len(surfaces[k][3]) >0:
            surfaces[k][3] =surfaces[k][3][filt]
    return surfaces



def interpolateSurface(x,y,z,d=None,debug=True):
    if debug:
        print("interpolating surfaces")
    if d:
        di = []
        rbfi = Rbf(x,y,z,function="thin_plate",smooth=20.0)
        xi,yi = generateMaskedMesh(x,y)
        zi = rbfi(xi,yi)
        for q in d:
            qrbfi = Rbf(x,y,z,q,function="thin_plate",smooth=20.0)
            di.append(qrbfi(xi,yi,zi))
        return xi,yi,zi,di
    else:
        rbfi = Rbf(x,y,z,function="thin_plate",smooth=10.0)
        xi,yi = generateMaskedMesh(x,y)
        zi = rbfi(xi,yi)
        return xi,yi,zi

def interpolateSurfaceGAM(x,y,z,d=None,debug=True):
    #print("######")
    X = np.zeros((len(x),2))
    X[:,0]=x
    X[:,1]=y
    xi,yi = generateMaskedMesh(x,y)
    if d:
        di=[]
        for q in [z] + d:
            gam = pygam.LinearGAM(pygam.te(0,1)).fit(X,q)
            Xgrid = np.zeros((yi.shape[0],2))
            Xgrid[:,0] = xi
            Xgrid[:,1] = yi
            di.append(gam.predict(Xgrid))
            #gam.summary()
        return xi,yi,di[0],di[1:]
    else:
        gam = pygam.LinearGAM(pygam.te(0,1)).fit(X,z)
        #gam.summary()
        #print(xi.ravel())
        Xgrid = np.zeros((yi.shape[0],2))
        Xgrid[:,0] = xi
        Xgrid[:,1] = yi
        zi = gam.predict(Xgrid)
        return xi,yi,zi

def homemadeXY(lon,lat):
    x=[]
    y=[]
    for i in range(len(lat)):
        theta = np.deg2rad(lon[i])
        r = ((90-lat[i]) *111*1000)
        x.append(r*np.cos(theta))
        y.append(r*np.sin(theta))
    return np.asarray(x),np.asarray(y)


def xyzToSurface(x,y,z,d,depth,debug = True):
    if debug:
        print("converting xyz back to lat lon a")
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum = 'WGS84',preserve_units=True)
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum = 'WGS84',preserve_units=True)
    surface = np.array([[],[],[[],[],[],[]],[]])
    lon, lat, a = pyproj.transform(ecef,lla,x,y,z,radians=False)
    surface[0]=lon
    surface[1]=lat
    surface[2]=[a]+d
    return {depth:surface}

def xyzToSurfacePolar(x,y,z,d,depth,debug = True):
    if debug:
        print("converting xyz back to lat lon a")
    surface = np.array([[],[],[[],[],[],[]],[]])
    lat = 90-(np.sqrt((x**2+y**2))/111000.0)
    lon = np.degrees(np.arctan2(y,x))
    #print("lat: ",lat, " lon: ",lon)
    surface[0]=lon
    surface[1]=lat
    surface[2]=[z]+d
    return {depth:surface}
     
        
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

def plotASpiral(profiles,center=None,x=None,y=None):
    fig,ax = plt.subplots(1,1) 
    center = getProfileById(profiles,"120389")
    x = getProfileById(profiles,"120387")
    y = getProfileById(profiles,"114688")
    dx = geodesic((x.lat,x.lon),(center.lat,center.lon)).m
    dy = geodesic((x.lat,x.lon),(center.lat,center.lon)).m
    u = 0
    v = 0
    f = center.f
    us = [0]
    vs = [0]
    for i in range(100,1700,100):
        dpdx = (x.densityAtPres(i)-center.densityAtPres(i))/dx
        dpdy = (y.densityAtPres(i)-center.densityAtPres(i))/dy
        v = (9.8/-f)*(dpdx)
        u = (9.8/f)*(dpdy)
        us.append(u)
        vs.append(v)
        plt.plot([0,u],[0,v])
        ax.annotate(i, (u, v))
    #plt.scatter(us,vs)
    #plt.plot(us,vs,c="r")
    plt.show()

def generateNeighborsList(x,y):
    print("calculating neighbors")
    ref = {}
    for center in range(len(x)):
        dist = ((x-x[center])**2 + (y-y[center])**2)
        ref[center]=np.argsort(dist)[1:5]
    return ref
                     
def componentDistance(surfaces,k,i1,i2):
    #x = surfaces[k][0][i1] - surfaces[k][0][i2]
    #if abs(x) > abs(360 - (x+360)%360):
        #x = np.sign(x)*(360-(x+360)%360)
    if  (surfaces[k][0][i1]+180 ) > (surfaces[k][0][i2]+180):
        x = surfaces[k][0][i1]+180 - (surfaces[k][0][i2]+180)
        if x>180:
            x = -(360-x)
    else:
        x = surfaces[k][0][i1]+180 - (surfaces[k][0][i2]+180)
        if x < -180:
            x= -(-360-x)
    x=x*np.cos(np.deg2rad(surfaces[k][1][i2]))*111000.0
    #print(surfaces[k][0][i1],surfaces[k][0][i2],x)

    y = (surfaces[k][1][i1]-surfaces[k][1][i2])*111000.0
    return x,y


def addPrimeToSurfaces(surfaces,neighbors,debug=False):
    for k in surfaces.keys():
        surfaces[k][2].append(np.zeros(len(surfaces[k][1])))
        surfaces[k][2].append(np.zeros(len(surfaces[k][1])))
    alldxs = []
    for k in neighbors.keys():
        print("adding primes to: ",k)
        for r in neighbors[k].keys():
            #alright so here k is our NS
            #r is the index of the point for which we are calculating these prime values
            #adjacent is the list of adjacent points
            adjacent = neighbors[k][r]
            dxsum = 0
            dysum = 0
            dxs = []
            dys = []
            dhs = []
            dists = 0
            for i in adjacent:
                dx,dy = componentDistance(surfaces,k,i,r)
                dxs.append(dx)
                dys.append(dy)
            dxindexs = [np.argmin(dxs),np.argmax(dxs)]
            dyindexs = [np.argmin(dys),np.argmax(dys)]
            dxfinal,b = componentDistance(surfaces,k,adjacent[dxindexs[1]],adjacent[dxindexs[0]])
            b,dyfinal = componentDistance(surfaces,k,adjacent[dyindexs[1]],adjacent[dyindexs[0]])
            #print(r,adjacent)
            dhx = surfaces[k][2][0][adjacent[dxindexs[1]]] - surfaces[k][2][0][adjacent[dxindexs[0]]]
            #dhx = surfaces[k][2][0][adjacent[dxindexs[1]]]
            dhy = surfaces[k][2][0][adjacent[dyindexs[1]]]-surfaces[k][2][0][adjacent[dyindexs[0]]]
            #dhy = surfaces[k][2][0][adjacent[dyindexs[1]]] - surfaces[k][2][0][adjacent[dyindexs[0]]]
            dhdtheta = dhx/dxfinal
            dhdr = dhy/dyfinal
            #surfaces[k][2][4][r] = dhdtheta *(1/((90-surfaces[k][2][0][r])*111000))*(1/np.tan(np.deg2rad(surfaces[k][2][0][r])))
            surfaces[k][2][4][r] = dhdtheta
            surfaces[k][2][5][r] = dhdr 
    return surfaces

def addStreamFunc(surfaces,profiles):
    neutraldepths={}
    for k in surfaces.keys():
        for i in range(len(surfaces[k][3])):
            if surfaces[k][3][i] not in neutraldepths:
                neutraldepths[surfaces[k][3][i]] =[[],[]]
            neutraldepths[surfaces[k][3][i]][0].append(k)
            neutraldepths[surfaces[k][3][i]][1].append(surfaces[k][2][0][i])
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
        if count > 100:
            break
        p = getProfileById(profiles,k)
        s.append(p.sals)
        t.append(p.temps)
        ip.append(p.pres)
        #p_ref.append([10.1235]*len(p.ipres))
        if len(np.unique(p.ipres))!=len(p.ipres):
            print("NOOOOOOOO")
        p_ref.append([0]*len(p.ipres))
        ns.append(neutraldepths[k][0])
        ns_p.append(neutraldepths[k][1])
        ks.append(k)
    
    results = geo_strf_isopycnal(s,t,ip,p_ref,ns,ns_p,ks)
    
    #with open('data/geoisopycnal.pickle', 'wb') as outfile:
        #pickle.dump(results, outfile)
    #print(results)

        
    

