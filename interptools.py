import bathtools
import nstools
import pygam
import numpy as np
import itertools
from geopy.distance import geodesic
from progress.bar import Bar
import pdb
import alphashape
from shapely.geometry import Point
from shapely import affinity
import matplotlib.pyplot as plt
from descartes import PolygonPatch

#add x and y to surfaces. x and y necesarry for interpolation
def addXYToSurfaces(surfaces,debug=True):
    if debug:
        print("converting surfaces to xyz")
    newsurfaces = {}
    for k in surfaces.keys():
        x,y = homemadeXY(surfaces[k]["lons"],surfaces[k]["lats"])
        surfaces[k]["x"]=x
        surfaces[k]["y"]=y
    return surfaces

def fillOutEmptyFields(surfaces):
    for k in surfaces.keys():
        datafields = ["u","v","hx","h","CKVB","hy","t","s","pv","pres",\
                     "curl","uabs","vabs","uprime","vprime","dsdx","dsdz","dsdy",\
                    "d2sdx2","d2sdy2","dtdx","dtdy","dpdx","dpdy","n^2",\
                    "dqnotdx","dqnotdy","d2thetads2","dalphadtheta",\
                    "alpha","beta","dalphads","dbetads","dalphadp",\
                    "dbetadp","psi","dqdz","dqdx","dqdy","toph","both",\
                    "d2qdz2","d2qdx2","d2qdy2","khp","khpdz","dpsidx","dpsidy"]
        for d in datafields:
            if d not in surfaces[k]["data"].keys():
                surfaces[k]["data"][d] = np.full(len(surfaces[k]["lons"]),np.nan)
    return surfaces


##sometimes points are too close together and the interpolation
## loses it so we just smooth em
def removeDiscontinuities(surface,radius=10,debug=True,inside={}):
    x=np.asarray(surface["x"])
    y=np.asarray(surface["y"])
    z=np.asarray(surface["data"]["pres"])
    final = np.zeros(x.shape)
    #print(x)
    for i in Bar("Removing discontinuities: ").iter(range(len(x))):
        if final[i] == False:
            r = (x- x[i])**2 + (y - y[i])**2
            inside = r<radius*(10**6)
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

def indexBoolMatrix(boolmatrix):
    indexcount = np.full(boolmatrix.shape,np.nan)
    count = 0
    for i in range(boolmatrix.shape[0]):
        for j in range(boolmatrix.shape[1]):
            if boolmatrix[i][j]:
                indexcount[i][j] = count
                count+=1
    return indexcount

def simpleBoundary(gridx,gridy,datax,datay,radius):
    radius=150
    mask = np.zeros(gridx.shape)
    for i in range(len(datax)):
        r = np.sqrt((gridx- datax[i])**2 + (gridy - datay[i])**2)
        inside = r<radius*1000
        mask = mask+inside
    return mask

def geoMask(gridx,gridy,datax,datay,radius):
    mask = simpleBoundary(gridx,gridy,datax,datay,radius)
    points = []
    for l in range(len(mask)):
        for j in range(len(mask[l])):
            if mask[l][j]:
                points.append((gridx[l][j],gridy[l][j]))
    shape = alphashape.alphashape(points,alpha=0)
    shape = affinity.scale(shape,xfact=0.99,yfact=0.99)
    newmask = np.zeros(gridx.shape)
    #plt.plot(*shape.exterior.xy)
    #plt.show()
    for l in range(len(newmask)):
        for j in range(len(newmask[l])):
            p =Point(gridx[l][j],gridy[l][j])
            if shape.contains(p):
                newmask[l][j] = 1
    return newmask


## neigbor in "x" direction
def findRowNeighbor(row,col,mask,indexcount):
    for i in range(1,3):
        if col+i < len(mask[row]) and mask[row][col+i]:
            if ~np.isnan(indexcount[row][col+i]):
                return int(indexcount[row][col+i])
    return False
## neigbor in "x" direction
def findColumnNeighbor(row,col,mask,indexcount):
    for i in range(1,3):
        if row+i < len(mask) and mask[row+i][col]:
            if ~np.isnan(indexcount[row+i][col]):
                return int(indexcount[row+i][col])
    return False
#find neighbor in corner
def findCornerNeighbor(row,col,mask,indexcount):
    for i in range(1,3):
        if col+i < len(mask[row]) and row+i < len(mask) and mask[row][col+i]:
            if ~np.isnan(indexcount[row+i][col+i]):
                return int(indexcount[row+i][col+i])
    return False

#generate a mesh and remove points in that mesh 
#which are too far away from locations with observations
def generateMaskedMesh(x,y,region,coord,radius=50):
    xi,yi = region.createMesh(x,y,coord)
    #Make sure grid points are within original data point

    mask = geoMask(xi,yi,x,y,radius)
    indexcount = indexBoolMatrix(mask)
    finalxi=[]
    finalyi=[]
    neighbors = []
    ids=[]
    for row in range(mask.shape[0]):
        for col in range(mask.shape[1]):
            if ~np.isnan(mask[row][col]) and ~np.isnan(indexcount[row][col]):
                finalxi.append(xi[row][col])
                finalyi.append(yi[row][col])
                ids.append(row*mask.shape[1]+col)
                rowneighbor = findRowNeighbor(row,col,mask,indexcount)
                columnneighbor = findColumnNeighbor(row,col,mask,indexcount)
                cornerneighbor = findCornerNeighbor(row,col,mask,indexcount)
                if rowneighbor and columnneighbor and cornerneighbor:
                    neighbors.append(tuple((int(indexcount[row][col]),\
                            rowneighbor,columnneighbor,cornerneighbor)))
    return np.asarray(finalxi),np.asarray(finalyi),neighbors,ids

def isGridPointIsolated(row,col,mask,radius):
    for r in range(max(row-radius,0),min(row+radius,mask.shape[0]-1)):
        for c in range(max(col-radius,0),min(col+radius,mask.shape[1]-1)):
            if r != row and c != col and mask[r][c]:
                return False
    return True


def bathVarCacheWrapper(lat,lon,region):
    if not hasattr(bathVarCacheWrapper,"bvardict"):
        bathVarCacheWrapper.bvardict={}
    if (lat,lon) not in bathVarCacheWrapper.bvardict.keys():
        bathVarCacheWrapper.bvardict[(lat,lon)] = np.var(bathtools.bathBox(lat,lon))
    return np.abs(bathVarCacheWrapper.bvardict[(lat,lon)])

def bathVarMask(gridx,gridy,region,mask=[]):
    if np.isnan(mask).any():
        mask = np.invert(np.zeros(gridx.shape))
    for row in Bar("Bath var Masking: ").iter(range(mask.shape[0])):
        for col in range(mask.shape[1]):
            if mask[row][col]:
                if row%2 == 0 and col%2 ==0:
                    mask[row][col]=True
                #else:
                    #mask[row][col]=False

                elif col!=mask.shape[1]-1 and row != mask.shape[0]-1:
                    lat,lon = xyToLatLon(gridx[row][col],gridy[row][col])
                    #lat2,lon2 = xyToLatLon(gridx[row][col+1],gridy[row][col+1])
                    #lat3,lon3 = xyToLatLon(gridx[row+1][col],gridy[row+1][col])
                    #xbathchange = abs(bathtools.searchBath(lat,lon,"nepb")-bathtools.searchBath(lat2,lon2,"nepb"))
                    #ybathchange = abs(bathtools.searchBath(lat,lon,"nepb")-bathtools.searchBath(lat3,lon3,"nepb"))
                    #if xbathchange<300 and ybathchange < 300:
                    bvar = bathVarCacheWrapper(lat,lon,region)
                    if bvar < 4000:
                        mask[row][col]=False

    return mask


#generate a mesh and remove points in that mesh 
#which are too far away from locations with observations
def smartMesh(x,y,region,coord,radius=500):
    xi,yi = region.createMesh(x,y,coord,spacingscale=2)
    #Make sure grid points are within original data point
    mask = geoMask(xi,yi,x,y,radius)
    mask = bathVarMask(xi,yi,region,mask)
    indexcount = indexBoolMatrix(mask)
    finalxi=[]
    finalyi=[]
    neighbors = []
    ids=[]
    for row in range(mask.shape[0]):
        for col in range(mask.shape[1]):
            if ~np.isnan(mask[row][col]) and ~np.isnan(indexcount[row][col]):
                finalxi.append(xi[row][col])
                finalyi.append(yi[row][col])
                ids.append(row*mask.shape[1]+col)
                rowneighbor = findRowNeighbor(row,col,mask,indexcount)
                columnneighbor = findColumnNeighbor(row,col,mask,indexcount)
                cornerneighbor = findCornerNeighbor(row,col,mask,indexcount)
                if rowneighbor and columnneighbor and cornerneighbor:
                    neighbors.append(tuple((int(indexcount[row][col]),\
                            rowneighbor,columnneighbor,cornerneighbor)))
    return np.asarray(finalxi),np.asarray(finalyi),neighbors,ids



## snap points from profiles onto
## closest surface point
def surfaceSnap(surface,xgrid,ygrid):
    interpdata = {}
    for d in surface["data"].keys():
        interpdata[d]=[np.nan]*len(xgrid)
    for l in range(len(xgrid)):
        idx = np.argmin((xgrid[l]-surface["x"])**2 + (ygrid[l]-surface["y"])**2)
        for d in surface["data"].keys():
            if l < len(surface["data"][d]):
                interpdata[d][l] = float(surface["data"][d][idx])

    for d in surface["data"].keys():
        interpdata[d]=np.asarray(interpdata[d])

    return interpdata

def surfacePrune(surfaces):
    for k in surfaces.keys():
        for d in surfaces[k]["data"].keys():
            if d not in ["pres","s","t","pv","psi","alpha","beta"]:
                del surfaces[k]["data"][d]
    return surfaces
 
##interpolate  a surface
## create the mesh, use gam interpolation
##also returns neighboring points for each points
## and the distance between those points
def interpolateSurface(surface,region,coord="xy",debug=True,interpmethod="gam",smart=True,splines=10):
    #print("######")
    interpsurf={}
    X = np.zeros((len(surface["x"]),2))
    X[:,0]=surface["x"]
    X[:,1]=surface["y"]
    if smart:
        xi,yi,neighbors,finalids = smartMesh(surface["x"],surface["y"],region,coord)
    else:
        xi,yi,neighbors,finalids = generateMaskedMesh(surface["x"],surface["y"],region,coord)

    interpdata={}
    interpsurf["x"] =xi
    interpsurf["y"] =yi
    interpsurf["ids"] =finalids
    if len(xi) != len(finalids):
        print("OH NOOOOOO")
    if interpmethod=="gam":
        for d in Bar("Interpolating: ").iter(surface["data"].keys()):
            notnan = ~np.isnan(surface["data"][d])
            if np.count_nonzero(notnan)>10:
                gam = pygam.GAM(pygam.te(0,1,n_splines=[splines,splines])).fit(X[notnan],np.asarray(surface["data"][d])[notnan])
                Xgrid = np.zeros((yi.shape[0],2))
                Xgrid[:,0] = xi
                Xgrid[:,1] = yi
                interpdata[d] = gam.predict(Xgrid)
            else:
                interpdata[d] = np.asarray([np.nan]*len(xi))
    elif interpmethod in ["linear","nearest"] :
        for d in Bar("Interpolating: ").iter(surface["data"].keys()):
            notnan = np.logical_and(~np.isnan(X[:,0]),~np.isnan(surface["data"][d]))
            if np.count_nonzero(notnan)>10:
                Xgrid = np.zeros((yi.shape[0],2))
                Xgrid[:,0] = xi
                Xgrid[:,1] = yi
                f = nstools.griddata(X[notnan],np.asarray(surface["data"][d])[notnan],Xgrid,method=interpmethod)
                if np.isinf(f).any():
                    print("oh no!")
                interpdata[d] = f
            else:
                interpdata[d] = np.asarray([np.nan]*len(xi))

    else:
        interpdata = surfaceSnap(surface,xi,yi)
 
    interpsurf["data"] = interpdata
    interpsurf["data"]["ids"] = finalids
    interpsurf = addLatLonToSurface(interpsurf)
    return interpsurf,neighbors

## interpolate all the surfaces vertically and store
## neighbors, and distances as well
def interpolateSurfaces(region,surfaces,coord="xy",debug=True,interpmethod="gam",smart=False,splines=10):
    surfaces = addXYToSurfaces(surfaces)
    interpolatedsurfaces = {}
    neighbors={}
    lookups={}
    for k in surfaces.keys():
        if (~np.isnan(surfaces[k]["data"]["pres"])).any():
            #surfaces[k] = removeDiscontinuities(surfaces[k],radius=0.1)
            interpolatedsurfaces[k],neighbors[k] = interpolateSurface(surfaces[k],region,coord=coord,interpmethod=interpmethod,smart=smart,splines=splines)
            lookups[k] = trueDistanceLookup(interpolatedsurfaces[k],neighbors[k])
    interpolatedsurfaces = fillOutEmptyFields(interpolatedsurfaces)
    return interpolatedsurfaces,neighbors,lookups

##lat lon to x y
def singleXY(coord):
    theta = np.deg2rad(coord[0])
    r = ((90-coord[1]) *111*1000)
    x = (r*np.cos(theta))
    y = (r*np.sin(theta))
    return x,y


#lat lon to x y for a bunch of lat lons
def homemadeXY(lon,lat):
    x=[]
    y=[]
    for i in range(len(lat)):
        theta = np.deg2rad(lon[i])
        r = ((90-lat[i]) *111*1000)
        x.append(r*np.cos(theta))
        y.append(r*np.sin(theta))
    return np.asarray(x),np.asarray(y)

#after interpolation everything is in x and y so we need to convert back 
# to latitude and longitude
def xyToLatLon(x,y):
    lat = 90-(np.sqrt((x**2+y**2))/111000.0)
    lon = np.degrees(np.arctan2(y,x))
    return lat,lon
    
def addLatLonToSurface(surface,debug = True):
    lat,lon = xyToLatLon(surface["x"],surface["y"])
    #print("lat: ",lat, " lon: ",lon)
    surface["lons"]=lon
    surface["lats"]=lat
    return surface
 
#what is the true distance between neighbors
def trueDistanceLookup(surface,neighbors):
    lookup = {}
    for square in Bar("distance calc: ").iter(neighbors):
        for edge in itertools.combinations(square,2):
            p = tuple(sorted(edge))
            if p not in lookup.keys():
                lookup[p] = geodesic((surface["lats"][p[0]],surface["lons"][p[0]]),(surface["lats"][p[1]],surface["lons"][p[1]])).m
                lookup[p[::-1]] = lookup[p]
                #latdist = abs(surface["lats"][p[0]] - surface["lats"][p[1]])*111.0*1000.0
                #londist = abs(surface["lons"][p[0]] - surface["lons"][p[1]]) *np.cos(np.deg2rad(surface["lats"][p[1]]))*111.0*1000.0
                #lookup[p] = latdist+londist 
    return lookup


