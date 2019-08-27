from nstools import *

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

##sometimes points are too close together and the interpolation
## loses it so we just smooth em
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

#create a mesh for interpolation
## here I put in my favorite x and y coordinates for the arctic
## they should be hardcoded because you want the grid points to 
## align vertically throughout the water column
def createMesh(n,xvals,yvals,custom=False):
    if custom:
        return np.meshgrid(np.linspace(np.min(xvals),np.max(xvals),n), np.linspace(np.min(yvals),np.max(yvals),n),indexing="xy")
    else:
        xmin=-1793163
        xmax=971927
        ymin=-1455096
        ymax=1200385
        return np.meshgrid(np.linspace(xmin,xmax,n), np.linspace(ymin,ymax,n),indexing="xy")

#generate a mesh and remove points in that mesh 
#which are too far away from locations with observations
def generateMaskedMesh(x,y,radius=200):
    xi,yi = createMesh(50,x,y)
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
        if final[0][i] > 4:
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

##interpolate  a surface
## create the mesh, use gam interpolation
##also returns neighboring points for each points
## and the distance between those points
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

## interpolate all the surfaces vertically and store
## neighbors, and distances as well
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
def addLatLonToSurface(surface,debug = True):
    lat = 90-(np.sqrt((surface["x"]**2+surface["y"]**2))/111000.0)
    lon = np.degrees(np.arctan2(surface["y"],surface["x"]))
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
    return lookup
