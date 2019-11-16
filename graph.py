from nstools import *
from interptools import *

import matplotlib
import os
import datetime
import git
import cmocean

##Graph a given quantity over a transect given a surfaces object
## savepath and show control whether the images should be saved and whether graphs should be displayed
def graphTransects(surfaces,quantindex,contour=False,profiles=None,deepestindex=None,show=True,maximize=True,savepath=None):
    quanttitlehash = {0:"Pressure Dbar",1:"Temperature C",2:"Salinity PSU",3:"PV"}
    quantfilehash = {0:"PRES",1:"TEMP",2:"SAL",3:"PV"}
    fig,ax = plt.subplots(1,1)
    for k in surfaces.keys():
        if len(surfaces[k][0])>3:
            #fig,ax = plt.subplots(1,1)
            #Plot the surface 

            lats = surfaces[k][1]
            for i in range(len(surfaces[k][0])):
                if surfaces[k][0][i] <0:
                    lats[i] = -(lats[i]-80)+100
            lats = (lats-70)*111
            sort = np.argsort(lats)
            lats = lats[sort]
            depths = surfaces[k][2][quantindex][sort]
            a = np.where(abs(np.gradient(depths,lats))<7.5)
            lats = lats[a]
            depths = depths[a]

            plt.plot(lats,depths,label=str(k))
            #ax.set_ylim(m-1*s,m+1*s)
            ax.set_xlim(0,40*111)
                      #map the reference profile
            fig.suptitle("Transect at 40/-140 degrees Lon"  )
            #if maximize:
                #fig.set_size_inches(16.5,12)
            if savepath:
                plt.savefig(savepath+quantfilehash[quantindex]+"/ns"+str(i)+".png")

    lats=[]
    depths=[]
    deep = []
    ax.set_ylim(-4000,0)
    for lat in np.linspace(70,90,10000):
        lats.append((lat-70)*111)
        depths.append(bathtools.searchBath(lat,30))
        lats.append(((-(lat-70)+110)-70)*111)
        depths.append(bathtools.searchBath(lat,-150))
        deep.append(-100000000)
        deep.append(-100000000)
    #print("bath",depths)
    lats = np.asarray(lats)
    depths = np.asarray(depths)
    deep = np.asarray(deep)
    asc_x =  np.argsort(lats)
    ax.fill_between(lats[asc_x],depths[asc_x],deep[asc_x])
    #plt.colorbar()
    if show:
        plt.legend()
        plt.show()
    plt.close()

##diagnostic tool to show that the calculated nearest points for each
## point is correct
def graphNeighbors(surfaces,neighbors):
    print("graphing")
    for k in neighbors.keys():
        for r in neighbors[k]:
            plt.scatter(surfaces[k]["lons"][r[0]],surfaces[k]["lats"][r[0]],c="green")
            plt.scatter(surfaces[k]["lons"][r[1]],surfaces[k]["lats"][r[1]],c="blue")
            plt.scatter(surfaces[k]["lons"][r[2]],surfaces[k]["lats"][r[2]],c="yellow")
            plt.scatter(surfaces[k]["lons"][r[3]],surfaces[k]["lats"][r[3]],c="red")
        plt.show()

## zoom given map and axis into the arctic.
def zoomGraph(m,ax,region):
    if region == "arctic":
        lllon = -136
        urlon = 78
        lllat = 55
        urlat = 63
    if region == "nepb":
        lllon = 170
        urlon = -100
        lllat = 0
        urlat = 60
    if region != "nepbmerc":
        xmin, ymin = m(lllon, lllat)
        xmax, ymax = m(urlon, urlat)

        ax = plt.gca()

        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])

def mapSetup(coords,region="arctic",newax=True):
    regions = {"arctic":[90,-60],
                "nepb":[40,-150],
                "nepbmerc":[20,60,-170,-120]}
    if newax:
        fig,ax = plt.subplots(1,1)
    else:
        fig=None
        ax=None
    if region != "auto" and region != "nepbmerc":
        mapy = Basemap(projection='ortho', lat_0=regions[region][0],lon_0=regions[region][1],area_thresh=10)
    elif region == "nepbmerc":
        mapy = Basemap(projection='merc', llcrnrlat=regions[region][0],urcrnrlat=regions[region][1],\
            llcrnrlon=regions[region][2],urcrnrlon=regions[region][3],lat_ts=40,area_thresh=10)

    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    zoomGraph(mapy,ax,region)
    return fig,ax,mapy

## analagous to graphSurfaces but you can provide a second surfaces object
## usually one which is not interpolated to test interpolation
def graphSurfacesComparison(surfaces,overlay,quantindex,contour=False,profiles=None,deepestindex=None,show=True,maximize=True,savepath=None,region="arctic"):
    newsurfaces = {}
    for k in surfaces.keys():
        tempSurf = nstools.emptySurface()
        tempSurf["lons"] = np.concatenate((surfaces[k]["lons"], overlay[k]["lons"]))
        tempSurf["lats"] = np.concatenate((surfaces[k]["lats"],overlay[k]["lats"]))
        for field in surfaces[k]["data"].keys():
            if field in overlay[k]["data"].keys():
                tempSurf["data"][field] = np.concatenate((surfaces[k]["data"][field],overlay[k]["data"][field]))
            else:
                print(field," not in overlay")
        if "ids" in surfaces.keys() and "ids" in overlay.keys():
            tempSurf["ids"] = np.concatenate((surfaces[k]["ids"],overlay[k]["ids"]))
        newsurfaces[k]=tempSurf
    print(newsurfaces.keys())
    graphSurfaces(newsurfaces,quantindex,contour,profiles,deepestindex,show,maximize,savepath,region=region)

## given a surfaces object, a quantity index, graph quantity
## if you really want you can supply a profiles object and a deepest index to display a point
## controls to save, graph or maximize
def graphSurfaces(surfaces,quantindex,contour=False,profiles=None,deepestindex=None,\
        show=True,maximize=True,savepath=None,idlabels=False,colorlimit=True,region="arctic"):
    quanttitlehash = {"pres":"Pressure Dbar","t":"Temperature C","s":"Salinity PSU","pv":"PV",\
                     "u":"relative U","v":"relative V","psi":"ISOPYCNAL STREAMFUNCTION","hx":"Neutral Gradient X",\
                    "hy":"Neutral Gradient Y","curl":"Curl","drdt":"Northward Velocity",\
                    "dthetadt":"Eastward Velocity","ids":"IDS","uabs":"Absolute U","vabs":"Absolute V",\
                    "uprime":"reference U velocity","vprime":"reference V velocity","h":"Thickness of ",\
                    "CKVB":"KV term with roughness","CKVO":"KV term without roughness","dsdx":"Salinity X gradient",\
                    "dsdy":"Salinity Y gradient","d2sdx2":"Salinity X curvature",\
                    "d2sdy2":"Salinity Y curvature","n^2":"N^2"}
    if savepath:
        try:
            os.makedirs(savepath+quantindex)
        except FileExistsError as e:
            print(e)
        writeInfoFile(savepath)
    for i in list(surfaces.keys()):
        if quantindex in surfaces[i]["data"].keys() and \
                len(surfaces[i]["lons"])>3 and\
                len(surfaces[i]["data"][quantindex])>3:
            fig,ax,mapy = mapSetup(surfaces,region=region)
            mapy.drawparallels(np.linspace(-90,90,19))
            mapy.drawmeridians(np.linspace(-180, 180, 37))
            x,y = mapy(surfaces[i]["lons"],surfaces[i]["lats"])
            d = np.asarray(surfaces[i]["data"][quantindex])
            ids = np.asarray(surfaces[i]["ids"])
            x = np.asarray(x)
            y = np.asarray(y)
            #Plot the surface 
            if contour:
                plt.tricontourf(x,y,np.asarray(surfaces[i]["data"][quantindex]),cmap=cmocean.cm.haline,levels=30)
            else:
                plt.scatter(x,y,c=d,cmap=cmocean.cm.haline)
            m = np.nanmedian(d)
            s = np.nanstd(d)
            print("################")
            if colorlimit:
                plt.clim(m-1*s,m+1*s)
                #plt.clim(i-400,i+400)
                mapy.colorbar()
            #map the reference profile
            if profiles and deepestindex:
                x,y = mapy(profiles[deepestindex].lon,profiles[deepestindex].lat)
                mapy.scatter(x,y,c="red")
            if idlabels:
                for j, eyed in enumerate(ids):
                    ax.annotate(eyed,(x[j],y[j]))
            if quantindex in quanttitlehash.keys():
                fig.suptitle(str(quanttitlehash[quantindex]) + " at NS: "+str(i))
            else:
                fig.suptitle(str(quantindex) + " at NS: "+str(i))

            if maximize:
                fig.set_size_inches(16.5,12)
            if show:
                plt.show()
            if savepath:
                plt.savefig(savepath+quantindex+"/ns"+str(i)+".png")

            plt.close()

def graphSurfacesOneContour(surfaces,surfacesnc,quantindex,contour=False,profiles=None,deepestindex=None,\
        show=True,maximize=True,savepath=None,idlabels=False,colorlimit=True,region="arctic"):
    quanttitlehash = {"pres":"Pressure Dbar","t":"Temperature C","s":"Salinity PSU","pv":"PV",\
                     "u":"relative U","v":"relative V","psi":"ISOPYCNAL STREAMFUNCTION","hx":"Neutral Gradient X",\
                    "hy":"Neutral Gradient Y","curl":"Curl","drdt":"Northward Velocity",\
                    "dthetadt":"Eastward Velocity","ids":"IDS","uabs":"Absolute U","vabs":"Absolute V",\
                    "uprime":"reference U velocity","vprime":"reference V velocity","h":"Thickness of ",\
                    "CKVB":"KV term with roughness","CKVO":"KV term without roughness","dsdx":"Salinity X gradient",\
                    "dsdy":"Salinity Y gradient","d2sdx2":"Salinity X curvature",\
                    "d2sdy2":"Salinity Y curvature","n^2":"N^2"}
    if savepath:
        try:
            os.makedirs(savepath+quantindex)
        except FileExistsError as e:
            print(e)
        writeInfoFile(savepath)
    for i in list(surfaces.keys()):
        if i in surfacesnc.keys() and len(surfacesnc[i]["lons"])>3 and len(surfaces[i]["lons"])>3 and len(surfaces[i]["data"][quantindex])>3:
            fig,ax,mapy = mapSetup(surfaces,region=region)
            mapy.drawparallels(np.linspace(-90,90,19))
            mapy.drawmeridians(np.linspace(-180, 180, 37))
            x,y = mapy(surfaces[i]["lons"],surfaces[i]["lats"])
            d = np.asarray(surfaces[i]["data"][quantindex])
            notnan = ~np.isnan(d)
            ids = np.asarray(surfaces[i]["ids"])
            x = np.asarray(x)[notnan]
            y = np.asarray(y)[notnan]
            d = d[notnan]
            xnc,ync = mapy(surfacesnc[i]["lons"],surfacesnc[i]["lats"])
            dnc = np.asarray(surfacesnc[i]["data"][quantindex])
            idsnc = np.asarray(surfacesnc[i]["ids"])
            xnc = np.asarray(xnc)
            ync = np.asarray(ync)

            #Plot the surface 
            mean = np.nanmean(np.concatenate((d, dnc)))
            stdev = np.nanstd(np.concatenate((d,dnc)))
            mi = mean -stdev
            ma = mean + stdev
            norm = matplotlib.colors.Normalize(vmin=mi,vmax=ma)
            plt.tricontourf(x,y,d,norm=norm,cmap=cmocean.cm.haline,levels=30)
            plt.scatter(xnc,ync,norm=norm,c=dnc,cmap=cmocean.cm.haline)
            m = np.nanmedian(d)
            s = np.nanstd(d)
            print("################")
            if colorlimit:
                #plt.clim(m-2*s,m+2*s)
                #plt.clim(i-400,i+400)
                mapy.colorbar()
            #map the reference profile
            if profiles and deepestindex:
                x,y = mapy(profiles[deepestindex].lon,profiles[deepestindex].lat)
                mapy.scatter(x,y,c="red")
            if idlabels:
                for j, eyed in enumerate(ids):
                    ax.annotate(eyed,(x[j],y[j]))
            if quantindex in quanttitlehash.keys():
                fig.suptitle(str(quanttitlehash[quantindex]) + " at NS: "+str(i))
            else:
                fig.suptitle(str(quantindex) + " at NS: "+str(i))

            if maximize:
                fig.set_size_inches(16.5,12)
            if show:
                plt.show()
            if savepath:
                plt.savefig(savepath+quantindex+"/ns"+str(i)+".png")

            plt.close()


## plots a given cruise and its reference cruise
def plotCruiseAndRef(cruises,refcruises,show=True,region="arctic"):
    fig,ax,mapy = mapSetup([],region=region)
    lats, lons, depths=[],[],[]
    for p in refcruises:
        lats.append(p.lat)
        lons.append(p.lon)
        depths.append(np.max(p.pres))
    #fig.suptitle(cruisename)
    x,y = mapy(lons,lats)
    plt.scatter(x,y)
    lats, lons, depths=[],[],[]
    for p in cruises:
        lats.append(p.lat)
        lons.append(p.lon)
        depths.append(np.max(p.pres))
    #fig.suptitle(cruisename)
    a,b,mapy = mapSetup([],newax=False,region=region) 
    x,y = mapy(lons,lats)
    plt.scatter(x,y,c=depths,cmap="plasma")
    mapy.colorbar()
    if show:
        plt.show()

##plot a certain cruise
def plotCruise(profiles,cruisename,fig=None,ax=None,show=True,region="arctic"):
    lats, lons, depths=[],[],[]
    for p in profiles:
        lats.append(p.lat)
        lons.append(p.lon)
        depths.append(np.max(p.pres))
    
    if not fig and not ax:
        mapSetup([],region=region)
    else:
        mapSetup([],region=region,newax=False)

    fig.suptitle(cruisename)

    x,y = mapy(lons,lats)
    plt.scatter(x,y,c=depths,cmap="plasma")
    mapy.colorbar()
    if show:
        plt.show()

##plot profiles with an option to supply one profile which should be highlighted
def plotProfiles(profiles,title,specialprofile=None,fig=None,ax=None,show=True,data="pres",depth=False,region="arctic"):
    lats, lons, depths=[],[],[]
    for p in profiles:
        lats.append(p.lat)
        lons.append(p.lon)
        if data == "pres":
            depths.append(np.max(p.pres))
        if data == "t" and depth:
            depths.append(p.atPres(depth)[0])

    if not fig and not ax:
        fig,ax,mapy=mapSetup([],region=region)
    else:
        fig,ax,mapy=mapSetup([],region=region,newax=False)


    fig.suptitle(title)
    x,y = mapy(lons,lats)
    plt.scatter(x,y,c=depths,cmap="plasma")
    mapy.colorbar()
    if specialprofile:
        x,y = mapy(specialprofile.lon,specialprofile.lat)
        plt.scatter(x,y,c="red")
    if show:
        plt.show()

##graph surfaces and then also the interpolated surfaces transect
## for a given quantity
def graphComparisonTransects(surfaces,interpsurfaces,profiles,quantindex,contour=False,deepestindex=None,show=True,maximize=True,savepath=None):
    quanttitlehash = {0:"Pressure Dbar",1:"Temperature C",2:"Salinity PSU",3:"PV"}
    quantfilehash = {0:"PRES",1:"TEMP",2:"SAL",3:"PV"}
    for k in surfaces.keys():
        if len(surfaces[k][0])>3:
            fig,ax = plt.subplots(1,1)

            #Plot the interpolated surface 

            lats = interpsurfaces[k][1]
            for i in range(len(interpsurfaces[k][0])):
                if interpsurfaces[k][0][i] <0:
                    lats[i] = -(lats[i]-80)+100
            lats = (lats-70)*111
            sort = np.argsort(lats)
            lats = lats[sort]
            depths = interpsurfaces[k][2][quantindex][sort]
            a = np.where(abs(np.gradient(depths,lats))<7.5)
            lats = lats[a]
            depths = depths[a]
            plt.plot(lats,depths,label=str(k))
            ax.set_xlim(0,40*111)
            fig.suptitle("Transect at 40/-140 degrees Lon at "+str(k) +" meters"  )

            #plot raw Neutral depth readings colorcoded by year
            years=[]
            lats = surfaces[k][1]
            for i in range(len(surfaces[k][1])):
                years.append(getProfileById(profiles,surfaces[k][3][i]).time.year)
                if surfaces[k][0][i] <0:
                    lats[i] = -(surfaces[k][1][i]-80)+100
            lats = (lats-70)*111

            plt.scatter(lats,surfaces[k][2][quantindex],c=years)
            plt.colorbar()

            if maximize:
                fig.set_size_inches(16.5,12)
            if savepath:
                plt.savefig(savepath+quantfilehash[quantindex]+"/ns"+str(k)+".png")
            if show:
                plt.show()

##graph surface with neighbors to test calculation of nearest points
def graphStaggeredSurface(surfaces,neighbors,debug=False):
    for k in surfaces.keys():
        surfaces[k]["data"]["uz"] = np.zeros(len(surfaces[k]["lons"]))
        surfaces[k]["data"]["vz"] = np.zeros(len(surfaces[k]["lons"]))
    alldxs = []
    
    for k in surfaces.keys():
        fig,ax = plt.subplots(1,1)
        mapy = Basemap(projection='ortho', lat_0=90,lon_0=-60)
        mapy.drawmapboundary(fill_color='aqua')
        mapy.fillcontinents(color='coral',lake_color='aqua')
        mapy.drawcoastlines()
        for i in neighbors[k]:
            i = np.asarray(i)
                
            x,y = mapy(surfaces[k]["lons"][i][[0,2,3,1,0]],surfaces[k]["lats"][i][[0,2,3,1,0]])
            plt.plot(x,y)
        plt.show()

    return surfaces

## graph a vector field given a surfaces object on a map
## any quantity can be supplied as a background field
def graphVectorField(surfaces,key1,key2,backgroundfield="pv",region="arctic",transform=True,savepath=False,show=True,metadata={}):

    if savepath:
        try:
            os.makedirs(savepath+key1+key2)
        except FileExistsError as e:
            print(e)
        writeInfoFile(savepath,metadata)
    for k in surfaces.keys():
        fig,ax,mapy = mapSetup([],region=region)
        urs=[]
        uthetas=[]
        lons = []
        lats = []
        for p in range(0,len(surfaces[k]["data"][key1])):
            u = surfaces[k]["data"][key1][p] 
            v = surfaces[k]["data"][key2][p]
            x = surfaces[k]["x"][p]
            y = surfaces[k]["y"][p]
            if transform:
                theta = np.deg2rad(surfaces[k]["lons"][p])
                #ur = u*np.cos(theta) + v*np.sin(theta)
                r = np.sqrt(x**2+y**2)
                ur = -(x*u +y*v)/r

                #utheta = v*np.cos(theta) - v*np.sin(theta)
                utheta = r*(x*v-y*u)/(r**2)

                urs.append(ur)
                uthetas.append(utheta)
            else:
                urs.append(v)
                uthetas.append(u)
            lons.append(surfaces[k]["lons"][p])
            lats.append(surfaces[k]["lats"][p])

        urs.append(0.01)
        uthetas.append(0)
        lons.append(-150)
        lats.append(40)

        urs.append(0)
        uthetas.append(0.01)
        lons.append(-150)
        lats.append(40)

        fig.suptitle(key1+"," + key2 + " NS: "+str(k))
        urs = np.asarray(urs)
        uthetas = np.asarray(uthetas)
        lons = np.asarray(lons)
        lats = np.asarray(lats)
        u,v,x,y = mapy.rotate_vector(uthetas,urs,lons,lats,returnxy=True)
        mag = np.sqrt(u**2+v**2)
        fig.set_size_inches(16.5,12)
        a = tuple([(abs(surfaces[k]["lats"]-90)>0.5) & (~np.isnan(surfaces[k]["data"][backgroundfield]))])
        if np.count_nonzero(a)>4:
            xpv,ypv = mapy(surfaces[k]["lons"][a],surfaces[k]["lats"][a])
            if backgroundfield != "f/h":
                bgfield = surfaces[k]["data"][backgroundfield][a]
            else:
                bgfield = gsw.f(surfaces[k]["lats"][a])/surfaces[k]["data"]["z"][a]
            plt.tricontourf(xpv,ypv,bgfield,levels=50,cmap="viridis")
            #plt.scatter(xpv,ypv,c=bgfield)
            m = np.nanmedian(bgfield)
            s = np.nanstd(bgfield)
            plt.clim(m-2*s,m+2*s)
            mapy.colorbar()
            mapy.quiver(x,y,u*2,v*2,mag,cmap="autumn",width = 0.004)
            if savepath:
                plt.savefig(savepath+key1+key2+"/ns"+str(k)+".png")

            if show:
                plt.show()
        plt.close()



def writeInfoFile(savepath,metadata=None):
    infofile = open(savepath+"info.txt","w")
    infofile.write("created: \n"+str(datetime.datetime.now())+"\n")
    infofile.write("#"*10+"\n")
    if metadata:
        infofile.write("hello. I am a small autogenerated note! the order of mixs is [kvo,kvb,kvh]. The order of scalecoeffs is [Ar,kvo,kvb,kvh]\n")
        infofile.write("metadata: \n" + str(metadata)+"\n")
        infofile.write("#"*10+"\n")
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    infofile.write("githash: \n" + str(sha)+"\n")
    infofile.write("#"*10+"\n")
    infofile.close()





#graphs vectors but not on map because sometimes thats hard and you want to check
def graphCartVectorField(surfaces,key1,key2,show=True,savepath=False):
    if savepath:
        try:
            os.makedirs(savepath+quantindex)
        except FileExistsError as e:
            print(e)
    for k in surfaces.keys():
        fig,ax = plt.subplots(1,1)
        us=[]
        vs=[]
        xs = []
        ys = []
        for p in range(0,len(surfaces[k]["data"]["uabs"]),2):
            u = surfaces[k]["data"][key1][p] 
            v = surfaces[k]["data"][key2][p]
            x = surfaces[k]["x"][p]
            y = surfaces[k]["y"][p]
            us.append(u)
            vs.append(v)
            xs.append(x)
            ys.append(y)

        fig.suptitle("NS: "+str(k))
        us = np.asarray(us)
        vs = np.asarray(vs)
        xs = np.asarray(xs)
        ys = np.asarray(ys)
        mag = np.sqrt(us**2 + vs**2)
        fig.set_size_inches(16.5,12)
        a = np.where(abs(surfaces[k]["lats"]-90)>0.5)
        plt.tricontourf(surfaces[k]["x"][a],surfaces[k]["y"][a],surfaces[k]["data"]["pv"][a])
        plt.quiver(xs,ys,us,vs,mag,cmap="spring",scale=1)
        if savepath:
            plt.savefig(savepath+key1+key2+"/ns"+str(k)+".png")
        if show:
            plt.show()

##simple file to plot some beta spirals
def twentyRandomSpirals(surfaces,reflevel=200):
    for index in np.random.choice(range(len(surfaces[reflevel]["x"])),40):
        eyed = int(surfaces[reflevel]["ids"][index])
        us = []
        uabs = []
        vs = []
        vabs = []
        fig,ax = plt.subplots(1,1)
        for k in sorted(list(surfaces.keys())):
            found = np.where(np.asarray(surfaces[k]["ids"])==eyed)
            if len(found)!=0 and len(found[0]) != 0:
                found = found[0][0]
                us.append(surfaces[k]["data"]["u"][found])
                vs.append(surfaces[k]["data"]["v"][found])
                uabs.append(surfaces[k]["data"]["uabs"][found])
                vabs.append(surfaces[k]["data"]["vabs"][found])

        ax.plot(us,vs,color="red",label="relative current")
        ax.plot(uabs,vabs,color="blue",label="absolute current")
        ax.legend()
        plt.show()
    
#function to plot a spiral
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

##given two start coords find nearest surface points
#coords of type ( lon,lat)
def getLinePoints(surfaces,startcoord,endcoord,level):
    startcoord = singleXY(startcoord)
    endcoord = singleXY(endcoord)
    slope=(endcoord[1]-startcoord[1])/(endcoord[0]-startcoord[0])
    res = 100
    ids = []
    progress = []
    for i in range(res):
        currx = ((endcoord[0]-startcoord[0])/res)*i+startcoord[0]
        curry = ((endcoord[0]-startcoord[0])/res)*i*slope +startcoord[1]
        #print("xdiff0: ",startcoord[0]-currx," ydiff0: ",startcoord[1]-curry)
        #print("xdif1: ",endcoord[0]-currx," ydif1: ",endcoord[1]-curry)
        smallestval =1000000000000000
        closest = -1
        for j in range(len(surfaces[level]["x"])):
            dist = np.sqrt((surfaces[level]["x"][j]-currx)**2+(surfaces[level]["y"][j]-curry)**2)
            if dist<smallestval and ~np.isnan(surfaces[level]["data"]["psinew"][j]):
                smallestval = dist
                closest = j
        if surfaces[level]["ids"][closest] not in ids:
            ids.append(surfaces[level]["ids"][closest])
            progress.append(i)
    return ids, (endcoord[0]-startcoord[0],endcoord[1]-startcoord[1]),progress

## graph how quantity changes along a surface and along given transect
#coords of type ( lon,lat)
def quantityLine(surfaces,startcoord,endcoord,quant,silldepth,factor=1):
    quants=[]
    coords = [[],[]]
    for k in surfaces.keys():
        ps,vec,progress = getLinePoints(surfaces,startcoord,endcoord,k)
        surfacequants = []

        if k<silldepth and ps and vec and progress:
            hs = []
            for i,eyed in enumerate(ps):
                curr = np.where(np.asarray(surfaces[k]["ids"]) == eyed)[0][0]
                surfacequants.append(surfaces[k]["data"][quant][curr]*factor)
                coords[0].append(surfaces[k]["lons"][curr])
                coords[1].append(surfaces[k]["lats"][curr])
            plt.plot(progress,surfacequants,label = "ns: "+str(k))
            quants.append(surfacequants)
            surfacetransports= np.asarray(surfacequants)
 
            #print(np.sum(surfacetransports))
    plt.legend()
    plt.show()
    mapy = Basemap(projection='ortho', lat_0=90,lon_0=-60)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()
    x,y = mapy(coords[0],coords[1])
    plt.scatter(x,y)
    mapy.colorbar()
    plt.show()

 
#graph transport across a line given a surfaces object
#coords of type ( lon,lat)
def transportLine(surfaces,startcoord,endcoord,silldepth,uway=False,show=False):
    transports=[]
    inflows =[]
    outflows = []
    coords = [[],[]]
    scales = []
    for k in surfaces.keys():
        ps,vec,progress = getLinePoints(surfaces,startcoord,endcoord,k)
        surfacetransports = []

        if k<silldepth and ps and vec and progress:
            hs = []
            for eyed in ps:
                curr = np.where(np.asarray(surfaces[k]["ids"]) == eyed)[0][0]
                hs.append(surfaces[k]["data"]["h"][curr])
            h = np.nanmean(hs)
            for i,eyed in enumerate(ps):
                if uway:
                    curr = np.where(np.asarray(surfaces[k]["ids"]) == eyed)[0][0]
                    if i != len(ps)-1:
                        nxt = np.where(np.asarray(surfaces[k]["ids"]) == ps[i+1])[0][0]
                        dist = geodesic((surfaces[k]["lats"][curr],surfaces[k]["lons"][curr]),(surfaces[k]["lats"][nxt],surfaces[k]["lons"][nxt])).m
                    uvec = (surfaces[k]["data"]["uabs"][curr],surfaces[k]["data"]["vabs"][curr])
                    mag = np.sqrt(vec[0]**2 +vec[1]**2)
                    proj = ((vec[0]*uvec[0],vec[1]*uvec[1])/(mag**2))*vec
                    #print("#####")
                    #print("set: ",uvec,vec,proj)
                    perp = uvec-proj
                    sign = np.sign(np.cross(vec,perp))
                    perp = np.sqrt(perp[0]**2 + perp[1]**2)
                    scales.append(sign)
                    transport=sign*perp*dist*h*(10**-6)
                    surfacetransports.append(transport)
                    coords[0].append(surfaces[k]["lons"][curr])
                    coords[1].append(surfaces[k]["lats"][curr])
                else:
                    if i != len(ps)-1:
                        curr = np.where(np.asarray(surfaces[k]["ids"]) == eyed)[0][0]
                        nxt = np.where(np.asarray(surfaces[k]["ids"]) == ps[i+1])[0][0]
                        dpsi = surfaces[k]["data"]["psinew"][nxt]-surfaces[k]["data"]["psinew"][curr]
                        f = gsw.f(surfaces[k]["lats"][curr])
                        transport=dpsi*h*(10**-6)*(1/f)
                        surfacetransports.append(transport)
                        coords[0].append(surfaces[k]["lons"][curr])
                        coords[1].append(surfaces[k]["lats"][curr])
                        scales.append(transport)
                    
            if show:
                if uway:
                    plt.plot(progress[:],surfacetransports,label = "ns: "+str(k))
                else: 
                    plt.plot(progress[:-1],surfacetransports,label = "ns: "+str(k))
            transports.append(surfacetransports)
            surfacetransports= np.asarray(surfacetransports)
            inflows.append(np.nansum(surfacetransports[surfacetransports>0]))
            outflows.append(np.nansum(surfacetransports[surfacetransports<0]))
 
            #print(np.sum(surfacetransports))
    inflow = np.nansum(inflows)
    outflow = np.nansum(outflows)
    print("inflow: ",np.nansum(inflows))
    print("outflow: ",np.nansum(outflows))
    if show:
        plt.legend()
        plt.show()
        mapy = Basemap(projection='ortho', lat_0=90,lon_0=-60)
        mapy.drawmapboundary(fill_color='aqua')
        mapy.fillcontinents(color='coral',lake_color='aqua')
        mapy.drawcoastlines()
        x,y = mapy(coords[0],coords[1])
        plt.scatter(x,y,c=scales)
        mapy.colorbar()
        plt.show()
    return inflow,outflow

##plot ts diagrams with neutral surfaces annotated
def tsNeutralExplore(profiles):
    fig,ax1 = plt.subplots()
    ns = {}
    for k in range(200,3900,200):
        ns[k] = [[],[]]
    for p in profiles[:]:
        i = p.presIndex(200)
        ax1.scatter(p.sals[i:],p.temps[i:],color="blue",s=0.2)
        for k in p.neutraldepth.keys():
            pres = p.neutraldepth[k]
            t,s = p.atPres(pres)
            ns[k][0].append(s)
            ns[k][1].append(t)
    flip = False
    for k in ns.keys():
        if flip:
            ax1.scatter(ns[k][0],ns[k][1],s=2,color="orange")
        else:
            ax1.scatter(ns[k][0],ns[k][1],s=2,color="red")
        flip = not flip
    fig.set_size_inches(16.5,12)
    fig.suptitle("Temperature and salinty with neutral surfaces overlayed in alternating orange and red")
    ax1.set_xlabel("Salinity (PSU)")
    ax1.set_ylabel("Temperature (C)")
    plt.show()


## graph a vector field given a surfaces object on a map
## any quantity can be supplied as a background field
def graphProfilesVectorField(profiles,depths=range(200,4000,200),savepath=False,show=True,region="arctic"):

    if savepath:
        try:
            os.makedirs(savepath+"eccouv")
        except FileExistsError as e:
            print(e)

    for k in depths[::-1]:
        urs=[]
        uthetas=[]
        lons = []
        lats = []
        bgfield = []
        fig,ax,mapy = mapSetup([],region=region)
        for p in profiles[::2]:
            if k in p.neutraldepth.keys() and 1000 in p.neutraldepth.keys():
                i = p.ipresIndex(p.neutraldepth[k])
                urs.append(p.knownv[i]-p.knownv[p.ipresIndex(1000)])
                uthetas.append(p.knownu[i]-p.knownu[p.ipresIndex(1000)])
                #urs.append(1)
                #uthetas.append(1)
                lons.append(p.lon)
                lats.append(p.lat)
                bgfield.append(p.atPres(p.ipres[i])[1])

        urs.append(0.1)
        uthetas.append(0)
        lons.append(89)
        lats.append(89)
        bgfield.append(np.nan)

        urs.append(0)
        uthetas.append(0.1)
        lons.append(89)
        lats.append(89)
        bgfield.append(np.nan)

        fig.suptitle("knownu,knownv NS: "+str(k))
        urs = np.asarray(urs)
        uthetas = np.asarray(uthetas)
        lons = np.asarray(lons)
        lats = np.asarray(lats)
        u,v,x,y = mapy.rotate_vector(uthetas,urs,lons,lats,returnxy=True)
        mag = np.sqrt(urs**2+uthetas**2)
        fig.set_size_inches(16.5,12)

        xpv,ypv = mapy(lons,lats)

        #plt.scatter(xpv,ypv,bgfield)
        #plt.clim(np.min(bgfield),np.max(bgfield))
        mapy.colorbar()
        mapy.quiver(x,y,u,v,mag,cmap="plasma",width = 0.002)
        if savepath:
            plt.savefig(savepath+"eccouv"+"/ns"+str(k)+".png")
        if show:
            plt.show()
        plt.close()

def distanceHist(distances):
    for k in distances.keys():
        plt.hist(distances[k].values())
        plt.show()

def saveAllQuants(surfaces,savepath,region="arctic"):
    for d in surfaces[list(surfaces.keys())[0]]["data"].keys():
        graphSurfaces(surfaces,d,show=False,savepath=savepath,region=region)

        

            

    
                

            
