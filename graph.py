from nstools import *

import os

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

def graphNeighbors(surfaces,neighbors):
    for k in neighbors.keys():
        for r in neighbors[k].keys():
            i = neighbors[k][r]
            plt.scatter(surfaces[k][0][i],surfaces[k][1][i])
            plt.scatter(surfaces[k][0][r],surfaces[k][1][r],c="red")
            plt.show()

def zoomGraph(m,ax):
    lllon = -136
    urlon = 78
    lllat = 55
    urlat = 63

    xmin, ymin = m(lllon, lllat)
    xmax, ymax = m(urlon, urlat)

    ax = plt.gca()

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

def graphSurfacesComparison(surfaces,overlay,quantindex,contour=False,profiles=None,deepestindex=None,show=True,maximize=True,savepath=None):
    newsurfaces = {}
    for k in surfaces.keys():
        tempSurf = emptySurface()
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
    graphSurfaces(newsurfaces,quantindex,contour,profiles,deepestindex,show,maximize,savepath)

def graphSurfaces(surfaces,quantindex,contour=False,profiles=None,deepestindex=None,show=True,maximize=True,savepath=None,idlabels=False):
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
    for i in surfaces.keys():
        if len(surfaces[i]["lons"])>3 and len(surfaces[i]["data"][quantindex])>3:
            fig,ax = plt.subplots(1,1)
            mapy = Basemap(projection='ortho', lat_0=90,lon_0=-60)
            mapy.drawmapboundary(fill_color='aqua')
            mapy.fillcontinents(color='coral',lake_color='aqua')
            mapy.drawcoastlines()
            x,y = mapy(surfaces[i]["lons"],surfaces[i]["lats"])
            d = np.asarray(surfaces[i]["data"][quantindex])
            ids = np.asarray(surfaces[i]["ids"])
            x = np.asarray(x)
            y = np.asarray(y)
            #Plot the surface 
            if contour:
                plt.tricontourf(x,y,np.asarray(surfaces[i]["data"][quantindex]),cmap="plasma")
            else:
                plt.scatter(x,y,c=d,cmap="plasma")
                m = np.nanmedian(d)
                s = np.nanstd(d)
                print("################")
                plt.clim(m-2*s,m+2*s)
                #plt.clim(i-400,i+400)
                mapy.colorbar()
            #map the reference profile
            if profiles and deepestindex:
                x,y = mapy(profiles[deepestindex].lon,profiles[deepestindex].lat)
                mapy.scatter(x,y,c="red")
            zoomGraph(mapy,ax)
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

def graphVectorField(surfaces,key1,key2,backgroundfield="pv",savepath=False,show=True):

    if savepath:
        try:
            os.makedirs(savepath+key1+key2)
        except FileExistsError as e:
            print(e)

    for k in surfaces.keys():
        fig,ax = plt.subplots(1,1)
        mapy = Basemap(projection='ortho', lat_0=90,lon_0=-60)
        mapy.drawmapboundary(fill_color='aqua')
        mapy.fillcontinents(color='coral',lake_color='aqua')
        mapy.drawcoastlines()
        urs=[]
        uthetas=[]
        lons = []
        lats = []
        for p in range(0,len(surfaces[k]["data"]["uabs"])):
            u = surfaces[k]["data"][key1][p] 
            v = surfaces[k]["data"][key2][p]
            x = surfaces[k]["x"][p]
            y = surfaces[k]["y"][p]
            theta = np.deg2rad(surfaces[k]["lons"][p])
            #ur = u*np.cos(theta) + v*np.sin(theta)
            r = np.sqrt(x**2+y**2)
            ur = -(x*u +y*v)/r

            #utheta = v*np.cos(theta) - v*np.sin(theta)
            utheta = r*(x*v-y*u)/(r**2)

            urs.append(ur)
            uthetas.append(utheta)
            lons.append(surfaces[k]["lons"][p])
            lats.append(surfaces[k]["lats"][p])

        urs.append(0.1)
        uthetas.append(0)
        lons.append(89)
        lats.append(89)

        urs.append(0)
        uthetas.append(0.1)
        lons.append(89)
        lats.append(89)

        fig.suptitle(key1+"," + key2 + " NS: "+str(k))
        urs = np.asarray(urs)
        uthetas = np.asarray(uthetas)
        lons = np.asarray(lons)
        lats = np.asarray(lats)
        u,v,x,y = mapy.rotate_vector(uthetas,urs,lons,lats,returnxy=True)
        mag = np.sqrt(u**2+v**2)
        zoomGraph(mapy,ax)
        fig.set_size_inches(16.5,12)
        a = np.where(abs(surfaces[k]["lats"]-90)>0.5)
        xpv,ypv = mapy(surfaces[k]["lons"][a],surfaces[k]["lats"][a])
        if backgroundfield != "f/h":
            bgfield = surfaces[k]["data"][backgroundfield][a]
        else:
            bgfield = gsw.f(surfaces[k]["lats"][a])/surfaces[k]["data"]["z"][a]
        plt.tricontourf(xpv,ypv,bgfield,levels=50)
        plt.clim(np.min(bgfield),np.max(bgfield))
        mapy.colorbar()
        mapy.quiver(x,y,u,v,mag,cmap="cool",width = 0.002)
        if savepath:
            plt.savefig(savepath+key1+key2+"/ns"+str(k)+".png")
        if show:
            plt.show()
        plt.close()




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

#def graphAndSaveAll(surfaces,savepath):
    #for k in surfaces[200]["data"].keys():
        #graphSurfaces(surfaces,k,show=False,savepath=savepath)
        #if contour:
            #plt.tricontourf(x,y,np.asarray(surfaces[i]["data"][quantindex]),cmap="plasma")
        #else:
            #plt.scatter(x,y,c=np.asarray(surfaces[i]["data"][quantindex]),cmap="plasma")
def framStraitTransport(surfaces):
    transectids=[1776,1901,2026,2151,2276,2401,2526,2651]
    plotLayerTransport(transectids,surfaces)

def barentsTransport(surfaces):
    transectids=[2651,2650,2899,2896,3019,3143,3267]
    plotLayerTransport(transectids,surfaces)

def plotLayerTransport(transectids,surfaces):
    inflows=[]
    outflows=[]
    for k in sorted(surfaces.keys())[::-1]:
        transports = []
        go=True
        indexs = []
        for eyed in transectids:
            if eyed not in surfaces[k]["ids"]:
                go = False
                print(eyed)
            else:
                indexs.append(np.where(np.asarray(surfaces[k]["ids"]) == eyed)[0][0])
        print(go)
        if k==200:
            h=200
        else:
            h = np.nanmean(surfaces[k]["data"]["h"][indexs])
        if go and k<2200 and  not np.isnan(h):
            for j in range(len(transectids)-1):
                curr = indexs[j]
                nxt = indexs[j+1]
                dpsi = surfaces[k]["data"]["psinew"][nxt]-surfaces[k]["data"]["psinew"][curr]
                print("h: ",h," dpsi: ",dpsi)
                transports.append(-dpsi*h*(10**-6)*(1/gsw.f(surfaces[k]["lats"][curr])))
        print(transports)
        if len(transports)>0:
            transports = np.asarray(transports)
            inflows.append(np.nansum(transports[transports>0]))
            outflows.append(np.nansum(transports[transports<0]))
            plt.plot(range(len(transports)),transports,label=("ns: " + str(k)))
            plt.xlabel("Grid points across Fram Strait E to W")
            plt.ylabel("Transport in Sverdrups")
    plt.title("Tranport across Fram Strait \n inflow: "+str(round(np.nansum(inflows),2))+ " outflow: "+str(round(np.nansum(outflows),2)))
    plt.legend()
    plt.show()

def getLinePoints(surfaces,startid,endid,level=1000):
    if startid in surfaces[level]["ids"] and endid in surfaces[level]["ids"]:
        start = np.where(np.asarray(surfaces[level]["ids"]) == startid)[0][0]
        end = np.where(np.asarray(surfaces[level]["ids"]) == endid)[0][0]
        startcoord = (surfaces[level]["x"][start],surfaces[level]["y"][start])
        endcoord = (surfaces[level]["x"][end],surfaces[level]["y"][end])
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
    return None, None,None


def refinedTransport(surfaces,startid,endid):
    transports=[]
    inflows =[]
    outflows = []
    coords = [[],[]]
    scales = []
    for k in surfaces.keys():
        ps,vec,progress = getLinePoints(surfaces,startid,endid,k)
        surfacetransports = []
        if k<2600 and ps and vec and progress:
            hs = []
            for eyed in ps:
                curr = np.where(np.asarray(surfaces[k]["ids"]) == eyed)[0][0]
                hs.append(surfaces[k]["data"]["h"][curr])
            h=np.nanmean(hs)
            for i,eyed in enumerate(ps[:-1]):
                curr = np.where(np.asarray(surfaces[k]["ids"]) == eyed)[0][0]
                nxt = np.where(np.asarray(surfaces[k]["ids"]) == ps[i+1])[0][0]
                dist = geodesic((surfaces[k]["lats"][curr],surfaces[k]["lons"][curr]),(surfaces[k]["lats"][nxt],surfaces[k]["lons"][nxt])).m
                uvec = (surfaces[k]["x"][nxt]-surfaces[k]["x"][curr],surfaces[k]["y"][nxt]-surfaces[k]["y"][curr])
                f = gsw.f(surfaces[k]["lats"][curr])
                uvec = (-1,2)
                vec=(10,0)
                mag = np.sqrt(vec[0]**2 +vec[1]**2)
                proj = ((vec[0]*uvec[0],vec[1]*uvec[1])/(mag**2))*vec
                print("#####")
                print("set: ",uvec,vec,proj)
                scale = (np.sqrt(proj[0]**2 + proj[1]**2)/np.sqrt(uvec[0]**2 + uvec[1]**2))
                print("scale: ",scale)
                scale = 1
                dpsi = surfaces[k]["data"]["psinew"][nxt]-surfaces[k]["data"]["psinew"][curr]
                transport=scale*-dpsi*h*(10**-6)*(1/f)
                surfacetransports.append(transport)
                coords[0].append(surfaces[k]["lons"][curr])
                coords[1].append(surfaces[k]["lats"][curr])
                scales.append(scale)
            plt.plot(progress[:-1],surfacetransports,label = "ns: "+str(k))
            transports.append(surfacetransports)
            surfacetransports= np.asarray(surfacetransports)
            inflows.append(np.nansum(surfacetransports[surfacetransports>0]))
            outflows.append(np.nansum(surfacetransports[surfacetransports<0]))
 
            #print(np.sum(surfacetransports))
    print("inflow: ",np.nansum(inflows))
    print("outflow: ",np.nansum(outflows))
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

            

    
                

            
