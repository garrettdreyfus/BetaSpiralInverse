import matplotlib, os, datetime, git , cmocean, julian
from nstools import *
from interptools import *
from progress.bar import Bar
from copy import copy

quanttitlehash = {"pres":"Pressure Dbar","t":"Temperature C","s":"Salinity PSU","pv":"PV",\
                 "u":"relative U","v":"relative V","psi":"ISOPYCNAL STREAMFUNCTION","hx":"Neutral Gradient X",\
                "hy":"Neutral Gradient Y","curl":"Curl","drdt":"Northward Velocity",\
                "dthetadt":"Eastward Velocity","ids":"IDS","uabs":"Absolute U","vabs":"Absolute V",\
                "uprime":"reference U velocity","vprime":"reference V velocity","h":"Thickness of ",\
                "CKVB":"KV term with roughness","CKVO":"KV term without roughness","dsdx":"Salinity X gradient",\
                "dsdy":"Salinity Y gradient","d2sdx2":"Salinity X curvature",\
                "d2sdy2":"Salinity Y curvature","n^2":"N^2"}

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

def graphProfileTransect(profiles,quant,lat=None,lon=None,contour=False,\
        deepestindex=None,show=True,maximize=True,presrange=[2000,np.inf],\
        cbarrange=None,savepath=None,coordrange=None):
    fig,ax = plt.subplots(1,1)
    xcoord = []
    ycoord = []
    quants = []
    for p in profiles:
        if (lat and abs(lat-p.lat)<1) or (lon and abs(lon-p.lon)<1):
            for d in range(len(p.ipres)):
                if p.ipres[d] > presrange[0] and p.ipres[d] < presrange[1]:
                    quants.append(getattr(p,quant)[d])
                    ycoord.append(-p.ipres[d])

                    if lat:
                        xcoord.append(p.lon)
                    if lon:
                        xcoord.append(p.lat)
    if coordrange:
        plt.gca().set_xlim(coordrange)
    plt.scatter(xcoord,ycoord,c=quants)
    if cbarrange:
        plt.clim(cbarrange)
    plt.colorbar()
    if lat:
        fig.suptitle("Transect at "+str(lat)+" degrees Latitude"  )
        ax.set_xlabel("Longitude")
    if lon:
        fig.suptitle("Transect at "+str(lon)+" degrees Longitude"  )
        ax.set_xlabel("Latitude")

    if maximize:
        fig.set_size_inches(16.5,12)
    if show:
        plt.legend()
        plt.show()
    if savepath:
        try:
            os.makedirs(savepath)
        except FileExistsError as e:
            print(e)
        if lat:
            plt.savefig(savepath+"/lat"+str(lat)+".png")
        if lon:
            plt.savefig(savepath+"/lon"+str(lon)+".png")

    plt.close()



##diagnostic tool to show that the calculated nearest points for each
## point is correct
def graphNeighbors(surfaces,neighbors):
    print("graphing")
    for k in neighbors.keys():
        for r in neighbors[k]:
            #plt.scatter(surfaces[k]["lons"][r[0]],surfaces[k]["lats"][r[0]],c="green")
            #plt.scatter(surfaces[k]["lons"][r[1]],surfaces[k]["lats"][r[1]],c="blue")
            #plt.scatter(surfaces[k]["lons"][r[2]],surfaces[k]["lats"][r[2]],c="yellow")
            #plt.scatter(surfaces[k]["lons"][r[3]],surfaces[k]["lats"][r[3]],c="red")
            s=[]
            s.append((surfaces[k]["x"][r[0]],surfaces[k]["y"][r[0]]))
            s.append((surfaces[k]["x"][r[1]],surfaces[k]["y"][r[1]]))
            s.append((surfaces[k]["x"][r[2]],surfaces[k]["y"][r[2]]))
            s.append((surfaces[k]["x"][r[3]],surfaces[k]["y"][r[3]]))
            plt.plot(s)

        plt.show()

## zoom given map and axis into the arctic.
def zoomGraph(m,ax,region):
    if region != "nepbmerc":
        xmin, ymin = m(region.mapbounds["lllon"], region.mapbounds["lllat"])
        xmax, ymax = m(region.mapbounds["urlon"], region.mapbounds["urlat"])

        ax = plt.gca()

        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])

def mapSetup(coords,region,newax=True):
    if newax:
        fig,ax = plt.subplots(1,1)
    else:
        fig=None
        ax=None
    mapy = Basemap(projection='ortho', lat_0=region.mapbounds["lat_0"],lon_0=region.mapbounds["lon_0"],area_thresh=10)
    mapy.drawmapboundary(fill_color='aqua')
    mapy.fillcontinents(color='coral',lake_color='aqua')
    mapy.drawcoastlines()

    parallels = np.arange(-90.,91,10.)
    mapy.drawparallels(parallels)
    meridians = np.arange(0.,351.,10.)
    mapy.drawmeridians(meridians)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(meridians[i]),xy=mapy(meridians[i],region.mapbounds["lllat"]),xycoords='data')
    for i in np.arange(len(parallels)):
        plt.annotate(np.str(parallels[i]),xy=mapy(region.mapbounds["urlon"]-10,parallels[i]),xycoords='data')
    zoomGraph(mapy,ax,region)
    return fig,ax,mapy

## analagous to graphSurfaces but you can provide a second surfaces object
## usually one which is not interpolated to test interpolation
def graphSurfacesComparison(surfaces,overlay,quantindex,contour=False,stds=1,profiles=None,deepestindex=None,show=True,maximize=True,savepath=None,region="arctic"):
    newsurfaces = {}
    for k in Bar("Graphing Surfaces: ").iter(surfaces.keys()):
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
    graphSurfaces(region,newsurfaces,quantindex,contour,profiles,deepestindex,show,maximize,savepath)

def threedTransect(profiles,ks):
    ax = plt.axes()
    for k in ks:
        x=[]
        pres = []
        for p in profiles:
            if k in p.neutraldepth.keys():
                x.append(p.lat)
                pres.append(-p.neutraldepth[k])
        a = np.argsort(x)
        x= np.asarray(x)[a]
        pres = np.asarray(pres)[a]
        ax.plot(x, pres);
    plt.show()


    
def minmaxBound(d,s):
    return np.nanmin(d),np.nanmax(d)

def stdevBound(d,stds):
    m = np.nanmedian(d)
    s = np.nanstd(d)
    return m-stds*s,m+stds*s



## given a surfaces object, a quantity index, graph quantity
## if you really want you can supply a profiles object and a deepest index to display a point
## controls to save, graph or maximize
def graphSurfaces(region,surfaces,quantindex,stds=2,contour=False,profiles=None,deepestindex=None,\
        show=True,maximize=True,savepath=None,idlabels=False,additional="",\
        colorlimit=True,select=range(0,10000),secondsurface=None,centerfunction=False,boundfunc=stdevBound):
    if savepath:
        try:
            os.makedirs(savepath+quantindex)
        except FileExistsError as e:
            print(e)
        writeInfoFile(savepath)

    for i in list(surfaces.keys()):
        if abs(i) in select\
            and quantindex in surfaces[i]["data"].keys() and \
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
            if centerfunction:
                m,s = centerfunction(surfaces[i],quantindex)
            else:
                m = np.nanmedian(d)
                s = np.nanstd(d)

            if contour:
                a = np.logical_and((abs(np.asarray(surfaces[i]["lats"])-90)>0.5) , (~np.isnan(surfaces[i]["data"][quantindex])))
                b,inds = np.unique(list(zip(x,y)),axis=0,return_index=True)
                c=np.full_like(a,False)
                c[inds] = True
                a=np.logical_and(a,c)
                if np.count_nonzero(a)>4:
                    plt.tripcolor(x[a],y[a],np.asarray(surfaces[i]["data"][quantindex])[a],cmap=cmocean.cm.haline,vmin=m-stds*s,vmax=m+stds*s)
            else:
                plt.scatter(x,y,c=d,cmap=cmocean.cm.haline,vmin=m-stds*s,vmax=m+stds*s)
            #map the reference profile
            if profiles and deepestindex:
                x,y = mapy(profiles[deepestindex].lon,profiles[deepestindex].lat)
                mapy.scatter(x,y,c="red")
            if secondsurface:
                print("SECONDSURFACE")
                x,y = mapy(secondsurface[i]["lons"],secondsurface[i]["lats"])
                plt.scatter(x,y,c=secondsurface[i]["data"][quantindex],vmin=m-stds*s,vmax=m+stds*s,cmap=cmocean.cm.haline)

            if colorlimit:
                plt.clim(boundfunc(d,stds))
                #plt.clim(i-400,i+400)
                mapy.colorbar()
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
                plt.savefig(savepath+quantindex+"/ns"+additional+str(i)+".png")

            plt.close()

def graphSurfacesOneContour(region,surfaces,surfacesnc,quantindex,contour=False,profiles=None,deepestindex=None,\
        show=True,maximize=True,savepath=None,idlabels=False,colorlimit=True):
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
def plotCruiseAndRef(region,cruises,refcruises,show=True):
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
def plotCruise(region,profiles,cruisename,fig=None,ax=None,show=True):
    lats, lons, depths=[],[],[]
    for p in profiles:
        if p.cruise == cruisename:
            lats.append(p.lat)
            lons.append(p.lon)
            depths.append(np.max(p.pres))
    
    if not fig and not ax:
        fig,ax,mapy = mapSetup([],region=region)
    else:
        mapSetup([],region=region,newax=False)

    fig.suptitle(cruisename)

    x,y = mapy(lons,lats)
    plt.scatter(x,y,c=depths,cmap="plasma")
    mapy.colorbar()
    if show:
        plt.show()

##plot profiles with an option to supply one profile which should be highlighted
def plotProfiles(region,profiles,title,specialprofile=None,fig=None,ax=None,show=True,data="pres",depth=False):
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
    clb = mapy.colorbar()
    clb.ax.set_title('Depth of CTD cast')
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
    
    for k in Bar("Graphing Surfaces: ").iter(surfaces.keys()):
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
def graphVectorField(region,surfaces,key1,key2,backgroundfield="pv",refarrow=0.01,select=range(0,10000),stdevs=2,\
        transform=True,savepath=False,show=True,metadata={},contour=True,scale=1,boundfunc=stdevBound):

    if savepath:
        try:
            os.makedirs(savepath+key1+key2)
        except FileExistsError as e:
            print(e)
        writeInfoFile(savepath,metadata)
    for k in Bar("Graphing Surfaces: ").iter(surfaces.keys()):
        if k in select:
            fig,ax,mapy = mapSetup([],region=region)
            urs=[]
            uthetas=[]
            lons = []
            lats = []
            for p in range(0,len(surfaces[k]["data"][key1])):
                u = surfaces[k]["data"][key1][p] 
                v = surfaces[k]["data"][key2][p]
                if ~np.isnan(u) and ~np.isnan(v) and np.sqrt(u**2 + v**2)<10:
                    if transform:
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
                    else:
                        urs.append(v)
                        uthetas.append(u)
                    lons.append(surfaces[k]["lons"][p])
                    lats.append(surfaces[k]["lats"][p])

            urs.append(refarrow)
            uthetas.append(0)
            lons.append(region.mapbounds["lon_0"])
            lats.append(region.mapbounds["lat_0"])

            urs.append(0)
            uthetas.append(refarrow)
            lons.append(region.mapbounds["lon_0"])
            lats.append(region.mapbounds["lat_0"])

            fig.suptitle(key1+"," + key2 + " NS: "+str(k))
            urs = np.asarray(urs)
            uthetas = np.asarray(uthetas)
            lons = np.asarray(lons)
            lats = np.asarray(lats)
            u,v,x,y = mapy.rotate_vector(uthetas,urs,lons,lats,returnxy=True)
            mag = np.sqrt(u**2+v**2)
            fig.set_size_inches(16.5,12)
            a = ~np.isnan(surfaces[k]["data"][backgroundfield])
            if np.count_nonzero(a)>4:
                lons = np.asarray(surfaces[k]["lons"])[a]
                lats = np.asarray(surfaces[k]["lats"])[a]
                xpv,ypv = mapy(lons,lats)
                bgfield = np.asarray(surfaces[k]["data"][backgroundfield])[a]
                if contour:
                    plt.tricontourf(xpv,ypv,bgfield,levels=50,cmap="viridis")
                else:
                    #print("no graph")
                    plt.scatter(xpv,ypv,c=bgfield,cmap="viridis")
                #plt.scatter(xpv,ypv,c=bgfield)
                plt.clim(boundfunc(bgfield,stdevs))
                plt.colorbar()
                mapy.quiver(x,y,u,v,mag,scale=scale,cmap="autumn",width = 0.004)
                plt.colorbar()
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
    for k in Bar("Graphing Surfaces: ").iter(surfaces.keys()):
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
    for k in Bar("Graphing Surfaces: ").iter(surfaces.keys()):
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


def meridionalTransport(surfaces,lat,startlon,endlon,startpres,endpres,show=False,label="",ax=None):
    transportsums = {}
    for k in surfaces.keys():
        if  startpres < int(k) < endpres:
            lons = []
            vabs = []
            height = []
            rhos = []
            for l in range(len(surfaces[k]["lats"])):
                if np.abs(surfaces[k]["lats"][l] - lat)<0.01 and  startlon< surfaces[k]["lons"][l] <endlon:
                    lons.append(surfaces[k]["lons"][l])
                    vabs.append(surfaces[k]["data"]["v"][l])
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
            #width = np.gradient(np.asarray(lons))*np.cos(np.deg2rad(lat))*111*10**3
            width =107438
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
    ax.plot(lons[s],np.nancumsum(finaltransport[s])*(10**-9),label=label)
    ax.legend()
    if show:
        plt.show()


def meridionalHeatMap(surfaces,lat,startlon,endlon,startpres,endpres,show=False,label="",ax=None):
    transportsums = {}
    lons = []
    vabs = []
    pres = []
    bottom=[]
    for k in surfaces.keys():
        if  startpres < int(k) < endpres:
            j=len(lons)
            height = []
            rhos = []
            for l in range(len(surfaces[k]["lats"])):
                if np.abs(surfaces[k]["lats"][l] - lat)<0.01 and  startlon< surfaces[k]["lons"][l] <endlon:
                    lons.append(surfaces[k]["lons"][l])
                    vabs.append(surfaces[k]["data"]["v"][l])
                    pres.append(-surfaces[k]["data"]["pres"][l])
                    height.append(surfaces[k]["data"]["h"][l])
                    rho = gsw.rho(surfaces[k]["data"]["s"][l],surfaces[k]["data"]["t"][l],surfaces[k]["data"]["pres"][l])
                    rhos.append(rho)
                    bottom.append(surfaces[k]["data"]["z"][l])
            plt.plot(lons[j:],pres[j:],c="gray",zorder=0)
    lons = np.asarray(lons)
    bottom = np.asarray(bottom)
    pres = np.asarray(pres)
    vabs = np.asarray(vabs)
    s= np.argsort(lons)
    bottom = bottom[s]
    lons = lons[s]
    pres = pres[s]
    vabs=vabs[s]
    plt.fill_between(lons,-5000,bottom,color="black")
    plt.xlabel("Longitude")
    plt.ylabel("Pressure (dbar)")
    plt.scatter(lons,pres,s=(np.abs(vabs))*2*10**5,c=vabs,vmax =0.007,vmin=-0.007,zorder=5,cmap="coolwarm")
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("North-South Velocity")
    plt.show()


def meridionalSurfaces(surfaces,lat,startlon,endlon,startpres,endpres,show=False,label="",ax=None):
    transportsums = {}
    lons = []
    vabs = []
    pres = []
    bottom=[]
    for k in surfaces.keys():
        if  startpres < int(k) < endpres:
            j=len(lons)
            height = []
            rhos = []
            for l in range(len(surfaces[k]["lats"])):
                if np.abs(surfaces[k]["lats"][l] - lat)<0.01 and  startlon< surfaces[k]["lons"][l] <endlon:
                    lons.append(surfaces[k]["lons"][l])
                    vabs.append(surfaces[k]["data"]["v"][l])
                    pres.append(-surfaces[k]["data"]["pres"][l])
                    height.append(surfaces[k]["data"]["h"][l])
                    rho = gsw.rho(surfaces[k]["data"]["s"][l],surfaces[k]["data"]["t"][l],surfaces[k]["data"]["pres"][l])
                    rhos.append(rho)
                    bottom.append(surfaces[k]["data"]["z"][l])
            plt.plot(lons[j:],pres[j:],zorder=0)
    lons = np.asarray(lons)
    bottom = np.asarray(bottom)
    pres = np.asarray(pres)
    vabs = np.asarray(vabs)
    s= np.argsort(lons)
    bottom = bottom[s]
    lons = lons[s]
    pres = pres[s]
    vabs=vabs[s]
    plt.fill_between(lons,-5000,bottom,color="black")
    plt.xlabel("Longitude")
    plt.ylabel("Pressure (dbar)")
    plt.show()




def HTtransports(inv):
    fig,(ax1,ax2,ax3) = plt.subplots(3,1)
    graph.meridionalTransport(inv,-30,-180,180,0,819,show=False,label="1-5",ax=ax1)
    #ax1.set_ylim(-60,30)
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Mass Transport in 10^9 kg*s^-1")
    ax1.set_title("Accumulated mass transport in layers 1-6")
    graph.meridionalTransport(inv,-30,-180,180,819,2150,show=False,label="6-9",ax=ax2)
    #ax2.set_ylim(-50,20)
    ax2.set_ylabel("Mass Transport in 10^9 kg*s^-1")
    ax2.set_xlabel("Longitude")
    ax2.set_title("Accumulated mass transport in layers 6-9")
    graph.meridionalTransport(inv,-30,-180,180,2150,100000,show=False,label="9-",ax=ax3)
    #ax3.set_ylim(-10,15)
    ax3.set_ylabel("Mass Transport in 10^9 kg*s^-1")
    ax3.set_xlabel("Longitude")
    ax3.set_title("Accumulated mass transport in layers 9 to the bottom")
    plt.show()
                
##plot ts diagrams with neutral surfaces annotated
def tsNeutralExploreProfiles(profiles):
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


def tsNeutralExploreSurfaces(surfaces):
    for k in surfaces.keys():
        plt.plot(surfaces[k]["data"]["s"],surfaces[k]["data"]["t"])
    plt.show()

## graph a vector field given a surfaces object on a map
## any quantity can be supplied as a background field
def graphProfilesVectorField(region,profiles,depths=range(200,4000,200),savepath=False,show=True):

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

        urs.append(0.01)
        uthetas.append(0)
        lons.append(region.mapbounds["lon_0"])
        lats.append(region.mapbounds["lat_0"])
        bgfield.append(np.nan)

        urs.append(0)
        uthetas.append(0.01)
        lons.append(region.mapbounds["lon_0"])
        lats.append(region.mapbounds["lat_0"])
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

def saveAllQuants(region,surfaces,savepath,select=range(0,100000),contour=False,secondsurface=None):
    for d in Bar("Saving Fields").iter(surfaces[list(surfaces.keys())[0]]["data"].keys()):
        if len(surfaces[list(surfaces.keys())[0]]["data"][d])>0 and (~np.isnan(surfaces[list(surfaces.keys())[0]]["data"][d])).any():
            graphSurfaces(region,surfaces,d,show=False,savepath=savepath,select=select,contour=contour,secondsurface=secondsurface)


def northSouthTransect(surfaces,quant,lat=None,lon=None,\
        show=False,savepath=None,maximize=True,\
        presrange=[1000,np.inf],cbarrange=None,coordrange=None):
    xs = []
    ps = []
    cs = []
    zs = []
    for k in surfaces.keys():
        for i in range(len(surfaces[k]["data"][quant])):
            if presrange[0]<surfaces[k]["data"]["pres"][i]<presrange[1]:
                if lat and abs(surfaces[k]["lats"][i]-lat)<0.5:
                    xs.append(surfaces[k]["lons"][i])
                    cs.append(surfaces[k]["data"][quant][i])
                    ps.append(surfaces[k]["data"]["pres"][i])
                    zs.append(surfaces[k]["data"]["z"][i])
                if lon and abs(surfaces[k]["lons"][i]-lon)<0.5:
                    xs.append(surfaces[k]["lats"][i])
                    cs.append(surfaces[k]["data"][quant][i])
                    ps.append(surfaces[k]["data"]["pres"][i])
                    zs.append(surfaces[k]["data"]["z"][i])
    m = np.nanmean(cs)
    s = np.nanstd(cs)
    #print(ps)
    srt = np.argsort(xs)
    plt.fill_between(np.asarray(xs)[srt],-np.asarray(zs)[srt],6000)
    plt.scatter(xs,ps, c=cs)
    if coordrange:
        plt.gca().set_xlim(coordrange)
    if cbarrange:
        plt.clim(cbarrange)
    else:
        plt.clim(m-2*s,m+2*s)
    if lat:
        plt.xlabel("Longitude in Degrees")
    if lon:
        plt.xlabel("Latitude in Degrees ")
    plt.ylabel("Depth in Dbar")
    plt.gca().invert_yaxis()
    #plt.gca().set_xlim([-38,-16])
    cbar = plt.colorbar()
    if lat:
        degline = " along " + str(lat) + " latitude"
    if lon:
        degline = " along " + str(lon) + " longitude"

    if quant in quanttitlehash.keys():
        cbar.set_label(quanttitlehash[quant])
        plt.title(quanttitlehash[quant] + degline)
    else:
        cbar.set_label(quant)
        plt.title(quant + degline)

    if show:
        plt.show()

    if maximize:
        plt.gcf().set_size_inches(8,6)

    if savepath:
        try:
            os.makedirs(savepath)
        except FileExistsError as e:
            print(e)
        if lat:
            plt.savefig(savepath+"/lat"+str(lat)+".png")
        if lon:
            plt.savefig(savepath+"/lon"+str(lon)+".png")
    plt.close()
 

def time_diagnostic(profiles,layer,lat,tol):
    ## Plot change of values with time over neutral surface
    pressures = []
    times = []
    pvs = []
    lons = []
    for p in profiles:
        if abs(p.lat-lat) <tol and layer in p.neutraldepth:
            d = p.neutraldepth[layer]
            pv,g = p.potentialVorticityAtHautala(d)
            if pv:
                pressures.append(d)
                times.append(julian.from_jd(p.time).year)
                pvs.append(pv)
                lons.append(p.lon.data)
    fig,(ax1,ax2) = plt.subplots(1,2)
    fig.suptitle(" {} decibar neutral surface from {} to {} latitude".format(layer,lat-tol,lat+tol))
    c = ax1.scatter(lons,pressures,c=times)
    ax1.set_xlabel("Longitude")
    ax2.set_xlabel("Longitude")
    ax1.set_ylabel("Pressure")
    ax2.set_ylabel("Potential Vorticity")
    plt.colorbar(c,ax=ax1)
    c = ax2.scatter(lons,np.log10(-np.asarray(pvs)),c=times)
    plt.colorbar(c,ax=ax2)
    plt.savefig("/home/garrett/Projects/arcticcirc-pics/20.png")
    plt.show()

def layerChooser(profile):
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    ax1.plot(profile.isals,-np.asarray(profile.ipres))
    ax2.axvline(x=1031.53,c="g")
    ax2.axvline(x=1031.98,c="g")
    ax2.axvline(x=1032.30,c="yellow")
    ax3.plot(profile.itemps,-np.asarray(profile.ipres))
    ax3.axvline(x=2,c="blue")
    ax1.plot(gsw.SA_from_SP([34.86]*len(profile.ipres),profile.ipres,profile.lon,profile.lat),-np.asarray(profile.ipres),c="blue")

    ax2.plot(gsw.rho(profile.isals,profile.itemps,1000),-np.asarray(profile.ipres))
    plt.show()
    #ax3.plot(profile.isals,profile.ipres)
