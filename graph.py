from nstools import *

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

def graphSurfaces(surfaces,quantindex,contour=False,profiles=None,deepestindex=None,show=True,maximize=True,savepath=None):
    quanttitlehash = {"pres":"Pressure Dbar","t":"Temperature C","s":"Salinity PSU","pv":"PV","uz":"Uz'","vz":"Vz'","psi":"ISOPYCNAL STREAMFUNCTION"}
    for i in surfaces.keys():
        if len(surfaces[i]["lons"])>3 and len(surfaces[i]["data"][quantindex])>3:
            fig,ax = plt.subplots(1,1)
            mapy = Basemap(projection='ortho', lat_0=90,lon_0=-60)
            mapy.drawmapboundary(fill_color='aqua')
            mapy.fillcontinents(color='coral',lake_color='aqua')
            mapy.drawcoastlines()
            x,y = mapy(surfaces[i]["lons"],surfaces[i]["lats"])
            #Plot the surface 
            if contour:
                plt.tricontourf(x,y,np.asarray(surfaces[i]["data"][quantindex]),cmap="plasma")
            else:
                plt.scatter(x,y,c=np.asarray(surfaces[i]["data"][quantindex]),cmap="plasma")
                m = np.median(np.asarray(surfaces[i]["data"][quantindex]))
                s = np.std(np.asarray(surfaces[i]["data"][quantindex]))
                plt.clim(m-2*s,m+2*s)
                #plt.clim(i-400,i+400)
                mapy.colorbar()
            #map the reference profile
            if profiles and deepestindex:
                x,y = mapy(profiles[deepestindex].lon,profiles[deepestindex].lat)
                mapy.scatter(x,y,c="red")
            zoomGraph(mapy,ax)
            fig.suptitle(str(quanttitlehash[quantindex]) + " at NS: "+str(i))
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
