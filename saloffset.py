import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import numpy as np
from functools import partial
import nstools
from geopy.distance import geodesic
from matplotlib.widgets import Button


def submit(profiles,cruisename,refprofile,fig,ax,s,text):
    val = eval(text)
    ax.clear()
    for p in profiles:
        a = np.asarray(p.ipres)
        indexs = np.where(abs(a) >= 1000)[0]
        if len(indexs)>0:
            ax.plot(p.isals[indexs[0]:indexs[-1]]+val,p.itemps[indexs[0]:indexs[-1]],linewidth=0.5)

    a = np.asarray(refprofile.ipres)
    indexs = np.where(abs(a) >= 1000)[0]
    ax.plot(refprofile.isals[indexs[0]:indexs[-1]],refprofile.itemps[indexs[0]:indexs[-1]],color="red",linewidth=1.0)
    ax.set_xlabel("salinity") 
    ax.set_ylabel("temperature") 
    plt.draw()
    s.answer = val

def passCruise(s,event):
    s.answer = None
    plt.close('all')


def selectorGraph(cruiseprofiles,cruisename,refprofile):
    fig, (ax1,ax2) = plt.subplots(1,2)
    nstools.plotCruise(cruiseprofiles,cruisename,fig=fig,ax=ax2,show=False)
    print(refprofile.time,cruiseprofiles[0].time)
    plt.subplots_adjust(bottom=0.2)
    submit(cruiseprofiles,cruisename,refprofile,fig,ax1,selectorGraph,"0.0")
    axbox = plt.axes([0.1, 0.05, 0.2, 0.04])
    axpass = plt.axes([0.81, 0.05, 0.1, 0.075])
    bpass = Button(axpass, 'Pass')
    bpass.on_clicked(partial(passCruise,selectorGraph))
    text_box = TextBox(axbox, 'Evaluate', initial="0.0")
    text_box.on_submit(partial(submit,cruiseprofiles,cruisename,refprofile,fig,ax1,selectorGraph))
    #mng = plt.get_current_fig_manager()
    #mng.frame.Maximize(True)
    plt.show()
    return selectorGraph.answer
    
def closestRefSearchAverage(refprofiles,profiles):
    mindistance = 100000000000000000
    minprofile = None
    for r in refprofiles:
        distancesum = 0 
        for p in profiles:
            distancesum += geodesic((p.lat,p.lon),(r.lat,r.lon)).meters
        if distancesum<mindistance:
            minprofile = r
            mindistancesum=distancesum
    return minprofile

def closestRefSearch(refprofiles,profiles):
    mindistance = 100000000000000000
    minprofile = None
    for r in refprofiles:
        for p in profiles:
            distance = geodesic((p.lat,p.lon),(r.lat,r.lon)).meters
            if distance<mindistance:
                minprofile = r
                mindistance = distance
    return minprofile



def singleSalinityOffsetRun(filename,cruisename,refcruisename):
    profiles,deepestindex = nstools.extractProfilesMonths(filename,range(13))
    cruiseprofiles = nstools.cruiseSearch(profiles,cruisename,1994)
    refprofiles = nstools.cruiseSearch(profiles,refcruisename)
    nstools.plotCruiseAndRef(cruiseprofiles,refprofiles,False)
    offset = selectorGraph(cruiseprofiles,refcruisename,closestRefSearch(refprofiles,cruiseprofiles))
    return offset

def runSalinityOffsetTool(filename,refcruisenames):
    profiles,deepestindex = nstools.extractProfilesMonths(filename,range(13))
    referenceprofiles = []
    for name in refcruisenames:
        referenceprofiles+= nstools.cruiseSearch(profiles,name)
    offsetDict = {}
    cruisenames = nstools.cruiseCount(profiles)
    while len(cruisenames) >0:
        for name in cruisenames.copy():
            cruiseprofiles = nstools.cruiseSearch(profiles,name)
            offset = selectorGraph(cruiseprofiles,name,closestRefSearch(referenceprofiles,cruiseprofiles))
            if offset:
                offsetDict[name] = offset
                referenceprofiles += cruiseprofiles
                cruisenames.remove(name)

    return offsetDict


runSalinityOffsetTool("data/3000mprofiles.json",["HUDSON_HUDSON2","LOUIS_S._ST._LAURENT_18SN940"])
#singleSalinityOffsetRun("data/2000mprofiles.json","LOUIS_S._ST._LAURENT_18SN940","HUDSON_HUDSON2")
#profiles,deepestindex = nstools.extractProfilesMonths("data/2000mprofiles.json",range(13))
#nstools.plotCruise(nstools.cruiseSearch(profiles,"LOUIS_S._ST._LAURENT_18SN940",1994),"name")

         
#for name in nstools.transArcticSearch(profiles):
    #print(name)
    #nstools.plotCruise(nstools.cruiseSearch(profiles,name),name)
        
    

#print("DONE WITH EXTRACTING PROFILES")
#surfaces = search(profiles,deepestindex)
#print("DONE FINDING SURFACES")
#with open('data/surfaces.json', 'w') as outfile:
    #json.dump(surfaces, outfile)
##json_file = open("data/surfaces.json") 
##surfaces = json.load(json_file)
#print("NOW GRAPHING")
#graphSurfaces(surfaces)
