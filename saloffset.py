import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import numpy as np
from functools import partial
import nstools
from geopy.distance import geodesic
import json
import glob
from matplotlib.widgets import Button


def submit(profiles,cruisename,refprofile,fig,ax,s,text):
    val = eval(text)
    ax.clear()
    ts = np.asarray([])
    sals= np.asarray([])
    lats =np.asarray([])
    for p in profiles:
        a = np.asarray(p.ipres)
        indexs = np.where(abs(a) >= 2000)[0]
        if len(indexs)>0:
            #sals=np.concatenate([(p.isals[indexs[0]:indexs[-1]]+val),sals])
            #ts=np.concatenate([ts,(p.itemps[indexs[0]:indexs[-1]])])
            #lats=np.concatenate([lats,([p.lat]*(indexs[-1]-indexs[0]))])
            ax.plot(p.isals[indexs[0]:indexs[-1]]+val,p.itemps[indexs[0]:indexs[-1]],linewidth=0.5)
    cm = ax.scatter(sals,ts,c=lats,s=0.1)

    a = np.asarray(refprofile.ipres)
    indexs = np.where(abs(a) >= 2000)[0]
    ax.scatter(refprofile.isals[indexs[0]:indexs[-1]],refprofile.itemps[indexs[0]:indexs[-1]],c=[refprofile.lat]*(indexs[-1]-indexs[0]),s=10,marker="o")
    ax.plot(refprofile.isals[indexs[0]:indexs[-1]],refprofile.itemps[indexs[0]:indexs[-1]])
    ax.set_xlabel("salinity") 
    ax.set_ylabel("temperature") 
    plt.draw()
    s.answer = val

def passCruise(s,event):
    s.answer = None
    plt.close('all')

def noMixingLine(s,event):
    s.answer = -999
    plt.close('all')




def selectorGraph(cruiseprofiles,cruisename,refprofile):
    fig, (ax1,ax2) = plt.subplots(1,2)
    selectorGraph.answer=0
    nstools.plotCruise(cruiseprofiles,cruisename,fig=fig,ax=ax2,show=False)
    plt.subplots_adjust(bottom=0.2)
    submit(cruiseprofiles,cruisename,refprofile,fig,ax1,selectorGraph,"0.0")
    axbox = plt.axes([0.1, 0.05, 0.2, 0.04])
    axpass = plt.axes([0.81, 0.05, 0.3, 0.075])
    axno = plt.axes([0.5, 0.05, 0.3, 0.075])
    bpass = Button(axpass, 'Mixing Line but Bad Ref')
    bno = Button(axno, 'No Clear Mixing Line')
    bpass.on_clicked(partial(passCruise,selectorGraph))
    bno.on_clicked(partial(noMixingLine,selectorGraph))
    text_box = TextBox(axbox, 'Evaluate', initial="0.0")
    text_box.on_submit(partial(submit,cruiseprofiles,cruisename,refprofile,fig,ax1,selectorGraph))
    #mng = plt.get_current_fig_manager()
    #mng.frame.Maximize(True)
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
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
    offset = selectorGraph(cruiseprofiles,refcruisename,closestRefSearchAverage(refprofiles,cruiseprofiles))
    return offset

def runSalinityOffsetTool(filenames,refcruisenames,box=False,
        lonleft=False,lonright=False,latbot=False,lattop=False):

    if box and lonleft and lonright and lattop and latbot:
        profiles,deepestindex = nstools.extractProfilesBox(filenames,lonleft,lonright,latbot,lattop)
    else:
        profiles,deepestindex = nstools.extractProfilesMonths(filenames,range(13))

    cruisenames = nstools.cruiseCount(profiles)
    print(cruisenames)
    referenceprofiles = []
    offsetDict = {}
    for name in refcruisenames:
        referenceprofiles+= nstools.cruiseSearch(profiles,name)
        offsetDict[name]=0.0
    while len(cruisenames) >0:
        for name in cruisenames.copy():
            cruiseprofiles = nstools.cruiseSearch(profiles,name)
            if len(cruiseprofiles)>10:
                refprofile = closestRefSearch(referenceprofiles,cruiseprofiles)
                titlestring = name + " " +str(cruiseprofiles[0].time)+"\n REF:" + refprofile.cruise + " " + str(refprofile.time)
                offset = selectorGraph(cruiseprofiles,titlestring,refprofile)
                if offset != None:
                    if abs(offset)<10:
                        offsetDict[name] = offset + offsetDict[refprofile.cruise]
                        referenceprofiles += cruiseprofiles
                        print("OFFSET ENTERED")
                    else:
                        print("BAD PROFILE FLAGGED")

                    cruisenames.remove(name)
                else:
                    print("PASS FOR NOW")
            else:
                cruisenames.remove(name)
        print("ONE PASS DONE")
        print(len(cruisenames))


    return offsetDict, profiles,deepestindex

def applyOffsets(profiles,offsetDict):
    for p in profiles:
        if p.cruise in offsetDict.keys():
            p.applyOffset(offsetDict[p.cruise])
    return profiles

