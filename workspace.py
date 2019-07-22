import nstools
import saloffset
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap
import numpy as np
import json
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyproj
import graph
import copy
from pprint import pprint
    

#offsets, profiles, deepestindex = saloffset.runSalinityOffsetTool(["data/1500mprofiles.json"],["ODEN_AGAVE"])

##profiles,deepestindex = nstools.extractProfilesBox(["data/1500mprofiles.json"],-180,180,65,90)
##profiles,deepestindex = nstools.removeNorwegianSea(profiles)
#nstools.plotCruise(profiles,"IPY 2007")
##print(offsets)
#profiles = nstools.filterCruises(profiles,offsets.keys())
#profiles = saloffset.applyOffsets(profiles,offsets)

#fileObject = open("1500NoNorwegian.pickle",'wb')  
## load the object from the file into var b
#b = pickle.dump([offsets,profiles,deepestindex],fileObject)  
#fileObject.close()


#USING CORRECTED PICKLED SURFACES TO MAP NEUTRAL SURFACES USING PEER SEARCH

#fileObject = open("1500NoNorwegian.pickle",'rb')  
#offsets,profiles,deepestindex = pickle.load(fileObject)
#fileObject.close()

#nstools.plotCruise(profiles,"IPY 2007")
#deepestindex = nstools.deepestProfile(profiles)
##profilechoice = random.choice(nstools.profileInBox(profiles,-180,180,60,90))
#profilechoice = profiles[deepestindex]
#surfaces = {}
#for d in range(200,4000,200)[::-1]:
    #print(d)
    #surfaces.update(nstools.peerSearch(profiles.copy(),deepestindex,d,profilechoice,1000))

#nstools.graphSurfaces(surfaces)

#fileObject = open(str(profilechoice.eyed)+".pickle",'wb')  
## load the object from the file into var b
#b = pickle.dump(surfaces,fileObject)  
#fileObject.close()

#####################################################################
#print(runSalinityOffsetTool(glob.glob("data/3000m2007profiles.json"),["Polarstern_ARK-XXIII_2"]))
#singleSalinityOffsetRun("data/2000mprofiles.json","LOUIS_S._ST._LAURENT_18SN940","HUDSON_HUDSON2")
#profiles,deepestindex = nstools.extractProfilesMonths("data/3000m2008profiles.json",range(13))
#nstools.plotCruise(nstools.cruiseSearch(profiles,"LOUIS_S._ST._LAURENT_18SN940",1994),"name")

#print("DONE WITH EXTRACTING PROFILES")
#surfaces = nstools.search(profiles,deepestindex)
#print("DONE FINDING SURFACES")
#with open('data/surfaces.json', 'w') as outfile:
    #json.dump(surfaces, outfile)
#json_file = open("data/surfaces.json") 
#surfaces = json.load(json_file)
#print("NOW GRAPHING")
#nstools.graphSurfaces(surfaces)
#######################################33333


#fileObject = open("data/286364.pickle",'rb')  
#surfaces = pickle.load(fileObject)
#fileObject.close()
fileObject = open("data/1500NoNorwegian.pickle",'rb')  
offsets,profiles,deepestindex = pickle.load(fileObject)
#nstools.findNeighboringPoints(profiles,profiles[2].lat,profiles[2].lon)
#nstools.plotASpiral(profiles)
fileObject.close()

#print("Everythings loaded up")
#surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
#with open('data/surfacesWithData.pickle', 'wb') as outfile:
    #pickle.dump(surfaces, outfile)
fileObject = open("data/surfacesWithData.pickle",'rb')  
surfaces = pickle.load(fileObject)
originalsurfaces = copy.deepcopy(surfaces)
surfaces = nstools.convertOldSurfaces(surfaces)
surfaces = nstools.addStreamFunc(surfaces,profiles)
#nstools.addStreamFuncFromFile(surfaces,profiles,"geoisopycnal.mat","np_vector.mat")
#graph.graphSurfaces(surfaces,"psi")

######################################################################################################
surfaces =nstools.addXYToSurfaces(surfaces)
#pprint(surfaces.keys())
#pprint(surfaces[3600])
###nstools.graphTransects(nstools.filterSurfacesByLine(originalsurfaces,40),0)
    ################
interpolatedsurfaces = {}
neighbors={}
lookups={}
#for k in surfaces.keys():
    #surfaces[k] = nstools.removeDiscontinuities(surfaces[k],radius=0.1)
    #interpolatedsurfaces[k],neighbors[k] = nstools.interpolateSurface(surfaces[k])
    #lookups[k] = nstools.trueDistanceLookup(interpolatedsurfaces[k],neighbors[k])

#surfaces = nstools.addHToSurfaces(interpolatedsurfaces)
#with open('data/lookupNeighborsSurfaces.pickle', 'wb') as outfile:
    #pickle.dump([interpolatedsurfaces,lookups,neighbors], outfile)

with open('data/lookupNeighborsSurfaces.pickle', 'rb') as outfile:
    [interpolatedsurfaces,lookups,neighbors]=pickle.load(outfile)

staggeredsurfaces = nstools.addPrimeToSurfacesCartesianTrueDistance(interpolatedsurfaces,neighbors,lookups)

#with open('data/ready4inverse.pickle', 'wb') as outfile:
    #pickle.dump(staggeredsurfaces,outfile)


#with open('data/ready4inverse.pickle', 'rb') as outfile:
    #staggeredsurfaces=pickle.load(outfile)

staggeredsurfaces = nstools.invert("salt",staggeredsurfaces)

#with open('data/svdinverted.pickle', 'wb') as outfile:
    #pickle.dump(staggeredsurfaces, outfile)

#with open('data/svdinverted.pickle', 'rb') as outfile:
    #staggeredsurfaces=pickle.load(outfile)

#graph.graphSurfaces(staggeredsurfaces,"u",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"v",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"uabs",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"vabs",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"uprime",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"vprime",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"hx",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"hy",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"h",show=False,savepath = "refpics/ref1000/")
#graph.graphSurfaces(staggeredsurfaces,"uabs")
##graph.graphSurfaces(staggeredsurfaces,"vabs")
graph.graphVectorField(staggeredsurfaces,"u","v",show=False,savepath="refpics/SaltInvertR1000/")
graph.graphVectorField(staggeredsurfaces,"uabs","vabs",show=False,savepath="refpics/SaltInvertR1000/")
#graph.graphVectorField(staggeredsurfaces,"uabs","vabs")
#graph.graphVectorField(staggeredsurfaces,"u","v")
#graph.graphVectorField(staggeredsurfaces,"uabs","vabs")
#graph.graphSurfaces(staggeredsurfaces,"v",show=False,savepath = "refpics/staggeredDerivatives/")
#graph.twentyRandomSpirals(staggeredsurfaces)

#interpolatedsurfaces = nstools.createStaggeredSurface(interpolatedsurfaces,neighbors)
#graph.graphSurfacesComparison(interpolatedsurfaces,surfaces,"uprime",show=False,savepath="refpics/interpPSI/")
