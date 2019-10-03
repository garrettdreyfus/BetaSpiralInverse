import nstools
import numpy as np
import inverttools
import eccotools
import interptools
import bathtools
import saloffset
import graph
import matplotlib.pyplot as plt
import parametertools as ptools
import pickle
import sensitivity
import random

with open('data/nepbctdprofiles.pickle', 'rb') as outfile:
    profiles = pickle.load(outfile)


#deepestindex=nstools.deepestProfile(profiles)
##54.062°N, 157.477°W,
##profilechoice = nstools.profileInBox(profiles,-180,-140,30,40,4000)[0]
profilechoice = nstools.getProfileById(profiles,178)

surfaces = nstools.runPeerSearch(profiles,100,4000,200,profilechoice,10**10)


with open('data/annotatednepbprofilessingleref.pickle', 'wb') as outfile:
    pickle.dump([surfaces,profiles],outfile)

with open('data/annotatednepbprofilessingleref.pickle', 'rb') as outfile:
    surfaces,profiles = pickle.load(outfile)

for l in profiles:
    l.knownns = l.knownns[()]
#graph.graphSurfaces(surfaces,"pres",region="nepb",savepath="refpics/surfaces/ecconepbns/",show=False)
##graph.tsNeutralExplore(profiles)
##graph.plotProfiles(profiles,"first ns search",specialprofile=profilechoice,region="nepb")
surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
graph.graphSurfaces(surfaces,"nsdiff",region="nepb")
#with open('data/nepbsurfaceswithdata.pickle', 'wb') as outfile:
    #pickle.dump(surfaces, outfile)


with open('data/nepbsurfaceswithdata.pickle', 'rb') as outfile:
    surfaces = pickle.load(outfile)



############################################
#graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewref/",region="nepb")
#surfaces,neighbors,distances = interptools.interpolateSurfaces(surfaces,fixedgrid="nepb")
#graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewrefinterp/",region="nepb")
#graph.graphSurfaces(surfaces,"pres",region="nepb")
#with open('data/readytoaddparamsnepb.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)
#with open('data/readytoaddparamsnepb.pickle', 'rb') as outfile:
    #surfaces,neighbors,distances = pickle.load(outfile)
#surfaces = nstools.addParametersToSurfaces(surfaces,neighbors,distances)

##graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquants/",region="nepb")
#with open('data/ready4inversenepb.pickle', 'wb') as outfile:
    #pickle.dump([surfaces,neighbors,distances], outfile)

#with open('data/ready4inversenepb.pickle', 'rb') as outfile:
    #surfaces,neighbors,distances = pickle.load(outfile)

############################################
#graph.saveAllQuants(surfaces,"refpics/surfaces/nepballquantsnewref/",region="nepb")
#graph.graphSurfaces(surfaces,"pres",region="nepb")
#params = {"reflevel":800,"upperbound":800,"lowerbound":1400,"mixs":{"kvo":False,"kvb":False,"kh":False}}
#out= inverttools.invert("coupled",surfaces,neighbors,distances,params=params)

##sensitivity.conditionErrorRefLevel("coupled",surfaces,neighbors,distances)
#inv = nstools.streamFuncToUV(out["surfaces"],neighbors,distances)

###graph.graphSurfaces(inv,"e")
#graph.graphVectorField(inv,"uabs","vabs","pv",metadata=out["metadata"],\
        #region="nepb",savepath="refpics/vectorfields/8001400nomix/",show=False)



##graph.graphSurfaces(surfaces,"pres",region="nepb")
##graph.graphSurfaces(surfaces,"t",region="nepb")

